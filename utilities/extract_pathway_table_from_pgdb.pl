#!/usr/bin/perl -w

# Extracts information about pathways and compounds from a Pathway/Genome
# Database (PGDB). Outputs a tsv file based on the following schema:
#
# (insert a schema here)
#
# Requires the PerlCyc module from the SOL Genomics Network
# (ftp://ftp.sgn.cornell.edu/programs/perlcyc/).


package perlcyc;
$VERSION="1.21";
use Socket;
use strict;
use Carp;


sub new {
    
    my $class = shift;
    my $organism =shift;
    my $debug = shift;

    if (!$organism) { 
	carp("Warning! No organism has been selected in the new call.\n");
    }
    # initialize object's variable hash
    my $vars = {};
    $vars->{_organism} = $organism;
    $vars->{_socket_name} = "/tmp/ptools-socket";
    $vars->{_debug} = $debug;
    $vars->{_socket} = undef;
    return bless $vars, $class;
}


sub makeSocket {
    my $self = shift;
    # get socket connection with pathway-tools
    socket(SOCK, PF_UNIX, SOCK_STREAM, 0) || die "socket: $!";
    connect (SOCK, sockaddr_un($self->{_socket_name})) || die "connect: $!";
    my $old_fh = select(SOCK); $|=1; select($old_fh);
    $self->{_socket} = \*SOCK;
}

sub cclose {
    my $self = shift;
    cclose($self->{_socket});
}

# the old lisp parser -- now deprecated
sub parse_lisp_list {
    my ($list_expr) = @_;
    my @results;
    while ($list_expr =~ m/
                           \|(.+?)\|      # match words delimited by pipes
                           |              # or
                           ([\w-]+)       # match a single word
	   /xg) {
	push @results, $1 || $2;
    }
    return @results;
}

## Given a lisp expression, returns a Perl list, where every inner
## Lisp list is converted into an equivalent Perl arrayref.
sub parselisp {
    my ($self, $lisp_string) = @_;
    my @tokens = tokenize($lisp_string);
    my $parsed_expr = parseExpr(\@tokens);

    if (wantarray && ref($parsed_expr)) {
	## Shallow flatten the list by one level if we're in list context.
	return @$parsed_expr;
    }
    else {
	return $parsed_expr;
    }
}

######################################################################
# Danny's Lisp parsing routines.
######################################################################

## Parses an lisp expression.
sub parseExpr {
    my ($token_ref) = @_;
    if (!(@$token_ref)) { return undef; }
    my $first = $token_ref->[0];
    if ($first eq '(') {
	shift @$token_ref;
	my @list_elements;
	while ($token_ref->[0] ne ')') {
	    push @list_elements, parseExpr($token_ref);
	}
	shift @$token_ref;
	return [@list_elements];
    }
    elsif ($first eq 'NIL') {
	shift (@$token_ref);
	return [];
    }
    else {
	return shift(@$token_ref);
    }
}
######################################################################

sub tokenize {
     my ($s) = @_;
     my $LPAREN = '\(';
     my $RPAREN = '\)';
     my $WSPACE = '\s+';
     my $STRING = '\"[^"]*?\"';
     my $PIPES = '\|[^\|]*?\|';
     my @tokens = grep {$_} split /($LPAREN|$RPAREN|$STRING|$PIPES)|$WSPACE/, $s;
     @tokens = map { $_ =~ s/\"([^\"]*)\"/$1/;  $_ } @tokens;  ## removes outer quotes from the string.
## note: we must not wipe out bars because they're used in frame ids.  Do not strip them out.
#     @tokens = map { $_ =~ s/\|([^|]*)\|/$1/;  $_ } @tokens;  ## removes outer pipes from the string.
     return @tokens;
}
#############
#
#########################################################################

sub send_query {
    # sends a query to server
    my $self = shift;
    my $query = shift;
    $self->makeSocket();
    $self-> debug_off();
    $self -> debug("Now sending query: $query");
    $self-> debug_off();
    my $s = $self->{_socket};
    print $s "$query\n";
}       

sub retrieve_results {
    # retrieves a result from the server after a query has been sent and parses it into an array.
    my $self = shift;
    my $s = $self->{_socket};
    my @results;
    while (<$s>)
    {
	push @results, $_;
    }
    return $self->parselisp(join("\n", @results));
}    

sub retrieve_results_string {
    my $self = shift;
    my $s = $self->{_socket};
    my @results;
    while (<$s>) {
      chomp;
      push @results, $_;
    }
    if ($results[0] eq "NIL") { $results[0]=""; }
    $results[0]=~s/\"(.*)\"/$1/;
    return ($results[0]);
}

sub retrieve_results_array {
    my $self = shift;
    my $s = $self->{_socket};
    my @results;
    while (<$s>) {
      chomp;
      push @results, $_;
    }
    return (@results);
}


sub wrap_query {
    my $self = shift;
    my $s = shift;
    return "(with-organism (:org-id\'".$self->{_organism}.") (mapcar #\'object-name ($s)))";
}

sub call_func {
    # used for functions that return a list. The list is converted to a list of object-names
    # and the resulting list is parsed into an array.
    #
    my $self = shift;
    my $function = shift;
##    print $self->wrap_query($function), "\n"; ## debug
    $self->send_query ($self->wrap_query($function));
    my @result = $self-> retrieve_results();

    $self -> debug_off();
    if ($self->debug) { 
      foreach my $r (@result) {
	print "READ TOKEN: $r\n";
      }
    }

    return @result;
}

## Similar to call_function, but forces scalar constext on retreive_results so that it maintains
## the inner list structure of 'result', array references and all.
sub call_func_with_structure {
    my $self = shift;
    my $function = shift;
    $self->send_query ($self->wrap_query($function));
    my $result = $self->retrieve_results();
    return $result;
}


# Caution! The following function is deprecated and doesn't work.
# Wrap your query instead in a multiple-value-list call.
sub call_func_that_returns_two_values {
  # used for functions that return two values (lisp lists). The lists 
  # are converted to two arrays, and two pointers to arrays are returned.
  my $self = shift;
  my $function = shift;
  $self->debug("Sending query...\n");
  $self -> send_query($self->wrap_query($function));
  $self->debug("QUERY SENT: $function\n");
  my ($res1, $res2) = $self -> retrieve_results_array();
  $self->debug("Results: $res1, $res2\n");
  my ($list1, $list2) = $self -> parse_lisp_list($res1);
  $self->debug("LISTS: $list1\n$list2\n\n");
  my (@res1) = $self -> parse_lisp_list($list1);
  my (@res2) = $self -> parse_lisp_list($list2);
  return (\@res1, \@res2);
}

sub call_func_that_returns_string {
    # use for functions that will return a string and not a list. 
    # this function doesn't call mapcar and doesn't parse the returned list.
    #
    my $self = shift;
    my $function = shift;
   $self-> send_query ("(with-organism (:org-id\'".$self->{_organism}.") (object-name ($function)))");
  
    return $self-> retrieve_results_string();
}

sub call_func_that_returns_boolean {
   # call this function for functions that return a boolean.
  my $self = shift;
  my $function = shift;
  my $s = $self -> call_func_that_returns_string($function);
  if ($s eq "") { return ""; }
  else {
    return "T";
  }
}

# GFP functions
#


sub protectFrameName {
    my $frame = shift;
    if ($frame =~ m/^\|.*\|$/) { return $frame; } ## if already pipe protected, don't do anything.
    if ($frame =~ m/\s/) {
	return "|$frame|";
    }
    return $frame;
}


sub variants_of_pathway{ # 
    my $self = shift;
    my $pathway=shift;
    return $self -> call_func("variants-of-pathway \'$pathway");
}

sub get_slot_values { # frame, slot_name -> slot_value
    my $self = shift;
    my $frame =shift;
    $frame = protectFrameName($frame);
    my $slot_name = shift;
    return $self->call_func("get-slot-values \'$frame \'$slot_name");
}

sub get_slot_value { # frame, slot_name -> slot_value
    my $self = shift;
    my $frame =shift;
    $frame = protectFrameName($frame);
    my $slot_name = shift;
    return $self->call_func_that_returns_string("get-slot-value \'$frame \'$slot_name");
}

sub get_class_slot_slotvalue {
  my $self = shift;
  my $class = shift;
  my $slot_name = shift;
  my $slot_value=shift;
  return $self -> call_func("get-class-slot-slotvalue \'$class \'$slot_name \'$slot_value");
}

sub get_class_all_instances { 
    my $self = shift;
    my $class =shift;
    return $self->call_func("get-class-all-instances \'$class");
}

sub get_class_all_subs {  # added 5/2008
    my $self = shift;
    my $class = shift;
    return $self->call_func("get-class-all-subs \'$class");
}

sub instance_all_instance_of_p {
    my $self = shift;
    my $instance = shift;
    my $class = shift;
    return $self->call_func_that_returns_boolean("instance-all-instance-of-p \'$instance \'$class");
}

sub member_slot_value_p {
    my $self = shift;
    my $frame =shift;
    $frame = protectFrameName($frame);
    my $slot = shift;
    my $value = shift;
    return $self->call_func_that_returns_boolean("member-slot-value-p \'$frame \'$slot \'$value");
}

sub current_kb {
    my $self = shift;
    #this may not make sense because we prefix everything with select-organism.
    #maybe we should just return $self->{_organism}?
    #return $self->call_func("current-kb");
    return $self->{_organism};
}

sub put_slot_values {
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  my $slot = shift;
  my $values = shift;
  return $self -> call_func("put-slot-values \'$frame \'$slot \'$values");
}

sub put_slot_value {
  my $self = shift;
  my $frame =shift;
  $frame = protectFrameName($frame);
  my $slot =shift;
  my $value = shift;
  return $self -> call_func("put-slot-value \'$frame \'$slot \'$value");
}

sub add_slot_value {
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  my $slot = shift;
  my $value = shift;
  return $self -> call_func("add-slot-value \'$frame \'$slot \'$value");
}

sub replace_slot_value {
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  my $slot  = shift;
  my $old_value = shift;
  my $new_value = shift;
  return $self -> call_func("replace-slot-value \'$frame \'$slot \'$old_value \'$new_value");
}

sub remove_slot_value {
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  my $slot = shift;
  return $self -> call_func("remove-slot-value \'$frame \'$slot");
}

sub coercible_to_frame_p { 
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  return $self -> call_func_that_returns_boolean("coercible-to-frame-p \'$frame");
}


sub class_all_type_of_p {
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  my $instance = shift;
  return $self -> call_func_that_returns_boolean("class-all-type-of-p \'$frame \'$instance");
}

sub get_instance_direct_types {
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  return $self -> call_func ("get-instance-direct-types \'$frame");
}


## Warning: specialized because we're getting the list of elements in an order
## we didn't expect. -Danny 11/11/03
sub get_instance_all_types {
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  my $lisp_query = "
  (let* ((results (with-organism (:org-id 'ARA) 
                                                    (get-instance-all-types '$frame))))
           (mapcar #'get-frame-pretty-name results))";
	$self->send_query($lisp_query);
  my @results = $self->retrieve_results();
  return split_up_types(@results);
}

sub split_up_types {
    my @results = @_;
    my $in_frames = 0;
    my @groups;

    my @titles;
    for my $e (@results) {
	if ($e ne 'THINGS') {
	    push @titles, $e;
	}
	else {
	    last;
	}
    }
    @titles = reverse @titles;

##    print "DEBUG dyoo: ", join(", ", @results), "\n";
    my @row;
    for my $e (@results, 'THINGS') {
	if ($e eq 'FRAMES') {
	    $in_frames = 1;
	}
	elsif ($e eq 'THINGS') {
	    shift @row; ## remove first element "Generalized-reactions"
	    if (@row) {
		push @groups, [@row, (shift @titles)];
	    }
	    $in_frames = 0;
	    @row = ();
	}
	elsif ($in_frames) {
	    push @row, $e;
	}
    }
    return @groups;
}



#  return $self -> call_func_with_structure("get-instance-all-types \'$frame");
#}

sub get_frame_slots {
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  return $self -> call_func ("get-frame-slots \'$frame");
}

sub put_instance_types {
  my $self = shift;
  my $frame = shift;
  $frame = protectFrameName($frame);
  my $new_types = shift;
  return $self -> call_func("put-instance-types \'$frame \'$new_types");
}

sub save_kb {
  my $self = shift;
  return $self -> call_func("save-kb");
}

sub revert_kb {
  my $self = shift;
  return $self -> call_func("revert-kb");
}

sub find_indexed_frame {
  my $self =shift;
  my $datum = shift;
  my $class = shift;
  return $self -> call_func("multiple-value-list (find-indexed-frame \'$datum \'$class)");
}

sub gfp_get_value_annot { 
    my $self = shift;
    my $frame1 = shift;
    my $term1 = shift;
    my $frame2 = shift;
    my $term2 = shift;

    return $self->call_func_that_returns_string("gkb:get_value_annot \'$frame1 \'$term1 \'$frame2 \'$term2");


}

sub gfp_get_value_annots { 
    my $self = shift;
    my $frame1 = shift;
    my $term1 = shift;
    my $frame2 = shift;
    my $term2 = shift;

    return $self->call_func("gkb:get_value_annots \'$frame1 \'$term1 \'$frame2 \'$term2");
}



# Pathway-Tools internal lisp functions
#

sub select_organism {
    # this just sets a variable in the object and doesn't access the 
    # select-organism function in lisp. The selected organism is prefixed to
    # every query sent to the socket server.
    #
    my $self = shift;
    my $organism = shift;
    $self -> {_organism} = $organism;
}

sub all_pathways {    
    my $self = shift;
    my $optional = shift || ":all";
    return $self-> call_func("all-pathways $optional");
}






## Converts a Perl boolean value to the equivalent Lisp value.
sub perl_boolean_to_lisp {
    my ($value) = @_;
    if ($value) {
	return "T";
    } else {
	return "NIL";
    }
}

## Converts a Perl value to NIL if necessary.
sub perl_undef_to_nil {
    my ($value) = @_;
    if (defined($value)) { return $value; }
    else {
	return "NIL";
    }
}

## Implements the GET-NAME-STRING wrapper.
sub get_name_string {
    my $self = shift;
    my $frame = shift;
    $frame = protectFrameName($frame);
    my %args = (
		rxn_eqn_as_name => 1,
		direction => undef,
		name_slot => undef,
		strip_html => undef,
		@_
		);
    return $self->call_func_that_returns_string("get-name-string \'$frame"
			    . ( " :rxn-eqn-as-name? " 
				. perl_boolean_to_lisp($args{rxn_eqn_as_name}))
			    . ( "  :direction " 
				. perl_undef_to_nil($args{direction}))
			    . ( "  :name-slot "
				. perl_undef_to_nil($args{name_slot}))
			    . ( "  :strip-html? "
				. perl_boolean_to_lisp($args{strip_html})));
}


sub all_orgs {
    my $self=shift;
    return $self -> call_func("all-orgs"); 
}

sub all_rxns {
    my $self = shift;
    return $self -> call_func("all-rxns");
}

sub genes_of_reaction {
    my $self = shift;
    my $reaction=shift;
    return $self -> call_func("genes-of-reaction \'$reaction");
}

sub substrates_of_reaction {
    my $self = shift;
    my $reaction=shift;
    return $self -> call_func("substrates-of-reaction \'$reaction");
}

#hypothetical function that may not exist
sub products_of_reaction {
  my $self = shift;
  my $reaction = shift;
  return $self -> call_func("products-of-reaction \'$reaction");
}

sub enzymes_of_reaction {
    my $self = shift;
    my $reaction=shift;
    return $self -> call_func("enzymes-of-reaction \'$reaction");
}

sub reaction_reactants_and_products {
    my $self = shift;
    my $reaction=shift;
    my $pathway=shift;
    return $self -> call_func("reaction-reactants-and-products \'$reaction \'$pathway");
}

sub get_predecessors {
    my $self = shift;
    my $reaction=shift;
    my $pathway=shift;
    return $self -> call_func("get-predecessors \'$reaction \'$pathway");
}

sub get_successors {
    my $self = shift;
    my $reaction=shift;
    my $pathway = shift;
    return $self -> call_func("get-successors \'$reaction \'$pathway");
}

sub get_reaction_list {
  my $self = shift;
  my $pathway=shift;
  return $self-> call_func("get-reaction-list \'$pathway");
}

sub genes_of_pathway {
    my $self = shift;
    my $pathway=shift;
    return $self -> call_func("genes-of-pathway \'$pathway");
}

sub enzymes_of_pathway {
    my $self = shift;
    my $pathway = shift;
    return $self -> call_func("enzymes-of-pathway \'$pathway");
}

sub compounds_of_pathway {
    my $self = shift;
    my $pathway = shift;
    return $self -> call_func("compounds-of-pathway \'$pathway");
}

sub substrates_of_pathway {
    my $self = shift;
    my $pathway = shift;
    return $self -> call_func("multiple-value-list (substrates-of-pathway \'$pathway)");
}

sub all_transcription_factors {
    my $self = shift;
    return $self -> call_func("all-transcription-factors");
}

sub transcription_factor_p {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func_that_returns_boolean("transcription-factor? \'$protein");
}

sub all_cofactors {
    my $self = shift;
    return $self -> call_func("all-cofactors");
}

sub all_modulators {
    my $self = shift;
    return $self -> call_func("all-modulators");
}

sub monomers_of_protein {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("monomers-of-protein \'$protein");
}

sub components_of_protein {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("multiple-value-list (components-of-protein \'$protein)");
}

sub genes_of_protein {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("genes-of-protein \'$protein");
}

sub reactions_of_enzyme {
    my $self = shift;
    my $enzyme = shift;
    return $self -> call_func("reactions-of-enzyme \'$enzyme");
}

sub enzyme_p {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func_that_returns_boolean("enzyme? \'$protein");
}

sub transport_p {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func_that_returns_boolean("transporter? \'$protein");
}

sub containers_of {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("containers-of \'$protein");
}

sub modified_forms {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("modified-forms \'$protein");
}

sub modified_containers {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("modified-containers \'$protein");
}

sub top_containers {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("top-containers \'$protein");
}

sub reactions_of_protein {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("reactions-of-protein \'$protein");
}

sub regulon_of_protein {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("regulon-of-protein \'$protein");
}

sub transcription_units_of_protein {
    my $self = shift;
    my $protein = shift;
    return $self -> call_func("transcription-units-of-protein \'$protein");
}

sub regulator_proteins_of_transcription_unit {
    my $self = shift;
    my $tu = shift;
    return $self -> call_func("regulator-proteins-of-transcription-unit \'$tu");
}

sub enzymes_of_gene {
    my $self = shift;
    my $gene = shift;
    return $self -> call_func("enzymes-of-gene \'$gene");
}

sub genes_in_same_operon {
    my $self = shift;
    my $gene = shift;
    return $self -> call_func("genes-in-same-operon \'$gene");
}

sub all_products_of_gene {
    my $self = shift;
    my $gene = shift;
    return $self -> call_func("all-products-of-gene \'$gene");
}

sub reactions_of_gene {
    my $self = shift;
    my $gene = shift;
    return $self -> call_func("reactions-of-gene \'$gene");
}

sub pathways_of_gene {
    my $self = shift;
    my $gene = shift;
    return $self -> call_func("pathways-of-gene \'$gene");
}

sub chromosome_of_gene {
    my $self = shift;
    my $gene = shift;
    return $self -> call_func_that_returns_string("chromosome-of-gene \'$gene");
}

sub transcription_units_of_gene {
    my $self = shift;
    my $gene = shift;
    return $self -> call_func("transcription-units-of-gene \'$gene");
}

sub transcription_unit_promoter {
    my $self = shift;
    my $tu = shift;
    return $self -> call_func_that_returns_string("transcription-unit-promoter \'$tu");
}

sub transcription_unit_genes {
    my $self = shift;
    my $tu = shift;
    return $self -> call_func("transcription-unit-genes \'$tu");
}

sub transcription_unit_binding_sites {
    my $self = shift;
    my $tu = shift;
    return $self -> call_func("transcription-unit-binding-sites \'$tu");
}

sub transcription_unit_transcription_factors {
    my $self = shift;
    my $bsite = shift;
    return $self -> call_func("transcription-unit-transcription-factors \'$bsite");
}

sub transcription_unit_terminators {
    my $self = shift;
    my $tu = shift;
    return $self -> call_func("transcription-unit-terminators \'$tu");
}

sub all_transported_chemicals {
    my $self = shift;
    return $self -> call_func("all-transported-chemicals");
}

sub reactions_of_compound {
    my $self = shift;
    my $compound = shift;
    return $self -> call_func("reactions-of-compound \'$compound");
}

sub full_enzyme_name {
    my $self = shift;
    my $enzyme = shift;
    return $self -> call_func_that_returns_string("full-enzyme-name \'$enzyme");
}

sub enzyme_activity_name {
    my $self = shift;
    my $enzyme = shift;
    my $reaction = shift;
    if (! $reaction) {
	return $self -> call_func_that_returns_string("enzyme-activity-name \'$enzyme");
    } else {
	return $self -> call_func_that_returns_string("enzyme-activity-name \'$enzyme \'$reaction");
    }
}

sub create_instance { 
    my $self = shift;
    my $name = shift;
    my @types = @_;
    my $types = undef;
    foreach my $t (@types) { 
	$types .= " \'$t";
    }
    if (@types > 1) { "(".$types.")"; }
    
    return $self -> call_func_that_returns_string("create-instance \'$name $types");
}

sub create_class {
    my $self = shift;
    my $name = shift;
    my @direct_supers = @_;
    my $supers = undef;
    foreach my $ds (@direct_supers) { 
	$supers .= " \'$ds";
    }
    if (@direct_supers > 1) { 
	$supers = "(".$supers.")";
    }
    
    return $self -> call_func_that_returns_string("create-class \'$name $supers");
}

# create-frame  name &key direct-types direct-supers doc template-slots template-facets own-slots own-facets kb error-p
sub create_frame {
    my $self = shift;
    my $name = shift;
    my @types = @_;
    my $types = undef;

    foreach my $t (@types) { 
	$types .= " :direct-types \'$t";
    }
    return $self -> call_func_that_returns_string("create-frame \'$name $types");
}

 
# added 5/2008. These functions are only supported 
# in Pathway Tools version 12 onwards.
#
sub pwys_of_organism_in_meta { 
    my $self = shift;
    my $gb_tax_id = shift;
    return $self->call_func("pwys-of-organism-in-meta \'$gb_tax_id");
}

sub enzymes_of_organism_in_meta { 
    my $self = shift;
    my $gb_tax_id = shift;
    return $self->call_func("enzymes-of-organism-in-meta \'$gb_tax_id");
}

sub lower_taxa_or_species_p { 
    my $self = shift;
    my $org_frame = shift;
    return $self->call_func_that_returns_boolean("lower-taxa-or-species-p \'$org_frame");
}


#### added per Suzanne's request, 5/2008
sub genes_regulating_gene { 
    my $self = shift;
    my $gene = shift;
    return $self->call_func("genes-regulating-gene \'$gene");
}

sub genes_regulated_by_gene { 
    my $self = shift;
    my $gene = shift;
    return $self->call_func("genes-regulated-by-gene \'$gene");
}

sub terminators_affecting_gene { 
    my $self = shift;
    my $gene = shift;
    return $self->call_func("terminators-affecting-gene \'$gene");
}

sub transcription_unit_mrna_binding_sites { 
    my $self = shift;
    my $transcription_unit = shift;
    return $self->call_func("transcription-unit-mrna-binding-sites \'$transcription_unit");
}

sub transcription_unit_activators { 
    my $self = shift;
    my $transcription_unit = shift;
    return $self->call_func("transcription-unit-activators \'$transcription_unit");
}

sub transcription_unit_inhibitors { 
    my $self = shift;
    my $transcription_unit = shift;
    return $self->call_func("transcription-unit-inhibitors \'$transcription_unit");

}

sub containing_tus { 
    my $self = shift;
    my $class = shift;
    return $self->call_func("containing-tus \'$class");
}

sub direct_activators { 
    my $self = shift;
    my $entity = shift;
    return $self->call_func("direct-activators \'$entity");
}

sub direct_inhibitors { 
    my $self = shift;
    my $entity = shift;
    return $self->call_func("direct-inhibitors \'$entity");
}


# Debugging functions
#

sub debug {
    my $self = shift;
    my $message = shift;
    if ($self -> {_debug}) { print "$message\n"; }
} 

sub debug_on {
    my $self = shift;
    $self->{_debug} = "TRUE";
}

sub debug_off {
    my $self = shift;
    $self->{_debug} = "";
}


use strict;
use Getopt::Long;
use Data::Dumper;
   
   my %replace=(gt => '>', lt => '<', quot => '"', amp => '&', alpha => 'alpha', beta => 'beta', gamma => 'gamma', delta => 'delta', epsilon => 'epsilon');
   
   my $DB_NAME;
   my $HELP;
   my $LIST;
   my $OUTFILE;
   my $TYPE;
   
   my $result = GetOptions(
       'help+' => \$HELP,
       'list' => \$LIST,
       'type' => \$TYPE,
       'output=s' => \$OUTFILE,
   );
   
   if (defined $HELP or !(defined $LIST or defined $OUTFILE)) {
       die "
   Usage: $0 [options] -o output.txt PGDB
     Extracts information from a PGDB, producing a tsv (tab delimited) file describing the base pathways
     of the database.
   
     Pathway Tools must be opened with the -api argument:
       /path/to/pathway-tools -api
   
     Options:
      -h, --help            show this help and exit
      -l, --list            list all PGDBs that can be extracted
     \n";
   }
   
   $DB_NAME = shift || "null";

   # run pathway tools in api mode
   #my $thr = threads->create('run_ptools_api_mode', 'argument');
   
   # Open the database
   
   # 1. Is the socket there?
   if (!-e "/tmp/ptools-socket") {die("Pathway-tools is not running in -api mode (/tmp/ptools-socket not found)\n");}
   
   # 2. Connect... hopefully.
   my $cyc = perlcyc->new($DB_NAME);
   
   # 3a. Check basic functionality:
   my @test1;
   while (! @test1) {
     eval {@test1 = $cyc->send_query("(* 222 3)")};
     if (! @test1 and $!=~m/Connection refused/) {
       warn("/tmp/ptools-socket: $!\nTrying again in 5 seconds.\n");
       sleep(5);
       next;
     }
     if (! @test1 and $!=~m/Socket operation on non-socket/) {
       die("/tmp/ptools-socket exists, but is not a socket.\n");}
     if (! @test1) {die("/tmp/ptools-socket: $!\n");}
   }
   
   # 3b. Get results of test query:
   my @test2;
   eval {@test2 = $cyc->retrieve_results();};
   if (! @test2) {die("Fail test2\n");}
   if ($test2[0] ne 666) {die("Pathway-tools is not running correctly (expected 222 * 3 = 666, but got $test2[0]\n");}
   
   # 4a. Get list of organisums in pathway-tools:
   my @test3;
   eval {@test3=$cyc->send_query("(mapcar #'object-name (all-orgs :all))");};
   if (! @test3) {die("Fail test3\n");}
   if ($test3[0] eq ":error") {die("test3 failed\n");}
   
   # 4b. Get reults of list query:
   my @test4;
   eval {@test4 = $cyc->retrieve_results();};
   if (! @test4) {die("Fail test4\n");}
   if ($test4[0] eq ":error") {die("test4 failed\n");}
   
   # 5. Fix up output to something listing just organisms:
   my @dblist=grep {!/^ECOBASE$/ && !/^[#@]/ && /BASE$/} @test4;
   map {s/BASE$//} @dblist;
   push @dblist,"ECOLI";  # Stupid exceptions to every rule
   
   # 6. Show all databases:
   if (defined $LIST) {print "Databases available:\n",map {lc($_)." \n"} sort @dblist;exit 0;}
   if (!grep {/^$DB_NAME$/i} @dblist) {die("There is no $DB_NAME in the list of databases.  Try '$0 -l' to list all available PGDB databases.\n");}
   
   ###############################################
   
   # Main loop:
   open(OUT,">$OUTFILE") or die("$OUTFILE: $!\n");
   my @my_base_pathways = $cyc->call_func("all-pathways :all T");
   
   print OUT "PWY_NAME\tPWY_COMMON_NAME\tNUM_REACTIONS\tNUM_COVERED_REACTIONS\tORF_COUNT\n";
   print  "PWY_NAME\tPWY_COMMON_NAME\tNUM_REACTIONS\tNUM_COVERED_REACTIONS\tORF_COUNT\n";
   my $pwy_count=0;
   foreach my $pathway (@my_base_pathways) {
       $pwy_count++;
       my @mygenes = $cyc->genes_of_pathway($pathway,"T");
       my @totalrxns = $cyc->get_slot_values($pathway, "REACTION-LIST");
     
       my $pathway_common_name = $cyc->get_slot_value($pathway,"common-name") || "?";
       my $num_reactions = scalar(@totalrxns);
       my $num_predicted_orfs = scalar(@mygenes);
     
       my $num_covered_rxns =0;
       my $num_genes =0;
       my %orf_strings = ();
       my $rxnOutputStr="";
       foreach my $reaction (@totalrxns) {
         my @rngenes = $cyc->genes_of_reaction($reaction,"T");
         my @rxngenes = ();
         map { push(@rxngenes, $cyc->get_slot_value($_,"common-name")) } @rngenes;
         if(scalar(@rxngenes) > 0) {
             $num_covered_rxns++;
         }
         my $rxn_name = $cyc->get_slot_value($reaction,"common-name") || "???";
         $num_genes += scalar(@rxngenes);
         map { $orf_strings{$_} =1 } @rxngenes;
         $rxnOutputStr .= "RXN:"."\t".$reaction."\t".$rxn_name."\t".scalar(@rxngenes);
         map{ $rxnOutputStr .= "\t".$_} @rxngenes;
         $rxnOutputStr.="\n";

       }
       my $outputstr =  "PATHWAY:"."\t".$pathway."\t".$pathway_common_name."\t".$num_reactions."\t".$num_predicted_orfs;
       map { $outputstr .= "\t".$_} keys %orf_strings;
       $outputstr.="\n";
       print  $pwy_count."\t".$pathway."\n";
       print OUT $outputstr;
       #print  $rxnOutputStr;
       print OUT $rxnOutputStr;
   }
   close(OUT);
   
   # For replacing entities like &xxxx; with real characters (but spelling out greek letters)
   sub cleanup {
     my $in=$_[0];
     my $search=join("|",keys %replace);
     $in =~ s/<\/?(i|SUB)>//ig;
     $in =~ s/&($search);/$replace{$1}/ieg; # Replace html character entities
     return $in;
   }
   

