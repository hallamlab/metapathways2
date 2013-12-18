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

sub close {
    my $self = shift;
    close($self->{_socket});
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

1;

__END__

=head1 NAME 

perlcyc - A Perl interface for Pathway Tools software.  Pathway Tools
software needs to run a special socket server program for this module
to work.

=head1 SYNOPSIS

perlcyc is a Perl interface for Pathway Tools software.

C<use perlcyc;>

C<my $cyc = perlcyc -E<gt> new("ARA");>
C<my @pathways = $cyc -E<gt> all_pathways();>

=head1 VERSION

Version 1.21 (May 2008).


=head1 HISTORY

Version History

=over 5

=item 0.1  March, 2002
     
[Lukas Mueller] initial version

=item 0.3  April 22, 2002
     
[Danny Yoo] Added better list parsing

=item 0.9  [Lukas Mueller] 

Added more functions

=item 1.0  August 28, 2002 [Lukas Mueller] 

Added pod documentation and eliminated some bugs.

=item 1.1  June 6, 2003
     
[Thomas Yan] Fixed some minor bugs.

=item 1.2  December 7, 2006 [Lukas Mueller] 

Added three functions: create-frame, create-class, and create-instance.

=item 1.21 May 7, 2008 [Lukas Mueller]

Added three functions that are new in PT v. 12: 
 genes_in_same_operon 
 pwys-of-organism-in-meta
 enzymes-of-organism-in-meta
 lower-taxa-or-species-p org-frame

 other new functions:
 get-class-all-subs

 genes-regulating-gene
 genes-regulated-by-gene
 terminators-affecting-gene
 transcription-unit-mrna-binding-sites
 transcription-unit-activators
 transcription-unit-inhibitors
 containing-tus
 direct-activators
 direct-inhibitors


=back

=head1 INSTALLATION

Installation is standard as for any Perl module. If you downloaded the
compressed tar file, uncompress and untar the file with the following
commands:

 gzip < perlcyc.tar.gz | tar xvf -

This will create a directory called perlcyc in your current
directory. To install the program, type

 make
 make install

The program should now be available in all your Perl programs. 'make
install' may require root access or access through sudo. For the
latter case, type

 sudo make install

You will be prompted for your password.

=head1 PATHWAY TOOL REQUIREMENTS

To use the Perl module, you also need the socket_server.lisp
program. In Pathway Tools version 8.0 or later, the server program can be started with the command line option "-api". On earlier versions, the server daemon needs to be loaded manually, as follows:
start Pathway-Tools with the -lisp option, at the prompt, type: (load "/path/to/socket_server.lisp"),
then start the socket_server by typing (start_external_access_daemon :verbose? t). The server is now ready to accept connections and queries.

=head1 DESCRIPTION

perlcyc.pm is a Perl module for accessing internal Pathway-Tools
functions. For the description of what the individual functions do,
please refer to the Pathway Tools documentation at
http://bioinformatics.ai.sri.com/ptools .  

In general, the Lisp
function name has to be converted to something compatible with Perl:
Dashes have to be replaced by underlines, and question marks with
underline p (_p).  

Note that the optional parameters of all functions
    are not supported in perlcyc, except for all_pathways() which can use the optional arguments :all T to get the base pathways only (no super-pathways).

=head3 Limitations

Perlcyc does not implement the GFP objects in Perl, rather it just sends snippets of code to Pathway-Tools through a socket connection. Only one such connection can be openend at any given time. Because the objects are not implemented in Perl, only object references are supported.


=head3 Available functions:

Object functions:

 new
 Parameters: The knowledge base name. Required!

GFP functions: 
More information on these functions can be found at:
http://www.ai.sri.com/~gfp/spec/paper/node63.html
  
 variants_of_pathway
 get_slot_values
 get_slot_value 
 get_class_all_instances 
 instance_all_instance_of_p
 member_slot_value_p 
 fequal 
 current_kb  
 put_slot_values
 put_slot_value
 add_slot_value
 replace_slot_value
 remove_slot_value
 coercible_to_frame_p
 class_all_type_of_p
 get_instance_direct_types
 get_instance_all_types
 get_frame_slots
 put_instance_types
 save_kb
 revert_kb
 
Pathway-tools functions:
More information on these functions can be found at:
http://bioinformatics.ai.sri.com/ptools/ptools-fns.html

 select_organism 
 all_pathways   
 all_orgs
 all_rxns
 genes_of_reaction
 substrates_of_reaction
 enzymes_of_reaction
 reaction_reactants_and_products
 get_predecessors
 get_successors
 genes_of_pathway
 enzymes_of_pathway
 compounds_of_pathway
 substrates_of_pathway
 transcription_factor_p
 all_cofactors
 all_modulators
 monomers_of_protein
 components_of_protein
 genes_of_protein 
 reactions_of_enzyme
 enzyme_p 
 transport_p
 containers_of 
 modified_forms 
 modified_containers
 top_containers 
 reactions_of_protein
 regulon_of_protein 
 transcription_units_of_protein
 regulator_proteins_of_transcription_unit
 enzymes_of_gene 
 all_products_of_gene
 reactions_of_gene 
 pathways_of_gene 
 chromosome_of_gene
 transcription_unit_of_gene 
 transcription_unit_promoter 
 transcription_unit_genes 
 transcription_unit_binding_sites 
 transcription_unit_transcription_factors
 transcription_unit_terminators 
 all_transported_chemicals 
 reactions_of_compound 
 full_enzyme_name 
 enzyme_activity_name 
 find_indexed_frame
 create-instance
 create-class
 create-frame

 pwys-of-organism-in-meta
 enzymes-of-organism-in-meta
 lower-taxa-or-species-p org-frame
 get-class-all-subs

added 5/2008 per Suzanne's request:
 genes-regulating-gene
 genes-regulated-by-gene
 terminators-affecting-gene
 transcription-unit-mrna-binding-sites
 transcription-unit-activators
 transcription-unit-inhibitors
 containing-tus
 direct-activators
 direct-inhibitors



not supported:
 get_frames_matching_value (why not?)

Internal functions:
 
 parselisp
 send_query
 retrieve_results
 wrap_query
 call_func
 debug 
 debug_on
 debug_off 

Deprecated functions
 parse_lisp_list

=head1 EXAMPLES

Change product type for all genes that are in a pathway to 'Enzyme'

 use perlcyc;

 my $cyc = perlcyc -> new ("ARA");
 my @pathways = $cyc -> all_pathways();

 foreach my $p (@pathways) {
   my @genes = $cyc -> genes_of_pathway($p);
   foreach my $g (@genes) {
     $cyc -> put_slot_value ($g, "Product-Types", "Enzyme");
   }
 }


Load a file containing two columns with accession and a comment into the comment field of the corresponding accession:

 use perlcyc;
 use strict;
 
 my $file = shift;

 my $added=0;
 my $recs =0;

 open (F, "<$file") || die "Can't open file\n";

 print STDERR "Connecting to AraCyc...\n";
 my $cyc = perlcyc -> new("ARA");

 print STDERR "Getting Gene Information...\n";
 my @genes = $cyc -> get_class_all_instances("|Genes|");

 my %genes;

 print STDERR "Getting common names...\n";
 foreach my $g (@genes) {
   my $cname = $cyc -> get_slot_value($g, "common-name"); 
   $genes{$cname}=$g; 
 } 

 print STDERR "Processing file...\n";
 while (<F>) { 
   my ($locus, $location, @rest) = split /\t/;
   $recs++;
   if (exists $genes{$locus}) { 
       my $product = $cyc -> get_slot_value($genes{$locus}, "product");
         if ($product) {
         $cyc -> add_slot_value($product, "comment", "\"\nTargetP location: $location\n\"");
          #print STDERR "Added to description of frame $product\n";
         $added++;
       }
     }
 }

 close (F);

 print STDERR "Done. Added $added descriptions. Total lines in file: $recs. \n";


Add a locus link to the TAIR locus page for each gene in the database

 use strict;
 use perlcyc;

 my $added =0;
 my $genesprocessed=0;

 print "Connecting to AraCyc...\n";
 my $cyc = perlcyc -> new ("ARA");

 print "Getting Gene Information...\n";
 my @genes = $cyc -> get_class_all_instances ("|Genes|");

 print "Adding TAIR links...\n";
 foreach my $g (@genes) {
   $genesprocessed++;
   my $common_name = $cyc -> get_slot_value($g, "common-name");
   if ($common_name && ($common_name ne "NIL")) {
     $cyc -> put_slot_value ($g, "dblinks", "(TAIR \"$common_name\")"); 
     $added++;
   }
   if ((!$genesprocessed ==0) && ($genesprocessed % 1000 == 0)) { print "$genesprocessed ";}
 }

 print "Done. Processed $genesprocessed genes and added $added links. Thanks!\n";
 $cyc -> close();

=head1 TROUBLESHOOTING

If your program terminates with the following error message:
C<connect: No such file or directory at perlcyc.pm line 166.>
then the lisp_server.lisp module in Pathway Tools is not running.
Refer to http://aracyc.stanford.edu for more information on how to 
run the server program.

Please send bug reports and comments to 
LAM87@cornell.edu

=head1 LICENSE

According to the MIT License, http://www.opensource.org/licenses/mit-license.php

Copyright (c) 2002 by Lukas Mueller, The Arabidopsis Information Resource (TAIR) and the Carnegie Institution of Washington.

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

=head1 AUTHOR & CURRENT MAINTAINER

Lukas Mueller (lam87@cornell.edu)

=head1 ACKNOWLEDGMENTS

Many thanks to Suzanne Paley, Danny Yoo and Thomas Yan.

=cut

1;
