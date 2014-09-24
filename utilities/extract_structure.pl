#!/usr/bin/perl -w

# Extracts information about pathways and compounds from a Pathway/Genome
# Database (PGDB). Outputs an XML file based on the following schema:
#
# (insert a schema here)
#
# Requires the PerlCyc module from the SOL Genomics Network
# (ftp://ftp.sgn.cornell.edu/programs/perlcyc/).

use strict;

use Getopt::Long;
use perlcyc;
use XML::Writer;
use IO::File;

use Data::Dumper;

my $DB_NAME;
my $HELP;
my $LIST;
my $FULL = undef;
my $OUTFILE;
sub recurse_pathway_classes ;

my $result = GetOptions(
    'help+' => \$HELP,
    'list' => \$LIST,
    'type' => \$FULL,
    'output=s' => \$OUTFILE,
);

if (defined $HELP or !(defined $LIST or defined $OUTFILE)) {
    die "Usage: $0 [options] -o output.xml PGDB
Extracts information from a PGDB, producing an XML file describing every aspect
of the database.

Pathway Tools must be opened with the -api argument:
    /path/to/pathway-tools -api

Options:
  -h, --help            show this help and exit
  -l, --list            list all PGDBs that can be extracted
  -t, --type            flag to produce the full structure
";
}

$DB_NAME = shift || "null";

# Open the database

# 1. Is the socket there?
if (!-e "/tmp/ptools-socket") {die("Pathway-tools is not running in -api mode (/tmp/ptools-socket not found)\n");}

# 2. Connect... hopefully.
my $cyc = perlcyc->new($DB_NAME);

# 3a. Check basic functionality:
my @test1;
eval {@test1 = $cyc->send_query("(* 222 3)")};
if (! @test1) {
  if ($!=~m/Connection refused/) {die("Pathway-tools is not running in -api mode (/tmp/ptools-socket exists, but is not active)\n");}
  if ($!=~m/Socket operation on non-socket/) {die("Pathway-tools is not running in -api mode (/tmp/ptools-socket exists, but is not a socket)\n");}
  die("/tmp/ptools-socket: $!");
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

# 7. Create a writer:
my $output = new IO::File(">$OUTFILE") or die("$OUTFILE: $!");
my $writer = new XML::Writer(DATA_MODE => 1, DATA_INDENT => 4,OUTPUT => $output);
$writer->startTag('pgdb', name => (uc $DB_NAME));

# Pull out reactions, and keep a hash of compounds to drag out later

my %compounds;
my %all_enzymes;

my @reactions = $cyc->call_func('all-rxns :all');
if ($reactions[0] eq ':error') { die("all-rxns said: @reactions\n"); }

my @enzymes = $cyc->call_func('all-enzymes');
if ($enzymes[0] eq ':error') { die("all-enzymes said: @enzymes\n"); }

my @pclass_pathways = $cyc->all_pathways();
if ($pclass_pathways[0] eq ':error') { die("all-pathways said: @pclass_pathways\n"); }

my @base_pathways = $cyc->call_func('all-pathways :all T');
if ($base_pathways[0] eq ':error') { die("all-pathways2 said: @base_pathways\n"); }

my @genes = $cyc->call_func("get-class-all-instances '|Genes|");
if ($genes[0] eq ':error') { die("get-class-all-instances genes said: @genes\n"); }

print STDERR "Number of reactions: ",scalar(@reactions),"\n";
print STDERR "Number of enzymes: ",scalar(@reactions),"\n";
print STDERR "Number of pclass_pathways: ",scalar(@pclass_pathways),"\n";
print STDERR "Number of base_pathways: ",scalar(@base_pathways),"\n";
print STDERR "Number of genes: ",scalar(@genes),"\n";


print STDERR "Fetching pathway classes\n";


# And keep track of what we need to traverse

my %pathways_to_traverse;

for (@pclass_pathways) {
    $pathways_to_traverse{$_}++;
}

# Perform parent traversal until we hit the |Pathways| frame and build a tree of pathway classes

my %pathway_classes;
my %pathway_class_parents;
my %pathway_class_children;


@pclass_pathways = $cyc->all_pathways();

#my @base_pathways = $cyc->call_func('all-pathways :all T');
my %tree;
for my $p (@base_pathways) {
    &recurse_pathway_classes($p, \%tree);
    my @totalrxns = $cyc->get_slot_values($p, "REACTION-LIST");
    if( !defined($tree{$p}) ) {
       %{$tree{$p}} = ();
    }
    for my $reaction (@totalrxns) {
        $tree{$p}{$reaction} = 2;
    }
}

if( defined($FULL) ) {
   &printTree('|Pathways|', '|Pathways|',  \%tree,  $cyc, 0);
}
else {
    for my $name (@base_pathways){
        &printTree($name, $name,  \%tree,  $cyc, 0);
    }
}


# Then output the pathway classes
exit(0);


for my $frame_id (keys %pathway_classes) {
    
    my $common_name = $cyc->get_slot_value($frame_id, 'Common-Name');

    if (not defined $common_name) {
        $common_name = $frame_id;
    }

    # A pathway class is a leaf if it has any children that are actual pathways

    my $leaf = 0;

    if (defined $pathway_class_children{$frame_id} and grep { /^[^\|]+$/ } (keys %{$pathway_class_children{$frame_id}})) {
        $leaf++;
    }


    if (defined $pathway_class_parents{$frame_id}) {
        for (keys %{$pathway_class_parents{$frame_id}}) {
            $writer->emptyTag('parent', 'frame-id' => $_);
        }
    }

    if (defined $pathway_class_children{$frame_id}) {
        for (keys %{$pathway_class_children{$frame_id}}) {
            $writer->emptyTag('child', 'frame-id' => $_);
        }
    }
}
exit(0);

# Pathways
print STDERR 'Fetching pathways';

# All pathways and base pathways

my @pw_pathways = $cyc->all_pathways();


for (@pw_pathways) {
    $pathways_to_traverse{$_}++;
}

my @pathways = keys %pathways_to_traverse;

my $num_pathways = scalar @pathways;
my $pathways_i = 0;


@pathways= qw( |Super-Pathways );
for my $pathway (@pathways) {

    $pathways_i++;

    my $pathway_name = $cyc->get_slot_value($pathway, 'Common-Name');
    $pathway_name =~ s/<.+?>//g;

    # Generate properties

    my %properties = ('frame-id' => $pathway, 'common-name' => $pathway_name);

    if (grep { $_ eq $pathway } @base_pathways) {
        $properties{base} = 'yes';
    }

    # Pathway classes

    my @ocelot_parents = $cyc->get_slot_values($pathway, 'Ocelot-GFP::Parents');
    $writer->startTag('pathway-classes');
    for my $ocp (@ocelot_parents) {
        if ($ocp =~ /^\|.+\|$/) {
            $writer->emptyTag('pathway-class', 'frame-id' => $ocp);
        }
    }

    # Superpathways
    my @superpathways = $cyc->get_slot_values($pathway, 'Super-Pathways');
    for my $sp (@superpathways) {
        $writer->emptyTag('superpathway', 'frame-id' => $sp);
    }

    # Subpathways

    my @subpathways = $cyc->get_slot_values($pathway, 'Sub-Pathways');
    my $pwy_common_name = $cyc->get_slot_value($pathway, 'Common-Name');

    for my $sp (@subpathways) {
        my $sp_common_name = $cyc->get_slot_value($sp, 'Common-Name');
        #print "P  ".$pathway."\t". $pathway_name."\n";
        #print $sp."\t".$sp_common_name."\n";
        #$writer->emptyTag('subpathway', 'frame-id' => $sp);
    }


}


sub recurse_pathway_classes {

    my $frame_id = shift;
    my $tree = shift;

    if (not $frame_id eq '|Pathways|') {


        my @parents = $cyc->get_slot_values($frame_id, 'ocelot-gfp::parents');

        # Assign child relationships for each parent whose frame type is :CLASS

        for my $p (@parents) {
            if( !defined($tree->{$p}) ) {
               %{$tree->{$p}} = ( ); 
            }
             $tree->{$p}{$frame_id} = 1;
        }

        for my $p (@parents) {
            &recurse_pathway_classes($p, $tree);
        }
    }

}

sub printTree{
   my $node = shift;
   my $pnode = shift;
   my $tree_ = shift;
   my $cyc = shift;
   my $indent = shift;
   
   my $common_name = $cyc->get_slot_value($node, 'Common-Name') || $node;

   my $outStr = "\t" x $indent . $node."\t".$common_name;
   if( !(ref $tree_->{$pnode}{$node} eq "HASH") && defined $tree_->{$pnode}{$node} &&   $tree_->{$pnode}{$node} == 2 ) {
         my @rxngenes = $cyc->genes_of_reaction($node,"T");
         $outStr .= " (" . scalar(@rxngenes) . ")";
   }


   $outStr =~ s/<[a-zA-Z]>//g;
   $outStr =~ s/<\/[a-zA-Z]>//g;

   $outStr =~ s/<em>//g;
   $outStr =~ s/<\/em>//g;
   $outStr =~ s/<small>//g;
   $outStr =~ s/<\/small>//g;
   $outStr =~ s/<sup>//ig;
   $outStr =~ s/<\/sup>//ig;
   $outStr =~ s/<sub>//ig;
   $outStr =~ s/<\/sub>//ig;

   print $outStr."\n";
   if( !defined($tree_->{$node} ) ) {
      return;
   }

   for my $cnode ( keys %{$tree_->{$node}} ) {
      &printTree($cnode, $node, $tree_, $cyc, $indent+1);
   }

}


1;
