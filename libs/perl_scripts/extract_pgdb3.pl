#!/usr/bin/perl -w

# Extracts information about pathways and compounds from a Pathway/Genome
# Database (PGDB). Outputs a CSV file based on the following schema:
#
# (insert a schema here)
#
# Requires the PerlCyc module from the SOL Genomics Network
# (ftp://ftp.sgn.cornell.edu/programs/perlcyc/).

use strict;

use lib './internal-tools/perllib';
use Getopt::Long;
use perlcyc;
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
    die "Usage: $0 [options] -o output.xml PGDB
Extracts information from a PGDB, producing a csv file describing the base pathways
of the database.

Pathway Tools must be opened with the -api argument:
    /path/to/pathway-tools -api

Options:
  -h, --help            show this help and exit
  -l, --list            list all PGDBs that can be extracted
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

foreach my $pathway (@my_base_pathways) {
  my @mygenes = $cyc->genes_of_pathway($pathway,"T");
  my @totalrxns = $cyc->get_slot_values($pathway, "REACTION-LIST");
  #print(Dumper(@mygenes)."\n");
  #print(Dumper(@my_base_pathways)."\n");
  #print(Dumper(@totalrxns)."\n"); # number of reactions in pathway
  
  #print OUT join("\t",$pathway,cleanup($cyc->get_slot_value($pathway,"common-name")),scalar(@mygenes), scalar(@totalrxns), (map {($cyc->get_slot_value($_, "common-name"))||"?"} @mygenes)),"\n";
  my $count = 1;
  foreach my $reaction (@totalrxns) {
     my @myrxngenes = $cyc->genes_of_reaction($reaction,"T");
     print $count."\t".$pathway."\t".$cyc->get_slot_value($pathway,"common-name")."\t".$reaction."\t".scalar(@myrxngenes)."\n";
     $count++;
  }
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
