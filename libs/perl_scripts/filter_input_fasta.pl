#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;

use Data::Dumper;
use CGI ':standard';
use LWP 5.64; # Loads all important LWP classes, and makes
                #  sure your version is reasonably recent.
no strict; # TODO fix this


use constant VERSION => 1.0;
use constant DEBUG => 1;


# returns the directory of the file
sub get_dir {
   my $file  = shift;
   my $dir; 
   my @fields = split("\/", $file);

   if( scalar( @fields ) > 1 ) {
       pop(@fields);
       $dir = join("\/",@fields);
       if( $dir eq "." ) {
          $dir="";
       }
   }
   else {
      $dir="";
   }
   if( length($dir) > 0 ) {
      $dir= $dir."\/";
   } 
   return $dir;
}


# get the current directory 
my $dir = get_dir(__FILE__);
require  "${dir}Utilities.pl";


{
#input fasta  file 
my $INPUT_FILE;

#output fasta  file 
my $OUTPUT_FILE;

#mapping file 
my $MAPPING_FILE;


# minimum length of the sequences after trimming  that are with
# at most with $MAX_AMGIG_X  X (ambiguous) characters
my $MIN_LENGTH;


#variable for the log file name
my $LOG_FILE;

my $result;
# If there are command line arguments:
if (defined $ARGV[0]) {
  # Process command line arguments
   $result=GetOptions(
    "i=s"=>\$INPUT_FILE,
    "o=s"=>\$OUTPUT_FILE,
    "m=s"=>\$MAPPING_FILE,
    "min_length=s"=>\$MIN_LENGTH,
    "log_file=s"=>\$LOG_FILE
  );
}



# Usage, or bad arguments:
if (!$result or !defined($INPUT_FILE) or !defined($MIN_LENGTH) 
    or !defined($OUTPUT_FILE)  or !defined($LOG_FILE) ) {
  die("Usage: $0 -i file.fna   --min_length N --log_file logfile.log -o outfasta  [ -m sequencefile ]
       This program will select sequences from with a  minimum length
       of  min_length after trimming  ambiguous characters
       The logs are recorded in the logfile. 
       If the  sequencename_map_file is  provided then the sequence names are translated
       and a mapping file sequencename_map_file is provided with the mapping\n"
     );  

}

my $sequences = read_fasta_sequences($INPUT_FILE);

my $sample_name = $INPUT_FILE;
$sample_name =~s/.*\///ig;
$sample_name =~s/\.fasta//ig;
$sample_name =~s/\.fna//ig;
$sample_name =~s/\.fas//ig;

#open the output file

open my $outfile, '>' , $OUTPUT_FILE or die( "cannot open output  file  $OUTPUT_FILE") ;
my $mappingfile;
if( defined($MAPPING_FILE) )  {
   open  $mappingfile, '>' , $MAPPING_FILE or die( "cannot open output  file  $OUTPUT_FILE") ;
}


my $seq;


use constant NUMSEQ => "Number of sequences :";   
use constant NUMSEQ_SHORTER => "Number of sequences shorter than ";
use constant AV_LENGTH=> "Average length of sequences:";
use constant MIN_LENGTH=> "Minimum length of sequences:";
use constant MAX_LENGTH=> "Maximum length of sequences:" ;

my %stats = ( 
              MIN_LENGTH() , { BEFORE => 10000000, AFTER=>1000000 }, 
              MAX_LENGTH() , { BEFORE => 0, AFTER=>0 }, 
              NUMSEQ() , { BEFORE => 0, AFTER=>0},   
              NUMSEQ_SHORTER() , { BEFORE => 0, AFTER=>0 },
              AV_LENGTH() , { BEFORE => 0, AFTER=>0 },
            );

my %length_distribution = ();

for(my $i=0; $i<=30; $i++) {
    $length_distribution{$i}= 0;
    $length_cumulative_distribution{$i}= 0;
}


my $seq_count=0;
while( my ($key, $seqvalue) =  each %{$sequences} ) {
    $seq = $seqvalue;

    my $index = int(length($seq) / 50);
    #print length($seq) ."\t".$index."\n";
    $length_distribution{$index}++;

    if( length($seq ) < $stats{MIN_LENGTH()}{BEFORE} ) {
          $stats{MIN_LENGTH()}{BEFORE} = length($seq); 
    }
    if( length($seq ) > $stats{MAX_LENGTH()}{BEFORE} ) {
          $stats{MAX_LENGTH()}{BEFORE} = length($seq); 
    }

    if( length($seq ) < $MIN_LENGTH )  {
       $stats{NUMSEQ_SHORTER()}{BEFORE}++; 
    }
    $stats{AV_LENGTH()}{BEFORE}  =  $stats{AV_LENGTH()}{BEFORE} + length($seq); 

    $seqvalue = filter_sequence($seq);
    $stats{NUMSEQ()}{BEFORE}++; 
    
    if( length($seqvalue ) >= $MIN_LENGTH )  {
       $stats{NUMSEQ()}{AFTER}++; 
       $stats{AV_LENGTH()}{AFTER}  =  $stats{AV_LENGTH()}{AFTER} + length($seqvalue); 
       if( !defined($MAPPING_FILE) )  {
          print $outfile $key."\n";
       }
       else {
            print $outfile ">${sample_name}_$seq_count\n";
            $key =~s/^>//g;
            print $mappingfile "${sample_name}_$seq_count\t$key\n";
            $seq_count++;
       }

       print $outfile $seqvalue."\n";
       
       if( length($seqvalue ) < $stats{MIN_LENGTH()}{AFTER} ) {
           $stats{MIN_LENGTH()}{AFTER} = length($seqvalue); 
       }
       if( length($seqvalue ) > $stats{MAX_LENGTH()}{AFTER} ) {
           $stats{MAX_LENGTH()}{AFTER} = length($seqvalue); 
       }
    }
}

if( $stats{NUMSEQ()}{BEFORE} > 0 ) { 
   $stats{AV_LENGTH()}{BEFORE}  = $stats{AV_LENGTH()}{BEFORE}/$stats{NUMSEQ()}{BEFORE}; 
}
else {
   $stats{AV_LENGTH()}{BEFORE}  = 0;
}

if( $stats{NUMSEQ()}{AFTER} > 0 ) { 
   $stats{AV_LENGTH()}{AFTER}  = $stats{AV_LENGTH()}{AFTER}/$stats{NUMSEQ()}{AFTER}; 
}
else {
   $stats{AV_LENGTH()}{AFTER}  = 0;
}


close($outfile);
if( defined( $MAPPING_FILE) )  {
  close($mappingfile);
}

open my $logfile, '>' , $LOG_FILE or die( "cannot open log  file  $LOG_FILE") ;

printf($logfile  "   %-40s\t%s\t%s\n", " ", "BEFORE","AFTER");
printf($logfile  "   %-40s\t%d\t%d\n", NUMSEQ(), $stats{NUMSEQ()}{BEFORE}, $stats{NUMSEQ()}{AFTER}); 
printf($logfile  "   %-40s\t%d\t%d\n", NUMSEQ_SHORTER().$MIN_LENGTH.":", $stats{NUMSEQ_SHORTER()}{BEFORE}, $stats{NUMSEQ_SHORTER()}{AFTER}); 
printf($logfile  "   %-40s\t%d\t%d\n", AV_LENGTH(), $stats{AV_LENGTH()}{BEFORE}, $stats{AV_LENGTH()}{AFTER}); 
printf($logfile  "   %-40s\t%d\t%d\n", MIN_LENGTH(), $stats{MIN_LENGTH()}{BEFORE}, $stats{MIN_LENGTH()}{AFTER}); 
printf($logfile  "   %-40s\t%d\t%d\n", MAX_LENGTH(), $stats{MAX_LENGTH()}{BEFORE}, $stats{MAX_LENGTH()}{AFTER}); 



printf($logfile  "\n\n"); 
printf($logfile  "   READ_LENGTH_RANGE\tFREQUENCY\t\tMIN_LENGTH\tCUMULATIVE_FREQUENCY\n"); 
printf($logfile  "   -----------------\t---------\t\t----------\t--------------------\n"); 

$length_cumulative_distribution{30} =  $length_distribution{30}; 
for(my $i=29; $i>=0; $i--){
     $length_cumulative_distribution{$i} = $length_cumulative_distribution{$i+1} + $length_distribution{$i}; 
}

for(my $i=0; $i<=30; $i++) {
     printf($logfile  "   %d-%d\t%9d\t\t\t%d\t\t%d\n", $i*50, ($i+1)*50, $length_distribution{$i}, ($i+1)*50,$length_cumulative_distribution{$i}); 
}

close($logfile);


}  #end of main



