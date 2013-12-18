
#########   Libraries of subroutines ##################

#removes the ambiguous characters from the sequences and 
# returns the longest clean sequence

sub filter_sequence($) {
   my $sequence = shift;
   $sequence =~ s/^[ATCG]/-/g;
   
   my @subsequences = split( "-", $sequence);
   #print Dumper(\@subsequences);
   
   my $max_length = 0;
   my $longest_sequence = "";
   foreach my $seq ( @subsequences ) {
      if( length($seq) > $max_length ) {
           $longest_sequence = $seq;
           $max_length = length($seq);
      }
   }
   return( $longest_sequence);

}



# Perl trim function to remove whitespace from the start and end of the string
 sub trim($)
 {
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}



sub standard_dev($) {
   my $array=shift;
   my @arraysq = qw();

   if( @${array} == 0 ) { die("Cannot compute standard deviation for an array of size 0\n"); }
   foreach my $val (@${array} ) {
      push( @arraysq, $val*$val);
   }
   return( sqrt( average(\@arraysq) - average($array)*average($array)) );
}

sub average($) {
   my $array=shift;
   
   if( @${array} == 0 ) { die("Cannot compute average for an array of size 0"); }
   return( sum(${array})/@${array} ); 
}


sub sum($) {
   my $array=shift;
   my $sum = 0;
   foreach my $val ( @${array}) {
      $sum = $sum + $val;
   }
   return($sum);
}

sub max($) {
   my $array=shift;
   my $max = $array->[0];
   foreach my $val ( @${array}) {
      if($val > $max) {
        $max =  $val;
      }
   }
   return($max);
}

sub min($) {
   my $array=shift;
   my $min = $array->[0];
   foreach my $val ( @${array}) {
      if($val < $min) {
        $min =  $val;
      }
   }
   return($min);
}
   

#
sub read_fasta_sequences
{
        my $filename =shift;
        my $line;
        my $first = 0;
        unless(open(DATAFILE, $filename)) 
        { 
	      print  "Cannot open file $filename \n";
          exit;
        }
        my @filedata = <DATAFILE>;
        
        my $seqname;
        my %sequences = ();
        foreach my $line (@filedata) {
              #  print $line;
                if ($line =~ /^>/)
                {
                        $seqname = trim($line);
                        $sequences{$seqname}="";
                        if ($first == 0)
                        {
                                $first = 1;
                        }
                        next;
                }
                else {
                    $sequences{$seqname} = $sequences{$seqname}.trim($line);
                }
                if ($first == 0)
                {
                        die"Not a standard FASTA file. Stop.\n";
                }
        }
        close(DATAFILE);
        return(\%sequences);
}
#
 sub get_file_data {
   my ($filename) = @_;
   use strict;
   use warnings;
   unless(open(DATAFILE, $filename)) 
   { 
	 print  "Cannot open file $filename \n";
     exit;	
   }
   my @filedata = <DATAFILE>;
   close(DATAFILE);

   return @filedata;
 }
 
 
 
#open_file  subroutine
#
#
sub open_file {
	my($filename) = @_;
	my $fh;

        unless(open(DATAFILE, $filename)) 
	{ 
		print  "Cannot open file $filename \n";
        	exit;	
        }
	return $fh;
}

sub binary_string {
      my ($array, $target) = @_;
  
      # $low is first element that is not too low;
      # $high is the first that is too high
      #
      my ( $low, $high ) = ( 0, scalar(@$array) );
   
      # Keep trying as long as there are elements that might work.
      #
      
      while ( $low < $high ) {
         # Try the middle element.
  
         use integer;
         my $cur = ($low+$high)/2;
         #   print $array->[$cur]."    ".$target."\n"; 
         if (quotemeta($array->[$cur]) lt quotemeta($target)) {
            $low  = $cur + 1;                     # too small, try higher
         } 
         else {
            $high = $cur;                         # not too small, try lower
         }
      }
      return $low;
}

sub array_unique($) {
   my $array=shift;
   my %unique_hash = ();
   foreach my $val ( @${array}) {
      $unique_hash{$val} = 1;
   }
   return( keys %unique_hash);
}

1;
