#! /usr/bin/perl -w

#####################################################
##
## mltreemap.pl
##
#####################################################

use strict;
use warnings;
use File::Basename;
#use Data::Dumper;

print "Start MLTreeMap v2.061\n";
my $mltreemap                           = mltreemap->new();
my $small_subroutines                   = small_subroutines->new();
my $tree_rerooter                       = reroot_RAxML_trees->new();

my $user_options                        = $mltreemap->read_user_input(\@ARGV);
my $blast_input_files                   = $mltreemap->split_input_file($user_options);
my ($cog_list,$text_of_analysis_type)   = $mltreemap->read_cog_file($user_options);
my $non_wag_cog_list                    = $mltreemap->get_non_wag_cogs;
my $tree_numbers_translation            = $mltreemap->read_species_translation_files($user_options,$cog_list);
my $blast_results_raw_files             = $mltreemap->blast($small_subroutines,$user_options,$blast_input_files);
my $blast_hits_purified                 = $mltreemap->read_blast_file($small_subroutines,$user_options,$blast_results_raw_files,$cog_list);
$small_subroutines->undef_hashes($blast_results_raw_files,$blast_input_files);
my ($contig_coordinates,$shortened_sequence_files)               = $mltreemap->produce_shortened_sequence_files_for_genewise($user_options, $blast_hits_purified);
my $genewise_outputfiles                = $mltreemap->start_genewise($user_options,$shortened_sequence_files,$blast_hits_purified);
my $genewise_summary_files              = $mltreemap->parse_genewise_results($small_subroutines,$user_options,$contig_coordinates,$genewise_outputfiles);
$mltreemap->get_rRNA_hit_sequences($user_options, $blast_hits_purified,$cog_list,$genewise_summary_files);
$small_subroutines->undef_hashes($blast_hits_purified,$contig_coordinates,$shortened_sequence_files,$genewise_outputfiles);
my $hmmalign_singlehit_files            = $mltreemap->prepare_and_run_hmmalign($user_options,$genewise_summary_files,$cog_list);
my ($concatenated_mfa_files,$nrs_of_sequences,$models_to_be_used)  = $mltreemap->concatenate_hmmalign_singlehits_files($user_options,$hmmalign_singlehit_files,$non_wag_cog_list);
my $gblocks_files                       = $mltreemap->start_gblocks($user_options,$concatenated_mfa_files,$nrs_of_sequences);
my $phy_files                           = $mltreemap->produce_phy_file($user_options,$gblocks_files,$nrs_of_sequences);
$small_subroutines->undef_hashes($genewise_summary_files,$hmmalign_singlehit_files,$concatenated_mfa_files,$nrs_of_sequences,$gblocks_files);
my $raxml_outfiles                      = $mltreemap->start_RAxML($small_subroutines,$user_options,$phy_files,$cog_list,$models_to_be_used);
my $final_RAxML_output_files            = $mltreemap->parse_RAxML_output($user_options,$tree_rerooter,$tree_numbers_translation,$raxml_outfiles,$text_of_analysis_type);
$mltreemap->concatenate_RAxML_output_files($user_options,$final_RAxML_output_files,$text_of_analysis_type);

my $final_RAxML_files  = BFile->new( 'file'=>$user_options->{-o}.'/final_RAxML_outputs.xml', 'permission' => 'w');
$final_RAxML_files->addFolder( $user_options->{-o}.'/final_RAxML_outputs/');
$final_RAxML_files->Close();

sub remove_dir;
sub create_consolidated_xml_file;

exit(0);
if( !exists $user_options->{-x} ) {
    print "Done\n";
    exit(0);
}

create_mltreemap_sequence_hits($user_options->{-o}.'/various_outputs/');
create_consolidated_xml_file($user_options->{-o}.'/various_outputs/', $user_options->{-o}.'/various_outputs.xml');
#remove_dir($user_options->{-o}.'/various_outputs/');

create_consolidated_xml_file($user_options->{-o}.'/final_RAxML_outputs/', $user_options->{-o}.'/final_RAxML_outputs.xml');
#remove_dir($user_options->{-o}.'/final_RAxML_outputs/');

print "Done.\n";

#attention: the identifiers of the COGs must be exactly 7 characters long.

sub create_consolidated_xml_file {
    my $foldername = shift;
    my $xmlname = shift;
    my $bfile  = BFile->new( 'file'=>$xmlname, 'permission' => 'w');
    $bfile->addFolder($foldername);
    $bfile->Close();
}

sub create_mltreemap_sequence_hits {
    my $foldername = shift;
    opendir(VDIR, $foldername) or croak('cannot open the folder'. $foldername);
    my @files = readdir(VDIR);
    foreach my $file (@files) {
       if( $file =~/'.*COG[0-9]*.*\.fa$/gi) {
          my (@matches) = ($file=~/[a-zA-Z]*_(.*)__[0-9]*__*(COG[0-9]*).*.fa$/);
       }
    }
    close(VDIR);
}


sub remove_dir { 
    my $foldername = shift;
    opendir(VDIR, $foldername) or croak('cannot open the folder'. $foldername);
    my @files_to_delete = readdir(VDIR);
    foreach my $file (@files_to_delete) {
       unlink($foldername.$file);
    }
    rmdir($foldername);
    close(VDIR);
}


#####################################################
#####################################################

# package mltreemap

#####################################################
#####################################################

package mltreemap; {
#use Data::Dumper;
#####################################################
# new
#####################################################

sub new {
    my $package = shift;
    my $reference = {};
    bless ($reference, $package);
    return ($reference);    
}

#####################################################
# read_user_input
#####################################################

sub read_user_input {    
    my $package = shift;
    my $argv = shift;
    
    my $error_message = "\nERROR: incorrect input.\nexample input:\n./mltreemap.pl -i example_input/test.txt -e binaries_folder\n\n";
    $error_message .= "mandatory input option:\n-i: your sequence input file. -e folder  containing blastp, blastx, hmmalign, genewise, gblocks and raxaml\n\n";
    $error_message .= "optional input options:\n-b: number of Bootstrap replicates (Default 0 i.e. no Bootstrapping).\n";
    $error_message .= "-c: usage a computer cluster (0 = no cluster (default), s = sun grid).\n";
    $error_message .= "-f: RAxML algorithm (v = Maximum likelihood (default), p = Maximum Parsimony).\n";
    $error_message .= "-g: minimal sequence length after Gblocks (default = 50).\n";
    $error_message .= "-l long input files will be split into files of n sequences each (default 2000)\n";
    $error_message .= "-m minimum memory on a sungrid_cluster in GB (default 0)\n";
    $error_message .= "-o: output directory (default: output/)\n";
    $error_message .= "-s: minimum bitscore for the blast hits (default 60)\n";
    $error_message .= "-t: phylogenetic reference tree (p = MLTreeMap reference tree (default), g = GEBA reference tree, i = fungi tree).\n";
    $error_message .= "-r: reference selection (a = all reference data, s = selection file).\n";
    $error_message .= "\n";
        
    my %mandatory_options = (-i => 0,  -e =>"");
    my %option_defaults = (-b => 0, -f => "v", -o => "output/", -m => 0, -t => "p", -s => 60, -g => 50, -c => 0, -l => 2000, -r => "a");
    
    my %user_options = %option_defaults;
    my $current_option = "";
    my $arg_count = 0;
    
    foreach my $argument (@$argv) {      
        if ( $arg_count % 2 ) {
            #is input
            $user_options{$current_option} = $argument;
        } else {
            #is option denominator
            die "$error_message" unless (exists ($option_defaults{$argument}) || exists ($mandatory_options{$argument}));
        }
        $current_option = $argument;        
        $arg_count++;   
    }
    
    #do some checks on the user input:
    
    #first, check if the mandatory options have been used:
    
    die "$error_message" unless $arg_count;
    foreach my $mandatory_option (sort {$a cmp $b} keys %mandatory_options) {
        die "$error_message" unless exists ($user_options{$mandatory_option});    
    }
    
    #second, do some (hard coded) checks on certain input variables:
    
    die "$error_message" unless ($user_options{"-b"} =~ /\A\d+\Z/);
    die "$error_message" unless (($user_options{"-f"} eq "v") || ($user_options{"-f"} eq "p"));
    die "$error_message" unless (($user_options{"-t"} eq "p") || ($user_options{"-t"} eq "g") || ($user_options{"-t"} eq "i"));
    die "$error_message" unless (($user_options{"-c"} eq "0") || ($user_options{"-c"} eq "s"));
    die "$error_message" unless ($user_options{"-s"} =~ /\A\d+\Z/);
    die "$error_message" unless ($user_options{"-m"} =~ /\A\d+\Z/);
    die "$error_message" unless ($user_options{"-g"} =~ /\A\d+\Z/);
    die "$error_message" unless (($user_options{"-l"} =~ /\A\d+\Z/) && ($user_options{"-l"} > 0));
    
    ############ Edited by Young C. Song                   ###############
    ############ new option: -r                            ###############
    ############ default value: a for 'all reference data' ###############
    die "$error_message" unless (($user_options{"-r"} eq "a") || ($user_options{"-r"} eq "s")) ;
    #third, set the prefix with the alignment data:
    
    $user_options{"reference_data_prefix"} = "";
    $user_options{"reference_data_prefix"} = "geba_" if ($user_options{"-t"} eq "g");
    $user_options{"reference_data_prefix"} = "fungi_" if ($user_options{"-t"} eq "i");
    
    #fourth, set the phylogenetic reference tree name.
    
    $user_options{"reference_tree"} = "MLTreeMap_reference.tree";
    $user_options{"reference_tree"} = "geba.tree" if ($user_options{"-t"} eq "g");
    $user_options{"reference_tree"} = "fungitr_tree.txt" if ($user_options{"-t"} eq "i");

    #create the output directories
        
    if (-e $user_options{-o}) {
#        print "WARNING: Your output directory \"$user_options{-o}\" allready exists!\n";
#        print "Overwrite [1], quit [2] or change directory [3]?\n";
        
#        my $alarm_timer = 30;
#        $SIG{ALRM} = sub { die "No input for more than $alarm_timer sek. Exit MLTreeMap.\n"};
#        alarm $alarm_timer;
#        my $answer = <STDIN>;
#        my $answer = "1";
#        $answer = "2" unless $answer; #quit if something goes wrong with <STDIN>.
#        chomp $answer;
#        alarm 0;
#        while (($answer ne "1") && ($answer ne "2") && ($answer ne "3")) {
#            print "invalid input. Please chose 1, 2 or 3.\n";
#            alarm $alarm_timer;
#            $answer = <STDIN>;
#            chomp $answer;
#            alarm 0;
#        }
#        if ($answer == 1) {
#            print "Do you really want to overwrite the old output directory?\n";
#            print "All data in it will be lost!\n";
#            print "Yes [y] or no [n]?\n";         
#            alarm $alarm_timer;
#            my $answer2 = <STDIN>;
#            chomp $answer2;
#            alarm 0;
#            while (($answer2 ne "y") && ($answer2 ne "n")) {
#                print "invalid input. Please chose y or n.\n";
#                alarm $alarm_timer;
#                $answer2 = <STDIN>;
#                chomp $answer2;
#                alarm 0;
#            }
#            if ($answer2 eq "y") {
#                my $rm_command = "rm -r $user_options{-o}";
#                system ($rm_command);
#            } else {
#                die "Exit MLTreeMap\n";    
#            }
#        } elsif ($answer == 2) {
#            die "Exit MLTreeMap\n";;
#        } elsif ($answer == 3) {
#            print "Please enter the new directory.\n";
#            alarm 600;
#            my $dir = <STDIN>;
#            chomp $dir;
#            alarm 0;
#            $user_options{-o} = $dir;
#        }
    }
    
    my @output_directories = "";
    $user_options{-o} .= "/" unless ($user_options{-o} =~ /\/\Z/);
    $output_directories[0] = $user_options{-o};
    $output_directories[1] = "$user_options{-o}various_outputs/";
    $user_options{-o_1} = $output_directories[1];
    $output_directories[2] = "$user_options{-o}final_RAxML_outputs/";
    $user_options{-o_2} = $output_directories[2];
    $output_directories[3] = "$user_options{-o}final_outputs/";
    $user_options{-o_3} = $output_directories[3];
    
    foreach my $output_directory (@output_directories) {
        my $mkdir_command = "mkdir $output_directory";     
        system ($mkdir_command) unless (-e $output_directory);
    }
    
    ############ Edited by Young C. Song                             ###############
    ############ Users will define the subset set of reference       ###############
    ############ data they wish to use in the file mentioned here.   ###############
    my $ref_select_file = "data/tree_data/reference_selection.txt";
    my $ref_selected = "";
    my %ref_select_hash = ();
    
    
    ############ Edited by Young C. Song:  added block of if statements ###############
    ############ If Users decide to select subset of BLAST databases,   ###############
    ############ notify the users of such decision.                     ###############
    if ($user_options{-r} eq "a") {
       ;
    }else{
        my $dirname = dirname(__FILE__);
        print "We will use the set of reference data as mentioned in $dirname/$ref_select_file.\n";
        print "You have selected following reference data sets:\n";
        open (SELECTION, "$dirname/$ref_select_file") or die "Can't open  $dirname/$ref_select_file\n";
        
        while (<SELECTION>) {
            chomp $_;
            my $line = $_;
            
            if ($line =~ /(\S+)\t(\S+)\t(\S+(\s+\S+)*)$/){
                $ref_selected = $1;
                $ref_select_hash{$ref_selected} = "1";
                print $ref_selected."\n";
            }#end if
        }#end while
    
        close(SELECTION);
    }#end if/else
    #done.
    
    return (\%user_options);    
}

#####################################################
# split_input_file
#####################################################

sub split_input_file {
    my $package = shift;
    my $user_options = shift;
    
    my $input_file = $$user_options{-i};
    my $output_directory_var = $$user_options{-o_1};
    my $split_size = $$user_options{-l};
    
    my $input_file_name = "";
    if ($input_file =~ /\A.+\/(.+)/) {
        $input_file_name = $1;    
    }
    
    my %sequence_names = ();
    
    my $prefix_for_qsub = "";
    $prefix_for_qsub = "a" if ($$user_options{-c} eq "s");
    
    my $output_basis_name = "$output_directory_var$prefix_for_qsub$input_file_name";
    my $blast_input_file = "$output_basis_name"."_0.txt";
    my $formated_input_file = "$output_basis_name"."_formated.txt";
    
    $$user_options{formated_input_file} = $formated_input_file;
    
    my $continue_run = 1;
    
    my %blast_input_files = ();
    
    die "ERROR, no input file specified!\n" unless $input_file;
    
   
    open (INPUT, "$input_file") or die "Can't open  $input_file\n";
    open (OUTPUT, "> $blast_input_file") or die "Can't create $blast_input_file";
    open (OUTPUT2, "> $formated_input_file") or die "Can't create $formated_input_file"; 
    $blast_input_files{$blast_input_file} = 1; 
  
    my $count1 = -1;
    my $count2 = 0;
    my $is_not_fasta = 1;  #i.e. if the first line does not begin with ">", this will produce a n error.
    my $is_first_character = 1;
    my $total_nr_of_nucleotides = 0;
    my %character_classes = ();
    $character_classes{atgc} = 0;
    $character_classes{xn} = 0;
    $character_classes{notdef} = 0;
    
    while (<INPUT>) {
        chomp $_;
        my $line = $_;       
        if (/\A>/) {
            $is_not_fasta = 0 if $is_first_character;
            $line =~ s/[^a-zA-Z0-9.]/_/g; #filenames are generated from the sequence names later on, so we have to limit the characters.
            $line =~ s/\A_/>/; #change the first "_" back to ">"
            if ($line =~ /\A(.{2,100})/) {
                #note, RAxML supports only filenames with a max length of 125. Let us be cautious and further restrict it to 100.
                $line = $1;
            } else {
                die "ERROR: Something is wrong with the sequence name $line.\n";  
            }
            die "ERROR: The sequence name $line is not unique!\n" if exists $sequence_names{$line};
            $sequence_names{$line} = 1;
            $count1++;
            if ($count1 == $split_size) {
                $count1 = 0;
                $count2++;
                close OUTPUT;
                $blast_input_file = "$output_basis_name"."_$count2.txt";
                open (OUTPUT, "> $blast_input_file") or die "Can't create $blast_input_file"; 
                $blast_input_files{$blast_input_file} = 1;   
            }
        } else {
            $line =~ s/[^a-zA-Z]//g;
            $line = lc $line;
            my @characters = split //, $line;
            foreach my $character (@characters) {
                $total_nr_of_nucleotides++;
                if      (($character eq "a") || ($character eq "t") || 
                         ($character eq "g") || ($character eq "c")) {
                    $character_classes{atgc}++; 
                } elsif (($character eq "x") || ($character eq "n")) {
                    $character_classes{xn}++          
                } else {
                    $character_classes{notdef}++
                }                
            }
        }
        $is_first_character = 0;   
        die "ERROR your file does not appear to be a proper FASTA file!\n" if $is_not_fasta;
        print OUTPUT "$line\n"; 
        print OUTPUT2 "$line\n";    
    }
    if ($total_nr_of_nucleotides) {           
        if ($character_classes{xn} == $total_nr_of_nucleotides) {
            die "ERROR your sequence(s) contain only X or N!\n";    
        }
        if (($character_classes{atgc} / $total_nr_of_nucleotides) < 0.5) {
            die "ERROR your sequence(s) most likely contain no DNA!\n";
        }
    } else {
        die "ERROR: your input file appears to be corrupted. No sequences were found!\n";    
    }
    close OUTPUT;
    close OUTPUT2;
    return (\%blast_input_files);   
}

#####################################################
# read_cog_file
#####################################################

sub read_cog_file {
    my $package = shift;
    my $user_options = shift;
    my %cog_list = ();
    my %text_of_analysis_type = ();
    my $alignment_set = $$user_options{"-t"};
  
    use File::Basename;
    my $dirname = dirname(__FILE__);
    open (COGIN, "$dirname/data/tree_data/cog_list.txt") or die "ERROR: can't open $dirname/data/tree_data/cog_list.txt!\n";
    my $kind_of_cog = "";
    while (<COGIN>) {
        chomp $_;
        if (/\A#(.+)/) {
            $kind_of_cog = $1;
            next;    
        }
        if ($kind_of_cog eq "phylogenetic_cogs") {
            $cog_list{$kind_of_cog}{$_} = $alignment_set;
            $cog_list{"all_cogs"}{$_} = $alignment_set;
            my $text_inset = "";
            $text_inset = " based on the GEBA reference" if ($alignment_set eq "g");
            $text_inset = " focusing only on fungi" if ($alignment_set eq "i");
            $text_of_analysis_type{$alignment_set} = "Phylogenetic analysis$text_inset:";
        } elsif ($kind_of_cog eq "phylogenetic_rRNA_cogs") {
            my ($cog, $denominator, $text) = split /\t/;
            $cog_list{$kind_of_cog}{$cog} = $denominator;
            $cog_list{"all_cogs"}{$cog} = $denominator;
            $text_of_analysis_type{$denominator} = "Phylogenetic analysis, $text:";
        } elsif ($kind_of_cog eq "functional_cogs") {
            my ($cog, $denominator, $text) = split /\t/;
            $cog_list{$kind_of_cog}{$cog} = $denominator;
            $cog_list{"all_cogs"}{$cog} = $denominator;
            $text_of_analysis_type{$denominator} = "Functional analysis, $text:";
        }
    }    
    close COGIN;
    return (\%cog_list,\%text_of_analysis_type);
}

#####################################################
# get_non_wag_cogs
#####################################################

sub get_non_wag_cogs {
    my $package = shift;
    my %non_wag_cog_list = ();
    use File::Basename;
    my $dirname = dirname(__FILE__);
    open (COGIN, "$dirname/data/tree_data/non_wag_cogs.txt") or die "ERROR: can't open $dirname/data/tree_data/non_wag_cogs.txt!\n";
    my $denominator = "";
    while (<COGIN>) {
        chomp $_;
        if (/\A#(.+)/) {
            $denominator = $1;    
        } else {
            my ($cog, $model) = split /\t/;
            $non_wag_cog_list{$denominator}{$cog} = $model;    
        }
    }    
    close COGIN;
    return (\%non_wag_cog_list);
}


#####################################################
# read_species_translation_files
#####################################################

sub read_species_translation_files {
    my $package = shift;
    my $user_options = shift;
    my $cog_list = shift;
    my %tree_numbers_translation = ();
    
    my %translation_files = ();
    my $phylogenetic_denominator = $$user_options{"-t"};
 
    if ($phylogenetic_denominator eq "g") {
        $translation_files{$phylogenetic_denominator} = "data/tree_data/tax_ids_geba_tree.txt";
    } elsif ($phylogenetic_denominator eq "i") {
        $translation_files{$phylogenetic_denominator} = "data/tree_data/tax_ids_fungitr.txt";
    } else {
        $translation_files{$phylogenetic_denominator} = "data/tree_data/tax_ids_nr.txt";
    }
    foreach my $functional_cog (sort {$a cmp $b} keys %{$$cog_list{"functional_cogs"}}) {
        my $denominator = $$cog_list{"functional_cogs"}{$functional_cog};
        my $filename = "tax_ids_"."$functional_cog".".txt";
        $translation_files{$denominator} = "data/tree_data/$filename";
    }
    foreach my $phylogenetic_rRNA_cog (sort {$a cmp $b} keys %{$$cog_list{"phylogenetic_rRNA_cogs"}}) {
        my $denominator = $$cog_list{"phylogenetic_rRNA_cogs"}{$phylogenetic_rRNA_cog};
        my $filename = "tax_ids_"."$phylogenetic_rRNA_cog".".txt";
        $translation_files{$denominator} = "data/tree_data/$filename";
    }
    
#    use File::Basename;
    my $dirname = dirname(__FILE__);

    foreach my $denominator (sort {$a cmp $b} keys %translation_files) {
        my $filename = $translation_files{$denominator};
        open (IN, "$dirname/$filename") or die "ERROR: Can't open $dirname/$filename!\n";
        while (<IN>) {
            chomp $_;
            my ($number, $translation) = split /\t/;
            $tree_numbers_translation{$denominator}{$number} = $translation;    
        }
    } 
    return (\%tree_numbers_translation);
}


#####################################################
# blast
#####################################################


#####################################################
# blast
#####################################################

sub blast {
#    use Data::Dumper;
    my $package = shift;
    my $small_subroutines = shift;
    my $user_options = shift;
    my $blast_input_files = shift;
    
    my $output_directory_var = $$user_options{-o_1};
    my $reference_data_prefix = $$user_options{"reference_data_prefix"};
    
    my %blast_results_raw_files = ();
    my %blast_results_raw_files_raw = ();
    
    my %sun_grid_jobs = ();
    
    my $ref_select_file = "data/tree_data/reference_selection.txt";
    my $ref_selected = "";
    my %ref_select_hash = ();
    
    print "run BLAST\n";
    my $dirname = dirname(__FILE__);
    my $binaries_dirname = $user_options->{-e};
    
    my $alignment_data_dir = "$dirname/data/$reference_data_prefix"."alignment_data/";
    

    my @allfiles;
    if ($$user_options{-r} eq "a") {
        opendir (PATH, "$alignment_data_dir") or die "Error, this directory \"$alignment_data_dir\" does not exist!\n";
        @allfiles =  readdir PATH;
        closedir (PATH);
    } 
    else{
        open (SELECTION, "$ref_select_file") or die "Can't open $ref_select_file\n";
        while (<SELECTION>) {
            chomp $_;
            my $line = $_;
            
             if ($line =~ /(\S+)\t(\S+)\t(\S+(\s+\S+)*)$/){
                $ref_selected = $1;
                my $ref_selected_file = $ref_selected.".fa";
                $ref_select_hash{$ref_selected_file} = "1";
                #print $ref_selected_file."****\n";
              }#end if
         }#end while
         close(SELECTION);
         @allfiles = keys %ref_select_hash;

     }
    my @files ;

    my %file_info = ();
    for my $file (@allfiles) {
       if(! ($file =~ /\.fa\Z/) ){ next; }

       if( $file =~ /rRNA/i){  
          $file_info{$file} = { 'prog'=>'blastn', 'marker'=>'rRNA_'};
       }
       else {
          $file_info{$file} = { 'prog'=>'blastx', 'marker'=>''};
       }
       $file_info{$file}{'dbfile'} = "$alignment_data_dir"."$file";
    }
    
#    print Dumper(\%file_info);
#    print Dumper($blast_input_files);


    foreach my $blast_input_file (sort {$a cmp $b} keys %$blast_input_files) {      
            my $blast_input_file_name = "";
            if ($blast_input_file =~ /\A.+\/(.+)\.txt\Z/) {
                $blast_input_file_name = $1; 
            } else {
                die "ERROR something is wrong with the directory of the BLAST input file!\n";    
            }


            my $blast_result_raw_file = "$output_directory_var/$blast_input_file_name"."rRNA_"."BLAST_results_raw.txt";
            $blast_results_raw_files_raw{$blast_result_raw_file} = $blast_input_file;
        
            my $blast_command = "";
            my $databases = "";
            foreach my $file (keys %file_info) {
                if($file_info{$file}{'prog'} eq 'blastn' ) {
                    $databases = $file_info{$file}{'dbfile'};
                }
                else {
                     next;
                }

                $blast_command = $binaries_dirname."/blastn"." -query ".$blast_input_file;
                $blast_command .= " -db ".$databases." -evalue 0.01";
                $blast_command .= " -max_target_seqs 20000"." -dbsize 1000000"." -outfmt 6";
                $blast_command .= " >> ".$blast_result_raw_file;
    
                if ($$user_options{-c} eq "s") {
                       my $qsub_id = $small_subroutines->do_sun_grid_computing($user_options,$blast_input_file_name,$blast_command);
                       $sun_grid_jobs{$qsub_id} = 1;
                } 
                else {
                    system($blast_command);  
                }  
            }


            $blast_result_raw_file = "$output_directory_var/$blast_input_file_name"."BLAST_results_raw.txt";
            $blast_results_raw_files_raw{$blast_result_raw_file} = $blast_input_file;
            $blast_command = "";
            $databases = "";
            foreach my $file (keys %file_info) {
                if($file_info{$file}{'prog'} eq 'blastx') {
                    $databases = " ".$file_info{$file}{dbfile};
                }
                else {
                    next;
                }

                $blast_command = $binaries_dirname."/blastx". " -query ".$blast_input_file;
                $blast_command .= " -db ".$databases." -matrix BLOSUM62"." -evalue 0.01";
                $blast_command .= " -max_target_seqs 20000"." -dbsize 1000000"." -outfmt 6";
                $blast_command .= " >> ".$blast_result_raw_file;

                 if ($$user_options{-c} eq "s") {
                    my $qsub_id = $small_subroutines->do_sun_grid_computing($user_options,$blast_input_file_name,$blast_command);
                    $sun_grid_jobs{$qsub_id} = 1;
                 } 
                 else {
                    system($blast_command);  
                 }  
            }
        
    }
    
    if ($$user_options{-c} eq "s") {
        $small_subroutines->halt_mltreemap_until_sungrid_is_done($user_options,\%sun_grid_jobs,\%blast_results_raw_files_raw);
    }
    
    foreach my $blast_result_raw_file (sort {$a cmp $b} keys %blast_results_raw_files_raw) {
        my $blast_input_file = $blast_results_raw_files_raw{$blast_result_raw_file};
        open (CHECK, "$blast_result_raw_file");    
        if (<CHECK>) {
            $blast_results_raw_files{$blast_result_raw_file} = 1;
        } else {
            my $remove_command0 = "rm $blast_result_raw_file";
            system ($remove_command0);     
        }     
        close CHECK; 
        if (-e $blast_input_file) {
            my $remove_command1 = "rm $blast_input_file";
            system ($remove_command1);
        }
    }
    
    #TODO if wished, reintegrate the special part for functional genes:
    #$blast_command .= "-F F" if $is_functional;
    #$blast_command .= "-F T" if !$is_functional;
    
    return (\%blast_results_raw_files);                  
}


#####################################################
# read_blast_file
#####################################################

sub read_blast_file {  
    my $package = shift;
    my $small_subroutines = shift;
    my $user_options = shift;
    my $blast_results_raw_files = shift;
    my $cog_list = shift;
    my %blast_hits_purified = ();
    my $bitscore_min = $$user_options{-s};
    my $output_directory_var = $$user_options{-o_1};

    foreach my $blast_results_raw_file (sort {$a cmp $b} keys %$blast_results_raw_files) {
        my %blast_hits_raw = ();
        my $blast_result_raw_identifier = 0;
        
#print "opening "."$blast_results_raw_file"."\n";
        open (IN, "$blast_results_raw_file") or die "ERROR: $blast_results_raw_file does not exist!!\n";   
        while (<IN>) {
            chomp $_;
            my ($contig, $detailed_cog, undef, undef, undef, undef, $seq_start, $seq_end, $ref_start, $ref_end, undef, $bitscore) = split /\t/;
            next if ($bitscore <= $bitscore_min);
            my $temp_var = "";
            my $direction = "forward";
            if ($ref_start > $ref_end) {
                $temp_var = $ref_start;
                $ref_start = $ref_end;
                $ref_end = $temp_var;
                $direction = "reverse";  
            }
            if ($seq_start > $seq_end) {
                $temp_var = $seq_start;
                $seq_start = $seq_end;
                $seq_end = $temp_var;
                die "ERROR: parsing error with the BLAST results. Please notify the authors\n" if ($direction eq "reverse");
                $direction = "reverse";   
            }
            my $cog = "";
            if ($detailed_cog =~ /(.{7})\Z/) {
                $cog = $1;    
            } else {
                die "ERROR: could not detect the COG of sequence $detailed_cog!\n";    
            }
                        
            $blast_hits_raw{$contig}{$blast_result_raw_identifier}{"bitscore"} = $bitscore;
            $blast_hits_raw{$contig}{$blast_result_raw_identifier}{"cog"} = $cog;
            $blast_hits_raw{$contig}{$blast_result_raw_identifier}{"seq_start"} = $seq_start;
            $blast_hits_raw{$contig}{$blast_result_raw_identifier}{"seq_end"} = $seq_end;
            $blast_hits_raw{$contig}{$blast_result_raw_identifier}{"direction"} = $direction;
            $blast_hits_raw{$contig}{$blast_result_raw_identifier}{"validity"} = 1;
            $blast_result_raw_identifier++;
        }
        close IN;
       
        #purify Blast hits
 
        foreach my $contig (sort {$a cmp $b} keys %blast_hits_raw) {
            my $identifier = 0;
            my %pairs_allready_seen = ();       
            foreach my $base_blast_result_raw_identifier (sort {$a <=> $b} keys %{$blast_hits_raw{$contig}}) {
            
                my $base_bitscore = $blast_hits_raw{$contig}{$base_blast_result_raw_identifier}{"bitscore"};
                my $base_cog = $blast_hits_raw{$contig}{$base_blast_result_raw_identifier}{"cog"};
                my $base_start = $blast_hits_raw{$contig}{$base_blast_result_raw_identifier}{"seq_start"};
                my $base_end = $blast_hits_raw{$contig}{$base_blast_result_raw_identifier}{"seq_end"};
                my $direction = $blast_hits_raw{$contig}{$base_blast_result_raw_identifier}{"direction"};
            
                next unless ($base_bitscore >= $bitscore_min);
                #if it's only a sidecog, its not valid (but it can still cancel out other hits, so continue.
                $blast_hits_raw{$contig}{$base_blast_result_raw_identifier}{"validity"} = 0 unless (exists $$cog_list{"all_cogs"}{$base_cog});                   
                            
                #ok, now we have the complete information about the blast hit (base). Now compare it to the others                           
                #BTW: there might be several opinions to do this. But this is how it works:
                #   A                  C
                #------------------- -----
                #        ---------------
                #              B
                # A kills B, B kills C. (Another approach would be to let C live... But i suspect such C's to be crap anyway.)
                        
                foreach my $check_blast_result_raw_identifier (sort {$a <=> $b} keys %{$blast_hits_raw{$contig}}) {
                    next if ($base_blast_result_raw_identifier == $check_blast_result_raw_identifier);
                    #activating the next lines would increase the speed of this step, but costs MUCH memory.
                    #if (exists $pairs_allready_seen{$check_blast_result_raw_identifier}) {
                    #    next if (exists $pairs_allready_seen{$check_blast_result_raw_identifier}{$base_blast_result_raw_identifier});
                    #}
                    #$pairs_allready_seen{$base_blast_result_raw_identifier}{$check_blast_result_raw_identifier} = 1;
                
                    my $check_bitscore = $blast_hits_raw{$contig}{$check_blast_result_raw_identifier}{"bitscore"};
                    my $check_cog = $blast_hits_raw{$contig}{$check_blast_result_raw_identifier}{"cog"};
                    my $check_start = $blast_hits_raw{$contig}{$check_blast_result_raw_identifier}{"seq_start"};
                    my $check_end = $blast_hits_raw{$contig}{$check_blast_result_raw_identifier}{"seq_end"};
                    my $direction = $blast_hits_raw{$contig}{$check_blast_result_raw_identifier}{"direction"};
                                                     
                    #ok, now we have the complete information about the blast hit, which is to be compared (check)...
                    my $base_length = $base_end - $base_start;
                    my $check_length = $check_end - $check_start;                         
                    my %info = ();    
                    $info{"base"}{"start"} = $base_start;
                    $info{"base"}{"end"} = $base_end;
                    $info{"check"}{"start"} = $check_start;
                    $info{"check"}{"end"} = $check_end;
                                        
                    my $overlap = $small_subroutines->calculate_overlap(\%info);
                    if ($overlap) {
                        if ((($overlap / $base_length) > 0.5) && ($base_bitscore < $check_bitscore)) {
                            $blast_hits_raw{$contig}{$base_blast_result_raw_identifier}{"validity"} = 0;
                        } elsif ((($overlap / $check_length) > 0.5) && ($check_bitscore < $base_bitscore)){
                            $blast_hits_raw{$contig}{$check_blast_result_raw_identifier}{"validity"} = 0;
                        } elsif ($base_start == $check_start && $base_end == $check_end){ #if both are the same keep only the one with the smaller identifier.
                            if ($check_blast_result_raw_identifier > $base_blast_result_raw_identifier) {
                                $blast_hits_raw{$contig}{$check_blast_result_raw_identifier}{"validity"} = 0;
                            } else {
                                $blast_hits_raw{$contig}{$base_blast_result_raw_identifier}{"validity"} = 0;
                            }
                        }
                    }
                }
                if ($blast_hits_raw{$contig}{$base_blast_result_raw_identifier}{"validity"}) {
                    $blast_hits_purified{$contig}{$identifier}{"bitscore"} = $base_bitscore;
                    $blast_hits_purified{$contig}{$identifier}{"cog"} = $base_cog;
                    $blast_hits_purified{$contig}{$identifier}{"start"} = $base_start;
                    $blast_hits_purified{$contig}{$identifier}{"end"} = $base_end;
                    $blast_hits_purified{$contig}{$identifier}{"direction"} = $direction;
                    $blast_hits_purified{$contig}{$identifier}{"is_allready_placed"} = 0;
                    $identifier++;                                                      
                }
            }
        }
    }
    #print the blast results for each contig:
    foreach my $contig (sort {$a cmp $b} keys %blast_hits_purified) {
        my $outfile = "$output_directory_var/$contig"."_blast_result_purified.txt";
        open (OUT, "> $outfile") or die "Can't create $outfile\n";
        my %sorting_hash = ();
        foreach my $identifier (sort {$a <=> $b} keys %{$blast_hits_purified{$contig}}) {
            my $bitscore = $blast_hits_purified{$contig}{$identifier}{"bitscore"};
            $sorting_hash{$bitscore}{$identifier} = 1;
        }
        foreach my $bitscore (sort {$b <=> $a} keys %sorting_hash) {
            foreach my $identifier (sort {$a <=> $b} keys %{$sorting_hash{$bitscore}}) {        
                my $cog      = $blast_hits_purified{$contig}{$identifier}{"cog"};
                my $start    = $blast_hits_purified{$contig}{$identifier}{"start"};
                my $end      = $blast_hits_purified{$contig}{$identifier}{"end"};
                my $direction= $blast_hits_purified{$contig}{$identifier}{"direction"};
                print OUT "$contig\t$start\t$end\t$direction\t$cog\t$bitscore\n";
            }
        }
        close OUT;
    }
    #done
    return (\%blast_hits_purified);
}

#####################################################
# produce_shortened_sequence_files_for_genewise
#####################################################

sub produce_shortened_sequence_files_for_genewise {
    my $package = shift;
    my $user_options = shift;
    my $blast_hits_purified = shift;
    my $output_directory_var = $$user_options{-o_1};
    my $input_file = $$user_options{formated_input_file};
    
    my $flanking_length = 1000; #recommended: 1000
    
    my %prae_contig_coordinates = ();
    my %contig_coordinates = ();
    my %shortened_sequence_files = ();
    
    #print "prepare Genewise\n";
    #get the information about the relevant parts of the sequences
    
    foreach my $contig (sort {$a cmp $b} keys %$blast_hits_purified) {
        foreach my $base_identifier (sort {$a cmp $b} keys %{$$blast_hits_purified{$contig}}) {
            next if ($$blast_hits_purified{$contig}{$base_identifier}{"cog"} =~ "rRNA"); # we deal with the rRNA hits at another place.
            next if $$blast_hits_purified{$contig}{$base_identifier}{"is_allready_placed"};
            $$blast_hits_purified{$contig}{$base_identifier}{"is_allready_placed"} = 1;
                       
            my $base_start = $$blast_hits_purified{$contig}{$base_identifier}{"start"} - $flanking_length;
            my $base_end = $$blast_hits_purified{$contig}{$base_identifier}{"end"} + $flanking_length;
            
            my $nr_of_blast_hits = keys (%{$$blast_hits_purified{$contig}});
            
            for (my $check_identifier = 0; $check_identifier < $nr_of_blast_hits; $check_identifier++) {
                next if ($$blast_hits_purified{$contig}{$check_identifier}{"cog"} =~ "rRNA"); # we deal with the rRNA hits at another place.
                next if $$blast_hits_purified{$contig}{$check_identifier}{"is_allready_placed"};
                my $check_start = $$blast_hits_purified{$contig}{$check_identifier}{"start"} - $flanking_length;
                my $check_end = $$blast_hits_purified{$contig}{$check_identifier}{"end"} + $flanking_length;
                
                #now check, if they overlap. Note: the small subroutine "calculate_overlap" can't be used here because here we merge_stuff.
                if (($base_start <= $check_start) && ($check_start <= $base_end) && ($base_end <= $check_end)) {
                    # Base     ----------       
                    # Check        ----------   
                    $base_end = $check_end;
                    $$blast_hits_purified{$contig}{$check_identifier}{"is_allready_placed"} = 1;
                    $check_identifier = -1;
                    next;
                }
                elsif (($base_start <= $check_start) && ($check_end <= $base_end)) {
                    # Base     ----------       
                    # Check        ---          
                    $$blast_hits_purified{$contig}{$check_identifier}{"is_allready_placed"} = 1;
                    $check_identifier = -1;
                    next;
                }
                elsif (($check_start <= $base_start) && ($base_start <= $check_end) && ($check_end <= $base_end)) {
                    # Base           ------     
                    # Check   ----------        
                    $base_start = $check_start;
                    $$blast_hits_purified{$contig}{$check_identifier}{"is_allready_placed"} = 1;
                    $check_identifier = -1;
                    next;
                }
                elsif (($check_start <= $base_start) && ($base_end <= $check_end)) {
                    # Base           ------     
                    # Check   ---------------- 
                    $base_start = $check_start;
                    $base_end = $check_end;  
                    $$blast_hits_purified{$contig}{$check_identifier}{"is_allready_placed"} = 1;
                    $check_identifier = -1;
                    next;          
                }
            }
            $prae_contig_coordinates{$contig}{$base_start}{$base_end} = 1;
        }
    }

    #now produce the genewise input files
    open (IN, "$input_file") or die "ERROR: Can't open $input_file!\n";
    my $contig_name = "";
    my $sequence = "";
    
    while (<IN>) {      
        chomp $_;
        $_ =~ s/\s/_/g;            
        if (/>(.+)/ || eof(IN)) {
            if (eof(IN)) {
                $sequence .= $_;
            }
            if (exists $prae_contig_coordinates{$contig_name}) {
                
                #ok, we have a sequence with a hit. Save the sequence file as well as a shortened version of it.
                
                my $sequence_length = length ($sequence);
                my @nucleotides = split //, $sequence;
                my $shortened_sequence = "";
                
                #start searching for the information to shorten the file.
                foreach my $start_B (sort {$a <=> $b} keys %{$prae_contig_coordinates{$contig_name}}) {
                    foreach my $end_B (sort {$a <=> $b} keys %{$prae_contig_coordinates{$contig_name}{$start_B}}) {
                        
                        #ok, now we have all information about the hit. Correct start and end if needed:                     
                        $start_B = 0 if ($start_B < 0);
                        $end_B = $sequence_length -1 if ($end_B >= $sequence_length);
                        
                        #Note: Genewise (GW) positions start with 1, Blast (B) positions with 0 -> thus differenciate between start_B and start_GW
                        my $shortened_start_GW = length ($shortened_sequence) + 1;
                        my $count = -1;
                        foreach my $nucleotide (@nucleotides) {
                            $count++;                           
                            next unless (($count >= $start_B) && ($count <= $end_B));
                            $shortened_sequence .= $nucleotide;
                        }
                                          
                        my $shortened_end_GW = length($shortened_sequence);
                        my $addition_factor = ($start_B + 1) - $shortened_start_GW; #$start_B + 1 == $start_GW
                        $contig_coordinates{$contig_name}{$shortened_start_GW}{$shortened_end_GW} = $addition_factor;
                    } 
                }
                #done
                
                open (OUT, "> $output_directory_var/$contig_name"."_sequence.txt")
                or die "ERROR: Can't create $output_directory_var/$contig_name"."_sequence.txt!\n"; 
                print OUT ">$contig_name\n$sequence";
                close OUT;
                
                my $prefix_for_qsub = "";
                $prefix_for_qsub = "a" if ($$user_options{-c} eq "s");
                               
                open (OUT2, "> $output_directory_var/$prefix_for_qsub$contig_name"."_sequence_shortened.txt")
                or die "ERROR: Can't create $output_directory_var/$prefix_for_qsub$contig_name"."_sequence_shortened.txt!\n"; 
                print OUT2 ">$contig_name\n$shortened_sequence";
                close OUT2;
                
                $shortened_sequence_files{"$output_directory_var/$prefix_for_qsub$contig_name"."_sequence_shortened.txt"} = $contig_name;
            }
                    
            $contig_name = $1;
            $sequence = "";
        } else {
            $sequence .= $_;    
        }      
    }
    close IN;
    return (\%contig_coordinates, \%shortened_sequence_files);         
}

#####################################################
# start_genewise
#####################################################

sub start_genewise {   
    my $package = shift;
    my $user_options = shift;
    my $shortened_sequence_files = shift;
    my $blast_hits_purified = shift;
    my $output_directory_var = $$user_options{-o_1};
    my $reference_data_prefix = $$user_options{"reference_data_prefix"};
    my %genewise_outputfiles = ();
    my %genewise_outputfiles_for_sun_grid_control = ();
    my %sun_grid_jobs = ();
    
    print "run Genewise\n";
    my $binaries_dirname = $user_options->{-e};
    
    foreach my $shortened_sequence_file (sort {$a cmp $b} keys %$shortened_sequence_files) {
        my $contig = $$shortened_sequence_files{$shortened_sequence_file};
        foreach my $identifier (sort {$a cmp $b} keys %{$$blast_hits_purified{$contig}}) {            
            my $cog = $$blast_hits_purified{$contig}{$identifier}{"cog"};
            
            my $prefix_for_qsub = "";
            $prefix_for_qsub = "a" if ($$user_options{-c} eq "s");
            
            my $genewise_outputfile = "$output_directory_var$prefix_for_qsub$contig"."_$cog"."_genewise.txt";
            $genewise_outputfiles{$contig}{$genewise_outputfile} = 1;
            $genewise_outputfiles_for_sun_grid_control{$genewise_outputfile} = 1;
            my $dirname = dirname(__FILE__);
            my $genewise_command = "$binaries_dirname/genewise $dirname/data/$reference_data_prefix"."hmm_data/$cog.hmm $shortened_sequence_file ";
            $genewise_command .= " -init local -quiet -gene $dirname/data/genewise_support_files/human.gf -matrix $dirname/data/genewise_support_files/blosum62.bla";
            $genewise_command .= " -codon $dirname/data/genewise_support_files/codon.table -hmmer -subs 0.01 -indel 0.01 -gap 11 -ext 1 -both -pep -sum > ";
            $genewise_command .= "$genewise_outputfile";
#          print "Genewise....".$genewise_command."\n";
            
            if ($$user_options{-c} eq "s") {
                my $qsub_id = $small_subroutines->do_sun_grid_computing($user_options,$shortened_sequence_file,$genewise_command);
                $sun_grid_jobs{$qsub_id} = 1;
            } else {
               system ($genewise_command);  
            }
        }
    }
    
    if ($$user_options{-c} eq "s") {
        $small_subroutines->halt_mltreemap_until_sungrid_is_done($user_options,\%sun_grid_jobs,\%genewise_outputfiles_for_sun_grid_control);
    }
     
    return (\%genewise_outputfiles);
}

#####################################################
# parse_genewise_results
#####################################################

sub parse_genewise_results {
    my $package = shift;
    my $small_subroutines = shift;
    my $user_options = shift;
    my $contig_coordinates = shift;
    my $genewise_outputfiles = shift;
    my $output_directory_var = $$user_options{-o_1};
        
    my %genewise_summary_files = ();
      
    foreach my $contig (sort {$a cmp $b} keys %$genewise_outputfiles) {
        
        my %genewise_results_raw = ();
        my %genewise_results = ();
        my $at_least_one_hit = 0;
        my $count = 0;
        
        #do the actual parsing      
        foreach my $genewise_outputfile (sort {$a cmp $b} keys %{$$genewise_outputfiles{$contig}}) {
            open (IN, "$genewise_outputfile") or die "Can't open $genewise_outputfile\n";
                
            my $header_count = 0;
            my $sequence_count = -1;
              
            while (<IN>) {
                chomp $_;
                if (/\A\d/) {
                    my ($bitscore, $query, undef, undef, undef, $start, $end) = split / +/, $_;
                    $at_least_one_hit = 1 if $query;
                    my $direction = "forward";
                    if ($start > $end) {
                        my $temp_var = $start;
                        $start = $end;
                        $end = $temp_var;
                        $direction = "reverse";   
                    }
                    
                    #correct the positions (Genewise has been run on a shortened sequence, thus calculate the true positions)
                    
                    foreach my $coords_start (sort {$a <=> $b} keys %{$$contig_coordinates{$contig}}) {
                        if ($start >= $coords_start) {
                            foreach my $coords_end (sort {$a <=> $b} keys %{$$contig_coordinates{$contig}{$coords_start}}) {
                                if  ($end <= $coords_end) {
                                    my $addition_factor = $$contig_coordinates{$contig}{$coords_start}{$coords_end};
                                    $start += $addition_factor;
                                    $end += $addition_factor;
                                    last;
                                }
                            }
                        }
                    }
                    
                    #done
                    
                    $genewise_results_raw{$contig}{$genewise_outputfile}{$header_count}{"start"} = $start;
                    $genewise_results_raw{$contig}{$genewise_outputfile}{$header_count}{"end"} = $end;
                    $genewise_results_raw{$contig}{$genewise_outputfile}{$header_count}{"cog"} = $query;
                    $genewise_results_raw{$contig}{$genewise_outputfile}{$header_count}{"bitscore"} = $bitscore;
                    $genewise_results_raw{$contig}{$genewise_outputfile}{$header_count}{"direction"} = $direction;
                    $header_count++;
                } elsif (/\A>/) {
                    $sequence_count++;           
                    $genewise_results_raw{$contig}{$genewise_outputfile}{$sequence_count}{"sequence"} = "";
                } elsif ((/\A\w/) && (! (/\ABits/)) && (! (/\AMaking/))) {
                    #/\AMaking/ is compared, because genewise makes some (unwanted) comments about the treatment of introns.
                    $genewise_results_raw{$contig}{$genewise_outputfile}{$sequence_count}{"sequence"} .= $_;
                }
            }    
            close IN;    
        }      
        #done
        #do the purifying step.
        next unless $at_least_one_hit;
        
        foreach my $base_genewise_outputfile (sort {$a cmp $b} keys %{$genewise_results_raw{$contig}}) {
            foreach my $base_count (sort {$a cmp $b} keys %{$genewise_results_raw{$contig}{$base_genewise_outputfile}}) {
                my $base_start = $genewise_results_raw{$contig}{$base_genewise_outputfile}{$base_count}{"start"};
                my $base_end = $genewise_results_raw{$contig}{$base_genewise_outputfile}{$base_count}{"end"};
                my $base_cog = $genewise_results_raw{$contig}{$base_genewise_outputfile}{$base_count}{"cog"};
                my $base_bitscore = $genewise_results_raw{$contig}{$base_genewise_outputfile}{$base_count}{"bitscore"};
                my $base_direction = $genewise_results_raw{$contig}{$base_genewise_outputfile}{$base_count}{"direction"};
                my $base_sequence = $genewise_results_raw{$contig}{$base_genewise_outputfile}{$base_count}{"sequence"};       
                unless ((defined $base_cog) && (defined $base_start) && (defined $base_end)) {
                    my $errorstring = "ERROR: the file \"$base_genewise_outputfile\" cannot be parsed!\n";
                    $errorstring .= "Please contact the outhors about it. As a quick solution to the problem, try to remove ";
                    $errorstring .= "the sequence, which produced this hit, from your input file.\n";
                    die "$errorstring";
            
                }
                my $base_length = $base_end - $base_start;
                my $is_valid = 1;
                
                foreach my $check_genewise_outputfile (sort {$a cmp $b} keys %{$genewise_results_raw{$contig}}) {
                    foreach my $check_count (sort {$a cmp $b} keys %{$genewise_results_raw{$contig}{$check_genewise_outputfile}}) {
                        next if ($base_count == $check_count);
                        my $check_start = $genewise_results_raw{$contig}{$check_genewise_outputfile}{$check_count}{"start"};
                        my $check_end = $genewise_results_raw{$contig}{$check_genewise_outputfile}{$check_count}{"end"};
                        my $check_cog = $genewise_results_raw{$contig}{$check_genewise_outputfile}{$check_count}{"cog"};
                        unless ((defined $check_cog) && (defined $check_start) && (defined $check_end)) {
                            my $errorstring = "ERROR: the file \"$check_genewise_outputfile\" cannot be parsed!\n";
                            $errorstring .= "Please contact the outhors about it. As a quick solution to the problem, try to remove ";
                            $errorstring .= "the sequence, which produced this hit, from your input file.\n";
                            die "$errorstring"; 
                        }
                        my $check_length = $check_end - $check_start;
                        my %info = ();    
                        $info{"base"}{"start"} = $base_start;
                        $info{"base"}{"end"} = $base_end;
                        $info{"check"}{"start"} = $check_start;
                        $info{"check"}{"end"} = $check_end;
                                        
                        my $overlap = $small_subroutines->calculate_overlap(\%info);
                        
                        #ok, now we have all needed information. purify.
                        
                        if (($overlap / $base_length) > 0.5) {
                            if ($base_cog eq $check_cog) {
                                if ($base_length < $check_length) {
                                    # the maior difference between the hits is the length. Keep the longer.
                                    $is_valid = 0;
                                }
                            } elsif ($base_length < ($check_length / 2)) {
                                #it's not the same cog. Thus only skip this one if it is <1/2 the lenght of the other...
                                $is_valid = 0;
                            }
                        }
                        
                        if ($is_valid && ($base_cog eq $check_cog)) {
                            #ok, there was no overlap. But i want also remove sidehits of the same COG 
                            if ($base_length < ($check_length * 0.7)) {
                                $is_valid = 0;
                            }                       
                        }        
                        #done                       
                    }
                }
                
                if ($is_valid) {
                    $genewise_results{$contig}{$count}{"start"} = $base_start;
                    $genewise_results{$contig}{$count}{"end"} = $base_end;
                    $genewise_results{$contig}{$count}{"cog"} = $base_cog;
                    $genewise_results{$contig}{$count}{"direction"} = $base_direction;
                    $genewise_results{$contig}{$count}{"sequence"} = $base_sequence;
                    $count++;
                }
            }
        }
        
        #ok, parsing and purification is done. Now write the summary file
        next unless ($count); 
        my $genewise_summary_file = "$output_directory_var$contig"."_genewise_result_summary.txt";
        open (OUT, "> $genewise_summary_file") or die "Can't open $genewise_summary_file\n";
        $genewise_summary_files{$contig}{$genewise_summary_file} = 1;
        foreach my $count (sort {$a <=> $b} keys %{$genewise_results{$contig}}) {
            my $start = $genewise_results{$contig}{$count}{"start"};
            my $end = $genewise_results{$contig}{$count}{"end"};
            my $cog = $genewise_results{$contig}{$count}{"cog"};
            my $direction = $genewise_results{$contig}{$count}{"direction"};
            my $sequence = $genewise_results{$contig}{$count}{"sequence"};
            print OUT "$cog\t$start\t$end\t$direction\t$sequence\n";
        }
        close OUT;      
    }
    return (\%genewise_summary_files);
}

#####################################################
# get_rRNA_hit_sequences
#####################################################

sub get_rRNA_hit_sequences {
    my $package = shift;
    my $user_options = shift;
    my $blast_hits_purified = shift;
    my $cog_list = shift;
    my $genewise_summary_files = shift;
    my $output_directory_var = $$user_options{-o_1};
    my $input_file = $$user_options{formated_input_file};
    
    my %contig_rRNA_coordinates = ();
    my %rRNA_hit_files = ();
    
    foreach my $contig (sort {$a cmp $b} keys %$blast_hits_purified) {
        
        my $prefix_for_qsub = "";
        $prefix_for_qsub = "a" if ($$user_options{-c} eq "s");
           
        #note: We skipped the Genewise step (we are dealing with rRNA) but we bring the rRNA files in the same structure as the Genewise summary files and bring them back into the ordinary pipeline.
        foreach my $identifier (sort {$a cmp $b} keys %{$$blast_hits_purified{$contig}}) {
            next unless ($$blast_hits_purified{$contig}{$identifier}{"cog"} =~ "rRNA"); # here we only deal with the rRNA hits.
            my $start = $$blast_hits_purified{$contig}{$identifier}{"start"};
            my $end = $$blast_hits_purified{$contig}{$identifier}{"end"};
            my $cog = $$blast_hits_purified{$contig}{$identifier}{"cog"};
            my $direction = $$blast_hits_purified{$contig}{$identifier}{"direction"};
            $contig_rRNA_coordinates{$contig}{$identifier}{"start"} = $start;
            $contig_rRNA_coordinates{$contig}{$identifier}{"end"} = $end;
            $contig_rRNA_coordinates{$contig}{$identifier}{"cog"} = $cog;
            $contig_rRNA_coordinates{$contig}{$identifier}{"direction"} = $direction;
            my $outfile_name = "$output_directory_var/$prefix_for_qsub"."$contig"."_rRNA_result_summary.txt";   
            $contig_rRNA_coordinates{$contig}{$identifier}{"outfile"} = $outfile_name;
            $$genewise_summary_files{$contig}{$outfile_name} = 1;             
            open (OUT, "> $outfile_name") or die "ERROR: Can't create $outfile_name.txt!\n"; 
            close OUT;
        }
    }
    
    #now produce the genewise input files
    open (IN, "$input_file") or die "ERROR: Can't open $input_file!\n";
    my $contig_name = "";
    my $sequence = "";
    
    while (<IN>) {      
        chomp $_;
        $_ =~ s/\s/_/g;            
        if (/>(.+)/ || eof(IN)) {
            if (eof(IN)) {
                $sequence .= $_;
            }
            if (exists $contig_rRNA_coordinates{$contig_name}) {
                
                #ok, we have a sequence with a hit. Save the sequence file as well as a shortened version of it.
                
                my $sequence_length = length ($sequence);
                my @nucleotides = split //, $sequence;
                my $shortened_sequence = "";
                
                #start searching for the information to shorten the file.
                
                foreach my $identifier (sort {$a <=> $b} keys %{$contig_rRNA_coordinates{$contig_name}}) {
                    my $start = $contig_rRNA_coordinates{$contig_name}{$identifier}{"start"};
                    my $end = $contig_rRNA_coordinates{$contig_name}{$identifier}{"end"};
                    my $cog = $contig_rRNA_coordinates{$contig_name}{$identifier}{"cog"};
                    my $direction = $contig_rRNA_coordinates{$contig_name}{$identifier}{"direction"};
                    my $outfile = $contig_rRNA_coordinates{$contig_name}{$identifier}{"outfile"};
                    my $denominator = $$cog_list{"all_cogs"}{$cog};
                    my $count = -1;
                    my $shortened_sequence = "";
                    foreach my $nucleotide (@nucleotides) {
                        $count++;                           
                        next unless (($count >= $start) && ($count <= $end));
                        $shortened_sequence .= $nucleotide;
                    }
                    if ($direction eq "reverse") {
                        #ok, our hit has been on the opposite strand of the reference.
                        #to get a proper alignment, we thus have to produce a negative strand version of the input
                        my @nucleotides2 = split //, $shortened_sequence;
                        $shortened_sequence = "";
                        @nucleotides2 = reverse (@nucleotides2);
                        foreach my $nucleotide (@nucleotides2) {
                            if ($nucleotide eq "t") {
                                $nucleotide = "a";
                            } elsif ($nucleotide eq "a") {
                                $nucleotide = "t";
                            } elsif ($nucleotide eq "c") {
                                $nucleotide = "g";
                            } elsif ($nucleotide eq "g") {
                                $nucleotide = "c";
                            }
                            $shortened_sequence .= $nucleotide;
                        }
                    }
                    
                    open (OUT, ">> $outfile") or die "ERROR: Can't create $outfile!\n";
                    print OUT "$cog\t$start\t$end\tn/a\t$shortened_sequence\n";
                    close OUT;
                    
                }
                
                #done
                
                open (OUT, "> $output_directory_var/$contig_name"."_sequence.txt")
                or die "ERROR: Can't create $output_directory_var/$contig_name"."_sequence.txt!\n"; 
                print OUT ">$contig_name\n$sequence";
                close OUT;
                
            }
                    
            $contig_name = $1;
            $sequence = "";
        } else {
            $sequence .= $_;    
        }      
    }
    close IN;
    return (\%contig_rRNA_coordinates, \%rRNA_hit_files);
    
}

#####################################################
# prepare_and_run_hmmalign
#####################################################

sub prepare_and_run_hmmalign {
    my $package = shift;
    my $user_options = shift; 
    my $genewise_summary_files = shift;
    my $cog_list = shift;

    my $output_directory_var = $$user_options{-o_1};
    my $reference_data_prefix = $$user_options{"reference_data_prefix"};
    
    my %hmmalign_singlehit_files = ();
    my %hmmalign_singlehit_files_for_sun_grid_control = ();
    my %sun_grid_jobs = ();
    
    my $binaries_dirname = $user_options->{-e};
    print "run hmmalign\n";
    
    foreach my $contig (sort {$a cmp $b} keys %$genewise_summary_files) {
        foreach my $genewise_summary_file (sort {$a cmp $b} keys %{$$genewise_summary_files{$contig}}) {
            open (IN, "$genewise_summary_file") or die "ERROR: Can't open $genewise_summary_file!\n";
            while (<IN>) {
                my ($cog, $start, $end, undef, $sequence) = split /\t/;
                my $denominator = $$cog_list{"all_cogs"}{$cog};
                my $f_contig = "$denominator"."_$contig";
                my $genewise_singlehit_file = "$output_directory_var"."$f_contig"."_$cog"."_$start"."_$end";
                $hmmalign_singlehit_files{$f_contig}{"$genewise_singlehit_file.mfa"} = 1;
                $hmmalign_singlehit_files_for_sun_grid_control{"$genewise_singlehit_file.mfa"} = 1;
                my $genewise_singlehit_file_fa = "$genewise_singlehit_file.fa";
                open (OUT, "> $genewise_singlehit_file_fa") or die "Can't create $genewise_singlehit_file.fa\n";
                print OUT ">query\n$sequence";
                close OUT;
                
                #run hmmalign 
                my $hmmalign_command = "";
                  
                my $dirname = dirname(__FILE__);
                $hmmalign_command = "$binaries_dirname/hmmalign -m --mapali $dirname/data/$reference_data_prefix"."alignment_data/$cog.fa ";
                $hmmalign_command .= "--outformat Clustal $dirname/data/$reference_data_prefix"."hmm_data/$cog.hmm ";
                $hmmalign_command .= "$genewise_singlehit_file_fa > $genewise_singlehit_file.mfa";
                
#print "HMMALIGN ..".$hmmalign_command."\n";
                
                if ($$user_options{-c} eq "s") {
                    my $qsub_id = $small_subroutines->do_sun_grid_computing($user_options,$genewise_singlehit_file_fa,$hmmalign_command);
                    $sun_grid_jobs{$qsub_id} = 1;
                } else {
                    system ($hmmalign_command);  
                }
            }
            close IN;
        }
    }
    
    if ($$user_options{-c} eq "s") {
        $small_subroutines->halt_mltreemap_until_sungrid_is_done($user_options,\%sun_grid_jobs,\%hmmalign_singlehit_files_for_sun_grid_control);
    }
    
    return (\%hmmalign_singlehit_files);
}

#####################################################
# concatenate_hmmalign_singlehits_files
#####################################################

sub concatenate_hmmalign_singlehits_files {
    my $package = shift;
    my $user_options = shift; 
    my $hmmalign_singlehit_files = shift;
    my $non_wag_cog_list = shift;
    
    my $output_directory_var = $$user_options{-o_1};
    my %concatenated_mfa_files = ();
    my %nrs_of_sequences = ();
    my %models_to_be_used = ();

    foreach my $f_contig (sort {$a cmp $b} keys %$hmmalign_singlehit_files) {
        #print $f_contig."\n";
        my %sequences = ();
        my $query_sequence = "";
        my $model_to_be_used = "";
        my $denominator = "";
        if ($f_contig =~ /^([^_]+)_([^_]+).*$/) {
            $denominator = $1;
            #print $denominator."\n";
        } else {
            die "ERROR the analysis type could not be parsed from $f_contig!\n";    
        }  
        foreach my $hmmalign_singlehit_file (sort {$a cmp $b} keys %{$$hmmalign_singlehit_files{$f_contig}}) {
            open (IN, "$hmmalign_singlehit_file") or die "ERROR: Can't open $hmmalign_singlehit_file!\n";
            my $reached_data_part = 0;
            #determine the best aa model
            my $cog = "";
            if ($hmmalign_singlehit_file =~ /\A.+_(.{7})_\d+_\d+\.mfa\Z/) {
                $cog = $1;    
            } else {
                die "ERROR the COG could not be parsed from $hmmalign_singlehit_file!\n";   
            }
            
            if ((exists $$non_wag_cog_list{$denominator}{$cog}) && ($model_to_be_used ne "PROTGAMMAWAG")) {
                $model_to_be_used = "$$non_wag_cog_list{$denominator}{$cog}";
            } else {
                $model_to_be_used = "PROTGAMMAWAG";
            }
            #done
            while (<IN>) {
                chomp $_;
                                      
                $reached_data_part = 1 if (/query/);
                next unless $reached_data_part;
    
                if (/\A(.+) (\S+)\Z/) {
                    my $name_long = $1;
                    my $sequence_part = $2;
                    my $sequence_name = "";
                    if ($name_long =~ /query/) {
                        $query_sequence .= $sequence_part;
                    }elsif ($name_long =~ /(\d+)_/) {
                        $sequence_name = $1;                   
                        if (exists $sequences{$sequence_name}) {
                            $sequences{$sequence_name} .= $sequence_part;    
                        } else {
                            $sequences{$sequence_name} = $sequence_part;
                        }
                    }
                }
            }
            close IN;
        }
        $models_to_be_used{$f_contig} = $model_to_be_used;
        $concatenated_mfa_files{$f_contig} = "$output_directory_var$f_contig.mfa";
        open (OUT, "> $output_directory_var$f_contig.mfa") or die "Can't create $output_directory_var$f_contig.mfa\n";
        #$query_sequence =~ s/\./X/g; #note: if we switch to trimAl, the sequences should be modified here
        #$query_sequence =~ s/\*/X/g;
        #$query_sequence =~ s/-/X/g;
        print OUT ">query\n$query_sequence\n";
        $nrs_of_sequences{$f_contig} = 1;
        foreach my $sequence_name (sort {$a <=> $b} keys %sequences) {
            $nrs_of_sequences{$f_contig}++;
            my $sequence = $sequences{$sequence_name};
            #$sequence =~ s/\./X/g;
            #$sequence =~ s/\*/X/g;
            #$sequence =~ s/-/X/g;
            print OUT ">$sequence_name\n$sequence\n";
        }
        close OUT;
    }  
    return(\%concatenated_mfa_files,\%nrs_of_sequences,\%models_to_be_used);
}

#####################################################
# start_gblocks
#####################################################

sub start_gblocks {
    my $package = shift;
    my $user_options = shift; 
    my $concatenated_mfa_files = shift;
    my $nrs_of_sequences = shift;
    my %gblocks_files = ();
    my %gblocks_files_for_sun_grid_control = ();
    my %sun_grid_jobs = ();
    
    print "run Gblocks\n";
    my $binaries_dirname = $user_options->{-e};
    
    foreach my $f_contig (sort {$a cmp $b} keys %$concatenated_mfa_files) {
        my $concatenated_mfa_file = $$concatenated_mfa_files{$f_contig};
        my $nr_of_sequences = $$nrs_of_sequences{$f_contig};
        my $min_flank_pos = int($nr_of_sequences * 0.55);
        my $gblocks_file = "$concatenated_mfa_file-gb";
        $gblocks_files{$f_contig} = $gblocks_file;
        $gblocks_files_for_sun_grid_control{$gblocks_file} = 1;
        my $dirname = dirname(__FILE__);
        my $gblocks_command = "$binaries_dirname/Gblocks ";
        $gblocks_command .= "$concatenated_mfa_file ";
        $gblocks_command .= "-t=p -s=y -u=n -p=t -b3=15 -b4=3 -b5=h -b2=$min_flank_pos  > /dev/null";
        if ($$user_options{-c} eq "s") {
            my $qsub_id = $small_subroutines->do_sun_grid_computing($user_options,$gblocks_file,$gblocks_command);
            $sun_grid_jobs{$qsub_id} = 1;
        } else {
#print "Gblocks ...".$gblocks_command."\n";
            system ($gblocks_command);  
        }        
    }
    
    if ($$user_options{-c} eq "s") {
        $small_subroutines->halt_mltreemap_until_sungrid_is_done($user_options,\%sun_grid_jobs,\%gblocks_files_for_sun_grid_control);
    }
    
    return (\%gblocks_files);
}

#####################################################
# produce_phy_file
#####################################################

sub produce_phy_file {  
    my $package = shift;
    my $user_options = shift;
    my $gblocks_files = shift;
    my $nrs_of_sequences = shift;
    
    my $output_directory_var = $$user_options{-o_1};
    my $minimum_length_after_gblocks = $$user_options{-g};
    
    my %phy_files = ();
    my %sequence_lengths = ();
    
    foreach my $f_contig (sort {$a cmp $b} keys %$gblocks_files) {
        
        my %sequences_for_phy = ();
        my $do_not_continue = 0;
        my %sequences_raw = ();
        my $gblocks_file = $$gblocks_files{$f_contig};
        my $seq_name = "";
        open (IN, "$gblocks_file") or die "ERROR: Can't open $gblocks_file!\n";
        while (<IN>) {
            chomp $_;
            if (/\A>(.+)/) {
                $seq_name = $1;
                my $error_msg = "ERROR: your reference alignment contains an element with the number -666.";
                $error_msg .= " Please change it, because this number is needed for internal purposes.\n";
                die $error_msg if ($seq_name eq "-666");
                $seq_name = -666 if ($seq_name eq "query");
            } else {
                s/ //g;
                if (exists $sequences_raw{$seq_name}) {
                    $sequences_raw{$seq_name} .= $_; 
                } else {
                    $sequences_raw{$seq_name} = $_;
                }
            }
        }
        foreach my $seq_name (sort {$a <=> $b} keys %sequences_raw) {      
            next if $do_not_continue;
            my $sequence = $sequences_raw{$seq_name};
            my $count = 0;
            $sequence_lengths{$f_contig} = length($sequence);
            $sequence =~ s/\./X/g;
            $sequence =~ s/\*/X/g;
            $sequence =~ s/-/X/g;
            if ($sequence =~ /\AX+\Z/) {
                $sequence =~ s/X/V/;   #RAxML cannot deal with sequences containing only X, thus change the first X to V. 
            }
            if ($seq_name == -666) {
                my $seq_dummy = $sequence;
                $seq_dummy =~ s/X//g;
                if (length($seq_dummy) < $minimum_length_after_gblocks) {
                    $do_not_continue = 1;
                    my $exit_file_name = "$output_directory_var$f_contig"."_exit_after_Gblocks.txt";
                    open (OUT, "> $exit_file_name");
                    print OUT "final alignment after gblocks is to short (<$minimum_length_after_gblocks aa) ";
                    print OUT "- unsufficient number of marker gene residues in query sequence.\n";
                    close OUT;
                    next;    
                } 
            }
            while ($sequence =~ /(.{1,50})/g) {
                $sequences_for_phy{$f_contig}{$count}{$seq_name} = $1;
                $count++;  
            }
        }
    
        # print the phy file.
        next if $do_not_continue;

        my $phy_file_name = "$output_directory_var$f_contig.phy";
        $phy_files{$f_contig} = $phy_file_name;
        open (OUT, "> $phy_file_name") or die "Can't open $phy_file_name\n";
        my $nr_of_sequences = $$nrs_of_sequences{$f_contig};
        print OUT " $nr_of_sequences  $sequence_lengths{$f_contig}\n";
        foreach my $count (sort {$a <=> $b} keys %{$sequences_for_phy{$f_contig}}) {
            foreach my $seq_name (sort {$a <=> $b} keys %{$sequences_for_phy{$f_contig}{$count}}) {           
                my $sequence_part = $sequences_for_phy{$f_contig}{$count}{$seq_name};
                unless ($count) {
                    my $print_seqname = $seq_name;
                    $print_seqname = "query" if ($seq_name eq "-666");
                    print OUT "$print_seqname";
                    my $length = length ($print_seqname);
                    for (my $c = $length; $c < 10; $c++) {
                        print OUT " ";    
                    } 
                }
                print OUT "$sequence_part\n";
            }
            print OUT "\n";
        }
        close OUT;
    }
    return (\%phy_files);
}

#####################################################
# start_RAxML
#####################################################

sub start_RAxML {
    my $package = shift;
    my $small_subroutines = shift;
    my $user_options = shift;
    my $phy_files = shift;
    my $cog_list = shift;
    my $models_to_be_used = shift;
    my $raxml_option = $$user_options{"-f"};
    my $output_directory_var = $$user_options{-o_1};
    
    my %sun_grid_jobs = ();
    my %expected_raxml_outfiles = ();
    my %raxml_outfiles = ();
    
    print "run RAxML\n";
    my $binaries_dirname = $user_options->{-e};
    
    my $bootstrap_replicates = $$user_options{"-b"};
    if (($bootstrap_replicates > 1) && ($raxml_option eq "p")) {
        print "ATTENTION: You intended to do $bootstrap_replicates bootstrap replicates. Unfortunately bootstrapping is ";
        print "disabled in the parsimony mode of MLTreeMap. The pipeline will continue without bootstrapping.\n";
        $bootstrap_replicates = 1;
    }
    
    my $dirname = dirname(__FILE__);
    foreach my $f_contig (sort {$a cmp $b} keys %$phy_files) {
        my $reference_tree_file = "$dirname/data/tree_data/$$user_options{\"reference_tree\"}";
        my $phy_file = $$phy_files{$f_contig};
        my $denominator = "";
        if ($f_contig =~ /^([^_]+)_([^_]+).*$/) {
            $denominator = $1;    
        }
        #check if it's a functional hit. If yes, change the reference tree.    
        unless (($denominator eq "p") || ($denominator eq "g") || ($denominator eq "i")) {
            foreach my $cog (sort {$a cmp $b} keys %{$$cog_list{"all_cogs"}}) {
                next unless ($$cog_list{"all_cogs"}{$cog} eq $denominator);
                $reference_tree_file = "$dirname/data/tree_data/$cog"."_tree.txt";
                last;  
            }
        }
        #done.
        
        #save the reference tree file for each denominator
        $$user_options{"reference_tree_file_of_denominator"}{$denominator} = $reference_tree_file;
        
        my @raxml_files = ("$output_directory_var"."RAxML_info.$f_contig","$output_directory_var"."RAxML_labelledTree.$f_contig",
        "$output_directory_var"."RAxML_classificaton.$f_contig");
        
        foreach my $raxml_file (@raxml_files) {
            my $remove_command = "rm $raxml_file";
#   print $remove_command."\n";
            system ($remove_command) if (-e $raxml_file);    
        }
        
        my $model_to_be_used = $$models_to_be_used{$f_contig};
        die "ERROR: no best aa model could be detected for the ML step!\n" unless $model_to_be_used;

        my $dirname = dirname(__FILE__);
        my $raxml_command = "$binaries_dirname/raxmlHPC -m $model_to_be_used ";
        $raxml_command .= "-x 12345 -# $bootstrap_replicates " if ($bootstrap_replicates > 1);
        $raxml_command .= "-s $phy_file -t $reference_tree_file -f $raxml_option -n $f_contig "; 
        $raxml_command .= "-w $output_directory_var > $output_directory_var$f_contig"."_RAxML.txt";
        
        #the following lines are meant for debugging. So they are outcommented by default.
        #my $raxml_command_file = "$output_directory_var"."RAxML_commands.txt";
        #if (-e $raxml_command_file) {
        #    open (OUT, ">> $raxml_command_file") or die "Can't append to $raxml_command_file!\n";
        #} else {
        #    open (OUT, "> $raxml_command_file") or die "Can't create $raxml_command_file!\n";
        #}   
        #print OUT "denominator: $denominator\n$raxml_command\n";
        #close OUT;
        #print "RAxML call: $raxml_command\n";
#        print $raxml_command."\n";



        if ($$user_options{-c} eq "s") {
            my $phy_file_name = "";
            if ($phy_file =~ /.+\/(.+)/) {
                $phy_file_name = $1;    
            }
            my $qsub_id = $small_subroutines->do_sun_grid_computing($user_options,$phy_file_name,$raxml_command);
            $expected_raxml_outfiles{"$output_directory_var"."RAxML_classification.$f_contig"} = 1 if ($raxml_option eq "v");
            $expected_raxml_outfiles{"$output_directory_var"."RAxML_originalLabelledTree.$f_contig"} = 1 if ($raxml_option eq "v");
            $expected_raxml_outfiles{"$output_directory_var"."RAxML_parsimonyTree.$f_contig"} = 1 if ($raxml_option eq "p");
            $sun_grid_jobs{$qsub_id} = 1;
        } else {
            system($raxml_command);
        }
        
    }
    
    if ($$user_options{-c} eq "s") {
        $small_subroutines->halt_mltreemap_until_sungrid_is_done($user_options,\%sun_grid_jobs,\%expected_raxml_outfiles);
    }
    
    foreach my $f_contig (sort {$a cmp $b} keys %$phy_files) {
        my $denominator = "";
        if ($f_contig =~ /^([^_]+)_([^_]+).*$/) {
            $denominator = $1;    
        }       
        my $move_command = "mv $output_directory_var"."RAxML_info.$f_contig ";
        $move_command .= "$output_directory_var"."$f_contig".".RAxML_info.txt";
        system($move_command);
        
        if ($raxml_option eq "v") {
            
            $raxml_outfiles{$denominator}{$f_contig}{"classification"} = "$output_directory_var$f_contig".".RAxML_classification.txt";
            $raxml_outfiles{$denominator}{$f_contig}{"labelled_tree"} = "$output_directory_var$f_contig".".originalRAxML_labelledTree.txt";        
            my $move_command1 = "mv $output_directory_var"."RAxML_classification.$f_contig ";
            $move_command1 .= "$raxml_outfiles{$denominator}{$f_contig}{\"classification\"}";    
            my $move_command2 = "mv $output_directory_var"."RAxML_originalLabelledTree.$f_contig ";
            $move_command2 .= "$raxml_outfiles{$denominator}{$f_contig}{\"labelled_tree\"}";
            my $remove_command = "rm $output_directory_var"."RAxML_labelledTree.$f_contig";
            system($move_command1);
            system($move_command2);
            system($remove_command);
            
        } elsif ($raxml_option eq "p") {
            
            $raxml_outfiles{$denominator}{$f_contig} = "$output_directory_var$f_contig".".RAxML_parsimonyTree.txt";
            my $move_command1 = "mv $output_directory_var"."RAxML_parsimonyTree.$f_contig ";
            $move_command1 .= "$raxml_outfiles{$denominator}{$f_contig}";
            system($move_command1);
            
        } else {
            die "ERROR: the chosen RAxML mode is invalid. This should have been noticed earlier by MLTreeMap. Please notify the authors\n";    
        }
    }
    
    return(\%raxml_outfiles);
}

#####################################################
# parse_RAxML_output
#####################################################

sub parse_RAxML_output {    
    my $package = shift;
    my $user_options = shift;
    my $tree_rerooter = shift;
    my $tree_numbers_translation = shift;
    my $raxml_outfiles = shift;
    my $text_of_analysis_type = shift;
    
    my $raxml_option = $$user_options{"-f"};
    
    print "finishing\n";
    
    my $output_directory_final_RAxML = $$user_options{-o_2};
    my %final_RAxML_output_files;
    
    foreach my $denominator (sort {$a cmp $b} keys %$raxml_outfiles) {
        
        my $description_text = "# $$text_of_analysis_type{$denominator}\n";
        my $reference_tree_file = $$user_options{"reference_tree_file_of_denominator"}{$denominator};
        my $terminal_children_strings_of_reference = $tree_rerooter->read_and_understand_the_reference_tree($reference_tree_file);


        my $content_of_previous_labelled_tree_file = "";
        my $rooted_labelled_trees = "";
        my $insertion_point_node_hash = "";
        
        my %final_assingment_target_strings = ();
        
        foreach my $f_contig (sort {$a cmp $b} keys %{$$raxml_outfiles{$denominator}}) {
            my $denominator = "";
            if ($f_contig =~ /^([^_]+)_([^_]+).*$/) {
                $denominator = $1;    
            }
            my $content_of_labelled_tree_file = "";
            my %assignments = ();
            my $nr_of_assignments = 0;
          
            #############################
            #depending on the RaxML algorithm, the results have to be queried differently
            #############################
            if ($raxml_option eq "v") {
                #i.e. ML analysis.
                my $classification_file = $$raxml_outfiles{$denominator}{$f_contig}{"classification"};
                my $labelled_tree_file = $$raxml_outfiles{$denominator}{$f_contig}{"labelled_tree"};
                
                open (IN, "$labelled_tree_file") or die "Can't open $labelled_tree_file!\n";
                while (<IN>) {
                    chomp $_;
                    $content_of_labelled_tree_file .= $_;    
                }
                if ($content_of_labelled_tree_file ne $content_of_previous_labelled_tree_file) {
                    ($rooted_labelled_trees,$insertion_point_node_hash) = $tree_rerooter->read_understand_and_reroot_the_labelled_tree
                                                                                          ($labelled_tree_file);
                    %final_assingment_target_strings = ();
                }
                my %new_assignments = ();
                my $at_least_one_new_assingment = 0;
                open (IN, "$classification_file") or die "ERROR: Can't open $classification_file!\n";
                while (<IN>) {
                    chomp $_;
                    my ($query, $insertion_point_l, $weight, $dummy) = split / /;
                    my $assignment = "";
                    if ($insertion_point_l =~ /I(\d+)/) {
                        $assignment = $1;
                        $nr_of_assignments++;
                    }
                    $assignments{$assignment} = $weight;
                    unless (exists $final_assingment_target_strings{$assignment}) {
                        $new_assignments{$assignment} = 1;
                        $at_least_one_new_assingment = 1;
                    }
                }
                if ($at_least_one_new_assingment) {
                    #only in case of new assignments (i.e. assignments, which in the whole dataset never before ocurred), we do the rerooting.
                    my $prae_assingment_target_strings = $tree_rerooter->identify_the_correct_terminal_children_of_each_assignment
                                                              ($terminal_children_strings_of_reference, $rooted_labelled_trees, 
                                                               $insertion_point_node_hash, \%new_assignments);
                    #add the new assignment strings to the hash, which contains them all.
                    foreach my $assignment (sort {$a cmp $b} keys %$prae_assingment_target_strings) {
                        my $assingment_target_string = $$prae_assingment_target_strings{$assignment};
                        $final_assingment_target_strings{$assignment} = $assingment_target_string;
                    }
                }
        
            } elsif ($raxml_option eq "p") {
                #i.e. MP analysis           
                my $mp_tree_file = $$raxml_outfiles{$denominator}{$f_contig};
                my $assignment = "mp_root";
                $assignments{$assignment} = 1;
                $nr_of_assignments = 1;        
                my $prae_assingment_target_strings = $tree_rerooter->get_correct_mp_assignment($terminal_children_strings_of_reference, 
                                                                                       $mp_tree_file,\%assignments);
                my $assingment_target_string = $$prae_assingment_target_strings{$assignment};
                $final_assingment_target_strings{$assignment} = $assingment_target_string;
            }
            #############################
            #done. Now print the information
            #############################
            
            
            my $final_RAxML_filename = "$output_directory_final_RAxML"."$f_contig"."_RAxML_parsed.txt";
            $final_RAxML_output_files{$denominator}{$final_RAxML_filename} = 1;
            
            
#            print "generating "."$final_RAxML_filename"."\n";
            
            open (OUT, "> $final_RAxML_filename") or die "ERROR: Can't create $final_RAxML_filename!\n";
            
            print OUT "$description_text\n";
            
            foreach my $assignment (sort {$a cmp $b} keys %assignments) {
            
                my $assingment_target_string = $final_assingment_target_strings{$assignment};
                my $weight = $assignments{$assignment};
                my $relative_weight = int((($weight / $nr_of_assignments) * 100)+0.5);
            
                my @assignment_terminal_targets = split / /, $assingment_target_string;
                my $nr_of_terminal_targets = @assignment_terminal_targets + 0;
                
                print OUT "Placement weight $relative_weight\%: Assignment of query to ";
                print OUT "the lowest common ancestor of " unless ($nr_of_terminal_targets == 1);
                my $count = "";
                for ($count = 1; $count <= $nr_of_terminal_targets; $count++) {          
                    my $assignment_terminal_target = $assignment_terminal_targets[$count-1];
                    my $is_last_element = 0;
                    $is_last_element = 1 if ($count == $nr_of_terminal_targets -1);
                    my $name_of_terminal_target = "";
                    $name_of_terminal_target = $$tree_numbers_translation{$denominator}{$assignment_terminal_target};
                    unless (defined $name_of_terminal_target) {
                        die "ERROR: $assignment_terminal_target could not be located in the tree with the denominator $denominator!\n";
                    }
                    print OUT "$name_of_terminal_target ($assignment_terminal_target)";
                    print OUT ", "    if ($count < $nr_of_terminal_targets - 1);
                    print OUT " and " if ($count == $nr_of_terminal_targets - 1);
                    print OUT "."     if ($count == $nr_of_terminal_targets);             
                }
            }
            close OUT;
            $content_of_previous_labelled_tree_file = $content_of_labelled_tree_file;
        }
    }
    return(\%final_RAxML_output_files);  
}

#####################################################
# concatenate_RAxML_output_files
#####################################################

sub concatenate_RAxML_output_files {  
    my $package = shift;
    my $user_options = shift;
    my $final_RAxML_output_files = shift;
    my $text_of_analysis_type = shift;
    
    my $output_directory_final = $$user_options{-o_3};
    
    foreach my $denominator (sort {$a cmp $b} keys %$final_RAxML_output_files) {
        my $nr_of_files = 0;
        my %assignments = ();
        
        my $description_text = "# $$text_of_analysis_type{$denominator}\n";
        my $final_output_file_name = "$output_directory_final"."$denominator"."_concatenated_RAxML_outputs.txt";
        
        foreach my $final_RAxML_output_file (sort {$a cmp $b} keys %{$$final_RAxML_output_files{$denominator}}) {
            $nr_of_files++;
            open (IN, "$final_RAxML_output_file") or die "ERROR: Can't open $final_RAxML_output_file!\n";
            while (<IN>) {
                chomp $_;
                if (/Placement weight (\d+)%: (.+)\Z/) {
                    my $weight = $1;
                    my $assignment = $2;
                    if (exists $assignments{$assignment}) {
                        $assignments{$assignment} += $weight;
                    } else {
                        $assignments{$assignment} = $weight;    
                    }
                } else {
                    next;    
                }
            }
            close IN;
        }
        my %assignments_with_relative_weights = ();
        foreach my $assignment (sort {$b cmp $a} keys %assignments) {
            my $weigth = $assignments{$assignment};
            my $relative_weight = (int((($weigth / $nr_of_files) * 10000)+0.5))/10000;
            $assignments_with_relative_weights{$relative_weight}{$assignment} = 1;
        }
        
        open (OUT, "> $final_output_file_name") or die "Can't create $final_output_file_name!\n";
        print OUT "$description_text\n";
        my $sum_of_relative_weights = 0;
        foreach my $relative_weight (sort {$b <=> $a} keys %assignments_with_relative_weights) {
            foreach my $assignment (sort {$b cmp $a} keys %{$assignments_with_relative_weights{$relative_weight}}) {
                $sum_of_relative_weights += $relative_weight;
#        print "Placement weight $relative_weight\%: $assignment\n";
                print OUT "Placement weight $relative_weight\%: $assignment\n";
            }
        }
        close OUT;
        print "$denominator"."_ sum of placement weights (should be 100): $sum_of_relative_weights\n";
    }
}

#####################################################
#####################################################
1;

#####################################################
#####################################################

# package small_subroutines

#####################################################
#####################################################

package small_subroutines;

#####################################################
# new
#####################################################

sub new {
    my $package = shift;
    my $reference = {};
    bless ($reference, $package);
    return ($reference);    
}

#####################################################
# calculate_overlap
#####################################################

sub calculate_overlap {   
    my $package = shift;  
    my $info = shift;
    
    my $base_start = $$info{"base"}{"start"};
    my $base_end = $$info{"base"}{"end"};
    my $check_start = $$info{"check"}{"start"};
    my $check_end = $$info{"check"}{"end"};
    
    my $overlap = 0;
    
    #now check, if they overlap:
    if (($base_start <= $check_start) && ($check_start <= $base_end) && ($base_end <= $check_end)) {
        # Base     ----------       
        # Check        ----------   
        $overlap = $base_end - $check_start;
    }
    elsif (($base_start <= $check_start) && ($check_end <= $base_end)) {
        # Base     ----------       
        # Check        ---          
        $overlap = $check_end - $check_start;
    }
    elsif (($check_start <= $base_start) && ($base_start <= $check_end) && ($check_end <= $base_end)) {
        # Base           ------     
        # Check   ----------        
        $overlap = $check_end - $base_start;
    }
    elsif (($check_start <= $base_start) && ($base_end <= $check_end)) {
        # Base           ------     
        # Check   ----------------   
        $overlap = $base_end - $base_start;            
    }
    
    return ($overlap);
}

#####################################################
# undef_hashes
#####################################################

sub undef_hashes {
    my $package = shift;

    foreach my $hash (shift) {
        undef %$hash;
    }
}

#####################################################
# do_sun_grid_computing
#####################################################

sub do_sun_grid_computing {
    my $package = shift;
    my $user_options = shift;
    my $input_file_name = shift;
    my $command = shift;
    
    my $output_directory_var = $$user_options{-o_1};
    my $minimum_memory_consumtion = $$user_options{-m};
    
    my $memory_sentence = "";
    $memory_sentence = "-l 'mem_free=$minimum_memory_consumtion"."G,h_vmem=$minimum_memory_consumtion"."G'" if $minimum_memory_consumtion;
        
    if ($input_file_name =~ /.*\/(.+)/) {
        $input_file_name = $1; #i.e. if a path to the file was given, extract the filename.    
    }   
    my $bash_script_file = "$output_directory_var"."$input_file_name".".sh";
    my $bash_script_file1 = "$output_directory_var"."$input_file_name"."_1.sh";
    my $report_file = "$output_directory_var"."$input_file_name"."_qsub_report.txt";
    my $blast_tmp_file_name = "/tmp/$input_file_name";
    
    my $results_raw_files_raw = "";
    my $is_blast_step = 0;
    $command =~ s/#/\\#/g; #the "#" symbol would otherwise be interpreted as the start of a comment by the bash...
    
    #ok if we are blasting somethings are different.
    my $rRNA_denominator = "";
    if ($command =~ /sub_binaries\/blastall.+> (.+)/) {
        $results_raw_files_raw = $1;
        if ($command =~ /rRNA_BLAST_results_raw.txt\Z/) {
            $bash_script_file = "$output_directory_var"."rRNA_$input_file_name".".sh";
            $bash_script_file1 = "$output_directory_var"."rRNA_$input_file_name"."_1.sh";
            $report_file = "$output_directory_var"."rRNA_$input_file_name"."_qsub_report.txt";
            $blast_tmp_file_name = "/tmp/rRNA_$input_file_name";
        }
        $is_blast_step = 1;
        if ($command =~ /\A(.+)>/) {
            $command = $1;
            $command .= "> $blast_tmp_file_name";  
        } else {
            die "Error, parsing problem with $command!\n";    
        }
    }
    #blast changes done
        
    open (OUT, "> $bash_script_file") or die "ERROR: Can't create $bash_script_file\n";    
    print OUT "#!/bin/bash\n";
    print OUT "qsub $memory_sentence -cwd -o $output_directory_var -e $output_directory_var $bash_script_file1 >$report_file\n";    
    close OUT;
    
    open (OUT, "> $bash_script_file1") or die "ERROR: Can't create $bash_script_file1\n";    
    print OUT "#!/bin/bash\n";
    print OUT "cd ";
    close OUT;
    system ("pwd >> $bash_script_file1"); #append the current dir to the script file
    open (OUT, ">> $bash_script_file1") or die "ERROR: Can't append to $bash_script_file1\n";
    print OUT "\n";
    print OUT "$command\n";
    print OUT "mv $blast_tmp_file_name $results_raw_files_raw\n" if ($is_blast_step);
    close OUT;
    my $command1 = "chmod +x $bash_script_file";
    my $command2 = "$bash_script_file";
    system ($command1);
    system ($command2);
    
    open (IN, "$report_file") or die "ERROR: $report_file can't be opened!\n";
    my $qsub_id = "";
    while (<IN>) {
        chomp $_;
        if (/\AYour job (\d+) /) {
            $qsub_id = $1;        
        } else {
            die "ERROR: your job has not been properly submitted to qsub!\n$_\n";    
        }
    }
    return $qsub_id;
}

#####################################################
# halt_mltreemap_until_sungrid_is_done
#####################################################

sub halt_mltreemap_until_sungrid_is_done {
    my $package = shift;
    my $user_options = shift;
    my $sun_grid_jobs = shift;
    my $result_files = shift;
    
    my $output_directory_var = $$user_options{-o_1};
    
    my $command = "";
    my $count = 0;
    my $count1 = 0;
    while (1) {
        #leave the loop only upon success or fail (or after a very long timeout).
        $count++;
        $count1++;
        my $sleep_timer1 = 60; #default 60
        sleep($sleep_timer1);
        my $leave_the_loop = 1;
        foreach my $result_file (sort {$a cmp $b} keys %$result_files) {
            $leave_the_loop = 0 unless (-e $result_file);
        }
        last if ($leave_the_loop);
        
        if ($count == 60) { #default 60
            $count = 0;
            #check if a job has crashed.
            my $qstat_file = "$output_directory_var"."qstat.txt";
            my $qstat_command = "qstat > $qstat_file";
            system($qstat_command);
            open (IN, "$qstat_file") or die "ERROR: Can't open $qstat_file!\n";
            while (<IN>) {
                chomp $_;
                if (/ (\d+) \S+ \S+ \S+      (\S+)/) {
                    my $job_id = $1;
                    my $state = $2;
                    if ((exists $$sun_grid_jobs{$job_id}) && ($state eq "Eqw")) {
                        die "ERROR: your job crashed on the cluster!\n";    
                    }  
                }    
            }
        }
        die "ERROR: your job takes suspiciously long on the cluster (>10 days). Please check.\n" if ($count1 == 14400); #default 14400 i.e. 10d when timer is 60 
    }
    print "Sun Grid jobs are done\n";
}

#####################################################
#####################################################
1;

#####################################################
#####################################################

# package reroot_RAxML_trees

#####################################################
#####################################################

package reroot_RAxML_trees;

#####################################################
# new
#####################################################

sub new {
    my $package = shift;
    my $reference = {};
    bless ($reference, $package);
    return ($reference);    
}

#####################################################
# read_and_understand_the_reference_tree
#####################################################

sub read_and_understand_the_reference_tree {
    my $package = shift;
    my $reference_tree_file = shift;
    
    my $reference_tree_elements                             = &read_the_reference_tree($reference_tree_file);
    my $reference_tree_info                                 = &create_tree_info_hash;
    &get_node_subtrees($reference_tree_elements,$reference_tree_info);

    &assign_parents_and_children($reference_tree_info);
    
    #build the terminal children strings of the reference tree
    my $terminal_children_strings_of_reference = &build_terminal_children_strings_of_reference_nodes($reference_tree_info);
    
    return($terminal_children_strings_of_reference);
}

#####################################################
# read_understand_and_reroot_the_labelled_tree
#####################################################

sub read_understand_and_reroot_the_labelled_tree {
    my $package = shift;
    my $labelled_tree_file = shift;

    
    my ($labelled_tree_elements,$insertion_point_node_hash) = &read_the_raxml_out_tree($labelled_tree_file);
    my $labelled_tree_info                                  = &create_tree_info_hash;

#print Dumper($labelled_tree_elements);
    &get_node_subtrees($labelled_tree_elements,$labelled_tree_info);
    &assign_parents_and_children($labelled_tree_info);
    &build_tree_info_quartets($labelled_tree_info);
    
    #root the labelled tree at every possible position
    my $rooted_labelled_trees                               = &build_newly_rooted_trees($labelled_tree_info);
    
    return($rooted_labelled_trees,$insertion_point_node_hash);
}

#####################################################
# identify_the_correct_terminal_children_of_each_assignment
#####################################################

sub identify_the_correct_terminal_children_of_each_assignment {
    my $package = shift;
    my $terminal_children_strings_of_reference = shift;
    my $rooted_labelled_trees = shift;
    my $insertion_point_node_hash = shift;
    my $assignments = shift;
    
    #build the terminal children string of the assignments
    my $terminal_children_strings_of_assignments            = &build_terminal_children_strings_of_assignments($rooted_labelled_trees,
                                                              $insertion_point_node_hash,$assignments);
                                                              
    #identify the correct terminal children of each assignment
    my $real_terminal_children_strings_of_assignments       = &compare_terminal_children_strings(
                                                              $terminal_children_strings_of_assignments,$terminal_children_strings_of_reference);
                                                              
    return($real_terminal_children_strings_of_assignments);
}

#####################################################
# get_correct_mp_assignment
#####################################################

sub get_correct_mp_assignment {  
    my $package = shift;
    my $terminal_children_strings_of_reference = shift;
    my $mp_tree_file = shift;
    my $assignments = shift;
          
    #read and the raxml mp tree return the two possible terminal children strings
    my $potential_terminal_children_strings                 = &read_the_raxml_mp_out_tree($mp_tree_file,$assignments);   
    
    #identify the correct terminal children of the MP tree assignment
    my $real_terminal_children_strings_of_assignments       = &compare_terminal_children_strings(
                                                              $potential_terminal_children_strings,$terminal_children_strings_of_reference);
    #done
    
    return($real_terminal_children_strings_of_assignments);
}

#####################################################
# read_the_reference_tree
#####################################################

sub read_the_reference_tree {   
    my $reference_tree_file = shift;
    
    my @reference_tree_elements = "";
    
    open (IN, "$reference_tree_file") or die "Can't open $reference_tree_file\n"; 
    my $tree_string = "";
    while (<IN>) {
        chomp $_;
        $tree_string .= $_;
    }
        
    #to avoid troubles with the pattern matching replace some characters
    $tree_string =~ s/\(/L/g;
    $tree_string =~ s/\)/R/g;
    
    #remove the branchlengths
    $tree_string =~ s/:\d+\.\d+//g;
    
    #assign internal node markers (negative numbers < -1)
    my $count = -2;
    while ($tree_string =~ s/R/Q$count/) {
        $count--;
    }
    $tree_string =~ s/Q-\d+;/Q;/; #remove the last number
        
    #revert the changed characters
    $tree_string =~ s/L/\(/g;
    $tree_string =~ s/Q/\)/g;
            
    #split the tree string
    
    my $reference_tree_elements = &split_tree_string($tree_string);
    
    return ($reference_tree_elements);
}

#####################################################
# read_the_raxml_out_tree
#####################################################

sub read_the_raxml_out_tree {       
    my $labelled_tree_file = shift;
    
    my @tree_elements = "";
    my %insertion_point_node_hash = ();
    
    open (IN, "$labelled_tree_file") or die "Can't open $labelled_tree_file\n"; 
    my $tree_string = "";
    while (<IN>) {
        chomp $_;
        $tree_string .= $_;
    }
    
    #arbitrarily root the tree
    my @tree_symbols_raw_1 = split //, $tree_string; 
    my $bracket_diff = 0;
    my $tree_string_neu = "(";
    my $comma_count = 0;
    foreach my $tree_symbol_raw_1 (@tree_symbols_raw_1) {
        if ($comma_count < 2) {
            $bracket_diff++ if ($tree_symbol_raw_1 eq "(");
            $bracket_diff-- if ($tree_symbol_raw_1 eq ")");
            $comma_count++ if (($tree_symbol_raw_1 eq ",") && ($bracket_diff == 1));
            $tree_string_neu .= "):1.0[I666999666]" if ($comma_count == 2);
        }
        $tree_string_neu .= $tree_symbol_raw_1;
    }
    $tree_string = $tree_string_neu;
        
    #to avoid troubles with the pattern matching replace some characters
    $tree_string =~ s/\(/L/g;
    $tree_string =~ s/\)/R/g;
    $tree_string =~ s/\[/Q/g;
    
    #remove the branchlengths (they are set to 1 anyway)
    $tree_string =~ s/:1.0//g;
    
    #shorten the tree around the terminal leaves.
    while ($tree_string =~ /((\D(\d+))QI(\d+)])/) {
        my $to_be_replaced = $1;
        my $replacement = $2;
        my $terminal_leaf = $3;
        my $insertion_point = $4;
        if ($terminal_leaf <= 0) {
            die "ERROR: your tree has terminal leaves with numbers <= 0. Please change them to positive values!\n";    
        }
        $insertion_point_node_hash{$insertion_point} = $terminal_leaf; 
        $tree_string =~ s/$to_be_replaced/$replacement/;
    }
    
    #replace the RAxML insertion points with internal node markers (negative numbers < -1)
    my $count = -2;
    while ($tree_string =~ s/QI(\d+)]/$count/) {
        $insertion_point_node_hash{$1} = $count;
        $count--;
    }
        
    #revert the changed characters
    $tree_string =~ s/L/\(/g;
    $tree_string =~ s/R/\)/g;
    $tree_string =~ s/Q/\[/g;
        
    #split the tree string
    my $tree_elements = &split_tree_string($tree_string);
      
    return ($tree_elements,\%insertion_point_node_hash);
}

#####################################################
# read_the_raxml_mp_out_tree
#####################################################

sub read_the_raxml_mp_out_tree {    
    my $mp_tree_file = shift;
    my $assignments = shift;
    
    #In this case, the query is allways at the root: (A,B,query);
    #           query
    #             |
    #           -----
    #           |   |
    #           A   B
    #Since the real root is either at A or B, it is also either
    #A or B, which contains all real terminal children of the assingment.
    #This subroutine just gets the terminal children of A and B.
 
    my %potential_terminal_children_strings = ();
    
    my $assignment = "";
    foreach my $assig (sort {$a cmp $b} keys %$assignments) {
        $assignment = $assig;
        last; #i.e. there is only one.
    }
    open (IN, "$mp_tree_file") or die "Can't open $mp_tree_file\n"; 
    my $tree_string = "";
    while (<IN>) {
        chomp $_;
        $tree_string .= $_;
    }
            
    #to avoid troubles with the pattern matching replace some characters
    $tree_string =~ s/\(/L/g;
    $tree_string =~ s/\)/R/g;
    
    #remove the query/root triforcation and die if the query is not at the root.
    unless ($tree_string =~ s/,queryR;\Z/R;/) {
        die "ERROR: the query is not at the root of $mp_tree_file!\n";
    }
    
    #remove the branchlengths
    $tree_string =~ s/:\d+\.\d+//g;
    
    #assign internal node markers (negative numbers < -1)
    my $count = -2;
    while ($tree_string =~ s/R/Q$count/) {
        $count--;
    }
    $tree_string =~ s/Q-\d+;/Q;/; #remove the last number
        
    #revert the changed characters
    $tree_string =~ s/L/\(/g;
    $tree_string =~ s/Q/\)/g;
    
    #split the tree into the two parts, separated by root
    my @tree_symbols = split //, $tree_string; 
    my $bracket_diff = 0;
    my $comma_count = 0;
    my @substrings = ("",",");
    
    foreach my $tree_symbol (@tree_symbols) {
        if ($comma_count < 1) {
            $bracket_diff++ if ($tree_symbol eq "(");
            $bracket_diff-- if ($tree_symbol eq ")");
            $comma_count++ if (($tree_symbol eq ",") && ($bracket_diff == 1));
            $substrings[0] .= $tree_symbol;
        } else {
            $substrings[1] .= $tree_symbol;
        }
    }
    
    foreach my $substring (@substrings) {
        my %terminal_children = ();
        while ($substring =~ /(\D)(\d+)/g) {
            next if ($1 eq "-");
            $terminal_children{$2} = 1;
        }
        my $potential_terminal_children_string = "";
        foreach my $potential_terminal_child (sort {$a <=> $b} keys %terminal_children) {
            $potential_terminal_children_string .= "$potential_terminal_child ";    
        }
        $potential_terminal_children_strings{$assignment}{$potential_terminal_children_string} = 1;
    }            
    return (\%potential_terminal_children_strings);
}

#####################################################
# split_tree_string
#####################################################

sub split_tree_string {
    my $tree_string = shift;
    my @tree_elements = "";
    
    my @tree_symbols_raw = split //, $tree_string;
    my $count = -1;
    my $previous_symbol = "";
    foreach my $tree_symbol_raw (@tree_symbols_raw) {
        if (($tree_symbol_raw =~ /\d/) && (($previous_symbol =~ /\d/) || ($previous_symbol eq "-"))) {
            $tree_elements[$count] .= $tree_symbol_raw;
        } else {
            $count++;
            $tree_elements[$count] = $tree_symbol_raw;
        }
        $previous_symbol = $tree_symbol_raw;
    }  
    return (\@tree_elements);   
}


#####################################################
# create_tree_info_hash
#####################################################

sub create_tree_info_hash {
    my %tree_info = ();
    $tree_info{"parent_of_node"} = undef;
    $tree_info{"children_of_node"} = undef;
    $tree_info{"subtree_of_node"} = undef;
    $tree_info{"quartets"} = undef;
    
    return (\%tree_info);
}

#####################################################
# get_node_subtrees
#####################################################

sub get_node_subtrees {
    my $tree_elements = shift;
    my $tree_info = shift;
    
    my $bracket_l_count = 0;
    my $bracket_r_count = 0;
    my %parents_of_node = ();
    
    my $tree_element_nr = -1;
    foreach my $tree_element (@$tree_elements) {
        $tree_element_nr++;
        if ($tree_element eq "(") {
            my $bracket_l_count = 1;
            my $bracket_r_count = 0;
            
            my $tree_sub_element_nr = $tree_element_nr;
            my $subtree_string = "(";
            while (1) {
                $tree_sub_element_nr++;
                my $tree_sub_element = $$tree_elements[$tree_sub_element_nr];
                $bracket_l_count++ if ($tree_sub_element eq "(");
                $bracket_r_count++ if ($tree_sub_element eq ")");
        
                if ($bracket_l_count == $bracket_r_count) {
                    my $nodename = $$tree_elements[$tree_sub_element_nr + 1];
                    $nodename = -1 if ($nodename eq ";");
                    $subtree_string .= ")$nodename";
                    $$tree_info{"subtree_of_node"}{$nodename} = $subtree_string;
                    last;
                } else {
                   $subtree_string .= $tree_sub_element;
                }
            }
        }
    }    
    foreach my $tree_element (@$tree_elements) {
        next unless ($tree_element =~ /\d+/);
        next if (exists $$tree_info{"subtree_of_node"}{$tree_element});
        $$tree_info{"subtree_of_node"}{$tree_element} = $tree_element;        
    }
}

#####################################################
# assign_parents_and_children
#####################################################
my $visit_count = 0;
sub assign_parents_and_children {
    my $tree_info = shift;

    
    foreach my $node (sort {$a <=> $b} keys %{$$tree_info{"subtree_of_node"}}) {
        
        #find the parent
        next if ($node == -1); #we are root. no parent to be found.
        my $subtree = $$tree_info{"subtree_of_node"}{$node};
                
        my $parent = undef;
        foreach my $potential_parent (sort {$a <=> $b} keys %{$$tree_info{"subtree_of_node"}}) {
            next if ($node == $potential_parent);
            my $potential_parent_subtree = $$tree_info{"subtree_of_node"}{$potential_parent};
            
            $subtree =~ s/\(/L/g;
            $subtree =~ s/\)/\#/g;
            $potential_parent_subtree =~ s/\(/L/g;
            $potential_parent_subtree =~ s/\)/#/g;
                       
            if (
                ($potential_parent_subtree =~ /\AL$subtree,.+#$potential_parent\Z/) || 
                ($potential_parent_subtree =~ /\AL.+,$subtree#$potential_parent\Z/)) {
                $parent = $potential_parent;
                last;    
            }
        }
        $$tree_info{"parent_of_node"}{$node} = $parent;
        $$tree_info{"children_of_node"}{$parent}{$node} = 1;      
    }   
}

#####################################################
# build_tree_info_quartets
#####################################################

sub build_tree_info_quartets {
    my $tree_info = shift;
    
    foreach my $node (sort {$a <=> $b} keys %{$$tree_info{"parent_of_node"}}) {           
        
        #a node has maximum 3 connections: 2 children and the parent. 
        #if the parent is root, the connection we want is the sister.
        
        #these quartets are universal, no matter where the tree is rooted.
        
        my $parent = $$tree_info{"parent_of_node"}{$node};
        if ($parent == -1) {
            foreach my $roots_child (sort {$a <=> $b} keys %{$$tree_info{"children_of_node"}{"-1"}}) {
                next if ($roots_child == $node);
                $parent = $roots_child;
            } 
        }
        $$tree_info{"quartets"}{$node}{$parent} = 1;
        if (exists $$tree_info{"children_of_node"}{$node}) {
            foreach my $child (sort {$a <=> $b} keys %{$$tree_info{"children_of_node"}{$node}}) {
                $$tree_info{"quartets"}{$node}{$child} = 1;
            }            
        }
    }
}

#####################################################
# build_newly_rooted_trees
#####################################################

sub build_newly_rooted_trees {
    my $tree_info = shift;
    my @rooted_trees = ();
    
    my $tree_number = 0;
    my %list_of_allready_used_attachments = ();
    
    foreach my $node (sort {$a <=> $b} keys %{$$tree_info{"quartets"}}) {
        next if (exists $list_of_allready_used_attachments{$node});
        foreach my $attachment (sort {$a <=> $b} keys %{$$tree_info{"quartets"}{$node}}) {
            #ok, now we would have a potential new root. lets build the tree...
            $list_of_allready_used_attachments{$attachment} = 1;
            my $tree_string = "";
            my $root = -1;
            my %node_infos = ();
            $node_infos{"previous_node"} = "";
            $node_infos{"node"} = ";";
            $node_infos{"open_attachments"}{$node} = 1;
            $node_infos{"open_attachments"}{$attachment} = 1;
            my $new_tree = &recursive_tree_builder($tree_info,\%node_infos,$tree_string);
            $rooted_trees[$tree_number] = $new_tree;
            $tree_number++;
        }
    }
    return(\@rooted_trees);
}

#####################################################
# recursive_tree_builder
#####################################################

sub recursive_tree_builder {
    my $tree_info = shift;
    my $node_infos = shift;
    my $tree_string = shift;
        
    #we assemlbe the tree by building modules as shown below:
    #(attachment1,attachment2)$node
    
    my $node = $$node_infos{"node"};
    
    my $count = 0;
    foreach my $attachment (sort {$a <=> $b} keys %{$$node_infos{"open_attachments"}}) { 
        $count++;
        $tree_string .= "(" if ($count == 1);
        
        my %node_infos2 = ();
        $node_infos2{"previous_node"} = $node;
        $node_infos2{"node"} = $attachment;
        
        my $count2 = 0;
        foreach my $attachment_of_used_attachment (sort {$a <=> $b} keys %{$$tree_info{"quartets"}{$attachment}}) { 
            next if (exists $$node_infos{"open_attachments"}{$attachment_of_used_attachment});
            next if $attachment_of_used_attachment eq $node;
            $count2++;
            $node_infos2{"open_attachments"}{$attachment_of_used_attachment} = 1;
        }
        if ($count2) {
            $tree_string = &recursive_tree_builder($tree_info,\%node_infos2,$tree_string);
        } else {     
            $tree_string .= "$attachment"; # i.e. this was a terminal leaf.
        }
        $tree_string .= "," if ($count == 1);
        $tree_string .= ")$node" if ($count ==2); 
    }
    return ($tree_string);
}

#####################################################
# build_terminal_children_strings_of_assignments
#####################################################

sub build_terminal_children_strings_of_assignments {
    my $rooted_trees = shift;
    my $insertion_point_node_hash = shift;
    my $assignments = shift;
     
    my %terminal_children_strings_of_assignments = ();
    
    foreach my $assignment (sort {$a <=> $b} keys %$assignments) {

        my $internal_node_of_assignment = $$insertion_point_node_hash{$assignment};
        
        foreach my $rooted_tree (@$rooted_trees) {
            my $rooted_tree_elements = split_tree_string($rooted_tree);
            my $rooted_tree_info = &create_tree_info_hash;
            &get_node_subtrees($rooted_tree_elements,$rooted_tree_info);
            my $assignment_subtree = $$rooted_tree_info{"subtree_of_node"}{$internal_node_of_assignment};        
            my %terminal_children = ();
            if ($assignment_subtree =~ /\A(\d+)\Z/) {               
                #i.e. it's a terminal leaf. Thus the subtree is only a number.
                $terminal_children{$1} = 1;
            } else {
                while ($assignment_subtree =~ /(\D)(\d+)/g) {
                    next if ($1 eq "-");
                    $terminal_children{$2} = 1;
                }
            }
            my $terminal_children_string_of_assingment = "";
            foreach my $terminal_child_of_assingment (sort {$a <=> $b} keys %terminal_children) {
                $terminal_children_string_of_assingment .= "$terminal_child_of_assingment ";    
            }
            $terminal_children_strings_of_assignments{$assignment}{$terminal_children_string_of_assingment} = 1;
        }
    }
    return(\%terminal_children_strings_of_assignments);
}

#####################################################
# build_terminal_children_strings_of_reference_nodes
#####################################################

sub build_terminal_children_strings_of_reference_nodes {
    my $reference_tree_info = shift;
    
    my %terminal_children_strings_of_reference = ();
        
    foreach my $node (sort {$a cmp $b} keys %{$$reference_tree_info{"subtree_of_node"}}) {
        my $reference_subtree = $$reference_tree_info{"subtree_of_node"}{$node};
        my %terminal_children = ();
        if ($reference_subtree =~ /\A(\d+)\Z/) {               
                #i.e. it's a terminal leaf. Thus the subtree is only a number.
                $terminal_children{$1} = 1;
        } else {
            while ($reference_subtree =~ /(.)(\d+)/g) {
                next if ($1 eq "-");
                $terminal_children{$2} = 1;
            }
        }
        my $terminal_children_string_of_reference = "";
        foreach my $terminal_child_of_reference (sort {$a <=> $b} keys %terminal_children) {
            $terminal_children_string_of_reference .= "$terminal_child_of_reference ";    
        }
        $terminal_children_strings_of_reference{$terminal_children_string_of_reference} = 1;
    }
    return(\%terminal_children_strings_of_reference);
}

#####################################################
# compare_terminal_children_strings
#####################################################

sub compare_terminal_children_strings {
    my $terminal_children_strings_of_assignments = shift;
    my $terminal_children_strings_of_reference = shift;
    
    my %real_terminal_children_strings_of_assignments = ();
    my $there_was_a_hit = 0;
    foreach my $assignment (sort {$a cmp $b} keys %$terminal_children_strings_of_assignments) {
        my $real_terminal_children_string = "";
        foreach my $terminal_children_string_of_assignment (sort {$a cmp $b} keys %{$$terminal_children_strings_of_assignments{$assignment}}) {
            if (exists $$terminal_children_strings_of_reference{$terminal_children_string_of_assignment}) {
                $real_terminal_children_string = $terminal_children_string_of_assignment;
                $real_terminal_children_strings_of_assignments{$assignment} = $real_terminal_children_string;
                $there_was_a_hit = 1;
                last;
            }
        }
        die "ERROR the RAxML output tree could not be rooted correctly!!!\n" unless ($real_terminal_children_string || ($assignment eq "mp_root"));
    }
    #note: the following error is only important for the MP analysis. The ML error is generated above.
    die "ERROR the RAxML output tree could not be rooted correctly!!!\n" unless ($there_was_a_hit);
    return(\%real_terminal_children_strings_of_assignments);
}

}

#####################################################
#####################################################
1;


package  BFile;

use strict;
use warnings;
use Carp;
#use Data::Dumper;

{
   my $count = 0;   
   my %files;
   my $dataAvail = 0;
   my $maxSize = 100000000;
   my $size = 0;



sub getContent {
   return \%files;
}



# Perl trim function to remove whitespace from the start and end of the string
sub trim($)
{
    my $string = shift;
    $string =~ s/^\s+//;
    $string =~ s/\s+$//;
    return $string;
}



sub ReadFile {

   my ($self, $filename) = @_;

   if( $dataAvail == 0)  {
       _readContent($self->{'_name'});
       $dataAvail = 1; 
   } 

   if( $self->{'_permission'} eq 'w' ) {
      croak( "Cannot read from a write only file " );
   }

   if( ! defined $filename) {
      croak( "File name to read is missing" );
   }

   if( !defined $files{$filename}) {
      return undef;
   }

   return $files{$filename}{'data'};
}



sub addFile {
   my ($self, $name, $data ) = @_;
   if( $self->{'_permission'} eq 'r' ) {
      croak( "Cannot add content to a read only file" );
   }

   if( ! defined $name) {
      croak( "File name missing" );
   }

   if( ! defined $data) {
      croak( "File datamissing" );
   }
   
   $files{$name} = {'data' => $data}; 
   $size +=  length($name);
   $size +=  length($data);
   $size += 2;
}


sub addFolder{
   my ($self, $foldername) = @_;
   opendir( DIR, $foldername) or croak("cannot open the folder $foldername");
   my $count = 0;
   while( my $file = readdir(DIR) ) {
       my $data = getData($foldername ."/". $file);
       addFile($self, $file, $data);
       if( $size > $maxSize ) {
           swapToDisk($self);
       }
       $count++;
   }
}

sub _readContent {
   my $filename = shift; 
   my %openers = ( "<file>"=> 0, "<name>" => 0, "<data>" => 0);

   my %closures=( "</file>"=>0, "</name>"=>0, "</data>"=>0);
   my %maps=( "</file>"=>"<file>", "</name>"=>"<name>", "</data>"=><"data"> );

   open  FILEHANDLE, "<", $filename  || croak('Cannot open file to write');

   my $nameStr="";
   my $dataStr="";
   my $lastop;
   while( <FILEHANDLE> ) {
     my $line =  trim($_);
     if ( ($line eq  "<file>" ) or ($line eq "<name>") or ($line eq  "<data>") or 
          ($line eq  "<\/file>" ) or ($line eq "</name>") or ($line eq  "</data>") ) {

         if ( ($line eq  "<file>" ) or ($line eq "<name>") or ($line eq  "<data>") ) {
          $openers{$line}++;
          $lastop = $line;
         }
         else {
           $openers{$maps{$line}}--;
           $lastop = "";
           if( $line eq "</file>") {
               $nameStr = trim($nameStr);
               $files{$nameStr} = {'data' => $dataStr}; 
               $nameStr = "";
               $dataStr = "";
           }
         }
     }
     elsif( $lastop eq "<name>" ) {
         $nameStr .= $_;
     }
     elsif( $lastop eq "<data>" ) {
         $dataStr .= $_;
     }

   }
   #print Dumper(\%files);


   close(FILEHANDLE)
}

sub Close{
   my $self = shift;
   open my $filehandle, '>', $self->{_name}  || croak('Cannot open file to write');
   for my $file (keys %files)  {
       print $filehandle "<file>\n";
       print $filehandle "  <name>\n";
       print $filehandle $file .  "\n";
       print $filehandle "  </name>\n";
       print $filehandle "  <data>\n";
       print $filehandle $files{$file}{'data'} ;
       print $filehandle "  </data>\n";
       print $filehandle "</file>\n";
   }
   close( $filehandle);
   %files = ();
   undef %files;

}

sub swapToDisk{
   my $self = shift;
   open my $filehandle, '>>', $self->{_name}  || croak('Cannot open file to write');
   for my $file (keys %files)  {
       print $filehandle "<file>\n";
       print $filehandle "  <name>\n";
       print $filehandle $file .  "\n";
       print $filehandle "  </name>\n";
       print $filehandle "  <data>\n";
       print $filehandle $files{$file}{'data'};
       print $filehandle "  </data>\n";
       print $filehandle "</file>\n";
   }
   close( $filehandle);

   %files = ();
   $size = 0;

   undef %files;

}


sub validPermission {
   my $permString = shift;
   my ($result)  = ( $permString =~ /(r|w)$/) ;
   if($result eq "r" or $result eq "w") {
      return $result;
   }
   else {
      return undef;
   }
}

sub new {
    my ($class, %args) = @_;
    $maxSize = 100000000;
    return bless {
                     _name =>  $args{'file'} ||  croak('You must provide a file\n'),
                     _permission=>  validPermission($args{'permission'})  || croak('You must provide a file permission w/r'),
                 }, $class;
 }


}

sub getData {
   my ($filename) = @_; 
   unless(open(DATAFILE, $filename)) 
   {   
     print  "Cannot open file $filename \n";
     exit;  
   }   
   my @filedata = <DATAFILE>;

   close(DATAFILE);

   return join('',@filedata);
 }
 

1;
