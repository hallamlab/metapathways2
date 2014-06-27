#! /usr/bin/perl -w

#####################################################
##
## mltreemap_imagemaker.pl
##
#####################################################

use strict;
use warnings;

use File::Basename;
#use lib 'lib';
#use TREEMAP_ml_svg_visualizer;

print "Start MLTreeMap imagemaker v2.061\n";
my $mltreemap_imagemaker                             = mltreemap_imagemaker->new();
my $user_options                                     = $mltreemap_imagemaker->read_user_input(\@ARGV);
my $input_files                                      = $mltreemap_imagemaker->get_input_files($user_options);
my ($concatenated_input_files,$text_of_denominator,$used_colors)  = $mltreemap_imagemaker->concatenate_files($user_options,$input_files);
$mltreemap_imagemaker->run_the_imagemaker($user_options,$concatenated_input_files,$text_of_denominator,$used_colors);
print "Done.\n";


#####################################################
#####################################################

# package mltreemap_imagemaker

#####################################################
#####################################################

package mltreemap_imagemaker;

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
    
    my $error_message = "\nERROR: incorrect input.\nexample input:\n./mltreemap_imagemaker.pl -i example_input/\n\n";
    $error_message .= "mandatory input option:\n-i: either a) a path to a directory with input files or b) to a specific input file.\n";
    $error_message .= "note that in case a) all input files of the same type (i.e. p, g, n, h, r, m, d) will be concatenated.\n\n";
    $error_message .= "optional input options:\n-o: output directory (default: output/)\n-d: use different colors for different datasets ";
    $error_message .= "(0: don't use this mode (default), 1: use this mode, 11: use this mode but concatenate percent values over datasets).";
    $error_message .= " See the user guide for a more detailed description.\n";
    $error_message .= "-t: display trees inclusive text and labels (1, default) or without text (0)\n";
    $error_message .= "-b: parameter for the size of the placement bubbles (default: 1)\n-r display 16s and 18s rRNA hits in different trees (2, default).";
    $error_message .= " or one tree of life (1)\n-h: high quality placement bubbles (1, default) or low quality bubbles (0). High quality bubbles work perfectly";
    $error_message .= "with inkscape but can cause trouble with Adobe Illustrator.\n";
        
    my %mandatory_options = (-i => 0);
    my %option_defaults = (-o => "output/", -d => 0, -b => 1, -r => 2, -h => 1, -t => 1, -m => "i",  -c => "custom_file.txt");
    
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
    
    #second, do some (hard coded) checks on certain input variables:
    my $custom_inputfile = $user_options{"-c"};
    print $custom_inputfile."***\n";
    
    die "$error_message" unless ($user_options{"-b"} =~ /\A\d+\.*\d*\Z/);
    die "$error_message" unless (($user_options{"-d"} eq 0) || ($user_options{"-d"} eq 1) || ($user_options{"-d"} eq 11)) || (-e $custom_inputfile);
    die "$error_message" unless (($user_options{"-r"} eq 1) || ($user_options{"-r"} eq 2));
    die "$error_message" unless (($user_options{"-t"} eq 1) || ($user_options{"-t"} eq 0));
    
    
    
    #add a "/" to the output directory, if not done so by the user.
    
    $user_options{-o} = "$user_options{-o}"."/" unless ($user_options{-o} =~ /\/\Z/);
    #check if the mandatory options have been used:
    
    die "$error_message" unless $arg_count;
    foreach my $mandatory_option (sort {$a cmp $b} keys %mandatory_options) {
        die "$error_message" unless exists ($user_options{$mandatory_option});    
    }
        
    #create the output directories
        
    if (-e $user_options{-o} && !( $user_options{-m} eq "a" ) ) {
        print "WARNING: Your output directory \"$user_options{-o}\" allready exists!\n";
        print "Overwrite [1], quit [2] or change directory [3]?\n";
        
        my $alarm_timer = 30;
        $SIG{ALRM} = sub { die "No input for more than $alarm_timer sek. Exit MLTreeMap.\n"};
        alarm $alarm_timer;
        my $answer = <STDIN>;
        $answer = "2" unless $answer; #quit if something goes wrong with <STDIN>.
        chomp $answer;
        alarm 0;
        while (($answer ne "1") && ($answer ne "2") && ($answer ne "3")) {
            print "invalid input. Please chose 1, 2 or 3.\n";
            alarm $alarm_timer;
            $answer = <STDIN>;
            chomp $answer;
            alarm 0;
        }
        if ($answer == 1) {
            print "Do you really want to overwrite the old output directory?\n";
            print "All data in it will be lost!\n";
            print "Yes [y] or no [n]?\n";         
            alarm $alarm_timer;
            my $answer2 = <STDIN>;
            chomp $answer2;
            alarm 0;
            while (($answer2 ne "y") && ($answer2 ne "n")) {
                print "invalid input. Please chose y or n.\n";
                alarm $alarm_timer;
                $answer2 = <STDIN>;
                chomp $answer2;
                alarm 0;
            }
            if ($answer2 eq "y") {
                my $rm_command = "rm -r $user_options{-o}";
                system ($rm_command);
            } else {
                die "Exit MLTreeMap\n";    
            }
        } elsif ($answer == 2) {
            die "Exit MLTreeMap\n";;
        } elsif ($answer == 3) {
            print "Please enter the new directory.\n";
            alarm 600;
            my $dir = <STDIN>;
            chomp $dir;
            alarm 0;
            $user_options{-o} = $dir;
        }
    }
    my $mkdir_command = "mkdir $user_options{-o}";     
    system ($mkdir_command) unless (-e $user_options{-o});
    
    
    if (-e $custom_inputfile) {
       my $cp_command = "cp $custom_inputfile $user_options{-o}";
       system ($cp_command);
    }
    #done.
    
    return (\%user_options);    
}

#####################################################
# get_input_files
#####################################################

sub get_input_files {
    my $package = shift;
    my $user_options = shift;
    
    my %input_files = ();
    my $input = $$user_options{-i};
    my $rRNA_display_option = $$user_options{-r};
    my $input2 = "/"."$input";    
    if ($input2 =~ /.*\/((.+)_concatenated_.*RAxML_.+.txt)\Z/) {
        if (-e $input) {
            my $filename = $1;
            print $filename."\n";
            my $denominator = $2;
            print $denominator."*********\n";
            $denominator = "A" if (($denominator eq "b" || $denominator eq "a") && ($rRNA_display_option == 1));
            $input_files{$denominator}{$input} = $filename;    
        } else {
            die "ERROR: you entered a file as input ($input), but it does not exist!\n";    
        }
    } else {
        print $input."^^^\n";
        $input .= "/" unless ($input =~ /\/\Z/);
        opendir (PATH, "$input") or die "Error, your input directory ($input) does not exist!\n";
        my @files = readdir PATH;
        closedir (PATH);
        
        my %files2 = ();
        
        foreach my $filename (@files) {
            if ($filename =~ /\A(.+)_concatenated_.*RAxML_.+.txt\Z/) {
                my $denominator = $1;
                #print $denominator."%%%%\n";
                my $filename_long = "$input$filename";
                $denominator = "A" if (($denominator eq "b" || $denominator eq "a") && ($rRNA_display_option == 1));
                $input_files{$denominator}{$filename_long} = $filename; 
            }    
        }       
    }
    return (\%input_files);    
}

#####################################################
#####################################################

use File::Basename;
sub concatenate_files {
    
    my $package = shift;
    my $user_options = shift;
    my $input_files = shift;
    
    my $output_dir = $$user_options{-o};
    my $rRNA_display_option = $$user_options{-r};
    
    my %concatenated_input_files = ();
    my %text_of_denominator = ();
    
    my @colors = "";
    my %used_colors = ();
    my $dirname = dirname(__FILE__);
    
    open (IN, "$dirname/tree_data/available_dataset_colors.txt") or die "ERROR: tree_data/available_dataset_colors.txt does not exist!\n";
    <IN>;
    my $count = 0;
    while (<IN>) {
        chomp $_;
        my ($red, $green, $blue) = split / /;
        next unless (defined $blue);
        my $color_tag = "rgb($red, $green, $blue)";
        $colors[$count] = $color_tag;
        $count++;
    }
        
    my $nr_of_colors = @colors + 0;
           
    foreach my $denominator (sort {$a cmp $b} keys %$input_files) {
        
        #print $denominator."\n";
        my $output_filename = "$denominator"."_concatenated_RAxML_outputs.txt";
        my $output_filename_long = "$output_dir"."$output_filename";
        open (OUTPUT, "> $output_filename_long") or die "ERROR: Can't create $output_filename_long\n";
        
        if ($$user_options{-d} > 0) {
            open (OUT2, "> $output_dir"."$denominator"."_color_legend.txt") or die "Error, can't create the color legend file\n";
        }
        
        $concatenated_input_files{$denominator}{$output_filename_long} = $output_filename;
        
        my %percentages_of_texts = ();
        my $nr_of_files = 0;

        foreach my $filename_long (sort {$a cmp $b} keys %{$$input_files{$denominator}}) {           
            $nr_of_files++;
            die "Error: this directory contains more dme = dirname(__FILE__);:1
atasets than can be displayed with different colors!\n" if ( ($nr_of_files > $nr_of_colors) && ($$user_options{-d} > 0));
            my $attachment = "colorcode_$colors[0]";
            $used_colors{"$colors[0]"} = 1;
            if ($$user_options{-d} > 0) {
                $attachment = "colorcode_$colors[$nr_of_files-1]";
                $used_colors{"$colors[$nr_of_files-1]"} = 1;
                print OUT2 "$colors[$nr_of_files-1]\t$filename_long\n";
            }
            open (FILE, "$filename_long") or die "can't open $filename_long";    
            
            while (<FILE>) {
                chomp $_;
            
                my $weight = "";
                my $text = "";
            
                if (/\APlacement weight (.+)%:(.+)/) {
                    $weight = $1;
                    $text = $2;
                    $text .= $attachment;
                } elsif (/\A# .+ analysis, (.+):/) {
                    $text_of_denominator{$denominator} = $1;
                    $text_of_denominator{$denominator} = "16s rRNA & 18s rRNA" if (($1 =~ /rRNA/) && ($rRNA_display_option == 1));
                    next;
                } elsif (/\A# Phylogenetic analysis(.*)/) {
                    $text_of_denominator{$denominator} = "MLTreeMap tree of life";
                    $text_of_denominator{$denominator} = "GEBA tree of life" if ($1 =~ /GEBA/);
                    $text_of_denominator{$denominator} = "Fungi" if ($1 =~ /fungi/);
                    next;
                } else {
                    next;    
                }
            
                if (exists $percentages_of_texts{$text}) {
                    $percentages_of_texts{$text} += $weight;
                } else {
                    $percentages_of_texts{$text} = $weight;    
                }     
            }
            close FILE;   
        }
        if ($$user_options{-d} > 0) {
            close OUT2;
        }

        my $check = 0;
        foreach my $text (sort {$a cmp $b} keys %percentages_of_texts) {
            my $weight = $percentages_of_texts{$text};
            my $relative_weight = $weight / $nr_of_files;
            $relative_weight = (int($relative_weight * 10000 + 0.5)) / 10000;
            print OUTPUT "Placement weight $relative_weight"."\%:$text\n";
            $check += $relative_weight;
        }    
        print "$denominator files: sum of percentages = $check\n";
        close OUTPUT;
    
    } 
    return(\%concatenated_input_files,\%text_of_denominator,\%used_colors);
}

#####################################################
#####################################################

sub run_the_imagemaker {
    my $package = shift;
    my $user_options = shift;
    my $concatenated_input_files = shift;
    my $text_of_denominator = shift;
    my $used_colors = shift;
    my $param_scale_bubble = $$user_options{-b};
    my $output_dir = $$user_options{-o};
    my $bubble_type = $$user_options{-h};
    my $color_mode = $$user_options{-d};
    my $text_mode = $$user_options{-t};
    my $custom_input = $$user_options{-c};
    my $input_filename = "";
    
    if (-e $custom_input) {
       foreach my $denominator (sort {$a cmp $b} keys %$concatenated_input_files) {
          #print $denominator."\n";
          $input_filename = $custom_input;
          my $custom_input_long = "$output_dir"."$input_filename";
          my $message = TREEMAP_ml_svg_visualizer::run_visualisation ($input_filename,$custom_input_long,$output_dir,$denominator,$param_scale_bubble,$text_of_denominator,$bubble_type,$used_colors,$color_mode,$text_mode);
       print "$message\n";
       }
    } else {
       foreach my $denominator (sort {$a cmp $b} keys %$concatenated_input_files) {
           foreach my $input_filename_long (sort {$a cmp $b} keys %{$$concatenated_input_files{$denominator}}) {
               $input_filename = $$concatenated_input_files{$denominator}{$input_filename_long};
           #    print $input_filename."\n";
           #    print $input_filename_long."\n";
               my $message = TREEMAP_ml_svg_visualizer::run_visualisation ($input_filename,$input_filename_long,$output_dir,$denominator,$param_scale_bubble,$text_of_denominator,$bubble_type,$used_colors,$color_mode,$text_mode);
               print "$message\n";
           }    
       }
    }
}

#####################################################
#####################################################
1;

package TREEMAP_ml_svg_visualizer;

#####################################################
#####################################################

# package TREEMAP_ml_svg_visualizer

#####################################################
#####################################################

use strict;
use warnings;

use SVG;
#use NEWICK_tree;
use Math::Trig;

#####################################################
## subroutine: run_visualisation ()
##
## note: "do_the_drawing_circular" is derived from the linear version. If you have troubles to understand
## the circular subroutines, have a look at their linear counterparts first.
#####################################################

sub run_visualisation {

    my ($input_filename,$input_filename_long,$output_dir,$denominator,$param_scale_bubble,$text_of_denominator,$bubble_type,$used_colors,$color_mode,$text_mode) = @_;
    
    $text_of_denominator = $$text_of_denominator{$denominator};
    print "$text_of_denominator\n";
    
    my $get_coordinates = get_coordinates->new();
    my $read_mltreemap_results = read_MLTreeMap_results->new();
    my $draw_tree = draw_tree->new();
      
    my $mltreemap_results           = $read_mltreemap_results->get_analysis_info($denominator);
    if( $mltreemap_results eq "NIL") { return "failure"; }
    
    $mltreemap_results              = $read_mltreemap_results->get_picture_dimensions($mltreemap_results,$text_mode);
    my $tree                        = $get_coordinates->read_tree_topology($mltreemap_results);

    $mltreemap_results              = $read_mltreemap_results->produce_terminal_children_of_strings_of_reference($tree,$mltreemap_results,$used_colors);
    $mltreemap_results              = $read_mltreemap_results->get_support_data($tree,$mltreemap_results);
    $mltreemap_results              = $read_mltreemap_results->read_RAxML_out($tree,$mltreemap_results,$input_filename_long);
    $draw_tree->draw_the_color_legend($denominator,$output_dir,$color_mode) if $color_mode > 0;
    $draw_tree->do_the_drawing($tree,$mltreemap_results,$param_scale_bubble,$input_filename,$output_dir,$text_of_denominator,$bubble_type,$used_colors,$color_mode,$text_mode);
    $draw_tree->do_the_drawing_circular($tree,$mltreemap_results,$param_scale_bubble,$input_filename,$output_dir,$text_of_denominator,$bubble_type,$used_colors,$color_mode,$text_mode);

        
    return "success";
}

#####################################################
#####################################################

#####################################################
#####################################################

# package read_MLTreeMap_results

#####################################################
#####################################################
1;

package read_MLTreeMap_results;

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
# produce_terminal_children_strings
#####################################################

sub get_analysis_info {
    my $package = shift;
    my $denominator = shift;
    
    my %mltreemap_results = ();
    my $tree_name = "";
    use File::Basename;
    my $dirname = dirname(__FILE__);
 
    open (IN, "$dirname/tree_data/drawing_info.txt") or die "Error, can't open drawing_info.txt!\n";
    <IN>;<IN>;<IN>;
    while (<IN>) {
        chomp $_;
	    #print "Drawing ".$_."\n";
        my ($dn, $tn) = split /\t/;
        #print "<<".$dn.">><<".$denominator.">>\n";
        $tree_name = $tn if ($dn eq "$denominator");       
    }
    if ($denominator eq "p") {
        $mltreemap_results{tree_file} = "MLTreeMap_reference.tree";
        $mltreemap_results{tax_ids_file} = "tax_ids_tree_of_life.txt";
        $mltreemap_results{descriptions_file} = "domain_and_color_descriptions.txt";   
    } elsif ($denominator eq "g") {
        $mltreemap_results{tree_file} = "geba.tree";
        $mltreemap_results{tax_ids_file} = "tax_ids_geba_tree.txt";
        $mltreemap_results{descriptions_file} = "domain_and_color_descriptions_geba_tree.txt";
    } else {
        #die "Error, $denominator is not recognized!\n" unless $tree_name;
        
        return("NIL") unless $tree_name;
	    #print $denominator."\n";
        $mltreemap_results{tree_file} = "$tree_name"."_tree.txt";
        $mltreemap_results{tax_ids_file} = "tax_ids_$tree_name.txt";
        $mltreemap_results{descriptions_file} = "domain_and_color_descriptions_$tree_name.txt";
    }
    
    return (\%mltreemap_results);
}

#####################################################
# get_picture_dimensions
#####################################################

sub get_picture_dimensions {
    my $package = shift;
    my $mltreemap_results = shift;
    my $text_mode = shift;
    my $dirname = dirname(__FILE__);
    
    my $tree_file = $$mltreemap_results{tree_file};
    open (IN, "$dirname/tree_data/$tree_file") or die "Error, can't open $tree_file!\n";
    my $species_count = 0;
    while (<IN>) {
        chomp $_;
        while (/\d+:/g) {
            $species_count++;    
        }    
    }
    close IN;
    
    my $image_width = 1000; #relations are as follows 267 species equal 3000 height, while width remains constant at 1000
    my $tree_height = ($species_count * ($image_width * 3)) / 267;
    my $y_offset = $image_width / 20;
    my $image_height = $tree_height + $y_offset * 2;
    $$mltreemap_results{image}{width} = $image_width;
    $$mltreemap_results{image}{tree_height} = $tree_height;
    $$mltreemap_results{image}{image_height} = $image_height;
    $$mltreemap_results{image}{y_offset} = $y_offset;
    
    #define tree and label margains
    
    $$mltreemap_results{x_coordinate_of_label_start} = 0.6 * $image_width;
    $$mltreemap_results{x_coordinate_of_label_end} = 0.9 * $image_width;
    $$mltreemap_results{x_coordinate_of_tree_end} = 0.7 * $image_width;
    
    unless ($text_mode) {
        #if we don't want the whole text to be displayed, we have to reset some values
        $$mltreemap_results{x_coordinate_of_label_start} = 0.7 * $image_width;
        $$mltreemap_results{x_coordinate_of_label_end} = 0.73 * $image_width;
    
    }
    
    return ($mltreemap_results);
}

#####################################################
# produce_terminal_children_strings
#####################################################

sub produce_terminal_children_of_strings_of_reference {
    my $package = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    my $used_colors = shift;
    my $species_count = 0;
    
    foreach my $n (sort {$a cmp $b} keys %{ $tree->{terminal_children_of_node}}) {
        foreach my $color (sort {$a cmp $b} keys %$used_colors) {
        $$mltreemap_results{counts_per_species}{$n}{$color} = 0; # put all counts to 0
        }
        $species_count++ if ($n > 0);
        my $terminal_children_string = "";
        foreach my $term_child (sort {$a cmp $b} keys %{ $tree->{terminal_children_of_node}{$n}}) {
           $terminal_children_string .= "@"."$term_child";
        }
        $$mltreemap_results{nodes_of_terminal_children_string_reference}{$terminal_children_string} = $n;
    }
    $$mltreemap_results{species_count} = $species_count;
    return ($mltreemap_results);
}
#####################################################
# read_RAxML_out
#####################################################

sub read_RAxML_out {
    my $package = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    my $input_filename_long = shift;
    
    my %nodes = ();
    
    my $input_filename = "$input_filename_long";
    
    print "input file: $input_filename\n";
    unless (-e $input_filename) { die "no such file '$input_filename'"; }
    open (RAXML_OUT, "$input_filename") or die "Can't open $input_filename\n";
    
    my %placements = ();
    my %terminal_children_of_placements = ();
    $$mltreemap_results{highest_count_per_species} = 0;
    
    my $placement = 0;
    while (<RAXML_OUT>) {
        chomp $_;
        $placement++;
        my $color = 0;
        if ($_ =~ /colorcode_(.+)\Z/) {
            $color = $1;
        }
        my $dirname = dirname(__FILE__);

        die "ERROR, color info could not be read!!!!\n" unless $color;
        
        while (/\((\d+)\)/g) {
            $terminal_children_of_placements{$placement}{$1} = 1;
        }
        my $terminal_children_string_of_placement = "";
        foreach my $terminal_child_of_placement (sort {$a cmp $b} keys %{ $terminal_children_of_placements{$placement}}) {
            $terminal_children_string_of_placement .= "@"."$terminal_child_of_placement";
        }
        if (exists $$mltreemap_results{nodes_of_terminal_children_string_reference}{$terminal_children_string_of_placement}) {
            my $node = $$mltreemap_results{nodes_of_terminal_children_string_reference}{$terminal_children_string_of_placement};
            if (/Placement weight (.+)%/)  {            
                my $bootstrap = $1;
                $nodes{$node}{$color} = $bootstrap;
                my $node_weight = $bootstrap / 100;               
                $mltreemap_results = &distribute_node_weight_by_topology ($tree, $node, $node_weight, $mltreemap_results,$color);
            }  
        } 
        else {
            #next; #Attention! this circumvents a problem with the unrooted RAxML tree!
            #die "ERROR: a subtree\n($terminal_children_string_of_placement)\n, as written in the parsed RAxML file, does not exist in the reference tree...\n";
        }
        
    }   
    close RAXML_OUT;
    
    foreach my $node (sort {$a cmp $b} keys %nodes) {
        my $is_first_color = 1;
        foreach my $color (sort {$a cmp $b} keys %{$nodes{$node}}) {
            $$mltreemap_results{counts_per_node}{$node}{"colors"}{$color} = $nodes{$node}{$color} / 100;
            $$mltreemap_results{counts_per_node}{$node}{"total_weight"} = 0 if $is_first_color;
            $is_first_color = 0;
            $$mltreemap_results{counts_per_node}{$node}{"total_weight"} += $nodes{$node}{$color} / 100; 
        }
    }
    return ($mltreemap_results);
}

#####################################################
# distribute_node_weight_by_topology ()
#####################################################

sub distribute_node_weight_by_topology {
    
    my ($tree, $node, $node_weight, $mltreemap_results, $color) = @_;

    my @terminal_children = keys %{$tree->{terminal_children_of_node}{$node}};
    my $nr_of_children = scalar @terminal_children;
    my $fractional_weight = $node_weight / $nr_of_children;
    foreach my $child (@terminal_children) {
        $$mltreemap_results{counts_per_species}{$child}{$color} += $fractional_weight;
        my $total_weight_of_child = 0;
        foreach my $color (sort {$a cmp $b} keys %{$$mltreemap_results{counts_per_species}{$child}}) {
            my $color_weight = $$mltreemap_results{counts_per_species}{$child}{$color};
            $total_weight_of_child += $color_weight; 
        }
        if ($total_weight_of_child > $$mltreemap_results{highest_count_per_species}) {
            $$mltreemap_results{highest_count_per_species} = $total_weight_of_child;
        }
    }
    return ($mltreemap_results);
}

#####################################################
# get_support_data
#####################################################

sub get_support_data {
    my $package = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
                          
    #get species names
    use File::Basename;
    my $dirname = dirname(__FILE__);
 
    open (IN, "$dirname/tree_data/$$mltreemap_results{tax_ids_file}") or die "cannot read 'tree_data/$$mltreemap_results{tax_ids_file}'!\n";
	
    while (<IN>) {
        chomp; next if /\A\#/;
        my ($species, @rest) = split;
        my $name = join " ", @rest;
        $$mltreemap_results{name_of_species}{$species} = $name;
    }
    close (IN);
    
    #get the color demarcations
    #open (IN, "tree_data/$$mltreemap_results{descriptions_file}") or die "Error, can't read tree_data/$$mltreemap_results{descriptions_file}\n";
    return ($mltreemap_results);
    while (<IN>) {
        chomp; 
        my ($start_taxon, $end_taxon, $background_red, $background_green, $background_blue, $subtree_name) = split;
        next if (! $start_taxon || /\A\#/);
        my $group_color = "rgb($background_red,$background_green,$background_blue)";
        my $y_coord_of_node_min = $tree->{y_coord_of_node_min}{$start_taxon};
        $$mltreemap_results{group_info}{$y_coord_of_node_min}{$start_taxon}{$end_taxon}{color} = $group_color;
        $$mltreemap_results{group_info}{$y_coord_of_node_min}{$start_taxon}{$end_taxon}{name} = $subtree_name;
    }
    close (IN);
    
    return ($mltreemap_results);
}


#####################################################
#####################################################
1;

#####################################################
#####################################################

# package draw_tree

#####################################################
#####################################################

package draw_tree;

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
# draw_the_color_legend
#####################################################

sub draw_the_color_legend {
    my $package = shift;
    my $denominator = shift;
    my $output_dir = shift;
    my $color_mode = shift;
    
    open (IN, "$output_dir$denominator"."_color_legend.txt") or die "Error, can't open $output_dir$denominator"."_color_legend.txt\n";
    my $nr_of_datasets = 0;
    my %color_info = ();
    my $namelength_max = 0;    
    while (<IN>) {
        my ($color, $filename) = split /\t/;
        die "Error, something is wrong with color_legend.txt!\n" unless ($color && $filename);
        $color_info{$color} = $filename;
        $namelength_max = length ($filename) if (length ($filename) > $namelength_max);
        $nr_of_datasets++;
    }
       
    close IN;
    
    my $output_file_name = "$output_dir$denominator"."_color_legend.svg";
    open (OUTPUT, "> $output_file_name") or die "Error, can't create $output_file_name!\n";
    my $legend_width = 7 * $namelength_max;
    my $y_offset = 12;
    my $legend_height = $y_offset + $y_offset * $nr_of_datasets;
    my $legend = SVG->new(width=>"$legend_width"."px",height=>"$legend_height"."px");
    
    my $y_pos = $y_offset;
    foreach my $color (sort {$a cmp $b} keys %color_info) {
        my $filename = $color_info{$color};
        $legend->rectangle(x => 10, y => $y_pos - $y_offset, width => $y_offset, height => $y_offset, fill => $color); 
        $legend->text ('x' => $y_offset + 10, 'y' => $y_pos, 'style' => {'font-family'=>'Verdana','font-size'=> 10})->cdata ($filename);
        $y_pos += $y_offset;
    } 
    print OUTPUT $legend->xmlify;
    close OUTPUT;
}

#####################################################
# do_the_drawing
#####################################################

sub do_the_drawing {
    my $package = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    my $param_scale_bubble = shift;
    my $input_file_name = shift;
    my $output_dir = shift;
    my $text_of_denominator = shift;
    my $bubble_type = shift;
    my $used_colors = shift;
    my $color_mode = shift;
    my $text_mode = shift;
    
    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $image_height = $$mltreemap_results{image}{image_height};
    my $y_offset = $$mltreemap_results{image}{y_offset};   
    my $species_count = $$mltreemap_results{species_count};
    
    my $output_file_name = "$output_dir"."$input_file_name"."_image_straight.svg";
    open (OUTPUT, "> $output_file_name") or die "Error, can't create $output_file_name!\n";
    
    my $image = SVG->new(width=>"$image_width"."px",height=>"$image_height"."px");
    $image->text ('x' => $image_width * 0.01, 'y' => $y_offset / 2, 'style' => {'font-family'=>'Verdana','font-size'=> $image_width / 50})->cdata ($text_of_denominator);
    
    my $placements = "";
      
    $image = &draw_group_colors($image,$tree,$mltreemap_results);
    $image = &draw_edges($image,$tree,$mltreemap_results);
    $image = &draw_guide_lines_and_leaf_names($image,$tree,$mltreemap_results,$text_mode);
    ($image,$placements) = &draw_placement_bubbles($image,$tree,$mltreemap_results,$param_scale_bubble,$bubble_type,0);
    $image = &draw_percents_and_placement_bars($image,$placements,$tree,$mltreemap_results,$used_colors,$color_mode,$text_mode);
    
    #print $image->xmlify;
    print OUTPUT $image->xmlify;
    close OUTPUT;
}


#####################################################
# do_the_drawing_circular
#####################################################

sub do_the_drawing_circular {
    my $package = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    my $param_scale_bubble = shift;
    my $input_file_name = shift;
    my $output_dir = shift;
    my $text_of_denominator = shift;
    my $bubble_type = shift;
    my $used_colors = shift;
    my $color_mode = shift;
    my $text_mode = shift;
    
    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $image_height = $$mltreemap_results{image}{image_height};
    my $y_offset = $$mltreemap_results{image}{y_offset};   
    my $species_count = $$mltreemap_results{species_count};
    
    my $image_diameter_circular = 2 * $image_width; # note for large trees, this has to be adjosted.
    $$mltreemap_results{image_circular}{diameter} = $image_diameter_circular;
    
    
    my $output_file_name = "$output_dir"."$input_file_name"."_image_circular.svg";
    open (OUTPUT, "> $output_file_name") or die "Error, can't create $output_file_name!\n";
    
    my $image = SVG->new(width=>"$image_diameter_circular"."px",height=>"$image_diameter_circular"."px");
    $image->text ('x' => $image_width * 0.01, 'y' => $y_offset / 2, 'style' => {'font-family'=>'Verdana','font-size'=> $image_width / 50})->cdata ($text_of_denominator);
    
    my $placements = "";
    ($tree,$mltreemap_results) = &prepare_coordinates_circular($tree,$mltreemap_results);
      
    $image = &draw_group_colors_circular($image,$tree,$mltreemap_results);
    $image = &draw_edges_circular($image,$tree,$mltreemap_results);
    $image = &draw_guide_lines_and_leaf_names_circular($image,$tree,$mltreemap_results,$bubble_type,$text_mode);
    ($image,$placements) = &draw_placement_bubbles($image,$tree,$mltreemap_results,$param_scale_bubble,$bubble_type,1);
    $image = &draw_percents_and_placement_bars_circular($image,$placements,$tree,$mltreemap_results,$bubble_type,$used_colors,$color_mode,$text_mode);
    
    #print $image->xmlify;
    print OUTPUT $image->xmlify;
    close OUTPUT;
}

#####################################################
# prepare_coordinates_circular
#####################################################

sub prepare_coordinates_circular {
    my $tree = shift;
    my $mltreemap_results = shift;
    
    foreach my $node (sort {$a cmp $b} keys %{$tree->{y_coord_of_node}}) {
            
        my $y_coord_of_node = $tree->{y_coord_of_node}{$node};
        my $x_coord_of_node = $tree->{x_coord_of_node}{$node};
        my $branch_length = $tree->{branch_length_of_node}{$node};
        my $x_coord_of_parent_node = $x_coord_of_node - $branch_length;
                
        #transform the coordinates of the node itself
        my ($px,$py) = &calculate_coordinates_circular($mltreemap_results,$x_coord_of_node,$y_coord_of_node);
        $tree->{x_coord_of_node_circular}{$node} = $px;
        $tree->{y_coord_of_node_circular}{$node} = $py;

        #transform the coordinates of the connecting point of the node to the parent
        ($px,$py) = &calculate_coordinates_circular($mltreemap_results,$x_coord_of_parent_node,$y_coord_of_node);
        $tree->{x_coord_of_node_connector_circular}{$node} = $px;
        $tree->{y_coord_of_node_connector_circular}{$node} = $py;
        #save the connecting point information with parent information
        unless ($node == -1) {
            my $parent_of_node = $tree->{parent_of_node}{$node};
            $tree->{connecting_points}{$parent_of_node}{$node}{x_coordinate_of_connector} = $px;
            $tree->{connecting_points}{$parent_of_node}{$node}{y_coordinate_of_connector} = $py;
        }
        #done
        
        if ($node > 0) {
            my $y_coord_of_node_min = $tree->{y_coord_of_node_min}{$node};
            my $y_coord_of_node_max = $tree->{y_coord_of_node_max}{$node};
            
            ($px,$py) = &calculate_coordinates_circular($mltreemap_results,$x_coord_of_node,$y_coord_of_node_min);
            $tree->{x_coord_of_node_min_circular}{$node} = $px;
            $tree->{y_coord_of_node_min_circular}{$node} = $py;
            
            ($px,$py) = &calculate_coordinates_circular($mltreemap_results,$x_coord_of_node,$y_coord_of_node_max);            
            $tree->{x_coord_of_node_max_circular}{$node} = $px;
            $tree->{y_coord_of_node_max_circular}{$node} = $py;  
        }       
    }         
    return ($tree,$mltreemap_results);
}

#####################################################
# calculate_coordinates_circular
# note: some reprogramming could merge this subroutine with the one that draws the pie-chart bubbles.
#####################################################

sub calculate_coordinates_circular {
    my $mltreemap_results = shift;
    my $x_coord_of_point = shift;
    my $y_coord_of_point = shift;
    
    my $image_diameter_circular = $$mltreemap_results{image_circular}{diameter};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $y_offset = $$mltreemap_results{image}{y_offset};
    my $pi = atan2 (1,1) * 4;
    
    #      la
    #  ----------P(px|py)
    #  |        /
    #  |alpha / 
    #  |    /   
    #lb|  /x-pos
    #  |/ 
    #  M
       
    my $alpha = (($y_coord_of_point - $y_offset) * 2 * $pi * 0.95) / $tree_height; #proportion: $tree_hight equals 2 $pi * 0.95.
    my $center_x = $image_diameter_circular / 2;
    my $center_y = $image_diameter_circular / 2;
    
    my $lb = cos($alpha) * $x_coord_of_point;
    my $la = sin($alpha) * $x_coord_of_point;
    my $px = $center_x + $la;
    my $py = $center_y - $lb;
    return ($px,$py);
}

#####################################################
# draw_group_colors_circular
#####################################################

sub draw_group_colors_circular {
    my $image = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    
    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $image_height = $$mltreemap_results{image}{image_height};
    my $y_offset = $$mltreemap_results{image}{y_offset};
    
    my $y_scaling_factor = $tree_height; #note: this works because y values range from 0 - 1.

    my $x_coordinate_of_label_start = $$mltreemap_results{x_coordinate_of_label_start};
    my $x_coordinate_of_label_end = $$mltreemap_results{x_coordinate_of_label_end};
    my $groups = $image->group (id=>"groups");
    my $fontsize = $image_width / 60;
    my $pi = atan2 (1,1) * 4;
    
    # #the group color band looks somewhat as follows:
    #
    #  C-----------D
    #   \         /
    #    \       /
    #     \     /
    #      A---B
    
    #A(xa,ya), B(xb,yb), C(xc,yc), D(xd,yd)
    
    my $is_first_label = 1;
    my $text_y = $y_offset /2;
    my $x_pos_of_text = $image_width * 2 - ($image_width / 4);
    
    foreach my $y_coord_of_node_min (sort {$a <=> $b} keys %{$$mltreemap_results{group_info}}) {
        foreach my $start_taxon (sort {$a <=> $b} keys %{$$mltreemap_results{group_info}{$y_coord_of_node_min}}) {
            foreach my $end_taxon (sort {$a <=> $b} keys %{$$mltreemap_results{group_info}{$y_coord_of_node_min}{$start_taxon}}) {
                my $color = $$mltreemap_results{group_info}{$y_coord_of_node_min}{$start_taxon}{$end_taxon}{color};
                my $name = $$mltreemap_results{group_info}{$y_coord_of_node_min}{$start_taxon}{$end_taxon}{name};
                $name = "" if ($name eq "#");
                      
                my $ya_linear = $y_coord_of_node_min; #eq yc_linear
                my $yb_linear = $tree->{y_coord_of_node_max}{$end_taxon}; #eq yd_linear
                my $xa_linear = $x_coordinate_of_label_start; #eq xb_linear
                my $xc_linear = $x_coordinate_of_label_end; #eq xd_linear
            
                &draw_trapezoid($mltreemap_results,$groups,$color,$xa_linear,$xc_linear,$ya_linear,$yb_linear);
            
                #prepare and write the group labels
                next unless $name;
            
                if ($is_first_label) {
                    my $text = "Group names (clockwise):";
                    $groups->text ('x' => $x_pos_of_text, 'y' => $text_y,'style' => {'font-family'=>'Verdana','font-size'=>$fontsize,'fill'=>"rgb(0,0,0)"})->cdata ($text);
                    $is_first_label = 0;    
                }
            
                $name =~ s/_/ /g;
                if ($color =~ /rgb\((.+),(.+),(.+)\)/) {
                    my $red = $1 - 80;
                    my $green = $2 - 80;
                    my $blue = $3 - 80;
                    $red = 0 if ($red < 0);
                    $green = 0 if ($green < 0);
                    $blue = 0 if ($blue < 0);
                    $color = "rgb($red, $green, $blue)"; #Make the colors darker... Otherwise the writing is almost invisible.    
                } else {
                    die "Parsing error with $color\n";    
                }
                $text_y += ($fontsize * 1.2);
                $groups->text ('x' => $x_pos_of_text, 'y' => $text_y,'style' => {'font-family'=>'Verdana','font-size'=>$fontsize,'fill'=>$color})->cdata ($name);
            }
        }
    }
    return ($image);
}

#####################################################
# draw_edges_circular
#####################################################

sub draw_edges_circular {
    my $image = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    
    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $y_offset = $$mltreemap_results{image}{y_offset};
    my $image_diameter_circular = $$mltreemap_results{image_circular}{diameter};
    
    my %allready_drawn_connector = ();
        
    my $edges = $image->group (id=>"edges");
    my $stroke_width = $image_width / 1000;
     
    
    foreach my $node (sort {$a cmp $b} keys %{$tree->{y_coord_of_node}}) {
        #next if ($node == -1);
        
        my $edge_color = "rgb(100,100,100)";
        
        unless ($node == -1) {
        
            #first get all needed information      
            my $x_coord_of_node = $tree->{x_coord_of_node_circular}{$node};
            my $y_coord_of_node = $tree->{y_coord_of_node_circular}{$node};
        
            my $x_coord_of_node_connector = $tree->{x_coord_of_node_connector_circular}{$node};
            my $y_coord_of_node_connector = $tree->{y_coord_of_node_connector_circular}{$node};
                
            #second, we draw the line from the node to the connection position of its parent.
            $edges->line(x1 => $x_coord_of_node, y1 => $y_coord_of_node, x2 => $x_coord_of_node_connector, y2 => $y_coord_of_node_connector, 'stroke-linecap' => 'round', 'stroke' => $edge_color, 'stroke-width' => $stroke_width, 'stroke-opacity' => 1);
        
        }
        
        #third, we draw the line between the connecting points (if not allready done)
        
        next if (exists $allready_drawn_connector{$node});
        $allready_drawn_connector{$node} = 1;
        
        my $sweep_flag = 1;
        my $large_arc_flag = 0;
        
        my $x1_coordinate_of_connector = "";
        my $y1_coordinate_of_connector = "";
        my $x2_coordinate_of_connector = "";
        my $y2_coordinate_of_connector = "";
        
        my $count = 0;
        foreach my $child (sort {$b <=> $a} keys %{$tree->{connecting_points}{$node}}) {
            if ($count == 0) {
                $x1_coordinate_of_connector = $tree->{connecting_points}{$node}{$child}{x_coordinate_of_connector};
                $y1_coordinate_of_connector = $tree->{connecting_points}{$node}{$child}{y_coordinate_of_connector};
            } else {
                $x2_coordinate_of_connector = $tree->{connecting_points}{$node}{$child}{x_coordinate_of_connector};
                $y2_coordinate_of_connector = $tree->{connecting_points}{$node}{$child}{y_coordinate_of_connector};
            }
            $count++;
        }
        #the points (x1/y1) (x2/y2) have to be sorted otherwise the arc function won't work properly.
        next unless ($x1_coordinate_of_connector);
        if ($y1_coordinate_of_connector > $y2_coordinate_of_connector) {
            my $temp_coordinate = $y1_coordinate_of_connector; 
            $y1_coordinate_of_connector = $y2_coordinate_of_connector;
            $y2_coordinate_of_connector = $temp_coordinate;
            $temp_coordinate = $x1_coordinate_of_connector; 
            $x1_coordinate_of_connector = $x2_coordinate_of_connector;
            $x2_coordinate_of_connector = $temp_coordinate;         
        }
        $sweep_flag = 0 if ($x1_coordinate_of_connector < ($image_diameter_circular / 2));
        #sorting done.
        
        my $radius_of_node = $tree->{x_coord_of_node}{$node};

        $edges->path("opacity"=>"1", 'stroke-width' => $stroke_width, 'stroke' => $edge_color, fill=>"none", d=>"M $x1_coordinate_of_connector,$y1_coordinate_of_connector A $radius_of_node,$radius_of_node 0 $large_arc_flag, $sweep_flag $x2_coordinate_of_connector,$y2_coordinate_of_connector");
    }
    
    return ($image);
    
}

#####################################################
# draw_guide_lines_and_leaf_names_circular
#####################################################

sub draw_guide_lines_and_leaf_names_circular {
    my $image = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    my $bubble_type = shift;
    my $text_mode = shift;
    
    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $y_offset = $$mltreemap_results{image}{y_offset};
    my $species_count = $$mltreemap_results{species_count};
        
    my $x_coordinate_of_label_start = $$mltreemap_results{x_coordinate_of_label_start};    
    my $guide_lines_and_labels = $image->group (id=>"guide_lines_and_labels");
    my $edge_color = "rgb(220,220,220)";
    my $stroke_width = $image_width / 1000;
    my $fontsize = ($image_width / 100);
    $fontsize *= (267 / $species_count) if ($species_count > 267); #the picture has been optimized for 267 species. If we have more, downsize the font.
    my $x_gap = $image_width / 100;
    
    foreach my $node (sort {$a cmp $b} keys %{$tree->{y_coord_of_node}}) {
        next unless ($node > 0);
        my $x_coord_of_node = $tree->{x_coord_of_node}{$node};
        my $y_coord_of_node = $tree->{y_coord_of_node}{$node};
        my $x_coord_of_text = $x_coordinate_of_label_start + $x_gap;
        my $y_coord_of_text = $y_coord_of_node;
        my $max_text_length = 37;
        
        if ($x_coord_of_node  + $x_gap < $x_coordinate_of_label_start) {
            my $x1_pos_linear = $x_coord_of_node + $x_gap / 2;
            my $x2_pos_linear = $x_coordinate_of_label_start - $x_gap / 2;
            my ($x1,$y1) = &calculate_coordinates_circular($mltreemap_results,$x1_pos_linear,$y_coord_of_node);
            my ($x2,$y2) = &calculate_coordinates_circular($mltreemap_results,$x2_pos_linear,$y_coord_of_node);
            $guide_lines_and_labels->line(x1 => $x1, y1 => $y1, x2 => $x2, y2 => $y2, 'stroke-linecap' => 'round', 'stroke' => $edge_color, 'stroke-width' => $stroke_width, 'stroke-opacity' => 1) if $text_mode;           
        } else {
            $x_coord_of_text = $x_coord_of_node + $x_gap if ($x_coord_of_text < $x_coord_of_node + $x_gap);
            $max_text_length = 32;
        }
        
        my $node_name = $$mltreemap_results{name_of_species}{$node};
        if ($node_name =~/\A(.{$max_text_length})./) {
            $node_name = "$1...";    
        }
        
        my $rot_angle = ((($y_coord_of_text - $y_offset) * 360 * 0.95) / $tree_height) - 90; ##proportion: $tree_hight equals 360 * 0.95.
        my ($x_text, $y_text) = &calculate_coordinates_circular($mltreemap_results,$x_coord_of_text,$y_coord_of_text);
        die "Error, $node has no name!\n" unless $node_name;
        
        if ($rot_angle + 90 <= 180) {
            $guide_lines_and_labels->text ('x' => $x_text, 'y' => $y_text, transform=>"rotate($rot_angle $x_text,$y_text)", 'style' => {'font-family'=>"Verdana",'font-size'=>$fontsize})->cdata ($node_name) if $text_mode;
        } else {
            $rot_angle += 180;
            $guide_lines_and_labels->tag("text", 'x' => $x_text, 'y' => $y_text, transform=>"rotate($rot_angle $x_text,$y_text)",style=>{'font-family'=>"Verdana",'font-size'=>$fontsize, ,'text-anchor'=>'end'})->cdata($node_name) if $text_mode;
        }
    }
    
    return ($image);
}

#####################################################
# draw_percents_and_placement_bars_circular
#####################################################

sub draw_percents_and_placement_bars_circular {
    my $image = shift;
    my $placements = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    my $bubble_type = shift;
    my $used_colors = shift;
    my $color_mode = shift;
    my $text_mode = shift;

    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $species_count = $$mltreemap_results{species_count};
    
    # #the placement bar looks somewhat as follows:
    #
    #  C-----------D
    #   \         /
    #    \       /
    #     \     /
    #      A---B
    
    #A(xa,ya), B(xb,yb), C(xc,yc), D(xd,yd)
        
    foreach my $node (sort {$a cmp $b} keys %{$tree->{y_coord_of_node}}) {
        next unless ($node > 0);
             
        #first get all needed information       
        my $y_coord_of_node = $tree->{y_coord_of_node}{$node};
        my $x_coordinate_of_label_end = $$mltreemap_results{x_coordinate_of_label_end};
        my $y_offset = $$mltreemap_results{image}{y_offset};
        my $highest_fraction_raw = $$mltreemap_results{highest_count_per_species} * 100;
        my $fraction_raw = 0;
        
        #prepare the text
        
        my $text = "";
        my $all_fractions_0 = 1;
        my $fraction_total = 0;
        
        foreach my $color (sort {$a cmp $b} keys %$used_colors) {
           #$text as generated here is only used if $color_mode == 1 (i.e. show different datasets in different colors AND percentages for all of them).
           $fraction_raw = $$mltreemap_results{counts_per_species}{$node}{$color} * 100 if (exists $$mltreemap_results{counts_per_species}{$node}{$color});
	       my $fraction = (int ($fraction_raw * 100 + 0.5)) / 100;
	       $fraction_total += $fraction;
	       $fraction = sprintf "%.2f", $fraction;
	    
	       $all_fractions_0 = 0 if ($fraction > 0);
	       $text .= "($fraction\%)";
	    
        }    
	    $fraction_total = sprintf "%.2f", $fraction_total;
	    $text = "(Total: $fraction_total\%)" if ($color_mode == 11);
	    
	    if ($text_mode) {
            my $fontsize = ($image_width / 100);
            $fontsize *= (267 / $species_count) if ($species_count > 267); #the picture has been optimized for 267 species. If we have more, downsize the font.
            my $y_coord_of_text = $y_coord_of_node;
            my $x_gap = $image_width / 200;
            my $x_coord_of_text = $x_coordinate_of_label_end - $x_gap;
        
            next unless ($all_fractions_0 == 0);
        
            my $rot_angle = ((($y_coord_of_text - $y_offset) * 360 * 0.95) / $tree_height) - 90; ##proportion: $tree_hight equals 360 * 0.95.
            my ($x_text, $y_text) = &calculate_coordinates_circular($mltreemap_results,$x_coord_of_text,$y_coord_of_text);
            if ($rot_angle + 90 <= 180) {
                $placements->text ('x' => $x_text, 'y' => $y_text, transform=>"rotate($rot_angle $x_text,$y_text)", 'style' => {'font-family'=>"Verdana",'font-size'=>$fontsize,'text-anchor'=>'end'})->cdata ($text);
            } else {
                $rot_angle += 180;
                $placements->tag("text", 'x' => $x_text, 'y' => $y_text, transform=>"rotate($rot_angle $x_text,$y_text)",style=>{'font-family'=>"Verdana",'font-size'=>$fontsize})->cdata($text);
            }
	    }
                
        my $only_one_color = 1;
        
        #prepare the placement bars
        my $x_offset = 0;    
        
        foreach my $color (sort {$a cmp $b} keys %$used_colors) {
        
            $fraction_raw = $$mltreemap_results{counts_per_species}{$node}{$color} * 100 if (exists $$mltreemap_results{counts_per_species}{$node}{$color});
            next unless ($fraction_raw > 0);
            my $y_min = $tree->{y_coord_of_node_min}{$node};
            my $y_max = $tree->{y_coord_of_node_max}{$node};
            my $height_tot = $y_max - $y_min;
            my $height = $height_tot * 0.9;
            my $ya_linear = $y_min + (($height_tot - $height) / 2); #eq yc_linear
            my $yb_linear = $ya_linear + $height; #eq yd_linear
            
            my $start_y = $y_min + (($height_tot - $height) / 2);
            my $x_gap2 = $image_width / 200;
            my $start_x = $x_coordinate_of_label_end + $x_gap2;      
            my $max_length = ($image_width - $image_width / 100) - ($start_x);
            $start_x += $x_offset;
            my $fractional_length = $max_length * ($fraction_raw / $highest_fraction_raw); #the highest fraction is max_length
            $x_offset += $fractional_length;
            
            my $xa_linear = $start_x; #eq xb_linear
            my $xc_linear = $start_x + $fractional_length; #eq xd_linear
            &draw_trapezoid($mltreemap_results,$placements,$color,$xa_linear,$xc_linear,$ya_linear,$yb_linear);
       }
    }
    return ($image);
}

#####################################################
# draw_trapezoid
#####################################################

sub draw_trapezoid {
    my ($mltreemap_results,$svg,$color,$xa_linear,$xc_linear,$ya_linear,$yb_linear) = @_;
    
    # #the trapezoid looks somewhat as follows:
    #
    #  C-----------D
    #   \         /
    #    \       /
    #     \     /
    #      A---B
    
    #A(xa,ya), B(xb,yb), C(xc,yc), D(xd,yd)
    
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $x_coordinate_of_label_start = $$mltreemap_results{x_coordinate_of_label_start};
    my $x_coordinate_of_label_end = $$mltreemap_results{x_coordinate_of_label_end};
    my $large_arc_flag = 0;
    my $sweep_flag = 1;
    my $sweep_flag2 = 0;
    if (($yb_linear - $ya_linear) >= (($tree_height / 0.95) / 2)) {
        #i.e. the group color band will cover more than 180°
        $large_arc_flag = 1;   
    }
            
    my ($xa,$ya) = &calculate_coordinates_circular($mltreemap_results,$xa_linear,$ya_linear);
    my ($xb,$yb) = &calculate_coordinates_circular($mltreemap_results,$xa_linear,$yb_linear);
    my ($xc,$yc) = &calculate_coordinates_circular($mltreemap_results,$xc_linear,$ya_linear);
    my ($xd,$yd) = &calculate_coordinates_circular($mltreemap_results,$xc_linear,$yb_linear);
    my $radius_AB = $x_coordinate_of_label_start;
    my $radius_CD = $x_coordinate_of_label_end;
            
    $svg->path("opacity"=>"1", fill=>"$color", d=>"M $xa,$ya A $radius_AB,$radius_AB 0 $large_arc_flag, $sweep_flag $xb,$yb L $xd,$yd A $radius_CD,$radius_CD 0 $large_arc_flag, $sweep_flag2, $xc,$yc z"); 
}

#####################################################
# draw_group_colors
#####################################################

sub draw_group_colors {
    my $image = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    
    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $y_scaling_factor = $tree_height; #note: this works because y values range from 0 - 1.

    my $x_coordinate_of_label_start = $$mltreemap_results{x_coordinate_of_label_start};
    my $x_coordinate_of_label_end = $$mltreemap_results{x_coordinate_of_label_end};
    my $groups = $image->group (id=>"groups");
    my $rx = $image_width / 400;
    my $ry = $rx;
    my $fontsize = $image_width / 125;
    
    foreach my $y_coord_of_node_min (sort {$a <=> $b} keys %{$$mltreemap_results{group_info}}) {
        foreach my $start_taxon (sort {$a <=> $b} keys %{$$mltreemap_results{group_info}{$y_coord_of_node_min}}) {
            foreach my $end_taxon (sort {$a <=> $b} keys %{$$mltreemap_results{group_info}{$y_coord_of_node_min}{$start_taxon}}) {
                my $color = $$mltreemap_results{group_info}{$y_coord_of_node_min}{$start_taxon}{$end_taxon}{color};
                my $name = $$mltreemap_results{group_info}{$y_coord_of_node_min}{$start_taxon}{$end_taxon}{name};
                $name = "" if ($name eq "#");
                my $start_y = $y_coord_of_node_min;
                my $end_y = $tree->{y_coord_of_node_max}{$end_taxon};
                my $width = $x_coordinate_of_label_end - $x_coordinate_of_label_start;
                my $height = $end_y - $start_y; 
                $groups->rectangle(x => $x_coordinate_of_label_start, y => $start_y, rx => $rx, ry => $ry, width => $width, height => $height, fill => $color);
                #prepare and write the group labels
                my $x_pos_of_text = $image_width - ($image_width / 100);
                $name =~ s/_/ /g;
                if ($color =~ /rgb\((.+),(.+),(.+)\)/) {
                    my $red = $1 - 80;
                    my $green = $2 - 80;
                    my $blue = $3 - 80;
                    $red = 0 if ($red < 0);
                    $green = 0 if ($green < 0);
                    $blue = 0 if ($blue < 0);
                    $color = "rgb($red, $green, $blue)"; #Make the colors darker... Otherwise the writing is almost invisible.    
                } else {
                    die "Parsing error with $color\n";    
                }
                $groups->text ('x' => $x_pos_of_text, 'y' => $start_y,'transform'=>"rotate(90 $x_pos_of_text,$start_y)",'style' => {'font-family'=>'Verdana','font-size'=>$fontsize,'fill'=>$color})->cdata ($name);          
            }
        }
    }
    return ($image);
}

#####################################################
# draw_edges
#####################################################

sub draw_edges {
    my $image = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    
    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $y_offset = $$mltreemap_results{image}{y_offset}; 
    
    my %allready_drawn_parents = ();
        
    my $edges = $image->group (id=>"edges");
    my $stroke_width = $image_width / 1000;
     
    # the tree looks schematically as follows and thus has horizontal and vertical lines:
    #     |-----
    # ----|
    #     |-------
    
    foreach my $node (sort {$a cmp $b} keys %{$tree->{y_coord_of_node}}) {
        next if ($node == -1);
        
        my $edge_color = "rgb(100,100,100)";
        
        #first get all needed information
        
        my $parent_of_node = $tree->{parent_of_node}{$node};
        my $branch_length = $tree->{branch_length_of_node}{$node};
        
        my $x_coord_of_node = $tree->{x_coord_of_node}{$node};
        my $y_coord_of_node = $tree->{y_coord_of_node}{$node};
        
        my $x_coord_of_parent_node = $tree->{x_coord_of_node}{$parent_of_node};
        my $y_coord_of_parent_node_center = $tree->{y_coord_of_node}{$parent_of_node} ;
        my $y_distance_parent_to_node = $y_coord_of_node - $y_coord_of_parent_node_center;
        
        #second, we draw the horizontal line from the node to the x- position to its parent.
                
        $edges->line(x1 => $x_coord_of_node, y1 => $y_coord_of_node, x2 => $x_coord_of_parent_node, y2 => $y_coord_of_node, 'stroke-linecap' => 'round', 'stroke' => $edge_color, 'stroke-width' => $stroke_width, 'stroke-opacity' => 1);
        
        #third, we draw the vertical line (if not allready done)
        
        next if (exists $allready_drawn_parents{$parent_of_node});
        $allready_drawn_parents{$parent_of_node} = 1;
        
        $edges->line(x1 => $x_coord_of_parent_node, y1 => $y_coord_of_node, x2 => $x_coord_of_parent_node, y2 => ($y_coord_of_node - 2 * $y_distance_parent_to_node), 'stroke-linecap' => 'round', 'stroke' => $edge_color, 'stroke-width' => $stroke_width, 'stroke-opacity' => 1);
    }
    
    return ($image);
    
}

#####################################################
# draw_guide_lines_and_leaf_names
#####################################################

sub draw_guide_lines_and_leaf_names {
    my $image = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    my $text_mode = shift;

    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    
    my $x_coordinate_of_label_start = $$mltreemap_results{x_coordinate_of_label_start};    
    my $guide_lines_and_labels = $image->group (id=>"guide_lines_and_labels");
    my $edge_color = "rgb(220,220,220)";
    my $stroke_width = $image_width / 1000;
    my $fontsize = $image_width / 125;
    my $y_offset2 = 0.3 * $fontsize;
    my $x_gap = $image_width / 100;
    foreach my $node (sort {$a cmp $b} keys %{$tree->{y_coord_of_node}}) {
        next unless ($node > 0);
        my $x_coord_of_node = $tree->{x_coord_of_node}{$node};
        my $y_coord_of_node = $tree->{y_coord_of_node}{$node};
        my $x_coord_of_text = $x_coordinate_of_label_start + $x_gap;
        my $max_text_length = 47;
        if ($x_coord_of_node  + $x_gap < $x_coordinate_of_label_start) {
            $guide_lines_and_labels->line(x1 => $x_coord_of_node + $x_gap / 2, y1 => $y_coord_of_node, x2 => $x_coordinate_of_label_start - $x_gap / 2, y2 => $y_coord_of_node, 'stroke-linecap' => 'round', 'stroke' => $edge_color, 'stroke-width' => $stroke_width, 'stroke-opacity' => 1) if $text_mode;           
        } else {
            $x_coord_of_text = $x_coord_of_node + $x_gap if ($x_coord_of_text < $x_coord_of_node + $x_gap);
            $max_text_length = 42;
        }
        my $node_name = $$mltreemap_results{name_of_species}{$node};
        if ($node_name =~/\A(.{$max_text_length})./) {
            $node_name = "$1...";    
        }
        die "Error, $node has no name!\n" unless $node_name;
        $guide_lines_and_labels->text ('x' => $x_coord_of_text, 'y' => $y_coord_of_node + $y_offset2, 'style' => {'font-family'=>'Verdana','font-size'=>$fontsize})->cdata ($node_name) if $text_mode;
    }
    
    return ($image);
}

#####################################################
# draw_placement_bubbles
# note: this subroutine could easily be merged into "draw_edges".
# the only reason why it exists is that I want the tree drawing to be separated from "results drawing".
#####################################################

sub draw_placement_bubbles {
    my $image = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    my $param_scale_bubble = shift;
    my $bubble_type = shift;
    my $is_circular_image = shift;
    
    my $pi = atan2 (1,1) * 4;
    my $placements = $image->group (id=>"placements");
    
    #prepare the code for the placement bubble, initial diameter 40, do not change, the bubble size is adjusted later.
    my $radius = 20;
    $placements = &create_aquabubblebody($placements, $radius);
    $placements = &create_aquabubblebrilliance($placements, $radius);
    # done
    my $placement_bubbles = $placements->group (id=>"placement_bubbles");
    
    #we want to print the big bubbles first so that they cannot completely cover smaller ones.
    
    my %total_weights = ();
    
    foreach my $node (sort {$a cmp $b} keys %{$tree->{y_coord_of_node}}) {
        next if ($node == -1);               
        my $total_node_weight = $$mltreemap_results{counts_per_node}{$node}{"total_weight"};
        next unless (defined $total_node_weight);
        $total_weights{$total_node_weight}{$node} = 1;
        
    }
    
    #ok, now draw the bubble
    
    foreach my $total_node_weight (sort {$b <=> $a} keys %total_weights) {  
        foreach my $node (sort {$a cmp $b} keys %{$total_weights{$total_node_weight}}) {
        
            #if $total node weight == 1, $radius_real == 20 * $param_scale_bubble
            my $radius_real = sqrt($total_node_weight * 400) * $param_scale_bubble;
            my $transform_factor = $radius_real / $radius;
            #$transform_factor = 1 if ($transform_factor == 0); #activate this to print empty trees

            my $branch_length_of_node = $tree->{branch_length_of_node}{$node};
            my $x_coord_of_node_placement_center = $tree->{x_coord_of_node}{$node} - $branch_length_of_node / 2;
            my $y_coord_of_node = $tree->{y_coord_of_node}{$node};            
        
            #transform the coordinates if we do the circular image
            if ($is_circular_image) {
                ($x_coord_of_node_placement_center,$y_coord_of_node) = &calculate_coordinates_circular($mltreemap_results,$x_coord_of_node_placement_center,$y_coord_of_node);  
            }
            #done
            $transform_factor = 1 unless $transform_factor;
            my $x_coord_of_node_placement = ($x_coord_of_node_placement_center - $radius_real) * (1/$transform_factor);
            my $y_coord_of_node_placement = ($y_coord_of_node - $radius_real) * (1/$transform_factor);
           
            #color the bubble
            my $previous_end_angle = $pi / 2;        
            foreach my $color (sort {$a cmp $b} keys %{$$mltreemap_results{counts_per_node}{$node}{"colors"}}) {
                my $fraction = $$mltreemap_results{counts_per_node}{$node}{"colors"}{$color} if (exists $$mltreemap_results{counts_per_node}{$node}{"colors"}{$color});
                next if $fraction < 0; 
                next if ($total_node_weight <= 0);           
                my $only_one_hit = 0;
                $only_one_hit = 1 if $fraction eq $total_node_weight;
            
                my $fract_of_bubble = $fraction / $total_node_weight;
       	
                my $start_angle = $previous_end_angle;
                if ($only_one_hit) {
                    $placement_bubbles->circle("opacity"=>"1", fill=>"$color",  "cx"=>$x_coord_of_node_placement_center, "cy"=>$y_coord_of_node, "r"=>$radius_real);
                } else {
                    ($placement_bubbles, $previous_end_angle) = &draw_pie_chart($placement_bubbles, $x_coord_of_node_placement_center, $y_coord_of_node, $radius_real, $start_angle, $fract_of_bubble, $color);
                }  
            }
            if ($bubble_type) {
            #draw the bubble 
                $placement_bubbles->tag("use", "xlink:href"=>"#body",x=>$x_coord_of_node_placement, y=>$y_coord_of_node_placement, "transform"=>"scale($transform_factor)");
                $placement_bubbles->tag("use", "xlink:href"=>"#brilliance",x=>$x_coord_of_node_placement, y=>$y_coord_of_node_placement, "transform"=>"scale($transform_factor)");
            }
        }
    }    

    return ($image,$placements);    
}

#####################################################
# draw_pie_chart
#####################################################

sub draw_pie_chart {
    my $placement_bubbles = shift;
    my $center_x = shift;
    my $center_y = shift;
    my $radius = shift;
    my $alpha_r = shift;
    my $fract_of_bubble = shift;
    my $color = shift;
    
    #get the coordinates
    
    #to do: now we go counter clockwise. To go clockwise would be more intuitive, but more difficult to code.
    
    # C (center_x|center_y) = center of circle; P (end_x|end_y), Q (start_x, start_y); alphatot = alpha_r + alpha
    # P-------- Q
    # |alpha  / |
    # |     /   | 
    # |alp/hatot|
    # | /alpha_r|
    # C----------
    #  
     
    my $pi = atan2 (1,1) * 4;
     
    my $alpha = ($fract_of_bubble) * 2 * $pi;
    my $alphatot = $alpha + $alpha_r;
        
    my $x_start = $center_x + cos($alpha_r) * $radius;
    my $y_start = $center_y - sin($alpha_r) * $radius;
        
    my $x_end = $center_x + cos($alphatot) * $radius;
    my $y_end = $center_y - sin($alphatot) * $radius;
    
    # next get the large_arc_flag
    
    my $large_arc_flag = 0;
    $large_arc_flag = 1 if ($alpha > $pi);
    
    # draw the pie
      
    $placement_bubbles->path("opacity"=>"1", fill=>"$color", d=>"M $center_x,$center_y L $x_start,$y_start A $radius,$radius 0 $large_arc_flag, 0 $x_end,$y_end z");
    
    return ($placement_bubbles,$alphatot);
}

#####################################################
# draw_percents_and_placement_bars
#####################################################

sub draw_percents_and_placement_bars {
    my $image = shift;
    my $placements = shift;
    my $tree = shift;
    my $mltreemap_results = shift;
    my $used_colors = shift;
    my $color_mode = shift;
    my $text_mode = shift;
    
    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
        
    foreach my $node (sort {$a cmp $b} keys %{$tree->{y_coord_of_node}}) {
        next unless ($node > 0);
             
        #first get all needed information       
        my $y_coord_of_node = $tree->{y_coord_of_node}{$node};
        my $x_coordinate_of_label_end = $$mltreemap_results{x_coordinate_of_label_end};
        my $highest_fraction_raw = $$mltreemap_results{highest_count_per_species} * 100;
        my $fraction_raw = 0;
                
        #prepare the text
        
        my $text = "";
        my $all_fractions_0 = 1;
        my $fraction_total = 0;
        
        foreach my $color (sort {$a cmp $b} keys %$used_colors) {
           #$text as generated here is only used if $color_mode == 1 (i.e. show different datasets in different colors AND percentages for all of them).
           $fraction_raw = $$mltreemap_results{counts_per_species}{$node}{$color} * 100 if (exists $$mltreemap_results{counts_per_species}{$node}{$color});
	       my $fraction = (int ($fraction_raw * 100 + 0.5)) / 100;
	       $fraction_total += $fraction;
	       $fraction = sprintf "%.2f", $fraction;
	    
	       $all_fractions_0 = 0 if ($fraction > 0);
	       $text .= "($fraction\%)";
	    
        }    
	    $fraction_total = sprintf "%.2f", $fraction_total;
	    $text = "(Total: $fraction_total\%)" if ($color_mode == 11);
        	    
	    if ($text_mode) {
	       my $fontsize = $image_width / 125;
            my $y_offset2 = 0.3 * $fontsize;
            my $x_gap = $image_width / 400;
        
            next unless ($all_fractions_0 == 0);
        
            $placements->text ('x' => $x_coordinate_of_label_end - $x_gap, 'y' => $y_coord_of_node + $y_offset2, 'style' => {'font-family'=>'Verdana','font-size'=>$fontsize},'text-anchor'=>'end')->cdata ("$text");
	    }
	    
        my $only_one_color = 1;
        
        #prepare the placement bars
        
        my $x_offset = 0;
        foreach my $color (sort {$a cmp $b} keys %{$$mltreemap_results{counts_per_species}{$node}}) {
        
            $fraction_raw = $$mltreemap_results{counts_per_species}{$node}{$color} * 100 if (exists $$mltreemap_results{counts_per_species}{$node}{$color});
            next unless $fraction_raw > 0;
            my $y_min = $tree->{y_coord_of_node_min}{$node};
            my $y_max = $tree->{y_coord_of_node_max}{$node};
            my $height_tot = $y_max - $y_min;
            my $height = $height_tot * 0.9;
            my $start_y = $y_min + (($height_tot - $height) / 2);
            my $x_gap2 = $image_width / 200;
            my $start_x = $x_coordinate_of_label_end + $x_gap2;      
            my $max_length = ($image_width - $image_width / 100) - ($start_x);
            $start_x += $x_offset;
            my $fractional_length = $max_length * ($fraction_raw / $highest_fraction_raw); #the highest fraction is max_length
            $x_offset += $fractional_length;
            $placements->rectangle(x => $start_x, y => $start_y, width => $fractional_length, height => $height, fill => $color);
       }
    }
    return ($image);
}

#####################################################
# create_aquabbublebody
#####################################################

sub create_aquabubblebody {
    
    my $placements = shift;
    my $radius = shift;
    my $svg = $placements->defs();
    
    #now the aquabubblebody:
    
    #partI
    
    my $bubble_gradient1 = $svg->tag("linearGradient", id=>"bubble_gradient1", x1=>"50%", y1=>"100%", x2=>"50%", y2=>"0%");

    $bubble_gradient1->tag("stop",  "offset"=>"0.2472", "style"=>"stop-color:#FAFAFA");
    $bubble_gradient1->tag("stop",  "offset"=>"0.3381", "style"=>"stop-color:#D0D0D0");
    $bubble_gradient1->tag("stop",  "offset"=>"0.4517", "style"=>"stop-color:#A2A2A2");
    $bubble_gradient1->tag("stop",  "offset"=>"0.5658", "style"=>"stop-color:#7C7C7C");
    $bubble_gradient1->tag("stop",  "offset"=>"0.6785", "style"=>"stop-color:#5E5E5E");
    $bubble_gradient1->tag("stop",  "offset"=>"0.7893", "style"=>"stop-color:#494949");
    $bubble_gradient1->tag("stop",  "offset"=>"0.8975", "style"=>"stop-color:#3C3C3C");
    $bubble_gradient1->tag("stop",  "offset"=>"1", "style"=>"stop-color:#383838");

    $svg->circle("id"=>"bubble_part1", "opacity"=>"0.4", "fill"=>"url(#bubble_gradient1)", "cx"=>"$radius","cy"=>"$radius","r"=>"$radius");   
    
    #partI done

    #partII
    
    my $bubble_gradient2 = $svg->tag("radialGradient", id=>"bubble_gradient2", cx=>"50%", cy=>"50%", r=>"50%");
    
    $bubble_gradient2->tag("stop",  "offset"=>"0", "style"=>"stop-color:#FFFFFF");
    $bubble_gradient2->tag("stop",  "offset"=>"0.3726", "style"=>"stop-color:#FDFDFD");
    $bubble_gradient2->tag("stop",  "offset"=>"0.5069", "style"=>"stop-color:#F6F6F6");
    $bubble_gradient2->tag("stop",  "offset"=>"0.6026", "style"=>"stop-color:#EBEBEB");
    $bubble_gradient2->tag("stop",  "offset"=>"0.68", "style"=>"stop-color:#DADADA");
    $bubble_gradient2->tag("stop",  "offset"=>"0.7463", "style"=>"stop-color:#C4C4C4");
    $bubble_gradient2->tag("stop",  "offset"=>"0.805", "style"=>"stop-color:#A8A8A8");
    $bubble_gradient2->tag("stop",  "offset"=>"0.8581", "style"=>"stop-color:#888888");
    $bubble_gradient2->tag("stop",  "offset"=>"0.9069", "style"=>"stop-color:#626262");
    $bubble_gradient2->tag("stop",  "offset"=>"0.9523", "style"=>"stop-color:#373737");
    $bubble_gradient2->tag("stop",  "offset"=>"0.9926", "style"=>"stop-color:#090909");
    $bubble_gradient2->tag("stop",  "offset"=>"1", "style"=>"stop-color:#000000");

    $svg->circle("id"=>"bubble_part2", "opacity"=>"0.1", fill=>"url(#bubble_gradient2)", "cx"=>"$radius", "cy"=>"$radius", "r"=>"$radius");
    
    my $group1_svg = $svg->group("id"=>"body");
    $group1_svg->tag("use", "xlink:href"=>"#bubble_part1");
    $group1_svg->tag("use", "xlink:href"=>"#bubble_part2");
    
    #part II done
    #aquabubblebody done
    return ($placements);
}


#####################################################
# create_aquabubblebrilliance
#####################################################

sub create_aquabubblebrilliance {
    
    my $placements = shift;
    my $radius = shift;
    my $svg = $placements->defs();
    
    my $cy = $radius / 2.5;
    my $rx = $radius / 2;
    my $ry = $radius / 4;
    
    #now the brillance effect:
    	
    #this way of generating a opacity gradient generates valid SVG code but cannot be interpreted by Adobe Illustrator CS3.    
    my $brilliance_gradient = $svg->tag("linearGradient", "id"=>"brilliance_gradient", x1=>"50%", y1=>"0%", x2=>"50%", y2=>"100%");
    $brilliance_gradient->tag("stop",  "offset"=>"0", "style"=>"stop-color:#FFFFFF; stop-opacity:1");
    $brilliance_gradient->tag("stop",  "offset"=>"0.1", "style"=>"stop-color:#FFFFFF; stop-opacity:0.99");
    $brilliance_gradient->tag("stop",  "offset"=>"1", "style"=>"stop-color:#FFFFFF; stop-opacity:0");

    $svg->ellipse("id"=>"bubble_brilliance", "opacity"=>"0.9", "fill"=>"url(#brilliance_gradient)", "cx"=>"$radius","cy"=>"$cy", "rx"=>"$rx","ry"=>"$ry");

    my $group1_svg = $svg->group("id"=>"brilliance");
    $group1_svg->tag("use", "xlink:href"=>"#bubble_brilliance");
    return ($placements);
}

#####################################################
#####################################################
1;


#####################################################
#####################################################

# package get_coordinates

#####################################################
#####################################################

package get_coordinates;

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
# read_tree_topology
#####################################################

sub read_tree_topology {
    my $package = shift;
    my $mltreemap_results = shift;

    print "reading tree topology ...\n";
    my $tree_file = "tree_data/$$mltreemap_results{tree_file}";
    my $tree = NEWICK_tree->new();
    $tree->read_tree_topology_from_file ($tree_file);
    $tree->optimize_species_display_order;

    #do some hard coded manipulations
    $tree->{branch_length_of_node}{-1} = 0.01;
    #done 

    &compute_node_positions ($tree);
    ($tree,$mltreemap_results) = &scale_node_positions($tree,$mltreemap_results);
    return ($tree);
}

########################################################
# scale_node_positions ()
# this subroutine scales the internal coordinates in order to fit them on the picture.
#########################################################

sub scale_node_positions {
    my $tree = shift;
    my $mltreemap_results = shift;
    
    my $image_width = $$mltreemap_results{image}{width};
    my $tree_height = $$mltreemap_results{image}{tree_height};
    my $y_offset = $$mltreemap_results{image}{y_offset};       
    my $y_scaling_factor = $tree_height; #note: this works because y values range from 0 - 1.
        
    #assign the highest possible tree x coordinate and calculate the x_scaling factor
    my $x_scaling_factor = 0;
    foreach my $x_val (sort {$b <=> $a} keys %{$tree->{nodes_of_x_coords}}) {
        $x_scaling_factor = $$mltreemap_results{x_coordinate_of_tree_end} / $x_val;
        die "Error, x scaling factor could not be calculated!\n" unless defined $x_scaling_factor;
        last;
    }
    die "Error, x scaling factor could not be determined!\n" unless $x_scaling_factor;
    #calculate the coordinates relative to the picture

    foreach my $node (sort {$a cmp $b} keys %{$tree->{y_coord_of_node}}) {
        my $y_coord_of_node = ($tree->{y_coord_of_node}{$node} * $y_scaling_factor) + $y_offset;
        my $x_coord_of_node = $tree->{x_coord_of_node}{$node} * $x_scaling_factor; 
        my $branch_length_of_node = $tree->{branch_length_of_node}{$node} * $x_scaling_factor;      
        $tree->{y_coord_of_node}{$node} = $y_coord_of_node;
        $tree->{x_coord_of_node}{$node} = $x_coord_of_node;
        $tree->{branch_length_of_node}{$node} = $branch_length_of_node;
        
        if ($node > 0) {
            my $y_coord_of_node_min = ($tree->{y_coord_of_node_min}{$node} * $y_scaling_factor) + $y_offset;
            my $y_coord_of_node_max = ($tree->{y_coord_of_node_max}{$node} * $y_scaling_factor) + $y_offset;
            $tree->{y_coord_of_node_min}{$node} = $y_coord_of_node_min;
            $tree->{y_coord_of_node_max}{$node} = $y_coord_of_node_max;
        }       
    }       
    return ($tree);
}

########################################################
# subroutine: compute_node_positions ()
#
# This routine performs the actual placement of the tree (branches, leafs) in terms of 
# x,y-coordinate with which they will later be put on the canvas. This routine does not
# perform the actual drawing. 
########################################################

sub compute_node_positions {
    
    my ($tree) = @_;
    
    my $y_position = 0;
        
    foreach my $node (sort { my $a_sort = 1; $a_sort = $tree->{species_display_order}{$a} if exists $tree->{species_display_order}{$a}; my $b_sort = 1; $b_sort = $tree->{species_display_order}{$b} if exists $tree->{species_display_order}{$b}; return $a_sort <=> $b_sort;} keys %{$tree->{parent_of_node}}) {
        next unless $node > 0;

        my $fraction = 1 / (scalar keys %{$tree->{terminal_nodes}});
        $tree->{y_coord_of_node_min}{$node} = $y_position;
        $tree->{y_coord_of_node}{$node} = $y_position + ($fraction * 0.5);
        $tree->{y_coord_of_node_max}{$node} = $y_position + $fraction;
        $y_position += $fraction;
    }

    &assign_y_coord_internal_node ($tree, -1);
    &assign_x_coord_to_node ($tree, -1, 0.01); #note 0.01 is the hardcoded branch length of node -1

    ## CVM WATCH: the vertical position of node -2 is meddled with here ... in order 
    ## to avoid the 'kink' at the root. 

}

########################################################
# subroutine: assign_y_coord_internal_node ()
#
# This recursive routine assigns the Y-position for all internal nodes. It also assigns
# positions needed for the 'taxon-illustration' bit. Call this routine only after
# you have already assigned Y-positions to the terminal (leaf) nodes.
########################################################

sub assign_y_coord_internal_node {
    
    my ($tree, $this_node) = @_;

    if (exists $tree->{y_coord_of_node}{$this_node}) {
        my $position = $tree->{y_coord_of_node}{$this_node};
        return $position;
    }

    return 0 if $this_node > 0;        ## for security. should never happen.

    my $min_position = 100000;
    my $max_position = 0;

    foreach my $child (keys %{$tree->{children_of_node}{$this_node}}) {
        my $position = &assign_y_coord_internal_node ($tree, $child);
        $min_position = $position if $position < $min_position;
        $max_position = $position if (($position > $max_position) and ($child < 1000000));
    }
    my $this_position = ($min_position + $max_position) / 2;

    $tree->{y_coord_of_node}{$this_node} = $this_position;

    return $this_position;
}

########################################################
# subroutine: assign_x_coord_internal_node ()
#
########################################################

sub assign_x_coord_to_node {

    my ($tree, $this_node, $current_x_position) = @_;

    $tree->{x_coord_of_node}{$this_node} = $current_x_position;
    $tree->{nodes_of_x_coords}{$current_x_position} = $this_node;
    
    if ($this_node < 0) {
	
        foreach my $child (keys %{$tree->{children_of_node}{$this_node}}) {
            my $branch_length = 0.02;
            if ((exists $tree->{branch_length_of_node}{$child}) and (defined $tree->{branch_length_of_node}{$child})) {
                $branch_length = $tree->{branch_length_of_node}{$child};
            } else {
                print "warning: no branch-length for node '$this_node'!\n";
            }
            &assign_x_coord_to_node ($tree, $child, $current_x_position + $branch_length);
        }
    }
}

#####################################################
#####################################################
1;

package NEWICK_tree;

use strict;
use warnings;

#########################################################################################
## constructor.
##
## this just makes the tree object and sets all values to 'empty'
#########################################################################################

sub new {

    my ($that) = @_;

    my $class = ref ($that) || $that;          
    my $self = {};                             
    bless $self, $class;

    $self->{node_id_counter} = -1;
    $self->{parent_of_node} = {};
    $self->{children_of_node} = {};
    $self->{terminal_children_of_node} = {};
    $self->{branch_length_of_node} = {};
    
    return $self;
}

####################################################################################################
## subroutine: read_tree_topology_from_file ()
##
## This reads the topology of the species-tree from a flat-file.
####################################################################################################

sub read_tree_topology_from_file {

    my ($self, $filename) = @_;

    my $input_concatenated = "";


    use File::Basename;
    my $dirname = dirname(__FILE__);
    open (FH, "$dirname/$filename") or die "cannot open inputfile '$filename'\n";

    while (<FH>) {
        chomp;
        s/\s+//g;
        $input_concatenated .= $_;
    }
    close FH;

    $self->parse_tree_topology_from_string ($input_concatenated);
}

####################################################################################################
## subroutine: parse_tree_topology_from_string ()
##
##
####################################################################################################

sub parse_tree_topology_from_string {

    my ($self, $tree_string) = @_;

    my @tree_input_symbols = split //, $tree_string;

    $self->species_tree_parse_node (\@tree_input_symbols);

    $self->connect_terminal_children_to_node (-1);

    $self->assign_cumulative_branchlength_to_node (-1, 0);

    $self->compute_average_cumulative_branch_length ();
}

##################################################################################################################
## subroutine: species_tree_parse_node ()
##
## This is a recursive routine to parse the species-tree topology out of an array of input symbols. The input
## symbols are brackets, letters, and numbers. The routine was initially written by Daniel, and is quite clever.
##################################################################################################################

sub species_tree_parse_node {

    my ($self, $inputarray) = @_;
    my $next_item;
    my $this_node = $self->{node_id_counter}--;
    my $nr = 0;
    

    my @children_of_this_node = ();

    ## first symbol must be open bracket.

    unless (($next_item = shift @$inputarray) eq "(") { die "#1 error parsing species tree\n"; }

    while (1) {  ## now loop over all subnodes of this node ...

        ## read next symbol: must be either open bracket or species_name.

        unless ($next_item = shift @$inputarray) { die "#2 error parsing species tree\n"; }

        if ($next_item eq "(") { ## open bracket ? call yourself recursively

            unshift @$inputarray, $next_item;
            my $sub_node = $self->species_tree_parse_node ($inputarray);
            $self->{parent_of_node}->{$sub_node} = $this_node;
            push @children_of_this_node, $sub_node;

        } elsif ($next_item =~ /\w+/) { ## word character ? this is a terminal node ('leaf' - species)

            my $node_summary = $next_item;
            while (($next_item = shift @$inputarray) =~ /[\w\.\:]+/) {
                $node_summary .= $next_item;
            }

	    my ($species, $branch_length) = (undef, undef);

	    if ($node_summary =~ /\A(\w+)\:([\d\.]+\z)/) {            ## are there branch-lengths ?
		($species, $branch_length) = ($1, $2);
	    } else {
		$species = $node_summary;
	    }
	    die "#4a error parsing species tree\n" unless defined $species;

            push @children_of_this_node, $species;
            $self->{parent_of_node}->{$species} = $this_node;
            $self->{branch_length_of_node}->{$species} = $branch_length if defined $branch_length;
	        $self->{terminal_nodes}->{$species} = 1;
            unless ($next_item) { die "#3 error parsing species tree\n"; }
            unshift @$inputarray, $next_item;
        
        }
        else { ## something else ? illegal.

            die "#4 error parsing species tree\n";
        }

        ## ok. one subnode has been parsed. next symbol must be either a comma
        ## (in which case we loop to the next subnode) or a closing bracket
        ## (in which case we are done).

        unless ($next_item = shift @$inputarray) { die "#5 error parsing species tree\n"; }

        next if ($next_item eq ",");
        next if ($next_item eq ":");

        if ($next_item eq ")") {

            ## code to get the branch_length here ...

            my $branch_summary = "";
            while ((@$inputarray) and (($next_item = shift @$inputarray) =~ /[\d\.\:]+/)) {
                $branch_summary .= $next_item;
            }
	    my ($bootstrap_support, $branch_length) = $branch_summary =~ /\A(\d*)\:([\d\.]+\z)/;
            $self->{branch_length_of_node}->{$this_node} = $branch_length if defined $branch_length;
	    $self->{bootstrap_support_of_node}->{$this_node} = $bootstrap_support if defined $bootstrap_support;
            unshift @$inputarray, $next_item;
            last;
        }

        ## anything else ? illegal.

        die "#6 error parsing species tree\n";
    }
    
    foreach my $child (@children_of_this_node) {
        $self->{children_of_node}->{$this_node}{$child} = 1;
            }
    
    return ($this_node);
}


####################################################################################################
## subroutine: connect_terminal_children_to_node ()
##
####################################################################################################

sub connect_terminal_children_to_node {

    my ($self, $node) = @_;

    if (($node =~ /\A[\-\d]+\z/) and ($node < 0)) {
	
	foreach my $child (keys %{$self->{children_of_node}->{$node}}) {
	    
	    if (($child =~ /\A[\-\d]+\z/) and ($child < 0)) {
		
		$self->connect_terminal_children_to_node ($child);
		foreach my $subchild (keys %{$self->{terminal_children_of_node}->{$child}}) {
		    $self->{terminal_children_of_node}->{$node}{$subchild} = 1;
		}
		
	    } else {
		
		$self->{terminal_children_of_node}->{$node}{$child} = 1;
		$self->{terminal_children_of_node}->{$child}{$child} = 1;
	    }
	}

    } else {       ## if the current node is a terminal node, attach it to itself.

	$self->{terminal_children_of_node}->{$node}{$node} = 1;
	$self->{terminal_children_of_node}->{$node}{$node} = 1;
    }
}

####################################################################################################
## subroutine: assign_cumulative_branchlength_to_node ()
##
####################################################################################################

sub assign_cumulative_branchlength_to_node {

    my ($self, $node, $previous_cumulative_length) = @_;

    my $branch_length_this_node = 0;
    $branch_length_this_node = $self->{branch_length_of_node}->{$node} if exists $self->{branch_length_of_node}->{$node};

    my $new_cumulative_length = $branch_length_this_node + $previous_cumulative_length;
    
    $self->{cumulative_branch_length_of_node}{$node} = $new_cumulative_length;

    if (($node =~ /\A[\-\d]+\z/) and ($node < 0)) {
	
	foreach my $child (keys %{$self->{children_of_node}->{$node}}) {

	    $self->assign_cumulative_branchlength_to_node ($child, $new_cumulative_length);
	}
    } 
}

####################################################################################################
## subroutine: compute_average_cumulative_branch_length ()
##
####################################################################################################
	    
sub compute_average_cumulative_branch_length {

    my ($self) = @_;

    my $valid_node_counter = 0;
    my $branch_length_sum = 0;

    foreach my $node (keys %{$self->{cumulative_branch_length_of_node}}) {
	
	next if (($node =~ /\A[\-\d]+\z/) and ($node < 0));       ## internal nodes not considered here.
	
	$valid_node_counter += 1;
	$branch_length_sum += $self->{cumulative_branch_length_of_node}->{$node};
    }

    $self->{average_cumulative_branch_length} = $branch_length_sum / $valid_node_counter;
}

####################################################################################################
## subroutine: connect_terminal_children_to_node ()
##
####################################################################################################

sub print_as_verbose_text {

    my ($self) = @_;

    foreach my $node (sort {$b <=> $a} keys %{$self->{children_of_node}}) {
    
	my $branch_length = "undefined";
        $branch_length = $self->{branch_length_of_node}->{$node} if exists $self->{branch_length_of_node}->{$node};
	my @children = sort {$a <=> $b} keys %{$self->{children_of_node}->{$node}};
	my $nr_children = scalar @children;
	my @terminal_children = sort {$a <=> $b} keys %{$self->{terminal_children_of_node}->{$node}};
	my $nr_terminal_children = scalar @terminal_children;
	
	print "node $node has branch-length $branch_length, $nr_children children [";
	print join ".", @children;
	print "], and $nr_terminal_children terminal children: [";
	print join ".", @terminal_children;
	print "]\n";
    }

    foreach my $node (sort {$a <=> $b} keys %{$self->{parent_of_node}}) {
	next unless $node > 0;
	my $branch_length = "undefined";
	$branch_length = $self->{branch_length_of_node}->{$node} if exists $self->{branch_length_of_node}->{$node};
	print "node $node has branch-length $branch_length\n";
    }
}

#######################################################################################################
## subroutine: get_species_to_species_tree_distance ()
##
## given two species, the cumulative branch-length needed to traverse from one to the other is 
## computed and provided by this routine. The routine will return -1 in case of errors.
#######################################################################################################

sub get_species_to_species_tree_distance {

    my ($self, $species1, $species2) = @_;

    return -1 unless defined $species1;
    return -1 unless defined $species2;

    return 0 if $species1 eq $species2;

    next unless exists $self->{parent_of_node}{$species1};
    next unless exists $self->{parent_of_node}{$species2};

    my %parent_traversion_count = ();

    my $node1 = $species1; $parent_traversion_count{$node1} = 1;
    while ((defined $node1) and (not ($node1 eq '-1'))) { $node1 = $self->{parent_of_node}{$node1}; $parent_traversion_count{$node1} += 1; }

    my $node2 = $species2; $parent_traversion_count{$node2} = 1;
    while ((defined $node2) and (not ($node2 eq '-1'))) { $node2 = $self->{parent_of_node}{$node2}; $parent_traversion_count{$node2} += 1; }

    my $cumulative_tree_distance = 0;

    $node1 = $species1; 
    while ($parent_traversion_count{$node1} < 2) { 
	$cumulative_tree_distance += $self->{branch_length_of_node}{$node1}; 
	$node1 = $self->{parent_of_node}{$node1};
    }

    $node2 = $species2; 
    while ($parent_traversion_count{$node2} < 2) { 
	$cumulative_tree_distance += $self->{branch_length_of_node}{$node2}; 
	$node2 = $self->{parent_of_node}{$node2};
    }

    return $cumulative_tree_distance;
}

#######################################################################################################
## subroutine: optimize_species_display_order ()
##
## this makes the tree optically nicer by swapping subtrees such that the internal nodes are maximally 
## 'staggered'. This only affects the displaying of the tree, not its topology or biological meaning. 
## If necessary, manual intervention into the display order is also possible, in the code below 
## (note: this will fuck up if the underlying tree is changed !!)
#######################################################################################################

sub optimize_species_display_order {

    my ($self) = @_;

    $self->{optimized_display_order_counter} = 1;

    $self->optimize_node_recursively (-1);
}

####################################################################################################
## subroutine: optimize_node_recursive ()
##
####################################################################################################

sub optimize_node_recursively {

    my ($self, $node) = @_;

    unless (($node =~ /\A[\-\d]+\z/) and ($node < 0)) {

	$self->{species_display_order}->{$node} = $self->{optimized_display_order_counter};
	$self->{optimized_display_order_counter} += 1;
	return;
    }

    ## ok, this node is an internal one - sort its children by 'terminal_size', and recurse down into them.

    my @children_this_node = keys %{$self->{children_of_node}->{$node}};
    my %nr_terminal_children_per_child = ();

    foreach my $child (@children_this_node) {
	my $nr_children = scalar keys %{$self->{terminal_children_of_node}->{$child}};
	$nr_terminal_children_per_child{$child} = $nr_children;
    }

    foreach my $child (sort { my $sort_value_a = $nr_terminal_children_per_child{$a};
			      my $sort_value_b = $nr_terminal_children_per_child{$b};
			      if ($sort_value_a == $sort_value_b) { $sort_value_a = $a; $sort_value_b = $b; }
			      return $sort_value_a <=> $sort_value_b;
			  } @children_this_node) {

	$self->optimize_node_recursively ($child);
    }
}

####################################################################################################
## subroutine: find_last_common_ancestor ()
##
####################################################################################################

sub find_last_common_ancestor {
    
    my ($self, $node1, $node2) = @_;
    return $node1 if $node1 == $node2;
    return $node1 if $self->is_in_parental_line ($node1, $node2);
    return $node2 if $self->is_in_parental_line ($node2, $node1);
    while (1) {
	return -1 unless exists $self->{parent_of_node}{$node1};
	$node1 = $self->{parent_of_node}{$node1};
	return $node1 if $self->is_in_parental_line ($node1, $node2);
	return -1 if $node1 == -1;
    }
}

####################################################################################################
## subroutine: is_in_parental_line ()
##
####################################################################################################

sub is_in_parental_line {

    my ($self, $putative_parent_node, $query_node) = @_;
    return 1 if $putative_parent_node == $query_node;
    while (1) {
	$query_node = $self->{parent_of_node}{$query_node};
	return 1 if $putative_parent_node == $query_node;
	last if $query_node == -1;
    }
    return 0;
}

    
1;
