#!/usr/bin/perl -w

use strict;
use Scalar::Util qw(looks_like_number);

############################################################
#
#       Tianyin Zhou 7/23/2014
#       For a given input file with extension '.txt.s', the
#       following feature combinations will be generated:
#       #Feature combination    -->     #File label
#       1mer  --> 1mer
#       1mer+2mer --> 1a2mer
#       1mer+2mer+3mer -->1a2a3mer
#       1mer+shape --> 1mer4shape
#       shape --> 4shape
#       1mer+MGW --> 1merM
#       1mer+Roll --> 1merR
#       1mer+HelT --> 1merT
#       1mer+ProT --> 1merP
#       1mer+MGW+Roll --> 1merMR
#       1mer+Roll+ProT --> 1merRP
#       1mer+2mer+shape --> 1a2mer4shape
# 
#       For each feature combinations the following three files 
#       will be generated: 
#       1. '.libsvm' -- feature combinations in libsvm format
#       2. '.matlab' -- feature combinations in matlab format
#       (The first column of the matrix is the response variable)
#       3. '.label'  -- The labels of all the features.
#
#############################################################


######################## User Options ########################################################

# Whether to include bilinear shape features.
my $use_bilinear = 0;
# If this is set to be 1, then the DNA shape features of terminal regions are predicted
# through averaging out all possible 2-mers on each side. For probes with very short 
# flanking regions such as those derived from SELEX-seq, setting this to be 1 will improve
# the performance of shape models.
my $encode_terminal = 0;

##############################################################################################

my $filename = $ARGV[0];

if (substr($filename, length($filename)-6, 6) ne '.txt.s') {
  die "Input file needs to have extention '.txt.s'\n";
}

# Suppose the filename has extension ".txt.s"
my $identifier = substr($filename, 0, length($filename)-6);

my $current_dir = substr($0, 0, length($0)-19); 
my $pred_bin = $current_dir."/DNAshape_v2.6/prediction";
#my $pred_bin = "/home/rcf-47/tianyinz/cmb/v2.5_fake1/prediction";



my $svm_add_intercept = 0;

# Numbers for normlaizing shape features
my $MGW_max_first = 6.2;
my $MGW_min_first = 2.85;
my $MGW_span_first = $MGW_max_first - $MGW_min_first;
my $MGW_max_second = 38.44;
my $MGW_min_second = 8.1225;
my $MGW_span_second = $MGW_max_second - $MGW_min_second;

my $ProT_max_first = -0.03;
my $ProT_min_first = -16.51;
my $ProT_span_first = $ProT_max_first - $ProT_min_first;
my $ProT_max_second = 272.5801;
my $ProT_min_second = 0.0396;
my $ProT_span_second = $ProT_max_second - $ProT_min_second;

my $Roll_max_first = 8.64;
my $Roll_min_first = -8.57;
my $Roll_span_first = $Roll_max_first - $Roll_min_first;
my $Roll_max_second = 43.3642;
my $Roll_min_second = -43.6135;
my $Roll_span_second = $Roll_max_second - $Roll_min_second;

my $HelT_max_first = 38.05;
my $HelT_min_first = 30.94;
my $HelT_span_first = $HelT_max_first - $HelT_min_first;
my $HelT_max_second = 1439.8188;
my $HelT_min_second = 1050.704;
my $HelT_span_second = $HelT_max_second - $HelT_min_second;

# The input file
my $file =$identifier.".txt.s";
my $seqfile = $identifier.".seq";
#my $outfile1 = $identifier."_1mer";
#my $outfile2 = $identifier."_1a2mer";
#my $outfile3 = $identifier."_1merMR";
#my $outfile4 = $identifier."_1merRP";
#my $outfile5 = $identifier."_1a2mer4shape";
#my $outfile6 = $identifier."_1a2a3mer";
#my $outfile7 = $identifier."_1merM";
#my $outfile8 = $identifier."_1merR";
#my $outfile9 = $identifier."_1merP";
#my $outfile10 =$identifier."_1merT";
#my $outfile11 = $identifier."_1mer4shape";
my $outfile12 = $identifier."_4shape";
#my $outfile13 = $identifier."_2mer";
#my $outfile14 = $identifier."_MGW";

my $minor_file = $seqfile.".MGW";
my $propel_file = $seqfile.".ProT";
my $roll_file = $seqfile.".Roll";
my $twist_file = $seqfile.".HelT";

open INF," <$file" or die "Cannot find $file\n";
open OUTF," >$seqfile" or die "Cannot write to $seqfile\n";

my $count = 0;
my @seq_list = ();
my @response_list = ();

my @four_nuc = qw/A T G C/;

while (<INF>){
  chomp;
  my @content = split /\s+/,$_;
  if (scalar(@content) == 2){
    $count++;

    if ($encode_terminal) {
      for (my $i = 0; $i < 4; $i++) {
        for (my $j = 0; $j <4; $j++) {
          print OUTF ">$count\n";
          print OUTF $four_nuc[$i].$four_nuc[$j].$content[1].$four_nuc[$i].$four_nuc[$j]."\n";
        }
      }
    } else {
      print OUTF ">$count\n";
      print OUTF $content[1]."\n";
    }

# Store the response variable. Assume that it is in the same column of the file. Otherwise change the index.
# You can also take the log of $content[0] if you want to transform the response variable. 
    push @response_list,$content[0];
# Store the sequence. Assume that it is in the second column of the file. Otherwise change the index. 
    push @seq_list,$content[1];  
  }
}
my $sequence_len=length($seq_list[0]);
close OUTF;
print $count."\n";

#predict all the parameters
system($pred_bin,"-i",$seqfile,"-width",1500);

#Read Minor groove widths
my @minor = ();
&load_shape_features($minor_file, \@minor);
print "Minor: ".scalar(@minor)."\n";

#Read Propeller twist
my @propel = ();
&load_shape_features($propel_file, \@propel);
print "Propel: ".scalar(@propel)."\n";

#Read Roll
my @roll = ();
&load_shape_features($roll_file, \@roll);
print "Roll: ".scalar(@roll)."\n";

#Read Helical Twist
my @twist = ();
&load_shape_features($twist_file, \@twist);
print "Twist: ".scalar(@twist)."\n";

unlink $seqfile;
unlink $minor_file;
unlink $propel_file;
unlink $roll_file;
unlink $twist_file;

if ((scalar(@propel)!=$count) or (scalar(@minor)!=$count)){
  die scalar(@propel)." ".scalar(@minor)." ".$count."\n";
}

my %complement = ('A'=>'T','T'=>'A','G'=>'C','C'=>'G');

#my @column_label1 = ();
#my @column_label2 = ();
#my @column_label3 = ();
#my @column_label4 = ();
#my @column_label5 = ();
#my @column_label6 = ();
#my @column_label7 = ();
#my @column_label8 = ();
#my @column_label9 = ();
#my @column_label10 = ();
#my @column_label11 = ();
my @column_label12 = ();
#my @column_label13 = ();
#my @column_label14 = ();

#my @output_matrix1 = ();
#my @output_matrix2 = ();
#my @output_matrix3 = ();
#my @output_matrix4 = ();
#my @output_matrix5 = ();
#my @output_matrix6 = ();
#my @output_matrix7 = ();
#my @output_matrix8 = ();
#my @output_matrix9 = ();
#my @output_matrix10 = ();
#my @output_matrix11 = ();
my @output_matrix12 = ();
#my @output_matrix13 = ();
#my @output_matrix14 = ();

#for (my $i=0; $i<$count; $i++){
for (my $i=0; $i<$count; $i++){
=pod
  my @output_array1 = ();   
  my @output_array2 = ();  
  my @output_array3 = (); 
  my @output_array4 = ();
  my @output_array5 = ();
  my @output_array6 = ();
  my @output_array7 = ();
  my @output_array8 = ();
  my @output_array9 = ();
  my @output_array10 = ();
  my @output_array11 = ();
=cut
  my @output_array12 = ();
#  my @output_array13 = ();
#  my @output_array14 = ();

=pod 
  push @output_array1, $response_list[$i];
  push @output_array2, $response_list[$i];
  push @output_array3, $response_list[$i];
  push @output_array4, $response_list[$i];
  push @output_array5, $response_list[$i];
  push @output_array6, $response_list[$i];
  push @output_array7, $response_list[$i];
  push @output_array8, $response_list[$i];
  push @output_array9, $response_list[$i];
  push @output_array10, $response_list[$i];
  push @output_array11, $response_list[$i];
=cut
  push @output_array12, $response_list[$i];
#  push @output_array13, $response_list[$i];
#  push @output_array14, $response_list[$i];

#sequence encoding
  for (my $j=0; $j<$sequence_len; $j++){
    my %letter_count = ('A'=>0,'T'=>0,'G'=>0,'C'=>0);
    my $nuc = substr($seq_list[$i],$j,1);
    $letter_count{$nuc}++;

    my @monomer_output = ();
    push @monomer_output,$letter_count{'A'};
    push @monomer_output,$letter_count{'T'};
    push @monomer_output,$letter_count{'G'};
    push @monomer_output,$letter_count{'C'};	

    my @monomer_label = ();
    push @monomer_label, 'A_at_pos_'.($j+1);
    push @monomer_label, 'T_at_pos_'.($j+1);
    push @monomer_label, 'G_at_pos_'.($j+1);
    push @monomer_label, 'C_at_pos_'.($j+1);
=pod
    push @output_array1,@monomer_output;
    push @output_array2,@monomer_output;
    push @output_array3,@monomer_output;
    push @output_array4,@monomer_output;
    push @output_array5,@monomer_output;
    push @output_array6,@monomer_output;
    push @output_array7,@monomer_output;
    push @output_array8,@monomer_output;
    push @output_array9,@monomer_output;
    push @output_array10,@monomer_output;
    push @output_array11,@monomer_output;
=cut
    if ($i==0) {
=pod
      push @column_label1,@monomer_label;
      push @column_label2,@monomer_label;
      push @column_label3,@monomer_label;
      push @column_label4,@monomer_label;
      push @column_label5,@monomer_label;
      push @column_label6,@monomer_label;
      push @column_label7,@monomer_label;
      push @column_label8,@monomer_label;
      push @column_label9,@monomer_label;
      push @column_label10,@monomer_label;
      push @column_label11,@monomer_label;
=cut
    }
  }

#2-mer encoding
  for (my $j=0; $j<($sequence_len-1); $j++){
    my @dimer_output = ();
    my @dimer_label = ();
    my @four_nuc = qw/A T G C/;
    my $dimer_seen = substr($seq_list[$i],$j,2);
    for (my $k1 = 0; $k1 < 4; $k1++) {
      for (my $k2 = 0; $k2 < 4; $k2++) {
        my $dimer = $four_nuc[$k1].$four_nuc[$k2];
        if ($dimer eq $dimer_seen) {
          push @dimer_output, 1;
        } else{
          push @dimer_output, 0;
        }
        push @dimer_label, $dimer."_at_pos_".($j+1);
      }
    }

    #push @output_array2,@dimer_output;
    #push @output_array5,@dimer_output;
    #push @output_array6,@dimer_output;
    #push @output_array13,@dimer_output;

    if ($i==0) {
      #push @column_label2,@dimer_label;
      #push @column_label5,@dimer_label;
      #push @column_label6,@dimer_label;
      #push @column_label13,@dimer_label;
    }
  }
  
#3-mer encoding
  for (my $j = 0; $j < ($sequence_len-2); $j++) {
    my @trimer_output = ();
    my @trimer_label = ();
    my @four_nuc = qw/A T G C/;
    my $trimer_seen = substr($seq_list[$i],$j,3);
    for (my $k1 = 0; $k1 < 4; $k1++) {
      for (my $k2 = 0; $k2 < 4; $k2++) {
        for (my $k3 = 0; $k3 < 4; $k3++) {
          my $trimer = $four_nuc[$k1].$four_nuc[$k2].$four_nuc[$k3];
          if ($trimer eq $trimer_seen) {
            push @trimer_output,1;
          } else {
            push @trimer_output,0;
          }
          push @trimer_label, $trimer."_at_pos_".($j+1);
        }
      }
    }

    #push @output_array6,@trimer_output;

    if ($i==0) {
      #push @column_label6, @trimer_label;
    }
  }

  my $skip_bp = 2;
  my $skip_step = 1;
  if ($encode_terminal) {
    $skip_bp = 0;
    $skip_step = 0;
  }
#minor encoding
  for (my $j = $skip_bp; $j < ($sequence_len-$skip_bp); $j++) {
    #my $minor_value = $minor[$i][$j];
    my $minor_value = 0;
    if (looks_like_number($minor[$i][$j])) {
      $minor_value = ($minor[$i][$j] - $MGW_min_first) / $MGW_span_first;
    } 

    my $minor_label = 'MGW_at_pos_'.($j+1);
    #push @output_array3,$minor_value;
    #push @output_array5,$minor_value;
    #push @output_array7,$minor_value;
    #push @output_array11,$minor_value;
    push @output_array12,$minor_value;
    #push @output_array14,$minor_value;
    if ($i==0) {
      #push @column_label3,$minor_label;
      #push @column_label5,$minor_label;
      #push @column_label7,$minor_label;
      #push @column_label11,$minor_label;
      push @column_label12,$minor_label;
      #push @column_label14,$minor_label;
    }
  }
#minor bilinear terms
  if ($use_bilinear) {
    for (my $j = $skip_bp; $j < ($sequence_len-$skip_bp-1); $j++) {
      my $minor_bilinear_term = 0;
      if (looks_like_number($minor[$i][$j]) and looks_like_number($minor[$i][$j+1])) {
        $minor_bilinear_term = $minor[$i][$j] * $minor[$i][$j+1];
      }
      $minor_bilinear_term = ($minor_bilinear_term - $MGW_min_second) / $MGW_span_second;
      my $minor_bilinear_label = 'product_of_MGW_between_pos_'.($j+1).'_and_'.($j+2);
      #push @output_array3, $minor_bilinear_term;
      #push @output_array5, $minor_bilinear_term;
      #push @output_array7, $minor_bilinear_term;
      #push @output_array11, $minor_bilinear_term;
      push @output_array12, $minor_bilinear_term;
      #push @output_array14, $minor_bilinear_term;

      if ($i==0) {
        #push @column_label3, $minor_bilinear_label;
        #push @column_label5, $minor_bilinear_label;
        #push @column_label7, $minor_bilinear_label;
        #push @column_label11, $minor_bilinear_label;
        push @column_label12, $minor_bilinear_label;
        #push @column_label14, $minor_bilinear_label;
      }
    }
  }

#roll encoding
  for (my $j = $skip_step; $j < ($sequence_len-$skip_step-1); $j++) {
    my $roll_value = 0;
    if (looks_like_number($roll[$i][$j])) {
      $roll_value = ($roll[$i][$j] - $Roll_min_first) / $Roll_span_first;
    }
    my $roll_label = 'Roll_at_step_'.($j+1);
    #push @output_array3,$roll_value;
    #push @output_array4,$roll_value;
    #push @output_array5,$roll_value;
    #push @output_array8,$roll_value;
    #push @output_array11,$roll_value;
    push @output_array12,$roll_value;
    if ($i==0) {
      #push @column_label3,$roll_label;
      #push @column_label4,$roll_label;
      #push @column_label5,$roll_label;
      #push @column_label8,$roll_label;
      #push @column_label11,$roll_label;
      push @column_label12,$roll_label;
    }
  }
#roll bilinear terms
  if ($use_bilinear) {
    for (my $j = $skip_step; $j < ($sequence_len-$skip_step-2); $j++) {
      my $roll_bilinear_term = 0;
      if (looks_like_number($roll[$i][$j]) and looks_like_number($roll[$i][$j+1])) {
        $roll_bilinear_term = $roll[$i][$j] * $roll[$i][$j+1];
      }
      $roll_bilinear_term = ($roll_bilinear_term - $Roll_min_second) / $Roll_span_second;
      my $roll_bilinear_label = 'product_of_Roll_between_step_'.($j+1).'_and_'.($j+2);
      #push @output_array3,$roll_bilinear_term;
      #push @output_array4,$roll_bilinear_term;
      #push @output_array5,$roll_bilinear_term;
      #push @output_array8,$roll_bilinear_term;
      #push @output_array11,$roll_bilinear_term;
      push @output_array12,$roll_bilinear_term;
      if ($i==0) {
        #push @column_label3,$roll_bilinear_label;
        #push @column_label4,$roll_bilinear_label;
        #push @column_label5,$roll_bilinear_label;
        #push @column_label8,$roll_bilinear_label;
        #push @column_label11,$roll_bilinear_label;
        push @column_label12,$roll_bilinear_label;
      }
    }
  }

#twist encoding
  for (my $j = $skip_step; $j < ($sequence_len-$skip_step-1); $j++) {
    my $twist_value = 0;
    if (looks_like_number($twist[$i][$j])) {
      $twist_value = ($twist[$i][$j] - $HelT_min_first) / $HelT_span_first;
    }
    my $twist_label = 'twist_at_step_'.($j+1);
    #push @output_array5,$twist_value;
    #push @output_array10,$twist_value;
    #push @output_array11,$twist_value;
    push @output_array12,$twist_value;
    if ($i==0) {
      #push @column_label5,$twist_label;
      #push @column_label10,$twist_label;
      #push @column_label11,$twist_label;
      push @column_label12,$twist_label;
    }
  }
#twist bilinear terms
  if ($use_bilinear) {
    for (my $j = $skip_step; $j < ($sequence_len-$skip_step-2); $j++) {
      my $twist_bilinear_term = 0;
      if (looks_like_number($twist[$i][$j]) and looks_like_number($twist[$i][$j+1])) {
        $twist_bilinear_term = $twist[$i][$j] * $twist[$i][$j+1];
      }
      $twist_bilinear_term = ($twist_bilinear_term - $HelT_min_second) / $HelT_span_second;
      my $twist_bilinear_label = 'product_of_twist_step_'.($j+1).'_and_'.($j+2);
      #push @output_array5,$twist_bilinear_term;
      #push @output_array10,$twist_bilinear_term;
      #push @output_array11,$twist_bilinear_term;
      push @output_array12,$twist_bilinear_term;
      if ($i==0) {
        #push @column_label5,$twist_bilinear_label;
        #push @column_label10,$twist_bilinear_label;
        #push @column_label11,$twist_bilinear_label;
        push @column_label12,$twist_bilinear_label;
      }
    }
  }

#propeller encoding
  for (my $j = $skip_bp; $j < ($sequence_len-$skip_bp); $j++) {    
    my $propel_value = 0;
    if (looks_like_number($propel[$i][$j])) {
      $propel_value = ($propel[$i][$j] - $ProT_min_first) / $ProT_span_first;
    }
    my $propel_label = 'ProT_at_pos_'.($j+1);
    #push @output_array4,$propel_value;
    #push @output_array5,$propel_value;
    #push @output_array9,$propel_value;
    #push @output_array11,$propel_value;
    push @output_array12,$propel_value;
    if ($i==0) {
      #push @column_label4,$propel_label;
      #push @column_label5,$propel_label;
      #push @column_label9,$propel_label;
      #push @column_label11,$propel_label;
      push @column_label12,$propel_label;
    }
  }
#propeller bilinear terms
  if ($use_bilinear) {
    for (my $j = $skip_bp; $j < ($sequence_len-$skip_bp-1); $j++) {    
      my $propel_bilinear_term = 0;
      if (looks_like_number($propel[$i][$j]) and looks_like_number($propel[$i][$j+1])) {
        $propel_bilinear_term = $propel[$i][$j] * $propel[$i][$j+1];
      }
      $propel_bilinear_term = ($propel_bilinear_term - $ProT_min_second) / $ProT_span_second;
      my $propel_bilinear_label = 'product_of_ProT_between_pos_'.($j+1).'_and_'.($j+2);
      #push @output_array4,$propel_bilinear_term;
      #push @output_array5,$propel_bilinear_term;
      #push @output_array9,$propel_bilinear_term;
      #push @output_array11,$propel_bilinear_term;
      push @output_array12,$propel_bilinear_term;
      if ($i==0) {
        #push @column_label4,$propel_bilinear_label;
        #push @column_label5,$propel_bilinear_label;
        #push @column_label9,$propel_bilinear_label;
        #push @column_label11,$propel_bilinear_label;
        push @column_label12,$propel_bilinear_label;
      }
    }
  }

  #push @output_matrix1, [@output_array1];
  #push @output_matrix2, [@output_array2];
  #push @output_matrix3, [@output_array3];
  #push @output_matrix4, [@output_array4];
  #push @output_matrix5, [@output_array5];
  #push @output_matrix6, [@output_array6];
  #push @output_matrix7, [@output_array7];
  #push @output_matrix8, [@output_array8];
  #push @output_matrix9, [@output_array9];
  #push @output_matrix10, [@output_array10];
  #push @output_matrix11, [@output_array11];
  push @output_matrix12, [@output_array12];
  #push @output_matrix13, [@output_array13];
  #push @output_matrix14, [@output_array14];
}

#&scaling_data(\@output_matrix1);
#&scaling_data(\@output_matrix2);
#&scaling_data(\@output_matrix3);
#&scaling_data(\@output_matrix4);
#&scaling_data(\@output_matrix5);
#&scaling_data(\@output_matrix6);
#&scaling_data(\@output_matrix7);
#&scaling_data(\@output_matrix8);
#&scaling_data(\@output_matrix9);
#&scaling_data(\@output_matrix10);
#&scaling_data(\@output_matrix11);
#&scaling_data(\@output_matrix12);
#&scaling_data(\@output_matrix13);
#&scaling_data(\@output_matrix14);

#&write_to_file($outfile1, \@output_matrix1, \@column_label1);
#&write_to_file($outfile2, \@output_matrix2, \@column_label2);
#&write_to_file($outfile3, \@output_matrix3, \@column_label3);
#&write_to_file($outfile4, \@output_matrix4, \@column_label4);
#&write_to_file($outfile5, \@output_matrix5, \@column_label5);
#&write_to_file($outfile6, \@output_matrix6, \@column_label6);
#&write_to_file($outfile7, \@output_matrix7, \@column_label7);
#&write_to_file($outfile8, \@output_matrix8, \@column_label8);
#&write_to_file($outfile9, \@output_matrix9, \@column_label9);
#&write_to_file($outfile10, \@output_matrix10, \@column_label10);
#&write_to_file($outfile11, \@output_matrix11, \@column_label11);
&write_to_file($outfile12, \@output_matrix12, \@column_label12);
#&write_to_file($outfile13, \@output_matrix13, \@column_label13);
#&write_to_file($outfile14, \@output_matrix14, \@column_label14);

# Arguments:
# output_file_base  matrx label_array
sub write_to_file {
  my $output_file_base = $_[0];
  my @matrix_to_output = @{$_[1]};
  my @column_label = @{$_[2]};
  
  my $matlab_outputfile = $output_file_base.".matlab";
  #my $libsvm_outputfile = $output_file_base.".libsvm";
  #my $label_outputfile = $output_file_base.".label";
  
  open MATLAB_OUT, " >$matlab_outputfile";
  #open LIBSVM_OUT, " >$libsvm_outputfile";
  #open LABEL_OUT, " >$label_outputfile";

  my $matrix_width = scalar(@{$matrix_to_output[0]});
  my @constant_column = (0)x$matrix_width;  # If the whole column contains constant value.
  #for (my $j = 1; $j < $matrix_width; $j++) {
  #  for (my $i = 1; $i < scalar(@matrix_to_output); $i++) {
  #    if ($matrix_to_output[$i][$j] != $matrix_to_output[0][$j]) {
  #      $constant_column[$j] = 0;
  #      last;
  #    }
  #  }
  #}

  for (my $i = 0; $i < scalar(@matrix_to_output); $i++) {
    print MATLAB_OUT $matrix_to_output[$i][0];
    #print LIBSVM_OUT $matrix_to_output[$i][0];
    my $column_count = 0;
    
    for (my $j = 1; $j < $matrix_width; $j++) {
      if (!$constant_column[$j]) {
        # If a column contains all constant value, this column will not be output.
        $column_count ++;
        # Check if the value is integer or not
        if ($matrix_to_output[$i][$j] =~ /\D/) {
          #printf MATLAB_OUT " %.3f", $matrix_to_output[$i][$j];
          #printf LIBSVM_OUT " %d:%.3f", $j, $matrix_to_output[$i][$j];
          printf MATLAB_OUT " %.4f", $matrix_to_output[$i][$j];
          #printf LIBSVM_OUT " %d:%.4f", $column_count, $matrix_to_output[$i][$j];
        } else{
          printf MATLAB_OUT " %d", $matrix_to_output[$i][$j];
          #printf LIBSVM_OUT " %d:%d", $column_count, $matrix_to_output[$i][$j];
        }
        # print column label
        if ($i==0) {
          #print LABEL_OUT $column_count.":".$column_label[$j-1]."\n";
        }
      }
    }
    if ($svm_add_intercept) {
      $column_count++;
      #print LIBSVM_OUT " ".$column_count.":1";
    }
    print MATLAB_OUT "\n";
    #print LIBSVM_OUT "\n";
  }
  close MATLAB_OUT;
  #close LIBSVM_OUT;
  #close LABEL_OUT;
}

# input:  output_matrix
# The first column is response value, which does not need to be scaled.
# Linear scaling on the attributes.
sub scaling_data {
  my @matrix = @{$_[0]};
  my $width = scalar(@{$matrix[0]});
  my @column_min = (0)x$width; 
  my @column_max = (0)x$width;
  my @column_kmer = (1)x$width;  # only apply scaling on shape features only.
  #my @column_kmer = (0)x$width;   # apply scaling on all the features
  for (my $i = 0; $i < scalar(@matrix); $i++) {
    for (my $j = 0; $j < $width; $j++) {	
      if ($i == 0) {
        $column_max[$j] = $matrix[$i][$j];
        $column_min[$j] = $matrix[$i][$j];
      } else {
        if ($matrix[$i][$j] > $column_max[$j]) {
          $column_max[$j] = $matrix[$i][$j];
        }
        if ($matrix[$i][$j] < $column_min[$j]) {
          $column_min[$j] = $matrix[$i][$j];
        }
    }
    }
  }

  # check if the column only contains value 0 or 1
  for (my $j = 0; $j < $width; $j++) {
    for (my $i = 0; $i < scalar(@matrix); $i++) {
      if ($matrix[$i][$j] != 0 && $matrix[$i][$j] != 1) {
        $column_kmer[$j] = 0;
        last;
      }
    }
  }
  # scaling data
  # The first column is the response value. Does not need to be normalized.
  #for (my $j = 1; $j < $width; $j++) {
  #  print $column_max[$j]." ".$column_min[$j]."\n";
  #}
=pod
  for (my $i = 0; $i < scalar(@matrix); $i++) {
    for (my $j = 1; $j < $width; $j++) {
      if ($column_max[$j] != $column_min[$j] && $column_kmer[$j] == 0) {
        $matrix[$i][$j] = ($matrix[$i][$j] - $column_min[$j]) / ($column_max[$j] - $column_min[$j]);
      }
    }
  }
=cut
}

#Load shape features from the output of DNAshape program
#Args:
#  files_to_open
#  array to store the shape features
sub load_shape_features {
  my $file_to_open = $_[0];
  my @shape_array = ();

  open INF, "< $file_to_open" or die "Cannot open $file_to_open\n";
  while (<INF>) {
    chomp;
    my $line = $_;
    if (substr($line, 0, 1) ne ">") {
      my @entry = split /,/, $line;
      if ($encode_terminal) {
        for (my $i = 0; $i < 15; $i++) {
          my $dummy_line = <INF>;
          my $real_line = <INF>;
          chomp $real_line;
          my @entry_to_add = split /,/, $real_line;
          for (my $j = 2; $j < (scalar(@entry) - 2); $j++) {
            $entry[$j] += $entry_to_add[$j];
          }
        }
        for (my $j = 2; $j < (scalar(@entry) - 2); $j++) {
          $entry[$j] = $entry[$j]/16;
        }
        push @shape_array, [@entry[2..(scalar(@entry)-3)]];
      } else {
        push @shape_array, [@entry];
      }
    }
  }
  @{$_[1]} = @shape_array;
  #print "Size of shape array: ".scalar(@shape_array)."\n";
}
