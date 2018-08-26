#!/usr/bin/perl -w 

use strict;
use Cwd 'abs_path';
use List::Util qw(shuffle);   #This module includes a method to shuffle arrays
use Scalar::Util qw(looks_like_number);

####################################################################
#
#
#           Cross Validation
#           Divide the simulations in the simulation pool into
#           10 bins. Each time take out 1 bin, build the query
#           table based on the remaining simulations, and use
#           that table to predict minor groove width of the 
#           sequences in the simulation taken out
#           Tianyin Zhou    02/11/2013
#
#
####################################################################

my $debug = 0;
my $minor_tolerance = 0.05;
my $propel_tolerance = 0.25;
my $roll_tolerance = 0.25;
my $twist_tolerance = 0.1;

my $abspath = abs_path($0)."\n";
my $prog_dir = substr($abspath,0,index($abspath,"Script"));
#print $prog_dir."\n";

my $original_list = $prog_dir."SimData/nonartifacts.txt";
my $ave_g_width_bin = $prog_dir."Script/ave_g_width.pl";
my $pred_bin = $prog_dir."prediction";

my $log_file = $prog_dir."Script/log.txt";

open LOG," >$log_file";

open INF," <$original_list" ;
my @orig_list = (<INF>);
foreach (@orig_list){
    chomp;
}
close INF;
print scalar(@orig_list)."\n";
#leave one out cross validation

my @minor_MC_overall = ();
my @minor_HT_overall = ();
my @propel_MC_overall = ();
my @propel_HT_overall = ();
my @roll_MC_overall = ();
my @roll_HT_overall = ();
my @twist_MC_overall = ();
my @twist_HT_overall = ();

for (my $i=0; $i<scalar(@orig_list); $i++){
    my $tmp_lib_file = $prog_dir."/Script/tmp_lib.txt";
    open OUTF,"> $tmp_lib_file";
    #build the lib file with one record removed
    for (my $j=0; $j<scalar(@orig_list); $j++){
	if ($j!=$i){
	    print OUTF $orig_list[$j]."\n";
	}
    }
    close OUTF;    
    #print $orig_list[$i]."\n";
    my $filename_len = length($orig_list[$i]);
    my $basename = substr($orig_list[$i],0,$filename_len-6);
    $basename = $prog_dir."SimData/".substr($basename,1,$filename_len);    
    my $seq_file = $basename.".seq";
    my $cclis_file = $basename."_cc.lis";
    my $minor_file = $basename.".minor";
    open INF," <$seq_file" or die "Cannot open $seq_file \n";
    my $seq = <INF>;
    chomp $seq;
    close INF;
    my $test_file = $prog_dir."Script/tmp_test.txt";
    open OUTF," > $test_file";
    print OUTF ">1\n";
    print OUTF $seq."\n";
    close OUTF;
    #print $seq."\n";
    #extract the minor groove out 
    my $command = $ave_g_width_bin." ".$minor_file;
    my @command_output = `$command`;
    my $output_line;
    my @minor_MC = ();
    my @minor_MC_index = ();
    foreach $output_line(@command_output){
	chomp $output_line;
	my @line_content = split /\s+/,$output_line;
	if (scalar(@line_content) == 3){
	    push @minor_MC_index,$line_content[1];
	    push @minor_MC,$line_content[2];
	}
    }    
    #extract the Propeller from block D
    my @propel_MC = ();
    my $line;
    my @line_content;
    open INF,"< $cclis_file";
    while (<INF>){
	if (/\|D\|/){
	    last;
	}
    }
    for (my $j=0; $j<7; $j++){
	$line = <INF>;
    }
    while (<INF>){
	$line = $_;
	chomp $line;
	@line_content = split /\s+/,$line;
	if (scalar(@line_content) == 13){
	    push @propel_MC,$line_content[9];
	}else{
	    last;
	}
    }
    #extract roll, twist from block H
    my @roll_MC = ();
    my @twist_MC = ();
    while (<INF>){
	if (/\|H\|/){
	    last;
	}
    }
    for (my $j=0; $j<7; $j++){
	$line = <INF>;
    }
    while (<INF>){
	$line = $_;
	chomp $line;
	@line_content = split /\s+/,$line;
	if (scalar(@line_content) == 12){
	    push @roll_MC,$line_content[9];
	    push @twist_MC,$line_content[10];
	}else{
	    last;
	}
    }

    close INF;

    #run the prediction
    system($pred_bin,"-i",$test_file,"-lib",$tmp_lib_file,"-width",50);
    #read the predicted minor
    my @minor_HT = ();
    open INF,"< $test_file.minor";
    $line = <INF>;  #skip one line
    chomp($line = <INF>);
    close INF;
    @line_content = split /,/,$line;    
    foreach (@minor_MC_index){
	my $current_index = $_-1;
	if (looks_like_number($line_content[$current_index])){
	    push @minor_HT,$line_content[$current_index];
	}
	else{
	    @minor_MC = ();
	    @minor_HT = ();
	    last;
	}
    }    
    
    #read the predicted propel
    my @propel_HT = ();
    open INF,"< $test_file.propel";
    while (<INF>){
	chomp;
	if (/NA,/){
	    $line = $_;
	    @line_content = split /,/,$line;
	    push @propel_HT,@line_content;
	}
    }
    close INF;
    #read the predicted roll
    my @roll_HT = ();
    open INF,"< $test_file.roll";
    while (<INF>){
	chomp;
	if (/NA,/){
	    $line = $_;
	    @line_content = split /,/,$line;
	    push @roll_HT,@line_content;
	}
    }
    #read the predicted twist
    my @twist_HT = ();
    open INF,"< $test_file.twist";
    while (<INF>){
	chomp;
	if (/NA,/){
	    $line = $_;
	    @line_content = split /,/,$line;
	    push @twist_HT,@line_content;
	}
    }
    my @minor_HT_rank = ();
    my @minor_MC_rank = ();
    &rank_array(\@minor_HT,\@minor_HT_rank,$minor_tolerance);
    &rank_array(\@minor_MC,\@minor_MC_rank,$minor_tolerance);

    push @minor_HT_overall,@minor_HT_rank;
    push @minor_MC_overall,@minor_MC_rank;   

    #remove the first two and last two
    pop @propel_HT;
    pop @propel_HT;
    shift @propel_HT;
    shift @propel_HT;
    pop @propel_MC;
    pop @propel_MC;
    shift @propel_MC;
    shift @propel_MC;    
    my @propel_HT_rank = ();
    my @propel_MC_rank = ();
    &rank_array(\@propel_HT,\@propel_HT_rank,$propel_tolerance);
    &rank_array(\@propel_MC,\@propel_MC_rank,$propel_tolerance);
    push @propel_HT_overall,@propel_HT_rank;
    push @propel_MC_overall,@propel_MC_rank;

    #remove the first and the last
    pop @roll_HT;
    shift @roll_HT;
    pop @roll_MC;
    shift @roll_MC;
    my @roll_HT_rank = ();
    my @roll_MC_rank = ();
    &rank_array(\@roll_HT,\@roll_HT_rank,$roll_tolerance);
    &rank_array(\@roll_MC,\@roll_MC_rank,$roll_tolerance);
    push @roll_HT_overall,@roll_HT_rank;
    push @roll_MC_overall,@roll_MC_rank;

    #remove the first and the last
    pop @twist_HT;
    shift @twist_HT;
    pop @twist_MC;
    shift @twist_MC;
    my @twist_HT_rank = ();
    my @twist_MC_rank = ();
    &rank_array(\@twist_HT,\@twist_HT_rank,$twist_tolerance);
    &rank_array(\@twist_MC,\@twist_MC_rank,$twist_tolerance);
    push @twist_HT_overall,@twist_HT_rank;
    push @twist_MC_overall,@twist_MC_rank;    

    if ((scalar(@minor_HT)>0)){
	if ((&calc_pearson_correlation(\@minor_HT,\@minor_MC))<0.3){
	    print LOG $orig_list[$i]."\n";	
	    for (my $j=0; $j<scalar(@minor_HT);$j++){
		print LOG $minor_HT[$j]." ";
	    }
	    print LOG "\n";
	    for (my $j=0; $j<scalar(@minor_MC);$j++){
		print LOG $minor_MC[$j]." ";
	    }
	    print LOG "\n";
	}
    }else{
	print LOG "NA detected $orig_list[$i]\n";
    }

    if ($debug){
	&print_list(\@minor_HT_overall);
	&print_list(\@minor_MC_overall);
	&print_list(\@propel_HT_overall);
	&print_list(\@propel_MC_overall);
	&print_list(\@roll_HT_overall);
	&print_list(\@roll_MC_overall);
	&print_list(\@twist_HT_overall);
	&print_list(\@twist_MC_overall);
    }     
       
}

print "Minor_HT length: ".scalar(@minor_HT_overall)."\n";
print "Minor_MC length: ".scalar(@minor_MC_overall)."\n";
print "Minor spearman: ".&calc_pearson_correlation(\@minor_HT_overall,\@minor_MC_overall)."\n";
print "Propel_HT length: ".scalar(@propel_HT_overall)."\n";
print "Propel_MC length: ".scalar(@propel_MC_overall)."\n";
print "Propel spearman: ".&calc_pearson_correlation(\@propel_HT_overall,\@propel_MC_overall)."\n";
print "Roll length: ".scalar(@roll_HT_overall)."\n";
print "Roll length: ".scalar(@roll_MC_overall)."\n";
print "Roll spearman: ".&calc_pearson_correlation(\@roll_HT_overall,\@roll_MC_overall)."\n";
print "Twist_HT length: ".scalar(@twist_HT_overall)."\n";
print "Twist_MC length: ".scalar(@twist_MC_overall)."\n";
print "Twist spearman: ".&calc_pearson_correlation(\@twist_HT_overall,\@twist_MC_overall)."\n";

close LOG;


sub print_list{
    my @to_print = @{$_[0]};
    for (my $i=0; $i<scalar(@to_print);$i++){
	print $to_print[$i]." ";
    }
    print "\n";
}


sub calc_pearson_correlation{
    my @arrayX = @{$_[0]};
    my @arrayY = @{$_[1]};
    if (scalar(@arrayX)!=scalar(@arrayY)){
	die "Cannot calculate pearson correlation";
    }
    my ($aveX,$aveY) = (0,0);
    for (my $i=0;$i<scalar(@arrayX);$i++){
	$aveX += $arrayX[$i];
	$aveY += $arrayY[$i];
    }
    $aveX = $aveX/scalar(@arrayX);
    $aveY = $aveY/scalar(@arrayY);
    my ($n1,$d1,$d2) = (0,0,0);

    for (my $i=0;$i<scalar(@arrayX);$i++){
	$n1 += ($arrayX[$i]-$aveX)*($arrayY[$i]-$aveY);
	$d1 += ($arrayX[$i]-$aveX)*($arrayX[$i]-$aveX);
	$d2 += ($arrayY[$i]-$aveY)*($arrayY[$i]-$aveY);
    }
    return $n1/sqrt($d1)/sqrt($d2);
}



sub rank_array{    
    my @before = @{$_[0]};
    my $tolerance = $_[2];

    #print "Before Rank:\n";
    #print join ",",@before;
    #print "\n";

    my @after = ();
    my $array_size = @before;
    my @pos = ();
    for (my $i=0;$i<$array_size; $i++){
	push @pos,$i;
    }
    
    #sorting 
    
    for (my $i=0;$i<$array_size-1;$i++){
	my $min = $before[$i];
	my $min_pos = $i;
	for (my $j=$i+1;$j<$array_size;$j++){
	    if ($before[$j]<$min){
		$min = $before[$j];
		$min_pos = $j;
	    }
	}
	if ($min_pos!=$i){
	    my $tmp = $before[$i];
	    $before[$i] = $min;
	    $before[$min_pos] = $tmp;

	    $tmp = $pos[$i];
	    $pos[$i] = $pos[$min_pos];      #rank-1 <->  position
	    $pos[$min_pos] = $tmp;
	}    
    }

    my %rank_hash = ();    #position <-> rank-1
    for (my $i=0; $i<$array_size;$i++){
	$rank_hash{$pos[$i]} = $i;
    }
    for (my $i=0; $i<$array_size; $i++){
	$after[$i] = $rank_hash{$i}+1;
    }

    # check for the same rank
    for (my $start=0; $start<$array_size-1; $start++){
	#print $start."\n";
	my $end = $start;
	while (($end+1<$array_size) and (abs($before[$end+1]-$before[$end])<$tolerance)){
	    $end++;
	}
	my $rank_value = ($start+1+$end+1)/2;
	#printf "End: %d\n",$end;
	#printf "Rank_value: %.2f \n",$rank_value;
	for (my $i=$start; $i<=$end; $i++){
	    $after[$pos[$i]] = $rank_value;
	}
	$start = $end;
    }    
    @{$_[1]} = @after;
}
