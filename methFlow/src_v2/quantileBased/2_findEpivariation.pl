use strict;
use Getopt::Long qw( :config posix_default bundling no_ignore_case );
use Scalar::Util qw(looks_like_number);
use List::Util qw( min max );
use POSIX;

####### INITIALIZING INPUT PARAMETERS WITH DEFAULT #######
my $file = '';
my $detPvalFile = '';
my $chr = 1;
my $start = 2;
my $end = 3;
my $window = 1000;
my $start_ind = '';
my $end_ind = -1;
my $b_thres = 0.1;
my $n_thres = 3;
my $quan_upper = 0.999;
my $quan_lower = 0.001;
my $mode = "BOTH";
my $outfile = "STDOUT";


####### GET INPUT FROM USER #######
help() if ( @ARGV < 1);
GetOptions('file|f=s' => \$file,
	   'detfile|d=s' => \$detPvalFile,
	   'chr|c=i' => \$chr,
	   'start|s=i' => \$start,
	   'end|e=i' => \$end,
	   'window|w=i' => \$window,
	   'start_ind|a=i' => \$start_ind,
	   'end_ind|b=i' => \$end_ind,
	   'beta_thres|t=f' => \$b_thres,
	   'p_thres|p=i' => \$n_thres,
           'quan_upper|9=f' => \$quan_upper,
           'quan_lower|1=f' => \$quan_lower,
	   'outfile|o=s' => \$outfile,
	   'mode|m=s' => \$mode,
	   'help|h' => \my $help,
	   'verbose|v' => \my $verbose);

help() if ( defined $help );
verbose() if(defined $verbose or !defined($file) or !defined($start_ind) or !defined($end_ind));


my $i = 0; my $j = 0;
my @data = ();
my $is_stdin = 0;
my $IN;
if($file eq "stdin" || $file eq "STDIN"){
    $IN = *STDIN;
    $is_stdin++;
}else{
    open $IN, "<", $file or die $!
}

####### READING THE NORMALIZED BETA VALUE FILE ########
my $header = <$IN>; chomp $header; my @h =  split(/\t/,$header);
while(<$IN>){
    chomp $_;
    my @s = split(/\s+/,$_);
    for($j = 0; $j<@s; $j++){
	$data[$i][$j] = $s[$j];
    }
    $i++;
}
close $IN unless $is_stdin;

my $nrow = $i;
my $ncol = $j;
undef $i;
undef $j;


my %detColHash; ### Hash to link column names to index (0 based) for Dection P value file
my %detRowHash; ### Hash to link row names to index (0 based) for Detection P value file
my @detp = ();

###### READING DETECTION P value file if provided ######
if($detPvalFile ne ""){
  open $IN, "<", $detPvalFile or die $!;
  my $detheader = <$IN>;
  chomp $detheader;
  my @deth =  split(/\t/,$detheader);
  @deth = map {s/.Detection Pval//g; $_; } @deth;
  @detColHash{@deth} = 0..$#deth;
  while(<$IN>){
      chomp $_;
      my @s = split(/\s+/,$_);
      for($j = 0; $j<@s; $j++){
          $detp[$i][$j] = $s[$j];
      }
      $detRowHash{$s[0]} = $i;
      $i++;
  }
  close $IN;
}


#### In case end_ind is not provided (default = -1), it is set to max column number #######
if($end_ind==-1){ $end_ind = $ncol-1; print STDERR "Col for end ind: $end_ind"; }

#### INITIALIZING OUTPUT VARIABLES TO BE USED FOR OUTPUT COLUMNS #######
my @sign = (0) x $nrow; ## Will be 0 or 1 if the probe in the window lies in significant DMR
my @sign_ind = ('NA') x $nrow; ## Will be list of comma separated sample IDs that carry the DMR
my @mean = ('NA') x $nrow; ## mean for each probe across all samples
my @direction_all = ('NA') x $nrow; ## Directionality of DMR for each carrier
my @num = (0) x $nrow;



my @uquan = (); ## Upper Quantile for each probe
my @lquan = (); ## Lower Quantile for each probe
##my @a=();
for(my $i = 0; $i<$nrow; $i++){
  my @a = @{$data[$i]}[($start_ind)..($end_ind)];
  $uquan[$i] = quantile(\@a,$quan_upper);
  $lquan[$i] = quantile(\@a,$quan_lower);

}

##### Looping througn each sample column and running sliding window analysis ####
for(my $k=$start_ind; $k<=$end_ind; $k++){
    my @direction = ('NA') x $nrow;
    for(my $i = 0; $i<$nrow; $i++){
        my $count = 0; my $count1=0; my $count2=0; my $total_p = 0; my $sum = 0;
        my $countq1=0; my $countq2=0; my @hyper; my @hypo;
	my $pvalcount=0; my @rle=();
        for(my $j = $i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$i][$chr]; $j++){
	    ## Count number of probes for each sample outside the quantile threshold and in same direction
	    ## array rle records if consecutive probes are are outside the quantiles or not. 
            if(looks_like_number($data[$j][$k])){
                if($data[$j][$k] >= $uquan[$j] + $b_thres){ $countq1++; }
                if($data[$j][$k] <= $lquan[$j] - $b_thres){ $countq2++; }
                $total_p++;
		if($data[$j][$k] >= $uquan[$j]){$rle[scalar(@rle)] = 2; push @hyper,$data[$j][2];}
		elsif($data[$j][$k] <= $lquan[$j]){$rle[scalar(@rle)] = 1; push @hypo,$data[$j][2];}
		else{$rle[scalar(@rle)] = 0;}
            }
	    if($detPvalFile ne ""){
  	       my $pi = $detRowHash{$data[$j][0]};
	       my $pj = $detColHash{$h[$k]};
	       if($detp[$pi][$pj] >0.01 && looks_like_number($detp[$pi][$pj]) ){ $pvalcount++;}
	    }
        }

	##flag tells whether the window passes the threshold to be cosidered as carrying Epivariation.
	##We required our Epivariation window to be at least 100 bp long
        my $flag='';
        if($mode eq "BOTH"){ $flag = (($countq1>=$n_thres && (max(@hyper)-min(@hyper)>=100))||($countq2>=$n_thres && (max(@hypo)-min(@hypo)>=100))); }
        elsif($mode eq "HYPER_ONLY") {$flag = ($countq1>=$n_thres && (max(@hyper)-min(@hyper)>=100)); }
        elsif($mode eq "HYPO_ONLY") {$flag = ($countq2>=$n_thres && (max(@hypo)-min(@hypo)>=100)); }
        else{ $flag = (($countq1>=$n_thres && (max(@hyper)-min(@hyper)>=100))||($countq2>=$n_thres && (max(@hypo)-min(@hypo)>=100))); }

	###Regject the window if there is even one probe with failed detection P value
	if($flag && $detPvalFile ne ""){$flag = $pvalcount==0;}

	###Check whether three consecutive probes are outide the quantile threshold range.
	if($flag){ 
		my $rleString = join '', @rle;
		if($countq1>=$n_thres){ $flag = $rleString =~ /222/}
		if($countq2>=$n_thres){ $flag = $rleString =~ /111/}
	}

	###if the window passes our criteria, add the sample name and direction of Epivariation to the list.
        if($flag){
            for(my $j =$i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$i][$chr]; $j++){
                $sign[$j] = 1;
                $h[$k] =~ s/.AVG_Beta$//;
                my $temp = ",",$h[$k].",";
                my @tt = split(",", $sign_ind[$j]);
                my $match_flag = 0;
                if($sign_ind[$j] eq "NA") {$sign_ind[$j] = $h[$k];}
                elsif(!($h[$k]~~@tt)){$sign_ind[$j] .= ",".$h[$k];}
            }
        }
        if(($countq1>=$n_thres)&&($countq2>=$n_thres) && $flag){
            for(my $j =$i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$i][$chr]; $j++){ $direction[$j] ="BOTH";}
        }elsif(($countq1>=$n_thres) && ($mode eq "BOTH" || $mode eq "HYPER_ONLY") && $flag){
            for(my $j =$i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$i][$chr]; $j++){ $direction[$j] ="HYPER";}
        }elsif(($countq2>=$n_thres) && ($mode eq "BOTH" || $mode eq "HYPO_ONLY") && $flag){
            for(my $j =$i; $j<$nrow && $data[$j][$start] <= $data[$i][$end]+$window && $data[$j][$chr] eq $data[$i][$chr]; $j++){ $direction[$j] ="HYPO";}
        }
    }
    for(my $i = 0; $i<@direction; $i++){
        if($direction_all[$i] eq "NA"){ $direction_all[$i] = $direction[$i]}
        else{
            if($direction[$i] ne "NA"){ $direction_all[$i] .= ",".$direction[$i]}
        }
    }

}


####OUTPUTTING THE RESULTS
my $OUT;
my $is_stdout = 0;
if($outfile eq "stdout" || $outfile eq "STDOUT"){
    $is_stdout++;
}else{
    open(STDOUT, ">".$outfile);
}

print STDOUT $header,"\tQuan",$quan_lower,"\tQuan",$quan_upper,"\tSign_individuals_t",$b_thres,"_n",$n_thres,"_w",$window/1000,"k_",$mode,"\tSign_direction_t",$b_thres,"_n",$n_thres,"_w",$window/1000,"k_",$mode
,"\tSign_window_t",$b_thres,"_n",$n_thres,"_w",$window/1000,"k_",$mode,"\n";

for(my $i = 0; $i <$nrow; $i++){
    for(my $j = 0; $j<$ncol; $j++){
        print STDOUT $data[$i][$j],"\t";
    }
    print STDOUT $lquan[$i],"\t",$uquan[$i],"\t",$sign_ind[$i],"\t",$direction_all[$i],"\t",$sign[$i],"\n";
}

close STDOUT unless $is_stdout;


###Function to calculate quantiles
sub quantile{
 my ($arr,$prob) = @_;
 my @sorted_arr = sort {$a <=> $b} @$arr;
 my $na_omitted = "NA";
 @sorted_arr = grep { $_ ne $na_omitted } @sorted_arr;
 my $k = ceil($prob*scalar(@sorted_arr));
 $sorted_arr[$k-1]
}


sub help{
    print STDERR "Usage: perl window_analysis.pl -f filename -d detPfile -c chr_col -s start_col -e end_col -w window -a ind_start -b ind_end  -9 upper_quantile -l lower_quantile -t beta_thres -p probe_thres -ooutfile -m mode\n";
    exit();
}

sub verbose{
print STDERR "Usage: perl window_analysis.pl -f filename -d detPfile -c chr_col -s start_col -e end_col -w window -a ind_start -b ind_end  -9 upper_quantile -l lower_quantile -t beta_thres -p probe_thres -ooutfile -m mode
-f, --filename=FILENAME\t\t\tInput filename sorted  by chromosome and probe start, can be stdin or STDIN while piping (REQUIRED)
-d, --detfile=DETP_FILENAME\t\t\tFile containing Detection P values, rows as probes, columns as samples (REQUIRED)
-a, --start_ind=IND_START_COL\t\tColumn from where sample beta value starts (REQUIRED)
-b, --end_ind=IND_END_COL\t\tColumn from where sample beta value ends (REQUIRED)
-o, --outfile=OUTFILE\t\t\tOutput file name with three extra columns for Significant indiduals, directionality and window, 
\t\t\t\t\tcan be STDOUT or stdout for printing commandline, can be used for piping and can be redirected to a file using '>' 
\t\t\t\t\t(default = STDOUT)
-c, --chr=CHR_COL\t\t\tColumn containing probe chromosome information (default = 1)
-s, --start=START_COL\t\t\tColumn containing probe start information (default = 2)
-e, --end=END_COL\t\t\tColumn containing probe end information (default = 3)
-w, --window=WINDOW_SIZE\t\tSize of the sliding window (default = 1000)
-1, --quan_lower=LOWER_QUANTILE\t\tThreshold for Lower quantile (for Hypo DMR) (default = 0.1%ile)
-9, --quan_upper=UPPER_QUANTILE\t\tThreshold for Upper quantile (for Hyper DMR) (default = 99.9%ile)
-t, --beta_thres=BETA_THRESHOLD\t\tBeta value threshold +/- of control min/max (default = 0.1)
-p, --p_thres=NO_OF_PROBES_THRESHOLD\tNumber of significant probes in window threshold (default = 2)
-m, --mode=MODE\t\t\t\tCan be BOTH, HYPER_ONLY or HYPO_ONLY. runs script in both or one direction (default=BOTH)
-h, --help\t\t\t\tUsage summary
-v, --verbose\t\t\t\tDetailed Usage Information\n";
exit();
}
