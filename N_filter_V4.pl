#!/usr/bin/perl -w

use strict;
use FindBin qw($Bin);
use File::Spec qw(catfile);
use File::Basename qw(dirname basename);
use Getopt::Long;
Getopt::Long::Configure qw(no_ignore_case);
use Switch;
use lib "$Bin";
use BISBasics;
use constant DEBUG => 0;
use Data::Dumper;

my $usage=<<END;

  perl $0 -fq1 XXXX.fq.gz [-fq2 XXXXX.fq.gz] -o output_dir
  return tite number as 1102,1107,1109

END
my ($fq1,$fq2,$output,$cycle_threshold,$tile_threshold)=('','','',0.5,100);
my %mess_tiles;
my %mess_cycles;

GetOptions(
	"fq1:s"=>\$fq1,
	"fq2:s"=>\$fq2,
	"-o:s"=>\$output
);

if(! $fq1 or !$output){
	print $usage;
	exit;
}



# Treating fq1 file
if( ! -e $fq1 ){
	print STDERR "File $fq1 don't exist\n";	
	exit;
}
my $fq_base_name=&basename($fq1);
my $fq_base_dir=&dirname($fq1);
#$fq_base_dir=~s/$fq_base_name//;
my $fq1_check=$fq_base_dir."/1.fqcheck";

if( ! -e $fq1_check ){
	print STDERR "File $fq1_check don't exist\n";	
	exit;
}

&_trim_N_before_filter($fq1_check, $fq1, \%mess_tiles,\%mess_cycles);

# Treating fq2 file

if($fq2 and -e $fq2){
	my $fq2_check=$fq_base_dir."/2.fqcheck";
	
	if( ! -e $fq2_check ){
		print STDERR "File $fq2_check don't exist\n";
		exit;
	}
	&_trim_N_before_filter($fq2_check, $fq2, \%mess_tiles,\%mess_cycles);	
}



my $n_tiles_str = join (",", keys %mess_tiles);
my $str='';
if($n_tiles_str){
	$str="--tile ".$n_tiles_str;
}
open F,">$output/N_tile.txt" or die "Can't open file $output/N_tile.txt\n";
print F $str;
close F;








sub _trim_N_before_filter
{
	my $fqcheck = shift;
	my $input_fq = shift;
	my $mess_tiles_ref = shift; # hash structure
	my $mess_cycles_ref = shift; # hash structure

	my $fqcheck_handle;

	if( !($fqcheck_handle = &bis_open("$fqcheck", 'r')) )
	{
		&bis_log('E', "can not open $fqcheck");
		exit(1);
	}

	##文件样式
	# the default quality shift value is: -64, 18894458 sequences, 925828442 total
	#length, Max length:49, average length:49.00
	#Standard deviations at 0.25:  total 0.00%, per base 0.01%
	#             A     C     G     T     N    
	#Total     26.0  23.6  23.5  26.9   0.0    
	#base   1  16.8  45.0  27.9  10.4   0.0    
	#base   2  20.1  15.3  25.3  39.2   0.0

	my @muddy_cycles;
	my %baseline;
	my $line_count = 0;
	while(<$fqcheck_handle>)
	{
		next if (!/^base/);
		chomp;
		my @tmp = split /\s+/;

		if ($tmp[6] >= $cycle_threshold){
		    push @muddy_cycles, $tmp[1];
		}elsif($tmp[1]>20){
		}
	}
	close $fqcheck_handle;



	## 如果没有超过cycle_threshold的base存在直接返回
	if (@muddy_cycles == 0){
		return;
	}

	print Dumper(\@muddy_cycles) if DEBUG;

	# existed N contamination,pick up tiles and cycles
	my $fq_handle;
	if( !($fq_handle = &bis_open("$input_fq", 'r')) )
	{
		&bis_log('E', "can not open $input_fq");
		exit(1);
	}

	my %composition_by_tile;
	my %percentage_by_tile;
	my $total_reads= 0;

	while(my $seq_id = <$fq_handle>){

		$total_reads++; 

                 my $tile_id=undef;
                # read id
                my @tmp = split /:/,$seq_id;
                if (@tmp==10){
                 $tile_id = $tmp[4];
                }else{
                 $tile_id = $tmp[2];}#edit by guanhaijiao,20160323

		# read sequence
		my $seq = <$fq_handle>;
		my @read_bases = split //,$seq;

		foreach my $m_cycle (@muddy_cycles){
			my $pos = $m_cycle -1; # cycle in fqcheck file start with 1-base
			my $base = $read_bases[$pos]; # get nucletide
			$composition_by_tile{$pos}{$tile_id}{$base}++;
		} 

		<$fq_handle>; # skip sequence comment str
		<$fq_handle>; # skip sequence quality str
	}
	close $fq_handle;

	foreach my $cycle (keys %composition_by_tile){
		foreach my $tile (keys %{$composition_by_tile{$cycle}}){
			my $a = $composition_by_tile{$cycle}{$tile}{"A"} || 0;
			my $c = $composition_by_tile{$cycle}{$tile}{"C"} || 0;
			my $t = $composition_by_tile{$cycle}{$tile}{"T"} || 0;
			my $g = $composition_by_tile{$cycle}{$tile}{"G"} || 0;
			my $n = $composition_by_tile{$cycle}{$tile}{"N"} || 0;
			
			my $sum = $a + $t + $c +$g + $n;

			$percentage_by_tile{$cycle}{$tile}{"A"} = sprintf("%.1f",$a/$sum*100);
			$percentage_by_tile{$cycle}{$tile}{"C"} = sprintf("%.1f",$c/$sum*100);
			$percentage_by_tile{$cycle}{$tile}{"T"} = sprintf("%.1f",$t/$sum*100);
			$percentage_by_tile{$cycle}{$tile}{"G"} = sprintf("%.1f",$g/$sum*100);
			$percentage_by_tile{$cycle}{$tile}{"N"} = sprintf("%.1f",$n/$sum*100);
		}
	}

	print Dumper(\%percentage_by_tile) if DEBUG;

	#TO-DO: move those parameters to config
	my $tile_minor_threshold = 1.0; # 1%
	my $tolerant_tile_ratio = 0.5; # 50% tile
	my $stable_cycle_start = 20; # start from 20
	my $tolerant_base_variance = 10; # 10%

	# %percentage_by_tile keys: cycle=>tile=>base
	my $s_tile_count = 0; # slightly messed
	foreach my $m_cycle (keys %percentage_by_tile){
		my $t_tile_count = 0;
		my %mess_swath;

		foreach my $tile (keys %{$percentage_by_tile{$m_cycle}}){
			$t_tile_count++;
			my $tile_a_percent = $percentage_by_tile{$m_cycle}{$tile}{"A"} || 0;
			my $tile_t_percent = $percentage_by_tile{$m_cycle}{$tile}{"T"} || 0;
			my $tile_c_percent = $percentage_by_tile{$m_cycle}{$tile}{"C"} || 0;
			my $tile_g_percent = $percentage_by_tile{$m_cycle}{$tile}{"G"} || 0;
			my $tile_n_percent = $percentage_by_tile{$m_cycle}{$tile}{"N"} || 0;

			if ($tile_n_percent >= $tile_threshold){
				${$mess_tiles_ref}{$tile}++;
				$mess_swath{int($tile/100)}++;
				print "swath: ",int($tile/100),"\n" if $mess_swath{int($tile/100)} >= 13; 

			}
			elsif($tile_n_percent >= $tile_minor_threshold){
				$s_tile_count++;
			}
		}



	}
}
