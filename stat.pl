#!/usr/bin/perl -w
# Author: liuchichuan@genomics.cn
# Date: 20171201

use strict;
use Getopt::Long;

my ($input, $output);
GetOptions(
       "in=s"=>\$input,
       "out=s"=>\$output
);


my @file = split(/\,/,$input);
my ($total, $remain, $delete, $percent);

foreach my $stat ( @file )
{
	open STAT, $stat or die "Can't open file $stat:$!\n";
	while (<STAT>)
	{
		chomp;
##########paired-end reads
		if ($_ =~ /(\S+)\s+\(\S+\)\s+were\s+paired/)	# 60686993 (100.00%) were paired; of these:
		{
			$total += $1;
        	}
		if ($_ =~ /(\S+)\s+\((\S+)%\)\s+aligned\s+concordantly\s+0\s+times/)	# 60004713 (98.88%) aligned concordantly 0 times
		{
			$remain += $1;
			last;
		}

##########single-end reads
		if ($_ =~ /(\S+)\s+\(\S+\)\s+were\s+unpaired/)	# 20000 (100.00%) were unpaired; of these:
		{
			$total += $1;
		}
		if ($_ =~ /(\S+)\s+\((\S+)%\)\s+aligned\s+0\s+times/ )	# 1247 (6.24%) aligned 0 times
		{
			$remain += $1;
			last;
		}	
	}
	close STAT;
}

$delete = $total - $remain;
$percent = sprintf("%.2f%%", $delete / $total * 100);
open OUT, ">$output"or die "Can't open file $output:$!\n";
	print OUT "total\t$total\ndelete\t$delete\t$percent\n";
close OUT;
