#!/usr/bin/perl
#Author: liuchichuan@genomics.cn
#Date: 2016-10-21

use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';

my ($list, $out);
GetOptions
(
        "list=s" => \$list,
        "out=s" => \$out
);

if (!$list) 
{
	die <<USAGE;
==============================================================================
Description: Union genes among samples
Usage: perl $0 [options]
Options:
        * -list     gene expression files,
                    format: a.gene.fpkm.xls,
                            b.gene.fpkm.xls,
                            c.gene.fpkm.xls,......
        * -out      output file
		    format: *.gene.fpkm.xls
E.g.:
        perl $0 -list sample.fpkm.xls,... -out file
==============================================================================
USAGE
}


my @files = split (/\,+/,$list);
my %expr;
my @name;
foreach my $f (@files) 
{
	$f=abs_path($f);
	my ($sample)=$f=~/.*\/(.*)\.gene\.fpkm\.xls/;
	push @name,$sample;
	open IN, $f or die $!;
	<IN>;
	while (<IN>) 
	{
		next if (/^\#/ || /^s+$/);
		chomp;
		my @col=split (/\t+/,$_);
		my $geneid=$col[0];
		if (!defined $expr{$geneid}) 
		{
			for (my $k=0;$k<@files;$k++) 
			{
				$expr{$geneid}[$k]=0;
			}
			$expr{$geneid}[(scalar @name)-1]=$col[4];
		}
		else	
		{
			$expr{$geneid}[(scalar @name)-1]=$col[4];
		}
	}
	close IN;
}


my $head;
my $line;

open UNION, ">$out" or die $!;
foreach my $s (@name) 
{
	$head .= "$s\t";
}
print UNION "gene\_id\t$head\n";

foreach my $g (sort keys %expr)
{
	for (my $i=0;$i<@name;$i++) 
	{
		$line .= "$expr{$g}[$i]\t";		
	}
	print UNION "$g\t$line\n";
	$line="";
}
close UNION; 
