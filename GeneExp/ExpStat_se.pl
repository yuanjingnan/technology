#!/usr/bin/perl
#Author: liaoxinhui@genomics.org.cn
#Date:  Thu Jul  9 04:02:19 CST 2015


use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';
use File::Basename;

my ($indir, $output);
GetOptions(
	"indir:s" => \$indir,
	"output:s" => \$output
);

if (! $indir || ! -d $indir || ! $output) {
	die "usage: perl $0 -indir indir -output output.xls\n";
}

my %stat = ();

foreach my $f (glob("$indir/*/*.gene.fpkm.xls")) {
	my ($s) = $f =~ /.*\/(.*)\.gene.fpkm.xls/;
	my $map_stat = "$indir/$s/$s.Bowtie2Gene.MapReadsStat.xls";	#lcc 
	my $trans_fpkm = "$indir/$s/$s.transcript.fpkm.xls";
	open F,$f or die $!;<F>;
	while (<F>) {
		next if (/^\s*$/ || /^#/);chomp;
		if (/^BGI_novel_.*/) {
			$stat{$s}{'novel_gene'} ++;
		}else {
			$stat{$s}{'know_gene'} ++;
		}
		$stat{$s}{'total_gene'} ++;
	}
	close F;
	open S, $map_stat or die $!;
	while (<S>) {
		next if (/^\s*$/ || /^#/);chomp;
        ########lcc
                if (/Total\s+Reads\s+(\d+)/) {
                        $stat{$s}{'total_reads'} = $1;
                }elsif (/Unique\s+Match\s+(\d+)\s+(\S+)/) {
                        $stat{$s}{'uniq_ratio'} = $2;
                }elsif (/Total\s+Mapped\s+Reads\s+(\d+)\s+(\S+)/) {
                        $stat{$s}{'map_ratio'} = $2;
                }
       ################
	}
	close F;
	open T, $trans_fpkm or die $!;<T>;
	while (<T>) {
		next if (/^\s*$/ || /^\s*#/);chomp;
		if (/^BGI_novel.*/) {
			$stat{$s}{'novel_tr'}++;
		}else {
			$stat{$s}{'know_tr'}++;
		}
		$stat{$s}{'total_tr'}++;
	}
	close T;
}

open O,">$output" or die $!;
print O "Sample\tTotal CleanReads\tTotal MappingRatio\tUniquely MappingRatio\tTotal GeneNumber\n";
foreach my $s (sort keys %stat) {
	my $out = $s;
	for('total_reads', 'map_ratio', 'uniq_ratio', 'total_gene') {
		$stat{$s}{$_} ||= 0;
		$out .= "\t$stat{$s}{$_}";
	}
	print O "$out\n";
}
close O;
