#!/usr/bin/perl

=pod
description: summarize all samples' gene expression
author: Du Shuai, dushuai@genomics.cn
        Zhang Fangxian, zhangfx@genomics.cn
created date: 20110127
modified date: 20110215
modified for GeneAndIsoformExp_RSEM:2014-08-28
=cut

use strict;
use warnings;
use File::Basename qw(dirname basename);

my $usage = << "USAGE";
description: summarize all samples' gene and isoform expression
usage: perl $0 indir output
USAGE

my ($indir, $output) = @ARGV;
die $usage if (!defined $indir or !defined $output);

my $header_gene = "gene_id\ttranscript_id(s)\tSymbol";
my (@samples_gene, %result_gene);

#########------- read files for all.gene.FPKM.xls -------#########

for my $file_gene (glob("$indir/*/*.gene.fpkm.xls.new")) {
	my $keyname_gene = basename($file_gene);
	$keyname_gene =~ s/\.gene\.fpkm\.xls\.new$//;
	push @samples_gene, $keyname_gene;

	&showLog("read file $file_gene");
	open GENE, "< $file_gene" or die $!;

	my $temp = <GENE>;
	chomp $temp;
	my @tabs = split /\t/, $temp;
	$header_gene .= "\t$keyname_gene\_FPKM";

	while (<GENE>) {
		chomp;
		@tabs = split /\t/, $_;
		@{$result_gene{"$tabs[0]\t$tabs[1]\t$tabs[5]"}{$keyname_gene}} = $tabs[4];
                
	}
	close GENE;
}

for my $gene_gene (keys %result_gene){
        if (!exists $result_gene{$gene_gene}{$samples_gene[0]}) {
                $result_gene{$gene_gene}{$samples_gene[0]}->[0] = 0;
        }
}
##########------- output all.gene.FPKM.xls --------#########

open OUT, "> $output" or die $!;
print OUT "$header_gene\n";
for my $gene_gene (sort {$result_gene{$b}{$samples_gene[0]}->[0] <=> $result_gene{$a}{$samples_gene[0]}->[0]} keys %result_gene) {
        print OUT $gene_gene;
        for (@samples_gene) {
                if (exists $result_gene{$gene_gene}{$_} && @{$result_gene{$gene_gene}{$_}} == 1) {

                          print OUT join("\t", "",@{$result_gene{$gene_gene}{$_}});

                } else {

                          print OUT "\t0" ;
                      }
                }

       print OUT "\n";
}
close OUT;
&showLog("done");
exit 0;

#######------ sub ------#######

sub showLog {
	my ($info) = @_;
	my @times = localtime; # sec, min, hour, day, month, year
	print STDERR sprintf("[%d-%02d-%02d %02d:%02d:%02d] %s\n", $times[5] + 1900, $times[4] + 1, $times[3], $times[2], $times[1], $times[0], $info);
}
