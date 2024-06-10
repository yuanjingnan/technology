#!/usr/bin/perl

#Author: liaoxinhui@genomics.org.cn
#Date: Mon May 26 10:15:33 CST 2014


use strict;
use warnings;
use FindBin '$Bin';
use Getopt::Long;

my ($Indir,$Outdir,$Algorithm,$Rscript,$Convert);
GetOptions(
	"indir:s" => \$Indir,
	"cor:s" => \$Algorithm,
	"Rpath:s" => \$Rscript,
	"convert:s" => \$Convert,
	"outdir:s" => \$Outdir
);
if (!defined $Indir) {
	print STDERR <<USAGE;
Description: Draw Hclust Tree for Samples
Usage: perl $0 [options]
Options:
	* -indir     input directory, containing *.Gene.rpkm.xls(gene expression level files)
	  -Rpath     path of Rscript, default: /opt/blc/genome/bin/Rscript
	  -convert   path of convert, default: /usr/bin/convert
	  -cor       algorithm for correlation, default: pearson
	  -outdir    output directory, default: current directory
E.g.:
	perl $0 -indir result/GeneExp -cor pearson -outdir ./

USAGE
	exit;
}

$Rscript ||= "unset R_LIBS_USER && \\\n unset R_LIBS && \\\n unset R_LIBS_SITE && \\\n export R_LIBS_USER=\"/ldfssz1/ST_CANCER/CGR/SHARE/tools/Anaconda/anaconda2/envs/r-3.5.1/lib/R/library:\$R_LIBS_USER\" && \\\n$Bin/../Bin_CentOS6/R/Rscript";
$Convert ||= "$Bin/../Bin_CentOS6/convert";
$Algorithm ||= "pearson";
$Outdir ||= "./";

my @files = glob("$Indir/*.gene.fpkm.xls");
my %exp = ();
my @samples = ();
foreach my $f (@files) {
	my ($sample) = $f =~ /.*\/(.*)\.gene\.fpkm\.xls/;
	push @samples,$sample;
	open my $fh_rpkm,"<",$f or die $!;<$fh_rpkm>;
	while (<$fh_rpkm>) {
		chomp;
		my @a = split /\s+/;
		$exp{$a[0]}{$sample} = $a[4]; # rpkm
	}
	close $fh_rpkm;
}

open my $fh_exp,">","$Outdir/all.exp.xls" or die $!;
print $fh_exp join("\t",("Sample",@samples)),"\n";
foreach my $g (keys %exp) {
	my $print = $g;
	foreach my $s (@samples) {
		$print .= ($exp{$g}{$s}) ? "\t$exp{$g}{$s}" : "\t0.01";
	}
	print $fh_exp "$print\n";
}
close $fh_exp;

my $rownames = join(",", map("\"$_\"",@samples));
open my $fh_rcode,">","$Outdir/hclust.R" or die $!;
print $fh_rcode <<CODE;
pdf("$Outdir/AllSamples.HclustTree.pdf",width=10,height=10,pointsize=15)
par(oma=c(8,0,8,0))
data <- read.table("$Outdir/all.exp.xls",head=TRUE)
len <- length(data)
exp <- data[,2:len]
S <- t(exp)
rownames(S)<-c($rownames)
AllSamples=dist(S,method="euclidean")
out.hclust=hclust(AllSamples,method="ward")
plot(out.hclust,sub="")
dev.off()
CODE

system("$Rscript $Outdir/hclust.R && rm $Outdir/hclust.R && rm $Outdir/all.exp.xls");
system("$Convert -density 300 -resize 40% $Outdir/AllSamples.HclustTree.pdf $Outdir/AllSamples.HclustTree.png");

exit;
