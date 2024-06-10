#!/usr/bin/perl

#Author: liaoxinhui@genomics.org.cn
#Date: Sat May 24 15:23:40 CST 2014


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
Description: Draw Correlation Heatmap between Samples
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

my %cor = ();
my $min_value = 1;
foreach my $i (0..$#samples) {
	foreach my $j ($i..$#samples) {
		my $cor;
		if ($j eq $i) {
			$cor = 1;
		}else {
			$cor = CalculateCorrelation($samples[$i],$samples[$j]);
		}
		$min_value = ($min_value > $cor) ? $cor : $min_value;
		$cor{$samples[$i]}{$samples[$j]} = $cor;
	}
}
open my $fh_stat,">","$Outdir/AllSamples.correlation.xls" or die $!;
print $fh_stat join("\t",("Sample",@samples)),"\n";
foreach my $k (@samples) {
	my $out = $k;
	foreach my $m (@samples) {
		if (exists $cor{$k}{$m}) {
			$out .= "\t$cor{$k}{$m}";
		}else {
			$out .= "\t$cor{$m}{$k}";
		}
	}
	print $fh_stat "$out\n";
}
close $fh_stat;

### draw heatmap
my $cellsize = (@samples > 6) ? @samples : 8;
#my $cellsize = 480 / @samples;
#my $fontsize = ($cellsize > 15) ? 15 : $cellsize;
my $legendName = $Algorithm."_value";
open my $fh_rcode,">","$Outdir/correlation-heatmap.R" or die $!;
print $fh_rcode <<CODE;
library(reshape2)
library(ggplot2)
library(RColorBrewer)
x <- read.table("$Outdir/AllSamples.correlation.xls", sep = "\t", head = T)
xx = as.matrix(x[,-1])
rownames(xx) = names(x)[-1]
xx = melt(xx)
names(xx)=c("Var1","Var2","$legendName");
pdf("$Outdir/AllSamples.CorrelationHeatmap.pdf",width=$cellsize,height=$cellsize)
ggplot(xx, aes(Var1, Var2, fill=$legendName))+
 #geom_tile(width=0.8, height=0.8)+
  geom_tile(color='black')+
  geom_text(label=round(xx\$$legendName, 3))+
  scale_fill_gradient(low='#DEEBF7',high='#08519C')+
  theme(axis.text = element_text(angle=30, hjust=1,size=11,vjust=0,color='black'),
  panel.background = element_rect(fill='transparent'),
  panel.grid=element_line(color='grey'),legend.title = element_text(size = 13))+ 
  labs(x="",y="")
dev.off()
CODE

system("$Rscript $Outdir/correlation-heatmap.R"); #&& rm $Outdir/correlation-heatmap.R");
# system("$Convert -density 300 -resize 40% $Outdir/AllSamples.CorrelationHeatmap.pdf $Outdir/AllSamples.CorrelationHeatmap.png");

exit;

sub CalculateCorrelation {
	my ($sa,$sb) = @_;
	my ($exp_a,$exp_b);
	foreach my $g (keys %exp) {
		next if (!exists $exp{$g}{$sa} && !exists $exp{$g}{$sb});
		$exp_a .= ($exp{$g}{$sa}) ? "$exp{$g}{$sa}," : "0.01,";
		$exp_b .= ($exp{$g}{$sb}) ? "$exp{$g}{$sb}," : "0.01,";
	}
	chop $exp_a;
	chop $exp_b;
	open my $fh_rcode,">$Outdir/$sa-vs-$sb.correlation.R" or die $!;
	print $fh_rcode "x <- c($exp_a)\n y <- c($exp_b)\ncor(x, y, method = \"$Algorithm\")\n";
	close $fh_rcode;
	my ($index,$value) = split /\s+/, `$Rscript $Outdir/$sa-vs-$sb.correlation.R && rm $Outdir/$sa-vs-$sb.correlation.R`;
	chomp $value;
	return $value;
}
