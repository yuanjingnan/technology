#!/usr/bin/perl -w 
use strict;
use Cwd 'abs_path';
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw/max min sum/;
=head1 Version
  Author: Wang Chongzhi, wangchongzhi@genomics.org.cn / wangcz85@163.com
  data : 2013.3.28,2013.4.24

=head1 Description
  generate R script for Principal Components Analysis to draw several figures displaying the analysis result.
  input: table of expression among samples. 

=head1 Usage
  perl princomp_draw.pl -expr <table file> [--o <output file>]
        --expr <str>   table format like the following:
########################################################################
GeneID  Sample_1  Sample_X sample_comtrol
ERCC1   1	2	3
ERCC2   3	6	9
......
########################################################################
        --o <str>          Out directory(default ".")
        --help             get help information and exit

=head1 Example
  eg 1: perl princomp_draw.pl --expr expression.xls
  eg 2: perl princomp_draw.pl --expr expression.xls -o workdir

=cut

my ($Indir, $Outdir, $Help, $Rbin, $Convert);
GetOptions(
	"i:s"=>\$Indir,
	"o:s"=>\$Outdir,
	"help"=>\$Help,
	"R"=>\$Rbin,
	"convert:s" => \$Convert
);
die `pod2text $0` if (!$Indir || $Help);
$Outdir ||= "./";
$Rbin = "unset R_LIBS_USER && \\\n unset R_LIBS && \\\n unset R_LIBS_SITE && \\\n export R_LIBS_USER=\"/ldfssz1/ST_CANCER/CGR/SHARE/tools/Anaconda/anaconda2/envs/r-3.5.1/lib/R/library\" && \\\n$Bin/../Bin_CentOS6/R" unless ($Rbin and -d $Rbin);
$Convert ||= "$Bin/../Bin_CentOS6/convert";

#######################generate expression file##############
my @samples = glob("$Indir/*.gene.fpkm.xls");
my $count = scalar @samples;
my %hash;
my @head;
for my $s (@samples){
	my ($name) = $s =~ /.*\/(.*).gene.fpkm.xls/;
	push @head, $name;
	open IN,$s or die $!;
	<IN>;
	while(<IN>){
		chomp;
		my @tmp = split(/\t/, $_);
		$hash{$tmp[0]}{$name} = $tmp[4];
	}
	close IN;
}

my $Expression_file = "$Outdir/$count"."Samples";
open OUT,">$Expression_file";
print OUT "Gene\t".join("\t", @head)."\n";
foreach my $g(sort keys %hash){
	print OUT "$g";
	for my $s(0 .. $#head){
		
		print OUT ($hash{$g}{$head[$s]}) ? "\t$hash{$g}{$head[$s]}" : "\t0.01";
	}
	print OUT "\n";
}
close OUT;

#######################generate R script#####################
open RSC,"> $Outdir/PCA.R" or die $!;
print RSC <<EOF;
expr <- read.table("$Expression_file",header = TRUE)
options(width = 160)
expr.dim <- dim(expr)
N_samp <- expr.dim[2] - 1
N_gene <- expr.dim[1]
expr.pca<- princomp(expr[,-1],cor = TRUE)
summary(expr.pca,loadings = TRUE)
pdf(file = "$Outdir/All.components.pdf", width=10, height=10)
library(plotrix, lib.loc = "/home/wangchongzhi/R/x86_64-unknown-linux-gnu-library/2.10")
pc.var <- expr.pca\$sdev ** 2
whichgap <- which.max(pc.var[-N_samp] - pc.var[-1])
gap.size <- 0.9 * (pc.var[whichgap] - pc.var[whichgap + 1])
gaps <- c(pc.var[whichgap + 1] + 0.05 * gap.size, pc.var[whichgap] - 0.05 * gap.size)
col <- color.gradient(c(0, 1), c(0, 1, 0), c(1, 0), N_samp)
gap.barplot(pc.var, gap=gaps, col = col, main = "Variances of Principal Components", ylim = c(0,expr.pca\$sdev[1] ** 2 - gap.size), xlab = "Components", ylab = "Variances", xlim = c(0,N_samp + 1), ytics=as.numeric(sprintf("%.3f", c(min(pc.var),max(pc.var),pc.var[whichgap+1],pc.var[whichgap]))))
box(lty = "solid", col = 'black')

load <- loadings(expr.pca)
pc <- data.frame(load[,1:3])
rnam <- rownames(pc)
rnam = gsub("rpkm", "", rnam, ignore.case = TRUE, perl = TRUE)
rnam = gsub("tpm", "", rnam, ignore.case = TRUE, perl = TRUE)
rnam = gsub("[_\\\\.]", "", rnam, ignore.case = TRUE, perl = TRUE)

#plot(pc[,-3], pch = "*", col = rainbow(N_samp), main = "loadings for PC1 and PC2 in samples")
#text(pc[,1],pc[,2], cex = 0.5, adj = c(-0.2,0.5))
#legend("topright", pch = 20, col = rainbow(N_samp),cex = 0.5,legend = rnam)

#plot(pc[,-2], pch = "*", col = rainbow(N_samp), main = "loadings for PC1 and PC3 in samples");
#text(pc[,1],pc[,3], cex = 0.5, adj = c(-0.2,0.5))
#legend("topright", pch = 20, col = rainbow(N_samp), cex = 0.5, legend = rnam)

#plot(pc[,-1], pch = "*", col = rainbow(N_samp), main = "loadings for PC2 and PC3 in samples");
#text(pc[,2],pc[,3], cex = 0.5, adj = c(-0.2,0.5))
#legend("topright", pch = 20, col = rainbow(N_samp), cex = 0.5, legend = rnam)

nouse = dev.off()

pdf(file = "$Outdir/All.scatter.pdf", width=10, height=10)
biplot(expr.pca, cex = 0.5, pc.biplot = FALSE)
title(main = "Genes Scatter Diagram")
nouse = dev.off()

pdf(file = "$Outdir/All.PCA-3D.pdf", width=10, height=10)
library(scatterplot3d, lib.loc = "/opt/blc/genome/biosoft/R/lib64/R/library")
pc.var.scaled <- 100 * pc.var/sum(pc.var)
pc.var.scaled = sprintf("%.2f", pc.var.scaled)
axis.names <- paste(names(pc.var), pc.var.scaled, sep = "(")
axis.names <- paste(axis.names, "%)", sep = "")
s3d <- scatterplot3d(pc, type = "h", angle = 40, highlight.3d = FALSE, scale.y = 0.7, pch = 19, main = "PCA 3D figure", las = 1, color = rainbow(N_samp), xlab = axis.names[1], ylab = axis.names[2], zlab = axis.names[3])
text(s3d\$xyz.convert(pc), labels = 1:N_samp, cex = 0.6)
legend("topright", pch = 20, col = rainbow(N_samp), cex = 0.6, legend = rnam)
nouse = dev.off()

pdf(file = "$Outdir/All.PCA-2D.pdf",width=10,height=10)
xlim_d=min(pc[,2])
xlim_u=max(pc[,2]) * 1.1
plot(pc[,2],pc[,3],pch=19,cex=1,col=rainbow(N_samp),xlab = axis.names[1], ylab = axis.names[2],main = "PCA Plot",las=1,xlim=c(xlim_d,xlim_u))
legend("topright", pch = 20, col = rainbow(N_samp), cex = 0.6, legend = rnam)
nouse = dev.off()
EOF
close RSC;

system("$Rbin/Rscript $Outdir/PCA.R && rm $Outdir/PCA.R && rm $Expression_file");
system("$Convert -density 300 -resize 40% $Outdir/All.components.pdf $Outdir/All.components.png");
system("$Convert -density 300 -resize 40% $Outdir/All.scatter.pdf $Outdir/All.scatter.png");
system("$Convert -density 300 -resize 40% $Outdir/All.PCA-3D.pdf $Outdir/All.PCA-3D.png");
system("$Convert -density 300 -resize 40% $Outdir/All.PCA-2D.pdf $Outdir/All.PCA-2D.png");

exit;
