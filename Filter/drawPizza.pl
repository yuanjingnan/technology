#!usr/bin/perl -w

#Author:  yumingqi@bgitechsolutions.com
#Date:  Aug 2014

use strict;
use Getopt::Long;
use File::Path;
use FindBin qw($Bin);

my ($infile, $outdir, $Rscript, $convert);
GetOptions(
	"infile=s" => \$infile,
	"outdir:s" => \$outdir,
	"Rscript:s" => \$Rscript,
	"convert:s" => \$convert
);
die "
Description:
	Draw pie chart of raw reads class. Output a pdf file and a png file.
Usage:
	*-infile   <file>    input raw reads class file: *.xls
	 -outdir   <path>    default: current
	 -Rscript  <app.>    Rscript path
	 -convert  <app.>    Convert path
e.g.
	perl $0 -infile .../.../sample.xls
Released:
	Aug 2014
Author:
	yumingqi\@bgitechsolutions.com\n
" unless defined $infile;
$Rscript ||= "unset R_LIBS_USER && \\\n unset R_LIBS && \\\n unset R_LIBS_SITE && \\\n export R_LIBS_USER=\"/ldfssz1/ST_CANCER/CGR/SHARE/tools/Anaconda/anaconda2/envs/r-3.5.1/lib/R/library\" && \\\n$Bin/../Bin_CentOS6/R/Rscript";
$convert ||= "$Bin/../Bin_CentOS6/convert";
# $convert ||= "convert";
$outdir = $outdir ? File::Spec->rel2abs($outdir) : "$Bin";
mkpath $outdir unless -e $outdir;
my ($sample, %hash, $content,$content2);
open IN, $infile or die "Cannot open file: $_\n";
($sample) = $infile =~ /.*\/(\S+)\.filter.stat.xls/;
my $sample_sub = $sample;
$sample_sub =~ s/\./_/g;
<IN>;
($hash{'Clean reads'}) = reverse split /\s+/,<IN>;
while (<IN>) {
	if (/^Discard .+ to (.+?)\t(\d+)/){
		my $a=$1; my $b=$2;
		$a=~s/low/Low/g if($a=~/low/);
		$hash{$a} = $b;
	}
}
close IN;
for ('N','Adapter','Low qual','Clean reads') {$content .= "'$_'=$hash{$_},"}
for ('N','Adapter','Low qual','Clean reads') {1 while $hash{$_}=~ s/(\d)(\d{3})((:?,\d\d\d)*)$/$1,$2$3/; $content2 .= "'$_'=\"$hash{$_}\",";}
chop $content2;
chop $content;
print "$content\n\n";
open RFILE, ">$outdir/$sample.R" or die "Cannot create and open Rscript file in $outdir.";
print RFILE <<CODE;
pdf("$outdir/$sample.RawReadsClass.pdf",width=8,height=6)
library(ggplot2)
values = c($content)
values2 = c($content2)
number = c(1:length(values))
colours = c("#FF8C00","#EEEE00","#90EE90","#87CEFA","#AB82FF","#FF3030")
Legend <- paste(number, ". ", labels(values), "  (", values2, "; ", round(values/sum(values)*100, 2), "%)", sep="")
data <- data.frame(Percentage = values)
pie <- ggplot(data, aes(x = "" ,y = Percentage, fill = Legend))
pie = pie + geom_bar(width = 3, stat = "identity")
pie = pie + coord_polar("y")
#pie = pie + scale_fill_manual(values = colours, name = "Reads of ($sample_sub) contains:")
pie = pie + scale_fill_manual(values = colours)
#pie = pie + labs(title = "Classification of Raw Reads ($sample_sub)")
pie = pie + theme(
	axis.title = element_blank(),
    axis.ticks = element_blank(),
    axis.text = element_blank(),
    panel.background = element_blank(),
    panel.grid = element_blank(),
    plot.title = element_text(face="bold")
)
pie
dev.off()
CODE
close RFILE;
`$Rscript $outdir/$sample.R && rm $outdir/$sample.R`;
`$convert -density 300 -resize 40% $outdir/$sample.RawReadsClass.pdf $outdir/$sample.RawReadsClass.png`;
