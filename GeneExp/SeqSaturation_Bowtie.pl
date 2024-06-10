#!/usr/bin/perl -w

#  Author:		yumingqi
#  Date:		2014-10-17
#  Updated:		2014-10-21	dealing with multi-match & using $prefix
#  Updated:		2014-10-29	debug: subroutine &points do not need argu., so remove the argu.
#  Updated:		2015-03-31	debug: accurate gene number

use strict;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use File::Path 'mkpath';

my ($bam, $seqType, $gene2tr, $prefix, $cutoff, $convert, $Rscript, $samtools);
GetOptions (
	"bam=s" => \$bam,
	"seqType=s" => \$seqType,
	"gene2tr=s" => \$gene2tr,
	"prefix=s" => \$prefix,
	"cutoff:i" => \$cutoff,
	"convert:s" => \$convert,
	"Rscript:s" => \$Rscript,
	"samtools:s" => \$samtools
);
die <<USAGE unless $bam and $seqType and $gene2tr and $prefix;
Description:
	draw a curve graph for sequencing saturation 
Usage:
	perl $0 [options]
Options:
	*-bam       <file>   bowtie bam/sam file
	*-seqType   <str>    PE or SE
	*-gene2tr   <file>   format: geneID <tab> transcriptID
	*-prefix    <str>    prefix = 'outdir/sampleName'
	 -cutoff    <int>    ignore gene if mapped reads less than cutoff, default: 1
	 -convert   <path>   convert path
	 -Rscript   <path>   Rscript path
	 -samtools  <path>   samtools path 
Released:
	2014-10-21
Author:
	yumingqi\@bgitechsolutions.com
USAGE
$convert ||= "$Bin/../Bin_CentOS6/convert";
$Rscript ||= "unset R_LIBS_USER && unset R_LIBS && \\\n unset R_LIBS_SITE && \\\n export R_LIBS_USER=\"/ldfssz1/ST_CANCER/CGR/SHARE/tools/Anaconda/anaconda2/envs/r-3.5.1/lib/R/library\" && $Bin/../Bin_CentOS6/R/Rscript";
$samtools ||= "$Bin/../Bin_CentOS6/samtools";
$cutoff ||= 1;
my $outdir = dirname $prefix;
mkpath $outdir unless -e $outdir;

my %g2tr = map {($a, $b) = split /\s+/; $b => $a} `cat $gene2tr`;
my %gene = map {$_ => 1} values %g2tr;
my $GeneNum = keys %gene;

my ($previousX, %existGene, @Xset, @Yset);
my $unit = 10 ** 5; # i.e. 0.1 M/unit
sub points { # generate a point and save to array
	my $x = ++ $previousX;
	my $count = grep $existGene{$_} >= $cutoff, keys %existGene; # filter by cutoff
	my $y = $count / $GeneNum * 100;
	push @Xset, $x;
	push @Yset, $y;
}
sub readBam { # split the bam into units and draw point for each unit
	my $pipe = shift;
	my ($readsNum, %hash);
	while (<$pipe>) {
		my ($readID, $placeholder, $transcript) = split /\s+/;
		unless (defined $hash{$readID}) {$readsNum ++; $hash{$readID} = 1}; # filter multi-match
		$existGene{$g2tr{$transcript}} ++ if $transcript ne '*';
		<$pipe> if $seqType eq 'PE';
		if ($readsNum == $unit) {&points(); $readsNum = 0} # recount reads number
	}
	&points(); # the last point, in fact meaningless
}

my $sampleName = basename $prefix;
open SAM, $bam =~ /\.sam$/ ? $bam : "$samtools view $bam |" or die "Cannot open $bam\n";
&readBam(*SAM);
close SAM;

$prefix .= '.SeqSaturation';
my ($matrix, $pdf, $png, $Rfile) = ("$prefix.tmp", "$prefix.pdf", "$prefix.png", "$prefix.R");
open MATRIX, ">$matrix" or die "Cannot create $matrix";
print MATRIX "$Xset[$_]\t$Yset[$_]\n" for 0 .. $#Xset;
close MATRIX;
my $code = "pdf('$pdf',width=8,height=6)
library(ggplot2)
data=read.table('$matrix')
ggplot(data, aes(data[,1], data[,2]))+
	geom_line(color='royalblue', size=1)+labs(
	#title='Sequencing saturation of $sampleName',
	x='Amount of $seqType reads ( x100k )',
	y='Gene identification ratio ( % )'
	)+theme(plot.title=element_text(face='bold',size=20),axis.text=element_text(color='black'),panel.background = element_rect(fill='transparent'),panel.grid=element_line(color='grey'),panel.border=element_rect(fill='transparent',color='black'),axis.title=element_text(size=15))
dev.off()";
# `echo "$code" >$Rfile && $Rscript $Rfile && rm $matrix $Rfile`;
`echo "$code" >$Rfile && $Rscript $Rfile && rm $Rfile`;
`gzip $matrix`;
`$convert -density 300 -resize 40% $pdf $png`;
