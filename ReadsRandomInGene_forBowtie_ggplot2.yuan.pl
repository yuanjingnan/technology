#!/usr/bin/perl -w

#Author: yumingqi
#Date: 2014-11-07

use strict;
use Getopt::Long;
use File::Basename;
use FindBin '$Bin';
use File::Path 'mkpath';

my ($bam, $len, $seqtype, $pr, $n, $Rscript, $convert, $samtools);
GetOptions (
	"bam=s" => \$bam,
	"len=s" => \$len,
	"seqType=s" => \$seqtype,
	"prefix=s" => \$pr,
	"n:i" => \$n,
	"Rscript:s" => \$Rscript,
	"convert:s" => \$convert,
	"samtools:s" => \$samtools
);
die <<USAGE unless $bam and $len and $seqtype and $pr;
Usage: Draw reads distrubution randomness in genes.
Options:
	*-bam        <str>    bam file(s), separated by comma ","
	*-len        <str>    gene length information, format: geneID <tab> length
	*-seqType    <str>    sequening type, PE or SE
	*-prefix     <str>    result prefix, the result file will be named as 'prefix.ReadsRandom.png'
	 -n          <int>    how many bars in the figure, the more the accurater, default: 200
	 -Rscript    <app>    Rscript path
	 -convert    <app>    convert path
	 -samtools   <app>    samtools path
Author: Yu Mingqi \@ 2014-11-07
USAGE
$n ||= 200;
$Rscript ||= "unset R_LIBS_USER && unset R_LIBS && \\\n unset R_LIBS_SITE && \\\n export R_LIBS_USER=\"/usr/local/envs/r_3_5_1/lib/R/library\" &&/usr/local/bin/Rscript";
$convert ||= "/usr/bin/convert";
$samtools ||= "/usr/local/RSEM-1.3.1/RSEM/1.3.1/bin/samtools-1.3/samtools";
my $name = basename $pr;
my $outdir = dirname $pr;
mkpath $outdir unless -e $outdir;
$pr = $pr . '.ReadsRandom';
$seqtype = 1 if $seqtype eq 'PE';
my %g2len = map {($a, $b) = split; $a => $b} `cat $len`;

my $oldID = 'for no warning';
my @counter = (0) x $n;
my @block;
sub Algorithm {
	my $ary = shift;
	my $readsNumber = @$ary;
	for (@$ary) {
		my ($geneID, $start, $seq) = (split)[2, 3, 9];
		my $glen = $g2len{$geneID};
		my $mlen = length $seq; # mapped read length
		my $end = $start + $mlen -1;
		$end = $glen if $end > $glen;
		my $increment = 1 / $mlen / $readsNumber; # for every mapped bp
		my $coefficient = $n / $glen; # every bp = ? windows length
		if ($coefficient <= 1) {
			$counter[int ($_ * $coefficient)] += $increment for $start .. $end
		} else {
			my $add = $increment / $coefficient; # increment of each window
				my $relative_start = ($start - 1) * $coefficient;
				my $relative_end = ($end - 1) * $coefficient;
				$counter[int $relative_start] += $add * (int $relative_start + 1 - $relative_start);
				$counter[int $relative_end] += $add * ($relative_end - int $relative_end);
				$counter[$_] += $add for int $start + 1 .. int $end - 1
		}
	}
}
sub CountLastBlock {
	my (@first, @second, @single);
	for (@block) { # filter & separate bam to blocks
		my ($flag, $geneID, $start) = (split)[1, 2, 3];
		next if $flag & 0x0004; # unmapped
		warn "Cannot find length info. of $geneID\n" and next unless my $glen = $g2len{$geneID};
		warn "Start $start > length of $geneID, check length info.\n" and next if $start > $glen;
		if ($seqtype) {
			if ($flag & 0x0040) { push @first, $_ } else { push @second, $_ }
		} else {
			push @single, $_;
		}
	}
	if ($seqtype) {
		&Algorithm(\@first);
		&Algorithm(\@second);
	} else {
		&Algorithm(\@single);
	}
}
my @files = split /,/, $bam;
$bam = "$pr.merge.bam" and `$samtools merge $bam @files` unless @files == 1;
open SAM, $bam =~ /\.sam$/ ? $bam : "$samtools view $bam |" or die "Cannot open $bam\n";
while (<SAM>) {
	my $readID = (split)[0];
	if ($readID ne $oldID) {
		&CountLastBlock();
		$oldID = $readID;
		@block = ();
	}
	push @block, $_;
}
close SAM;
&CountLastBlock(); # the last readID
`rm $bam` unless @files == 1;

my $print;
$print .= "$_\t$counter[$_]\n" for 0 .. $#counter;
my $R = <<CODE;
library(ggplot2)
pdf('$pr.pdf', width=8, height=8)
d=read.table(text='$print')
d\\\$V1=d\\\$V1/$n
ggplot(d, aes(V1,V2))+geom_line(colour='royalblue', size=1)+
	labs(x='Relative Position in Genes (5\\\' -> 3\\\') ($n windows)', y='Reads Number of Each Window')+
	theme(axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),plot.title=element_text(face='bold', size=15),axis.text=element_text(color='black'),panel.background = element_rect(fill='transparent'),panel.grid=element_line(color='grey'),panel.border=element_rect(fill='transparent',color='black'),axis.title=element_text(size=15))
dev.off()
CODE
`echo "$R" > $pr.R && $Rscript $pr.R && $convert -density 1000 -resize 10% $pr.pdf $pr.png`;
`rm $pr.R`;
