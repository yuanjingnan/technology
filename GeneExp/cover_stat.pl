#!/usr/bin/perl
#Author: liaoxinhui@genomics.org.cn
#Date:  Tue Nov 10 15:39:36 CST 2015


use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';
use File::Basename;

die "usage: perl $0 [depth.bed] [transcript_len] [outprefix] [outdir]\n" unless (@ARGV == 4);

my $depth = shift;
my $trlen = shift;
my $prefix = shift;
my $outdir = shift;

my %len;
open LEN,$trlen or die $!;
while (<LEN>) {
	chomp; my @a = split /\s+/;
	$len{$a[0]} = $a[1];
}
close LEN;

my %tr;
open DEP,$depth or die $!;
while (<DEP>) {
	chomp;
	my @a = split /\s+/;
	next if ($a[2] == 0);
	$tr{$a[0]}++;
}
close DEP;

my $l = "";
my %stat;
foreach my $t (keys %tr) {
		$tr{$t} ||= 0;
		my $ptg = sprintf("%.2f", $tr{$t}/$len{$t}*100);
		my $transform = $ptg/10;
		$l .= "$transform,";
		if ($ptg <= 10) {
			$stat{'0-10'}++;
		}elsif ($ptg <= 20) {
			$stat{'10-20'}++;
		}elsif ($ptg <= 30) {
			$stat{'20-30'}++;
		}elsif ($ptg <= 40) {
			$stat{'30-40'}++;
		}elsif ($ptg <= 50) {
			$stat{'40-50'}++;
		}elsif ($ptg <= 60) {
			$stat{'50-60'}++;
		}elsif ($ptg <= 70) {
			$stat{'60-70'}++;
		}elsif ($ptg <= 80) {
			$stat{'70-80'}++;
		}elsif ($ptg <= 90) {
			$stat{'80-90'}++;
		}else {
			$stat{'90-100'}++;
		}
}
chop $l;

my $total_trans = scalar(keys %tr);
my $c = "";
foreach my $ptg ('0-10', '10-20', '20-30', '30-40', '40-50', '50-60', '60-70', '70-80', '80-90', '90-100') {
	$stat{$ptg} ||= 0;
	$stat{$ptg} = sprintf("%.2f", $stat{$ptg}/$total_trans*100);
	$c .= "$stat{$ptg},";
}
chop $c;

open R,">$outdir/draw_coverage.R" or die $!;
print R <<CODE;
data<-c($c)
data<-data
data <- matrix(data,nrow=1)
colnames(data) <- c('0-10','10-20','20-30','30-40','40-50','50-60','60-70','70-80','80-90','90-100')
den <- density(c($l))
pdf("$outdir/$prefix.ReadsCoverage.pdf",width=10,height=8)
par(omi=c(0.165,0.3,0,1.1))
barplot(data,space=0.3,col="#377EB8",xlab="Percent covered",axes=F,ylab="",cex.lab=1.5,border=0,width=0.5)
axis(2,col="#377EB8",col.axis="#377EB8",las=2)
mtext(side=2,"Percentage of transcripts",line=3,cex=1.5,col="#377EB8")
par(new=T)
par(omi=c(0,0.3,0,1.1))
plot(den,yaxt="n",ylab="",xaxt="n",axes=F,xlab="",main="")
lines(den,lty=1,col="grey40",lwd=2)
axis(4,col="grey40",col.axis="grey40",las=2)
mtext(side=4,"Density",line=5,cex=1.5,col="grey40")
dev.off()
CODE
close R;
system("unset R_LIBS_USER && \\\n unset R_LIBS && \\\n unset R_LIBS_SITE && \\\n export R_LIBS_USER=\"/ldfssz1/ST_CANCER/CGR/SHARE/tools/Anaconda/anaconda2/envs/r-3.5.1/lib/R/library:\$R_LIBS_USER\" && \\\n$Bin/../Bin_CentOS6/R/Rscript $outdir/draw_coverage.R");
system("$Bin/../Bin_CentOS6/convert -density 300 -resize 25% $outdir/$prefix.ReadsCoverage.pdf $outdir/$prefix.ReadsCoverage.png");
