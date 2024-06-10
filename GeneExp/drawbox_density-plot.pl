#!/usr/bin/perl -w 
use strict;
use Cwd 'abs_path';
use Getopt::Long;
use FindBin qw($Bin $Script);
use File::Basename qw(basename dirname);
use List::Util qw/max min sum/;

my ($Indir, $Outdir, $Help, $Rbin, $Convert);
GetOptions(
	"i:s"=>\$Indir,
	"o:s"=>\$Outdir,
	"help"=>\$Help,
	"R"=>\$Rbin,
	"convert:s" => \$Convert
);
die "erro:\nuse: perl $0 -i fpkmdir -o outdir\n" if (!$Indir || $Help);
$Outdir ||= "./";
$Rbin = "unset R_LIBS_USER && \\\n unset R_LIBS && \\\n unset R_LIBS_SITE && \\\n export R_LIBS_USER=\"/ldfssz1/ST_CANCER/CGR/SHARE/tools/Anaconda/anaconda2/envs/r-3.5.1/lib/R/library:\$R_LIBS_USER\" && \\\n$Bin/../Bin_CentOS6/R" unless ($Rbin and -d $Rbin);
$Convert ||= "$Bin/../Bin_CentOS6/convert";

#######################generate expression file##############
my @samples = glob("$Indir/*.gene.fpkm.xls");
my $count = scalar (@samples);
my %hash;
my @sample;
open OUT,">$Outdir/boxplot.matrix";
print OUT "log10FPKM\tSample\n";
for my $s (sort @samples){
	my ($name) = $s =~ /.*\/(.*).gene.fpkm.xls/;
	push @sample, $name;
	open IN,$s or die $!;
	<IN>;
	while(<IN>){
		chomp();
		my @a=split /\t+/,$_;
		print OUT log($a[4])/log(10)."\t$name\n"
	}
	close IN;
}
close OUT;

my $colors='"#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999","#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F","#E5C494","#B3B3B3","#8DD3C7","#FFFFB3","#BEBADA","#FB8072","#80B1D3","#FDB462","#B3DE69","#FCCDE5","#D9D9D9","#BC80BD","#CCEBC5","#FBB4AE","#B3CDE3","#CCEBC5","#DECBE4","#FED9A6","#FFFFCC","#E5D8BD","#FDDAEC","#F2F2F2","#B3E2CD","#FDCDAC","#CBD5E8","#F4CAE4","#E6F5C9","#FFF2AE","#F1E2CC","#CCCCCC"';
my @colors;
my $sample_count=scalar(@sample);
my $xsize=($sample_count>5)?"96/len":17;
my $iflegend=($sample_count>20)?"+guides(fill=FALSE)":"";
#######################generate R script#####################
open RSC,"> $Outdir/box.R" or die $!;
print RSC <<EOF;
pdf("$Outdir/Boxplot.pdf",width=10,height=10)
library(ggplot2)
library(RColorBrewer)
options(stringsAsFactors=F)
data <- read.table("$Outdir/boxplot.matrix",header=T,check.names=F)
#attach(data)
data\$Sample<-as.character(data\$Sample)
len <- length(unique(data\$Sample))
a<-c(as.character(data\$Sample))
a = as.factor(ifelse(is.na(a), "NA", a))
# colours<-c(rep(c($colors),20)) # changed by wumeizhen
cols<-brewer.pal(9, "YlOrRd")
pal<-colorRampPalette(cols)
colours<-pal($sample_count)
ggplot(data, aes(x=Sample, y=log10FPKM)) + geom_boxplot(aes(fill=factor(a))) +
theme_bw()+
#theme(panel.grid.major=element_line(colour=NA))+
scale_fill_manual(values=colours,name="Sample",na.value=colours[length(unique(data\$Sample))]) +theme(axis.text.x = element_text(colour="grey20",size=11,angle=30, hjust =1)) +
labs(x="")+ #去掉x轴的sample
guides(fill = guide_legend(title = "Sample",keywidth = 2, keyheight = 2))+
#theme(axis.title.y = element_text(face="bold",size = 18))+
theme(axis.title.y = element_text(size = 18))+
theme(axis.text.y = element_text(size=16,hjust =1))+
#theme(axis.text.x = element_text(face="bold",size=$xsize,angle=30, hjust =1),panel.grid.major=element_line(colour=NA),panel.grid =element_blank())+
theme(axis.text.x = element_text(size=$xsize,angle=30, hjust =1))+
theme(legend.title=element_text(size=15))+
theme(legend.text=element_text(size=15)) $iflegend
dev.off()
EOF
close RSC;

system("$Rbin/Rscript $Outdir/box.R && rm $Outdir/box.R");
#system("$Rbin/Rscript $Outdir/box.R");
system("$Convert -density 300 -resize 40% $Outdir/Boxplot.pdf $Outdir/Boxplot.png");

open RSC,"> $Outdir/Density.R" or die $!;
print RSC <<EOF;
options(stringsAsFactors=F)
library(ggplot2)
library("grid")
data <- read.table("$Outdir/boxplot.matrix",header=T,check.names=F)
a<-c(as.character(data\$Sample))
a <- as.factor(ifelse(is.na(a), "NA", a))


attach(data)
pdf("$Outdir/Density.pdf",width=10,height=10)
#ggplot(data,aes(x=log10FPKM, fill=Sample,colours="black")) + geom_density(alpha=0.4) + guides(fill = guide_legend(keywidth = 0.8, keyheight = 0.8)) + theme(panel.background=element_rect(fill=NA,colour="grey"), panel.grid=element_line(color='grey'), panel.border=element_rect(fill='transparent',color='black'), legend.title=element_blank(), legend.text = element_text(size = 8), plot.title=element_text(face='bold', size=20), axis.title=element_text(size=20), axis.text.x=element_text(color='black', size=15), axis.text.y=element_text(color='black', size=15)) + labs(x='log10FPKM', y='Density',title = "Expression Density Distribution")
ggplot(data,aes(x=log10FPKM, fill=factor(a),colours="black")) + geom_density(alpha=0.4)+theme(panel.background=element_rect(fill=NA,colour="grey"), panel.grid=element_line(color='grey'), panel.border=element_rect(fill='transparent',color='black'),plot.title=element_text(face='bold', size=20),axis.title=element_text(size=20), axis.text.x=element_text(color='black', size=15),axis.text.y=element_text(color='black', size=15))+labs(x='log10FPKM', y='Density',title = "Expression Density Distribution")+theme(legend.key.size= unit(0.4, "cm"))
dev.off()
EOF
close RSC;
system("$Rbin/Rscript $Outdir/Density.R && rm $Outdir/Density.R && rm $Outdir/boxplot.matrix");
#system("$Rbin/Rscript $Outdir/Density.R ");
system("$Convert -density 300 -resize 40% $Outdir/Density.pdf $Outdir/Density.png");
exit;
