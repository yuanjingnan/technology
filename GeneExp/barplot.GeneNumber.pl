#!/usr/bin/perl
use strict;
use Getopt::Long;
use FindBin '$Bin';
#use File::Path 'mkpath';

my ($summary,$outdir,$Rscript, $convert);
GetOptions(
	"summary:s" => \$summary,
	"outdir:s" => \$outdir,
	"Rscript:s" => \$Rscript,
	"convert:s" => \$convert
);
die <<USAGE unless ( $summary && $outdir );
Draw barplot for samples' known gene number and novel gene number
	*-summary	<file>	GeneExpressionSummary.xls
	*-outdir	<path>	output directory path
	-Rscript	<path>	Rscript path
	-convert	<path>	convert path
perl $0 -summary GeneExpressionSummary.xls -outdir output_directory
USAGE

$Rscript ||= "unset R_LIBS_USER && \\\n unset R_LIBS && \\\n unset R_LIBS_SITE && \\\n export R_LIBS_USER=\"/ldfssz1/ST_CANCER/CGR/SHARE/tools/Anaconda/anaconda2/envs/r-3.5.1/lib/R/library:\$R_LIBS_USER\" && \\\n$Bin/../Bin_CentOS6/R/Rscript";
$convert ||= "$Bin/../Bin_CentOS6/convert";

my $sampleNum=`cat $summary | cut -f 1 | sort -u | wc -l`;
my $width=6;
my $ylim=1.45;
if ($sampleNum>8)
{
	$width=$sampleNum;
}


open R, ">$outdir/knowngene_novelgene_number_barplot.R" or dir $!;
print R <<RCODE;
pdf("$outdir/knowngene_novelgene_number_barplot.pdf",pointsize=14,height=6,width=$width)
par(mar=c(4,5,5,3)+0.1,xpd=TRUE)
data <- read.table("$summary",sep="\t",skip=1,row.names=1)
novel_col="forestgreen"
known_col="dodgerblue3"
ylim <- max(data\$V6+data\$V7) * $ylim
barplot(matrix(rbind(data\$V6,data\$V7),nrow=2),beside=F,col=c(known_col,novel_col),space=1,ylim=c(0,ylim),cex.axis=1,main="Number of Known Genes and Novel Genes",cex.main=1.2,ylab="Gene Number",cex.lab=1.2,)
end_point <- (length(rownames(data)) - 1) * 2
text(seq(1.5,end_point+1.5,by=2),max(data\$V6+data\$V7),pos=3,labels=data\$V6,cex=1.1,col=known_col)
text(seq(1.5,end_point+1.5,by=2),max(data\$V6+data\$V7)*1.06,pos=3,labels=data\$V7,cex=1.1,col=novel_col)
text(seq(1.8,end_point+1.8,by=2),0,srt=25,adj=c(1,2),xpd=TRUE,labels=rownames(data),cex=1.2)
legend('topright', c("Novel genes","Known genes"),col=c(novel_col,known_col),pch=15, pt.cex=1.2,bty = "n", cex=1.2)
box()
dev.off()
RCODE
close R;
`$Rscript $outdir/knowngene_novelgene_number_barplot.R && rm $outdir/knowngene_novelgene_number_barplot.R && $convert -density 300 $outdir/knowngene_novelgene_number_barplot.pdf $outdir/knowngene_novelgene_number_barplot.png`;





