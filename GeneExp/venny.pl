#!/usr/bin/perl

#Author: liaoxinhui@genomics.org.cn
#Date:2013-04-11 

use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';

my ($infile, $header, $name, $font_family, $name_font_size, $num_font_size, $colour, $desc, $outdir, $imgname, $help);

GetOptions(
        "infile:s" => \$infile,
	"header!" =>\$header,
        "name:s" => \$name,
	"desc:s" => \$desc,
	"font_family:s" => \$font_family,
	"name_font_size:f" => \$name_font_size,
	"num_font_size:f" => \$num_font_size,
        "outdir:s" => \$outdir,
	"imgname:s" => \$imgname,
        "help|?" => \$help
);
if (!defined $infile || defined $help) {
	die <<USAGE

Description: Generate tables & graph for Venny (2D,3D,4D,5D)
Usage: perl $0 [options]
Options:
	* -infile			input files with IDs in column 1, separated by ","
	* -name				name of each category, separated by ","
	  -desc				description file for each ID, format: GeneID, description1, description2... separated by "\\t" (header needed)
	  -header			only active when files with header
	  -font_family			font family, default is Verdana
	  -name_font_size		font size of name, default is 18
	  -num_font_size		font size of numbers, default is 13
	  -outdir			directory of output files, default is current directory
	  -imgname          output image name, default is "Venn-[input.number]D"
	  -help 			help information
E.g.:
	perl $0 -infile A.xls,B.xls.C.xls,D.xls -name A,B,C,D -color

USAGE
}
my @files = split /,+/,$infile;
my @names = split /,+/,$name;
die "\nWarning: Files' number not equal to Names' number!\nWarning: Please double check the parameters!\n\n" if (@files != @names);
$font_family ||= "Verdana";
$name_font_size ||= 18;
$num_font_size ||= 13;
$outdir ||= "./";
`mkdir -p $outdir` unless (-d $outdir);
$imgname ||= "Venn-".scalar(@names)."D";
#read ID files
my %all_ids;
foreach my $index (0..$#files) {
	open INFILE,$files[$index];
	<INFILE> if defined $header;
	my %check_dup = ();
	while (<INFILE>) {
		chomp;my @temp = split /\s+/;
		next if ($check_dup{$temp[0]});
		$check_dup{$temp[0]} = 1;
		$all_ids{$temp[0]} .= "$names[$index]_";
	}
	close INFILE;
}
#read description files
my %desc; my $desc_header;
if (defined $desc) {
	open DESC,$desc;
	my $head = <DESC>;
	my @head = split /\t+/,$head;shift @head; $desc_header = join "\t",@head;
	while (<DESC>) {
		chomp;my @temp = split /\t+/,$_;
		my $id = shift @temp;
		$desc{$id} = join "\t",@temp;
	}
	close DESC;
}
#grouping & print
my %all_groups;
foreach my $id (keys %all_ids) {
	chop $all_ids{$id};
	my $id_desc = $desc{$id} ? $desc{$id} : "NULL";
	$all_groups{$all_ids{$id}}{'ids'} .= (defined $desc) ? "$id\t$id_desc\n" : "$id\n";
	$all_groups{$all_ids{$id}}{'num'} += 1;
}
foreach my $group (keys %all_groups) {
	open OUTFILE,">$outdir/$group.List.xls";
	print OUTFILE (defined $desc) ? "GeneID\t$desc_header" : "GeneID\n";
	print OUTFILE $all_groups{$group}{'ids'};
	close OUTFILE;
}
#draw graph
&draw_venn2(\%all_groups,\@names) if @names == 2;
&draw_venn3(\%all_groups,\@names) if @names == 3;
&draw_venn4(\%all_groups,\@names) if @names == 4;
&draw_venn5(\%all_groups,\@names) if @names == 5;

system("$Bin/../Bin_CentOS6/R/Rscript $outdir/$imgname.R 2>/dev/null");
system("$Bin/../Bin_CentOS6/convert $outdir/$imgname.pdf $outdir/$imgname.png");

exit 0;
#sub programs
sub draw_venn2 {
	my ($groups, $names) = @_;
	my ($g12,$g1,$g2) = (
		$$groups{"$$names[0]_$$names[1]"} ? $$groups{"$$names[0]_$$names[1]"}{'num'} : 0,
		$$groups{"$$names[0]"} ? $$groups{"$$names[0]"}{'num'} : 0,
		$$groups{"$$names[1]"} ? $$groups{"$$names[1]"}{'num'} : 0
	);
	my ($t1,$t2);
	$t1=$g1+$g12;
	$t2=$g2+$g12;
	open CMD,">$outdir/$imgname.R";
	print CMD <<CMD;
library(plotrix)
pdf(file="$outdir/$imgname.pdf",width=8.8,height=8)
color <- c("#E41A1C","#377EB8")
color_transparent <- adjustcolor(color, alpha.f = 0.2) 
color_transparent1 <- adjustcolor(color, alpha.f = 1)
########################################
p_x   <- c(-13,0,13)
p_y   <- c(0,0,0)
p_lab <- c($g1,$g12,$g2)

title_x <- c(-16,16)
title_y <- c(-16,-16)
title_lab <- c("$$names[0]\\n($t1)","$$names[1]\\n($t2)")
########################################
par(mar=c(7,10,7,8)+0.1,xpd=TRUE)
plot(c(-18,18), c(-18,18), type="n",xaxt = "n", xlab="",ylab="",yaxt = "n", axes=F,main="")
draw.ellipse(c(-5,5), c(0,0), c(14,14), c(14,14),border=color_transparent1,
angle=c(45,45), lty=1,col = color_transparent,lwd = 2)
text (p_x,p_y,p_lab,cex=1,col="grey20")
text (title_x,title_y,title_lab,cex=1.2,col=color_transparent1)
dev.off()
CMD
}
#system("$Bin/../software/Rscript $outdir/$imgname.R 2>/dev/null");
#system("$Bin/../software/convert $outdir/$imgname.pdf $outdir/$imgname.png");


sub draw_venn3 {
	my ($groups, $names) = @_;
	my ($g123,$g12,$g13,$g23,$g1,$g2,$g3) = (
		$$groups{"$$names[0]_$$names[1]_$$names[2]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]"} ? $$groups{"$$names[0]_$$names[1]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]"} ? $$groups{"$$names[0]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]"} ? $$groups{"$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]"} ? $$groups{"$$names[0]"}{'num'} : 0,
		$$groups{"$$names[1]"} ? $$groups{"$$names[1]"}{'num'} : 0,
		$$groups{"$$names[2]"} ? $$groups{"$$names[2]"}{'num'} : 0
	);
	my ($t1,$t2,$t3);
	$t1=$g1+$g12+$g13+$g123;
	$t2=$g2+$g12+$g23+$g123;
	$t3=$g3+$g23+$g13+$g123;
	open CMD,">$outdir/$imgname.R";
print CMD <<CMD;
library(plotrix)
pdf(file="$outdir/$imgname.pdf",width=8.8,height=8)
color <- c("#E41A1C","#377EB8","#FDB462")
color_transparent <- adjustcolor(color, alpha.f = 0.2) 
color_transparent1 <- adjustcolor(color, alpha.f = 1)
########################################
p_x   <- c(0,13,-13,  10.4,-10.4,0, 0)
p_y   <- c(13,-9,-9,  4,4,-13.5, -1)
p_lab <- c($g1,$g2,$g3,$g12,$g13,$g23,$g123)

title_x <- c(0,17,-17)
title_y <- c(19.8,-18,-17)
title_lab <- c("$$names[0]\\n($t1)","$$names[1]\\n($t2)","$$names[2]\\n($t3)")
########################################
par(mar=c(7,10,7,8)+0.1,xpd=TRUE)
plot(c(-18,18), c(-18,18), type="n",,xaxt = "n", xlab="",ylab="",yaxt = "n", axes=F,main="")
draw.ellipse(c(0,4,-4), c(3.02,-3.912,-3.912), c(14,14,14), c(14,14,14),border=color_transparent1,
angle=c(60,120,0), lty=1,col = color_transparent,lwd = 2)
text (p_x,p_y,p_lab,cex=1,col="grey20")
text (title_x,title_y,title_lab,cex=1.2,col=color_transparent1)
dev.off()
CMD
#system("$Bin/../software/Rscript $outdir/$imgname.R 2>/dev/null");
#system("$Bin/../software/convert $outdir/$imgname.pdf $outdir/$imgname.png");
}

sub draw_venn4 {
	my ($groups, $names) = @_;
	my ($g1234,$g123,$g124,$g134,$g234,$g12,$g13,$g14,$g23,$g24,$g34,$g1,$g2,$g3,$g4) = (
		$$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[2]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[3]"} ? $$groups{"$$names[0]_$$names[1]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]_$$names[3]"} ? $$groups{"$$names[0]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]_$$names[3]"} ? $$groups{"$$names[1]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]"} ? $$groups{"$$names[0]_$$names[1]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]"} ? $$groups{"$$names[0]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[3]"} ? $$groups{"$$names[0]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]"} ? $$groups{"$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[3]"} ? $$groups{"$$names[1]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[2]_$$names[3]"} ? $$groups{"$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]"} ? $$groups{"$$names[0]"}{'num'} : 0,
		$$groups{"$$names[1]"} ? $$groups{"$$names[1]"}{'num'} : 0,
		$$groups{"$$names[2]"} ? $$groups{"$$names[2]"}{'num'} : 0,
		$$groups{"$$names[3]"} ? $$groups{"$$names[3]"}{'num'} : 0
	);
	my ($t1,$t2,$t3,$t4);
	$t1=$g1+$g12+$g13+$g14+$g123+$g124+$g134+$g1234;
	$t2=$g2+$g12+$g23+$g24+$g123+$g124+$g234+$g1234;
	$t3=$g3+$g13+$g23+$g34+$g123+$g134+$g234+$g1234;
	$t4=$g4+$g14+$g24+$g34+$g124+$g134+$g234+$g1234;
	open CMD,">$outdir/$imgname.R";
print CMD <<CMD;
library(plotrix)
pdf(file="$outdir/$imgname.pdf",width=8.8,height=8)
color <- c("#E41A1C","#377EB8","#FDB462","#4DAF4A") 
color_transparent <- adjustcolor(color, alpha.f = 0.2) 
color_transparent1 <- adjustcolor(color, alpha.f = 1)

p_x   <- c(-15,15,-7.3,7.3,  0,-10,-9,10,9,0,  4.9,-4.9,-6.5,6.5 ,0)
p_y   <- c(2,2,11,11,  -12,6.3,-4,6.3,-4,8.8,  -7.4,-7.4,2.6,2.6 ,-2)
p_lab <- c ($g1,$g2,$g3,$g4,$g12,$g13,$g14,$g23,$g24,$g34,$g123,$g124,$g134,$g234,$g1234)

title_x <- c(-18,18,-15,15)
title_y <- c(-16,-16,16,16)
title_lab <- c("$$names[0]\\n($t1)","$$names[1]\\n($t2)","$$names[2]\\n($t3)","$$names[3]\\n($t4)")
#########################################
par(mar=c(7,10,7,8)+0.1,xpd=TRUE)
plot(c(-18,18), c(-18,18), type="n",xaxt = "n", xlab="",ylab="",yaxt = "n", axes=F,main="")
draw.ellipse(c(-6,6,0,0), c(-4,-4,2,2), c(16,16,14,14), c(9.89,9.89,8.65,8.65),border=color_transparent1,
angle=c(-40,40,-45,45), lty=1,col = color_transparent,lwd = 2)
text (p_x,p_y,p_lab,cex=1,col="grey20")
text (title_x,title_y,title_lab,cex=1.2,col=color_transparent1)
dev.off()
CMD
#system("$Bin/../software/Rscript $outdir/$imgname.R 2>/dev/null");
#system("$Bin/../software/convert $outdir/$imgname.pdf $outdir/$imgname.png");
}

sub draw_venn5 {
	my ($groups, $names) = @_;
	my ($g12345,$g1234,$g1235,$g1245,$g1345,$g2345,$g123,$g124,$g125,$g134,$g135,$g145,$g234,$g235,$g245,$g345,$g12,$g13,$g14,$g15,$g23,$g24,$g25,$g34,$g35,$g45,$g1,$g2,$g3,$g4,$g5) =(
		$$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]_$$names[4]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[4]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[3]_$$names[4]"} ? $$groups{"$$names[0]_$$names[1]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]_$$names[3]_$$names[4]"} ? $$groups{"$$names[0]_$$names[2]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]_$$names[3]_$$names[4]"} ? $$groups{"$$names[1]_$$names[2]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[2]"} ? $$groups{"$$names[0]_$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[3]"} ? $$groups{"$$names[0]_$$names[1]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]_$$names[4]"} ? $$groups{"$$names[0]_$$names[1]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]_$$names[3]"} ? $$groups{"$$names[0]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]_$$names[4]"} ? $$groups{"$$names[0]_$$names[2]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[3]_$$names[4]"} ? $$groups{"$$names[0]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]_$$names[3]"} ? $$groups{"$$names[1]_$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]_$$names[4]"} ? $$groups{"$$names[1]_$$names[2]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[3]_$$names[4]"} ? $$groups{"$$names[1]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[2]_$$names[3]_$$names[4]"} ? $$groups{"$$names[2]_$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[1]"} ? $$groups{"$$names[0]_$$names[1]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[2]"} ? $$groups{"$$names[0]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[3]"} ? $$groups{"$$names[0]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[0]_$$names[4]"} ? $$groups{"$$names[0]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[2]"} ? $$groups{"$$names[1]_$$names[2]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[3]"} ? $$groups{"$$names[1]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[1]_$$names[4]"} ? $$groups{"$$names[1]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[2]_$$names[3]"} ? $$groups{"$$names[2]_$$names[3]"}{'num'} : 0,
		$$groups{"$$names[2]_$$names[4]"} ? $$groups{"$$names[2]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[3]_$$names[4]"} ? $$groups{"$$names[3]_$$names[4]"}{'num'} : 0,
		$$groups{"$$names[0]"} ? $$groups{"$$names[0]"}{'num'} : 0,
		$$groups{"$$names[1]"} ? $$groups{"$$names[1]"}{'num'} : 0,
		$$groups{"$$names[2]"} ? $$groups{"$$names[2]"}{'num'} : 0,
		$$groups{"$$names[3]"} ? $$groups{"$$names[3]"}{'num'} : 0,
		$$groups{"$$names[4]"} ? $$groups{"$$names[4]"}{'num'} : 0,
	);
	my ($t1,$t2,$t3,$t4,$t5);
	$t1=$g1+$g12+$g13+$g14+$g15+$g123+$g124+$g125+$g134+$g135+$g145+$g12345;
	$t2=$g2+$g12+$g23+$g24+$g25+$g123+$g124+$g125+$g234+$g235+$g245+$g12345;
	$t3=$g3+$g13+$g23+$g34+$g35+$g123+$g134+$g135+$g234+$g235+$g345+$g12345;
	$t4=$g4+$g14+$g24+$g34+$g45+$g124+$g134+$g145+$g234+$g145+$g345+$g12345;
	$t5=$g5+$g15+$g25+$g35+$g45+$g125+$g135+$g145+$g235+$g245+$g345+$g12345;
	open CMD,">$outdir/$imgname.R";
print CMD <<CMD;
library(plotrix)
pdf(file="$outdir/$imgname.pdf",width=8.8,height=8)
par(mar=c(6,5,4,11)+0.1,xpd=TRUE)
color <- c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FDB462") 
color_transparent <- adjustcolor(color, alpha.f = 0.2) 
color_transparent1 <- adjustcolor(color, alpha.f = 1)
########################################################
p_x   <- c(-4.5,14,13,-6,-18,  7.2,7.2,0,-9,13.2,-7.5,12.05,1.5,-12,-12.5,  9,3.7,9.1,4,-10.5,-4,-3.8,12.2,-10,-11.5,  2.2,9.9,3.7,-7.8,-8.5,  0)
p_y   <- c(17,9.5,-10.5,-15.5,1,  11,-10.8,12.7,10,-3.4,-10.8,4,-13.5,3.5,-5,  -7.2,12.1,6.5,-11.7,7.4,10.3,-11,-0.1,-7.2,0.3,  -10.1,-1,9,6.3,-5.8,  0)
p_lab <- c ($g1,$g2,$g3,$g4,$g5,$g12,$g13,$g14,$g15,$g23,$g24,$g25,$g34,$g35,$g45,$g123,$g124,$g125,$g134,$g135,$g145,$g234,$g235,$g245,$g345,$g1234,$g1235,$g1245,$g1345,$g2345,$g12345)

############################################
title_x <-c(-16,17,22,-14,-23)
title_y <-c(21,18.5,-18,-23,-9)
title_lab <- c("$$names[0]\\n($t1)","$$names[1]\\n($t2)","$$names[2]\\n($t3)","$$names[3]\\n($t4)","$$names[4]\\n($t5)")
########################################################
par(mar=c(7,10,7,8)+0.1,xpd=TRUE)
plot(c(-18,18), c(-18,18), type="n",xaxt = "n", xlab="",ylab="",yaxt = "n", axes=F,main="")
draw.ellipse(c(0,4,2.4,-2.4,-4), c(4.35,1.1,-3.55,-3.55,1.1), c(18,18,18,18,18), c(10.5,10.5,10.5,10.5,10.5),border=color_transparent1,
 angle=c(290,218,146,74,2), lty=1,col = color_transparent,lwd = 2)
text (p_x,p_y,p_lab,cex=1,col="grey20")
text (title_x,title_y,title_lab,cex=1.2,col=color_transparent1)
dev.off()
CMD
#system("$Bin/../software/Rscript $outdir/$imgname.R 2>/dev/null");
#system("$Bin/../software/convert $outdir/$imgname.pdf $outdir/$imgname.png");
}
