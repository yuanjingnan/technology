#!/usr/bin/perl
#Author: liaoxinhui@genomics.org.cn
#Date:  Mon Oct 12 16:57:37 CST 2015

use strict;
use warnings;
use Getopt::Long;
use FindBin '$Bin';
use Cwd 'abs_path';
use File::Basename;
use lib "$Bin/../self_perllib";
use CommonMethods; 

my ($conf,$list);
GetOptions (
	"conf:s" => \$conf,
	"list:s" => \$list
);

if (!$conf || !$list) {
	print STDERR <<USAGE;
=============================================================================
Descriptions: Generate shell scripts for gene expression analysis
Usage:
	perl $0 [options]
Options:
	* -conf          input software parameters and file paths info
	* -list          input clean data file list
E.g.:
	perl $0 -conf para.conf -list clean.list
=============================================================================
USAGE
	exit;
}
my %para = ();
$para{'SeqType'}="PE";
my %sample = ();
open CONF,$conf or die $!;
while (<CONF>) {
	next if (/^\s*$/ || /^\s*#/ || !/=/);
	chomp; s/^([^#]+)#.*$/$1/;
	my ($key,$value) = split /=/,$_,2;
	$key =~ s/^\s*(.*?)\s*$/$1/;
	$value =~ s/^\s*(.*?)\s*$/$1/;
	next if ($value =~ /^\s*$/);
	$para{$key} = $value;
}
close CONF;

my ($Option,$seqType);
my $ifPE = "";
my @a;
open LIST,$list or die $!;
while (<LIST>) {
	next if (/^\s*$/ || /^\s*#/);
	chomp; my @t = split /\s+/;
	@a = split /,+/, $t[1];
	if(@a > 1){
		$sample{$t[0]} = [$a[0], $a[1], $t[3]]; # fq1, fq2, depends
		$Option = $para{'Bowtie2_Options'};
		$seqType ='PE';
		$ifPE = "--paired-end";
	}else{
		$sample{$t[0]} = [$t[1], $t[3]];
		$Option = $para{'Bowtie2_Options'};
		$seqType ='SE';
	}
}
close LIST;

for('SeqType','Gene_Fasta','Gene2Tr','RSEM_Options',
	'Shell_Dir','Result_Dir','Cp_Dir',
	'Depend_List','DownStream_List') {
	 unless ($para{$_}) {ThrowError("ERROR: Can't find parameter [$_]!\n",*STDERR)}
}

open DEPL,">$para{'Depend_List'}" or die $!;
open DOWL,">$para{'DownStream_List'}" or die $!;
open RES,">$para{'Result_List'}" or die $!;

my @depends = ();
my @depends2 = ();
my %fpkmlist = ();
my $shdir = Mkdir("$para{'Shell_Dir'}");
my $outdir = Mkdir("$para{'Result_Dir'}");
my $fpkmdir = Mkdir("$para{'Cp_Dir'}/GeneExpression");
my $mapdir = Mkdir("$para{'GeneMap_Dir'}");
my $randomdir = Mkdir("$para{'Cp_Dir'}/ReadsRandom");
my $coverdir = Mkdir("$para{'Cp_Dir'}/ReadsCoverage");
my $saturationdir =Mkdir("$para{'Cp_Dir'}/SeqSaturation");

###software
my $perl="$Bin/../Bin_CentOS6/perl";
my $python="$Bin/../Bin_CentOS6/python";

my $g2id='';
my $refer = "";
my $gene2tr="";
my $Transcriptbed="";
my $Transcripttxt="";
my $exp_shell ="";
my @p = split(/,/,$para{'Gene_Fasta'});
if (@p>1){
$exp_shell .= "if [ ! -d $outdir/rsem-build ];then mkdir -p $outdir/rsem-build;fi && \\\n";
$exp_shell .= "cat " . join(" ", split /,+/, $para{'Gene_Fasta'}) . " >$outdir/rsem-build/refMrna.fa && \\\n";
$exp_shell .= "cat " . join(" ", split /,+/, $para{'Gene2Tr'}) . " >$outdir/gene2tr.txt && \\\n";

$exp_shell .= "$Bin/../Bin_CentOS6/rsem/rsem-prepare-reference $outdir/rsem-build/refMrna.fa $outdir/rsem-build/refMrna.fa --bowtie2 --bowtie2-path $Bin/../Bin_CentOS6/bowtie2/bowtie2 --transcript-to-gene-map $outdir/gene2tr.txt && \\\n";
# $exp_shell .= "$Bin/../software/rsem/rsem-prepare-reference $outdir/rsem-build/refMrna.fa $outdir/rsem-build/refMrna.fa --bowtie2 --bowtie2-path $Bin/../Bin_CentOS6/bowtie2/bowtie2 --transcript-to-gene-map $outdir/gene2tr.txt && \\\n";

$exp_shell .= "$perl $Bin/fastaDeal.pl -attr id:len $outdir/rsem-build/refMrna.fa >$outdir/TranscriptLength.txt && \\\n";
$exp_shell .= "awk '{print \$1\"\\t1\\t\"\$2}' $outdir/TranscriptLength.txt >$outdir/TranscriptLength.bed ";
Mkshell("$shdir/build_index.sh", $exp_shell);
print DEPL "$shdir/build_index.sh:0.5G\n" if(!$para{'ERCC'});
print DEPL "$para{'ERCC'}:3G\t$shdir/build_index.sh:0.5G\n" if($para{'ERCC'});
print RES "GeneExp_RSEM\t$shdir/build_index.sh\t$outdir/TranscriptLength.bed\t0\n";
$refer="$outdir/rsem-build/refMrna.fa";
$gene2tr="$outdir/gene2tr.txt";
$Transcripttxt = "$outdir/TranscriptLength.txt";
$Transcriptbed = "$outdir/TranscriptLength.bed";
}else{
$refer=$para{'Gene_Fasta'};
$refer=~ s/,//;
$gene2tr = $para{'Gene2Tr'};
$gene2tr=~ s/,//;
$Transcriptbed = dirname($refer)."/TranscriptLength.bed";
$Transcripttxt = dirname($refer)."/TranscriptLength.txt";
$g2id = $para{'gene2id'};
}
# # # replace $outdir/rsem-build/refMrna.fa with $refer
$para{'Bowtie2_Options'}=~ /-p (?<processer>\d)/;
my $cpu_num1= "$+{processer}cpu";
$para{'RSEM_Options'}=~ /-p (?<processer2>\d)/;
my $cpu_num2= "8cpu";
if (defined $+{processer2}){
    my $cpu_num2= "$+{processer2}cpu";
}else{
    $para{'RSEM_Options'}.= ' -p 8 ';
    }


`mkdir -p  $para{'Result_Dir'}/../../QC/ERCC` unless -d "$para{'Result_Dir'}/../../QC/ERCC";
foreach my $s (keys %sample) {
	my $resultdir = Mkdir("$para{'Result_Dir'}/$s");
	# $exp_shell = "source /ldfssz1/ST_CANCER/CGR/SHARE/CancerPipeline/mRNAseq_V1/RNA_RNAref_2018a/software/SourceMe2.sh && \\\n";
	# $exp_shell = "export PATH=$Bin/../software/:\$PATH && \\\n";
	if(@a > 1){
		$exp_shell = "$Bin/../Bin_CentOS6/bowtie2/bowtie2 $para{'Bowtie2_Options'} -x $refer -1 $sample{$s}[0] -2 $sample{$s}[1] 2>$resultdir/$s.Map2GeneStat.xls | $Bin/../Bin_CentOS6/samtools view -S -b -o $resultdir/$s.bam -  ";
	}else{
		$exp_shell = "$Bin/../Bin_CentOS6/bowtie2/bowtie2 $para{'Bowtie2_Options'} -x $refer -U $sample{$s}[0] 2>$resultdir/$s.Map2GeneStat.xls | $Bin/../Bin_CentOS6/samtools view -S -b -o $resultdir/$s.bam - ";
	}
	Mkshell("$shdir/GeneExp1_$s.sh", $exp_shell);
	# $exp_shell = "source /zfssz5/BC_PUB/Software/03.Soft_ALL/SourceMe2.sh && \\\n";
	# $exp_shell .= "export PATH=$Bin/../software/:\$PATH && \\\n";

	$exp_shell = "$Bin/../Bin_CentOS6/rsem/rsem-calculate-expression $para{'RSEM_Options'} $ifPE --sort-bam-by-coordinate --bam $resultdir/$s.bam $refer $resultdir/$s && \\\n";
	# --append-names
        $exp_shell .= "awk '{if(\$7!=0.00)print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$5\"\\t\"\$7}' $resultdir/$s.genes.results |grep -v '^ERCC' >$resultdir/$s.gene.fpkm.xls && \\\n";
	$exp_shell .= "awk '{if(\$7!=0.00)print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$5\"\\t\"\$7}' $resultdir/$s.genes.results |grep -e '^gene_id' -e '^ERCC' >$para{'Result_Dir'}/../../QC/ERCC/$s.ERCC.gene.fpkm.xls && \\\n" if($para{'ERCC'});
	$exp_shell .= "awk '{if(\$7!=0.00)print \$1\"\\t\"\$2\"\\t\"\$3\"\\t\"\$5\"\\t\"\$7}' $resultdir/$s.isoforms.results |grep -v '^ERCC' >$resultdir/$s.transcript.fpkm.xls && \\\n";
	$exp_shell .= "cp $resultdir/$s.gene.fpkm.xls $resultdir/$s.transcript.fpkm.xls $fpkmdir && \\\n";
	# Mkshell("$shdir/GeneExp2_$s.sh", $exp_shell);
	$exp_shell .= "$perl $Bin/BowtieMapStat.pl -bam $resultdir/$s.bam -strand -key $resultdir/$s.Bowtie2Gene -seqType $seqType -samtools $Bin/../Bin_CentOS6/samtools -gene2tr $gene2tr && \\\n"; #lcc 
	$exp_shell .= "$Bin/BowtieErccMapStat.pl -bam $resultdir/$s.bam -key $resultdir/$s.Ercc.Bowtie2Gene -seqType $seqType -samtools $Bin/../Bin_CentOS6/samtools -gene2tr $gene2tr && \\\n" if($para{'ERCC'});
	$exp_shell .= "cp $resultdir/$s.Ercc.Bowtie2Gene.MapReadsStat.xls $para{'Result_Dir'}/../../QC/ERCC && \\\n" if($para{'ERCC'});
	$exp_shell .= "cp $resultdir/$s.Bowtie2Gene.MapReadsStat.xls $mapdir && \\\n";
	# Mkshell("$shdir/Mapstat_$s.sh", $exp_shell);
	print DEPL "$sample{$s}[-1]\t$shdir/GeneExp1_$s.sh:15G:$cpu_num1\n";
	# print DEPL "$shdir/GeneExp1_$s.sh:5G:8cpu\t$shdir/GeneExp2_$s.sh:5G:8cpu\n";
	# print DEPL "$shdir/GeneExp2_$s.sh:5G:8cpu\t$shdir/Mapstat_$s.sh:1.5G\n";
	# push @depends, "$shdir/GeneExp2_$s.sh:5G:8cpu";
	# push @depends2, "$shdir/Mapstat_$s.sh:1.5G";
	$fpkmlist{$s} = "$resultdir/$s.gene.fpkm.xls";
	if (@p>1){
	print DEPL "$shdir/build_index.sh:0.5G\t$shdir/GeneExp1_$s.sh:15G:$cpu_num1\n";
	}

        # print RES "GeneExp_RSEM\t$shdir/GeneExp2_$s.sh\t$fpkmdir/$s.gene.fpkm.xls\t1\n";
        # print RES "GeneExp_RSEM\t$shdir/GeneExp2_$s.sh\t$fpkmdir/$s.transcript.fpkm.xls\t1\n";
        # print RES "GeneExp_RSEM\t$shdir/GeneExp2_$s.sh\t$mapdir/$s.Bowtie2Gene.MapReadsStat.xls\t1\n";

	$exp_shell .= "$perl $Bin/ReadsRandomInGene_forBowtie_ggplot2.pl -len $Transcripttxt -bam $resultdir/$s.bam -seqType $seqType -prefix $resultdir/$s && \\\n";
	$exp_shell .= "cp $resultdir/$s.ReadsRandom.pdf $resultdir/$s.ReadsRandom.png $randomdir && \\\n";
	# Mkshell("$shdir/ReadsRandom_$s.sh", $exp_shell);
	# print DEPL "$shdir/GeneExp2_$s.sh:5G:8cpu\t$shdir/ReadsRandom_$s.sh:0.5G\n";
	# print RES "GeneExp_RSEM\t$shdir/ReadsRandom_$s.sh\t$randomdir/$s.ReadsRandom.png\t1\n";
	# $exp_shell .= "rm $resultdir/$s.bam && \\\n"; 
	$exp_shell .= "$Bin/../Bin_CentOS6/samtools depth -b $Transcriptbed $resultdir/$s.transcript.sorted.bam >$resultdir/$s.depth.txt && \\\n";
	$exp_shell .= "$perl $Bin/cover_stat.pl $resultdir/$s.depth.txt $Transcripttxt $s $resultdir && \\\n";
	$exp_shell .= "gzip $resultdir/$s.depth.txt && \\\n";
	$exp_shell .= "cp $resultdir/$s.ReadsCoverage.pdf $resultdir/$s.ReadsCoverage.png $coverdir && \\\n";
	$exp_shell .= "rm $resultdir/$s.transcript.sorted.bam $resultdir/$s.transcript.sorted.bam.bai && \\\n";
	# Mkshell("$shdir/ReadsCoverage_$s.sh", $exp_shell);
	# print DEPL "$shdir/GeneExp2_$s.sh:5G:8cpu\t$shdir/ReadsCoverage_$s.sh:0.5G\n";
	# print RES "GeneExp_RSEM\t$shdir/ReadsCoverage_$s.sh\t$coverdir/$s.ReadsCoverage.png\t1\n";
	$exp_shell .= "$Bin/SeqSaturation_Bowtie.pl -bam $resultdir/$s.transcript.bam -seqType $seqType -gene2tr $gene2tr -prefix $resultdir/$s -cutoff 1 && \\\n";
	$exp_shell .= "mv $resultdir/$s.SeqSaturation* $saturationdir && \\\n";
        $exp_shell .= "rm $resultdir/$s.transcript.bam $resultdir/$s.bam";  #the first bam from geneexp1,should leave to last
	# Mkshell("$shdir/SeqSaturation_$s.sh", $exp_shell);
             
	# print DEPL "$shdir/GeneExp2_$s.sh:5G:8cpu\t$shdir/SeqSaturation_$s.sh:1G\n";
	# print RES "GeneExp_RSEM\t$shdir/SeqSaturation_$s.sh\t$coverdir/$s.SeqSaturation.png\t1\n";
        Mkshell("$shdir/GeneExp2_$s.sh", $exp_shell);
        print DEPL "$shdir/GeneExp1_$s.sh:15G:$cpu_num1\t$shdir/GeneExp2_$s.sh:15G:$cpu_num2\n";
        print DOWL "$s\t$resultdir/$s.gene.fpkm.xls\t$shdir/GeneExp2_$s.sh:15G:$cpu_num2\n";
        print RES "GeneExp_RSEM\t$shdir/GeneExp2_$s.sh\t$fpkmdir/$s.gene.fpkm.xls\t1\n";
        print RES "GeneExp_RSEM\t$shdir/GeneExp2_$s.sh\t$fpkmdir/$s.transcript.fpkm.xls\t1\n";
        print RES "GeneExp_RSEM\t$shdir/GeneExp2_$s.sh\t$mapdir/$s.Bowtie2Gene.MapReadsStat.xls\t1\n";
        print RES "GeneExp_RSEM\t$shdir/GeneExp2_$s.sh\t$randomdir/$s.ReadsRandom.png\t1\n";
        print RES "GeneExp_RSEM\t$shdir/GeneExp2_$s.sh\t$coverdir/$s.ReadsCoverage.png\t1\n";
        print RES "GeneExp_RSEM\t$shdir/GeneExp2_$s.sh\t$coverdir/$s.SeqSaturation.png\t1\n";
        push @depends, "$shdir/GeneExp2_$s.sh:15G:$cpu_num2";
	
}

#$exp_shell = "# remove bam file\n";
#foreach my $s (keys %sample)
#{
#$exp_shell .= "rm $outdir/$s/$s.transcript.sorted.bam $outdir/$s/$s.bam $outdir/$s/$s.transcript.bam $outdir/$s/$s.transcript.sorted.bam.bai\n";
#}
# add fpmklist from $para{'Result_Dir'}
my @temp_fpkmlist=`ls $para{'Result_Dir'}`;
foreach (@temp_fpkmlist){
	# print $_;
	if ($_ =~ /PCA/){
		next;
	}
	else{
		chomp($_);
		$fpkmlist{$_}="$para{'Result_Dir'}/$_/$_.gene.fpkm.xls" unless ($fpkmlist{$_});
		# print $fpkmlist{$_};
	}
}



$exp_shell = "# all samples fpkm statistic\n";###lcc
$exp_shell .= "rm $outdir/*gene.fpkm.xls \n";
$exp_shell .= "$python $Bin/trans_genexp.py -i $outdir -s $g2id && \\\n";
$exp_shell .= "$perl $Bin/AllGeneStat.pl $outdir $outdir/AllSamples.GeneExpression.FPKM.xls && \\\n";
$exp_shell .= "$perl $Bin/AllTranscriptStat.pl $outdir $outdir/AllSamples.TranscriptExpression.FPKM.xls && \\\n"; ###lc 20170110 for TFcluster
$exp_shell .= "$perl $Bin/ExpStat_se.pl -indir $outdir -output $outdir/GeneExpressionSummary.xls && \\\n" if(@a < 2);
$exp_shell .= "$perl $Bin/ExpStat.pl -indir $outdir -output $outdir/GeneExpressionSummary.xls && \\\n" if(@a > 1);
$exp_shell .= "cp $outdir/AllSamples.GeneExpression.FPKM.xls $fpkmdir && \\\n";
if(@a > 1 && $para{'Gene_Fasta'}=~/,/){
$exp_shell .= "$perl $Bin/stackedplot.SampleGene.pl -genelist $outdir/AllSamples.GeneExpression.FPKM.xls -outdir $outdir && \\\n"; ###lcc 20170108
$exp_shell .= "$perl $Bin/barplot.GeneNumber.pl -summary $outdir/GeneExpressionSummary.xls -outdir $outdir && \\\n";     ###lcc 20170105
$exp_shell .= "cp $outdir/knowngene_novelgene_number_barplot.pdf $outdir/knowngene_novelgene_number_barplot.png $fpkmdir && \\\n";       ###lc
$exp_shell .= "cp $outdir/Expressed_gene_percentage.pdf $outdir/Expressed_gene_percentage.png $fpkmdir && \\\n"; ###lcc 20170108
}

$exp_shell .= "cp $outdir/GeneExpressionSummary.xls $para{'Cp_Dir'} && \\\n";
if (keys %sample >= 2) {
	my $cordir = Mkdir("$para{'Cp_Dir'}/CorrelationHeatmap");
	$exp_shell .= "# correlation heatmap\n";
	$exp_shell .= "$perl $Bin/drawCorrelationHeatmap.pl -indir $fpkmdir -outdir $cordir && \\\n";
	$exp_shell .= "mv $cordir/correlation-heatmap.R $outdir && \\\n";
	$exp_shell .= "$perl $Bin/CorrelationQC.pl $cordir/AllSamples.correlation.xls $para{'Group'} >$outdir/CorrelationQC.xls && \\\n" if($para{'Group'}); #licong 20170907 add for QC
	print RES "GeneExp_RSEM\t$shdir/statistic.sh\t$outdir/CorrelationQC.xls\t0\n";
}
if (keys %sample >= 3) {
	my $hclustdir = Mkdir("$para{'Cp_Dir'}/HclusterTree");
	$exp_shell .= "# hcluster tree\n";
	$exp_shell .= "$perl $Bin/drawHclustTree.pl -indir $fpkmdir -outdir $hclustdir && \\\n";
}

#######lcc
my %linkfpkm;
if ($para{'Venn_Group'})
{
        $exp_shell .= "#union genes among samples\n";
        my (%Venn_g, %check_s);
        my @Venn_group = split (/\s+/, $para{'Venn_Group'});##HBRR:HBRR1,HBRR2,HBRR3 UHRR:UHRR1,UHRR2,UHRR3
        foreach my $V (@Venn_group)#HBRR:HBRR1,HBRR2,HBRR3 UHRR:UHRR1,UHRR2,UHRR3
        {
                my ($Venn_name, $Venn_sample)=split(/:+/,$V);####HBRR   HBRR1,HBRR2,HBRR3
		ThrowError("ERROR: Venn_Group [$Venn_name] was named repeatedly,please check!\n",*STDERR) if ($Venn_g{$Venn_name});$Venn_g{$Venn_name}=1;
                my $linkin = "";
				%check_s=();
                foreach my $Vs (split /,+/,$Venn_sample)###HBRR1 HBRR2 HBRR3
                {
                        ThrowError("ERROR: can't find sample [$Vs] in Venn_Group [$Venn_name]\n", *STDERR) unless ($fpkmlist{$Vs});
                        ThrowError("ERROR: sample [$Vs] exists more than Once in Venn_Group [$Venn_name]\n",*STDERR) if ($check_s{$Vs}); $check_s{$Vs} = 1;
                        $linkin .= "$fpkmlist{$Vs},";
                }
                chop $linkin;
                #my $linkdir = Mkdir ("$para{'Result_Dir'}/$Venn_name");
		#$linkfpkm{$Venn_name}="$linkdir/$Venn_name.gene.fpkm.xls";
		$linkfpkm{$Venn_name}="$para{'Result_Dir'}/$Venn_name.gene.fpkm.xls";#20170105
		$exp_shell .= "$perl $Bin/Union_genes_fpkm.pl -list $linkin -out $linkfpkm{$Venn_name}  && \\\n";
        }
}
#########

if ($para{'Venn_Plan'}) {
	$exp_shell .= "#venn graph for expressed gene among samples or groups\n";
	my %repeat = ();
	my @plan = split /\s+/,$para{'Venn_Plan'};
	foreach my $p (@plan) {
		my @samples = split /,+/, $p;
		my $repeat = join("#",sort @samples);
		print STDERR "WARN: Venn_Plan [$p] be seted repeatedly, ignoring!\n" if ($repeat{$repeat}); $repeat{$repeat} = 1;
		my $in = "";
		my $name = "";
		my $img = "";
		ThrowError("Error: each venn plan should more than 1 sample and less than 6 samples!\n",*STDERR) unless (@samples > 1 && @samples < 6);
		foreach my $s (@samples) {
			ThrowError("ERROR: no such sample named [$s]\n",*STDERR) unless ($fpkmlist{$s} || $linkfpkm{$s});###lcc
			if (defined $fpkmlist{$s})#lcc
			{
				$in .= "$fpkmlist{$s},";
				$name .= "$s,";
				$img .= "$s-";
			}
			elsif (defined $linkfpkm{$s})#lcc
			{
				$in .= "$linkfpkm{$s},";
				$name .= "$s,";
				$img .= "$s-";
			}			
		}
		chop $in; chop $name; chop $img;
		my $venndir = Mkdir("$para{'Cp_Dir'}/VennDiagram/$img");
		$exp_shell .= "$perl $Bin/venny.pl -infile $in -name $name -header -outdir $venndir -imgname $img.venn && \\\n";
		$exp_shell .= "mv $venndir/$img.venn.R $outdir && \\\n";
	}
}

if ($para{'PCA_Plan'}) {
	my $pca_shell = "# pca\n";
	my %group = ();
	my $pca_flag = 0;
	if ($para{'PCA_Group'}) {
		my @group = split /\s+/,$para{'PCA_Group'};
		foreach my $g (@group) {
			my ($name,$samples) = split /:+/, $g;
			ThrowError("ERROR: PCA_Group [$name] be named repeatedly,please check!\n",*STDERR) if ($group{$name});
			my %check = ();
			foreach my $s (split /,+/,$samples) {
				ThrowError("ERROR: can't find sample [$s] in PCA_Group [$name]\n", *STDERR) unless ($fpkmlist{$s});
				ThrowError("ERROR: sample [$s] exists more than 1 times in PCA_Group [$name]\n",*STDERR) if ($check{$s});$check{$s} = 1;
			}
			$group{$name} = $samples;
		}
	}
	my $group_index = 1;
	foreach my $p (split /\s+/, $para{'PCA_Plan'}) {
		open TMP,">$para{'Result_Dir'}/PCA_Group.$group_index.txt" or die $!;
		open LST,">$para{'Result_Dir'}/PCA_GeneLsit.$group_index.txt" or die $!;
		my @samples = split /,+/, $p;
		my $num = 0; my %check = ();
		foreach my $s (@samples) {
			if ($fpkmlist{$s}) {
				$num++;
				ThrowError("ERROR: sample [$s] exists more than 1 times in PCA_Plan [$p]\n", *STDERR) if ($check{$s}); $check{$s} = 1;
				print TMP "$s\t$s\n";
				print LST "$s\t$fpkmlist{$s}\n";
			}elsif ($group{$s}) {
				my @samp = split /,+/,$group{$s};
				$num += scalar(@samp);
				print TMP "$s\t$group{$s}\n";
				for(@samp) {
					ThrowError("ERROR: sample [$_] exists more than 1 times in PCA_Plan [$p]\n", *STDERR) if ($check{$_}); $check{$_} = 1;
					print LST "$_\t$fpkmlist{$_}\n";
				}
			}else {
				ThrowError("ERROR: no such sample or PCA_Group [$s]!\n",*STDERR);
			}
		}
		ThrowError("ERROR: PCA should not less than 5 samples,please check PCA_Plan [$p]\n", *STDERR) unless ($num >= 5);
#		my $pca_processdir = Mkdir("$para{'Result_Dir'}/PCA/".join("-",@samples));
#		my $pca_resultdir = Mkdir("$para{'Cp_Dir'}/PCA/".join("-",@samples));
		my $pca_processdir = Mkdir("$para{'Result_Dir'}/PCA/plan_$group_index");
		my $pca_resultdir = Mkdir("$para{'Cp_Dir'}/PCA/plan_$group_index");
		$pca_shell .= "$perl $Bin/PCA_deal.pl -list $para{'Result_Dir'}/PCA_GeneLsit.$group_index.txt -IDcolumn 1 -Exprcolumn 5 -output $para{'Result_Dir'}/PCA_GeneExp.$group_index.txt && \\\n";
		if($para{'PCA_Group'}){
			$pca_shell .= "$perl $Bin/princomp_draw.pl -expr $para{'Result_Dir'}/PCA_GeneExp.$group_index.txt -g $para{'Result_Dir'}/PCA_Group.$group_index.txt -o $pca_processdir && \\\n";
		}else{
			$pca_shell .= "$perl $Bin/princomp_draw.pl -expr $para{'Result_Dir'}/PCA_GeneExp.$group_index.txt -o $pca_processdir && \\\n";
		}
		$pca_shell .= "cp $pca_processdir/*pdf $pca_processdir/*png $pca_processdir/PCA_result.xls  $pca_resultdir && \\\n";
		$group_index ++;$pca_flag = 1;
	}
	$exp_shell .= $pca_shell if ($pca_flag == 1);
}
$exp_shell .="$perl $Bin/drawbox_density-plot.pl -i $fpkmdir -o $fpkmdir && \\\n";
$exp_shell .="$perl $Bin/drawstacked_bar.pl -indir $fpkmdir -outdir $fpkmdir && \\\n";
$exp_shell .= "echo done ";
Mkshell("$shdir/statistic.sh", $exp_shell);
foreach my $dep (@depends) {
	print DEPL "$dep\t$shdir/statistic.sh:0.1G\n";
}
foreach my $dep2 (@depends2) {
	print DEPL "$dep2\t$shdir/statistic.sh:0.1G\n";
}
print RES "GeneExp_RSEM\t$shdir/statistic.sh\t$outdir/GeneExpressionSummary.xls\t1\n";

print STDERR "DONE";
