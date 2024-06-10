#!/use/bin/perl -w
=head1 Command-line Option
        -bam      ? bowtie bam result(or sam), separated by comma ","
        -sam      ? bowtie sam result, separated by comma ","
	-strand     stat strand information
        -key      * prefix of output
        -seqType  * sequencing type,PE or SE
        -gene2tr  ? (geneID \t transcriptID), necessary for gene
        -samtools ? samtools, necessaru for bam format
        -help       help information
=cut
#modified: 20140809 linruichai, 20120829, 20120724
#modified: 20150701 lwy, 20161124
use strict;
use File::Basename;
use Getopt::Long;
use PerlIO::gzip;
my($bam,$sam,$samtools,$gene2tr,$key,$seqType,$strand,$help);
GetOptions(
        "bam=s"     =>\$bam,
        "sam=s"     =>\$sam,
        "key=s"     =>\$key,
        "gene2tr=s" =>\$gene2tr,
        "seqType=s" =>\$seqType,
	"strand"    =>\$strand,
        "help"      =>\$help,
        "samtools=s"=>\$samtools,
);
die `pod2text $0` if((!defined $bam && !$sam) || !defined $key || defined $help);
die `pod2text $0` if($bam=~/\.bam$/ && !defined $samtools);
my (%reads,$length,$sampleName,%gene2tr);
my ($total_reads,$total_mapped,$perfect,$mismatch,$unique,$multiple,$unmap) = (0,0,0,0,0,0,0);
my @strand_counts = (0, 0, 0, 0, 0, 0);
#my ($uniq_p, $uniq_m, $mult_p, $mult_m, $mult_mult) = ($strand_counts[0],[1],[2],[3],[4]);

%gene2tr=&Gene2Tr($gene2tr) if(defined $gene2tr);
my $infile = ($sam)?$sam:$bam;
my @bam = split /,/, $infile;
if(@bam==1){
        $bam=$bam[0];
}else{
        system("$samtools merge $key.merge.bam @bam");
        $bam="$key.merge.bam";
}
open OUT,"> $key.MapReadsStat.xls";
print OUT "Mapping Type\tReads Number\tPercentage\n";
### reading bam file...
if($bam=~/\.sam$/){
        open BAM,$bam;
}else{
        open(BAM,($bam =~ /\.bam$/)?"$samtools view $bam|":$bam) or die $!;
}
my $oldReadID = "";
my $tmpLine = "";

while(<BAM>){
	chomp;
        next if(/^\@/);
        my $line=$_;
        my @tmp=split /\t/,$line;
        $length=length ($tmp[9]) if(!$length);
        $reads{$tmp[0]} = 1;
        if($tmp[2] eq "*"){
                $unmap++;
                next;
        }
        if($tmp[0] ne $oldReadID && $tmpLine ne ""){
		$oldReadID = $tmp[0];
                if($seqType eq "SE"){  #for SE result
			$total_mapped++;
			my @read1 = split /\n/,$tmpLine;
			&get_count(\@read1, 1);
		}else{  #for PE result,与SE处理方式相同，所read id相同的行放一起处理，只是会成对出现，所以还要区分read1 read2
			$total_mapped++; #read1加一次
                        $total_mapped++; #read2也要加一次
			my @reads = split /\n/,$tmpLine;
                        my @read1 = ();
                        my @read2 = ();
			for(my $i=0;$i<@reads;$i++){ #分别将read1 read2存入不同数组
                                push @read1,$reads[$i] if ($i%2!=0);
                                push @read2,$reads[$i] if ($i%2==0);
                        }
			#对read1进行处理，方式与SE一样
			&get_count(\@read1, 1);
			#对read2进行处理，方式与SE一样
			&get_count(\@read2, 2);
		}
		$tmpLine = ""; #清空
                $tmpLine .= "$_\n"; #将新起一行的ID存起来
	}else{
                $oldReadID = $tmp[0];
                $tmpLine .= "$_\n"; #将相同readID的比对结果存放起来
        }
}
close BAM;
## for last tmpLine
if($tmpLine ne ""){
	my @reads = split /\n/,$tmpLine;
	if($seqType eq "SE"){
		$total_mapped++;
		&get_count(\@reads, 1);
	}else{
		$total_mapped += 2;
		my (@read1, @read2);
		for(my $i=0;$i<@reads;$i++){
			push @read1,$reads[$i] if ($i%2!=0);
			push @read2,$reads[$i] if ($i%2==0);
		}
		&get_count(\@read1, 1);
		&get_count(\@read2, 2);
	}
}

$total_reads = (scalar keys %reads)*2 if ($seqType eq "PE");
$total_reads = (scalar keys %reads) if ($seqType eq "SE");

print OUT"Total Reads\t$total_reads\t100.00\nTotal BasePairs\t".$total_reads*$length."\t100.00\n\nTotal Mapped Reads\t$total_mapped\t".&per($total_mapped)."\n\n"."Perfect Match\t$perfect\t".&per($perfect)."\nMismatch\t$mismatch\t".&per($mismatch)."\n\nUnique Match\t$unique\t".&per($unique)."\nMulti-position Match\t$multiple\t".&per($multiple)."\n";
#print OUT "uniq +: $strand_counts[0]*2\tuniq -: $strand_counts[1]*2\nmult +: $strand_counts[2]*2\tmult -: $strand_counts[3]*2\tmult mult: $strand_counts[4]*2\n$strand_counts[5]\n";
if($strand){
	my $positive = ($strand_counts[0] + $strand_counts[2]) * 2;
	my $negative = ($strand_counts[1] + $strand_counts[3]) * 2;
	my $obscure = $strand_counts[4] * 2;
	print OUT "\nSense-strand Match\t$positive\t".&per($positive)."\nTemplate-strand Match\t$negative\t".&per($negative)."\nObscure-strand Match\t$obscure\t".&per($obscure)."\n";
}
print OUT "\nTotal Unmapped Reads\t$unmap\t".&per($unmap)."\n";
close OUT;
`rm -rf $key.merge.bam`;

################# sub
sub per{
        my $member=shift;
        my $per=sprintf("%.2f",$member/$total_reads*100);
		#return $per."%";
		return $per;
}
sub Gene2Tr{
        my $file=shift;
        my %hash;
        open Gene2Tr,"<", $file||die "Can't open $file: $!";
        while(<Gene2Tr>){
                chomp;
                my @tmp=split /\t/,$_;
                $hash{$tmp[1]}=$tmp[0];
        }
        close Gene2Tr;
        return %hash;
}

sub get_count{
	my ($read_a, $s_flag) = @_;
	
	my @gene = ();
	my @MD = ();
	my @flags = ();

	for(my $i=0; $i<scalar(@{$read_a}); $i++){ #处理所有相同readID的行
		my $read_l = $$read_a[$i];
        	my @column = split /\t/,$read_l;
                if(defined $gene2tr && exists $gene2tr{$column[2]}){
                	push @gene,$gene2tr{$column[2]};
                }else{
                        push @gene,$column[2];
                }
                my @MD_Z = grep /MD:Z/, @column;
                my $MD_Z =($MD_Z[0]=~tr/[ATCG]/[ATCG]/);
		push @MD, $MD_Z;
                push @flags, $column[1];
	}

	#unique or multiple
        if (@gene == 1){
                $unique++;
        	get_strand_flag(\@flags, $seqType, 'uniq', \@strand_counts) if(defined $strand && $s_flag == 1);
        }else{
                if (!$gene2tr){
                	$multiple++;
                        get_strand_flag(\@flags, $seqType, 'mult', \@strand_counts) if(defined $strand && $s_flag == 1);
                }else{
                        my %count;
                        my @uniq_gene = grep { ++$count{ $_ } < 2; } @gene; #去掉重复基因ID
                        if(@uniq_gene == 1){
                        	$unique++;
                                get_strand_flag(\@flags, $seqType, 'uniq', \@strand_counts) if(defined $strand && $s_flag == 1);
                        }else{
                                $multiple++;
                                get_strand_flag(\@flags, $seqType, 'mult', \@strand_counts) if(defined $strand && $s_flag == 1);
                        }
                }
	}

	#perfect match or mismatch
        my $perfect_flag = 0;
        for(@MD){
        	if($_ == 0){ $perfect_flag = 1;last;}
        }
        $perfect++ if($perfect_flag == 1);
        $mismatch++ if ($perfect_flag == 0);
}

sub get_strand_flag{
	my ($arr, $pese, $type, $arr2) = @_;
	my $out_f;
        my $abs_arr0 = abs($$arr[0]);
        if($type eq 'uniq' && $pese eq "SE"){
                $out_f = ($abs_arr0 & 0x10)?1:0;
                if($out_f == 0){
                        $$arr2[0] ++;
                }else{
                        $$arr2[1] ++;
                }
        }elsif($type eq 'mult' && $pese eq "SE"){
                my $pre_flag = ($abs_arr0 & 0x10)?1:0;
                foreach my $a(@{$arr}){
                        if($pre_flag == 0 && ($abs_arr0 & 0x10)){
                                $out_f = 2; last;
                        }elsif($pre_flag == 1 && ($abs_arr0 | 0x10)){
                                $out_f = 2; last;
                        }
                }
                $out_f = $pre_flag unless($out_f && $out_f == 2);
                if($out_f == 0){
                        $$arr2[2] ++;
                }elsif($out_f ==1){
                        $$arr2[3] ++;
                }else{
                        $$arr2[4] ++;
                }
        }elsif($type eq 'uniq' && $pese eq "PE"){
                if($abs_arr0 & 34){ ## 2+32=34;
                        $$arr2[0] ++;
                #}elsif($abs_arr0 & 18){ ## 2+16=18;
                }elsif(($abs_arr0 & 0x2) && ($abs_arr0 & 0x16)){
                        $$arr2[1] ++;
                }
        }elsif($type eq 'mult' && $pese eq "PE"){
                my $pre_flag;
                foreach my $a(@{$arr}){
                        $a = abs($a);
                        if(!defined $pre_flag){
                                if($a & 34){
                                        $pre_flag = 0;
                                }elsif(($a & 0x2) && ($a & 0x16)){
                                        $pre_flag = 1;
                                }
                        }else{
                                if(($pre_flag == 0) && ($a & 0x10)){
                                        $out_f = 2; last;
                                }elsif(($pre_flag == 1) && ($a & 0x20)){
                                        $out_f = 2; last;
                                }
                        }
                }
                if(!defined $pre_flag){
                        $out_f = 3;
                }else{
                        $out_f = $pre_flag unless(defined $out_f && $out_f == 2);
                }
                if($out_f == 0){
                        $$arr2[2] ++;
                }elsif($out_f ==1){
                        $$arr2[3] ++;
                }elsif($out_f == 2){
                        $$arr2[4] ++;
                }else{
                        $$arr2[5] ++;
                }
        }
}
