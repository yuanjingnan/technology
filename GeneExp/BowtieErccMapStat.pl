#!/usr/bin/perl -w
=head1 Command-line Option
	-bam      * bowtie bam result(or sam), separated by comma ","
	-key      * prefix of output
	-seqType  * sequencing type,PE or SE
	-gene2tr  ? (geneID \t transcriptID), necessary for gene
	-samtools ? samtools, necessaru for bam format
	-help       help information
=cut
#modified:20140809 linruichai, 20120829, 20120724
use strict;
use File::Basename;
use Getopt::Long;
use PerlIO::gzip;
my($bam,$samtools,$gene2tr,$key,$seqType,$help);
GetOptions(
	"bam=s"=>\$bam,
	"key=s"=>\$key,
	"gene2tr=s"=>\$gene2tr,
	"seqType=s"=>\$seqType,
	"help"=>\$help,
	"samtools=s"=>\$samtools,
);
die `pod2text $0` if(!defined $bam || !defined $key || defined $help);
die `pod2text $0` if($bam=~/\.bam$/ && !defined $samtools);
my (%reads,$total_reads,$total_mapped,$perfect,$mismatch,$unique,$multiple,$unmap,$length,$sampleName,%gene2tr);

%gene2tr=&Gene2Tr($gene2tr) if(defined $gene2tr);
my @bam=split /,/,$bam;
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
while(<BAM>)
{
	chomp;
	next if(/^\@/);
	my $line=$_;
	my @tmp=split /\t/,$line;
	$length=length ($tmp[9]) if(!$length);
	$reads{$tmp[0]} = 1;
	if($tmp[2]!~/ERCC/)
	{
	#	$unmap++;
		next;
	}
	elsif($tmp[0] ne $oldReadID && $tmpLine ne "")
	{
		$oldReadID = $tmp[0];
		if($seqType eq "SE")
		{ #for SE result
			$total_mapped++;
			my @read1 = split /\n/,$tmpLine;
			my @gene = ();
			my @MD = ();
			for(@read1)
			{ #处理所有相同readID的行
				my @column = split /\t/,$_;
				my @MD_Z = grep /MD:Z/,@column;
				my $MD_Z =($MD_Z[0]=~tr/[ATCG]/[ATCG]/);
				if(defined $gene2tr && exists $gene2tr{$column[2]})
				{
					push @gene,$gene2tr{$column[2]};
				}
				elsif((defined $gene2tr && !exists $gene2tr{$column[2]}) || undef $gene2tr)
				{
					push @gene,$column[2];
				}
					push @MD, $MD_Z;
			}
			
			#unique or multiple	
			if (@gene == 1)
			{
				$unique++;				
			}
			else
			{
				if (!$gene2tr)
				{
					$multiple++;
				}
				else
				{
					my %count;
					my @uniq_gene = grep { ++$count{ $_ } < 2; } @gene; #去掉重复基因ID
					if(@uniq_gene == 1){ $unique++;}else{$multiple++;}
				}
			}
	
			#perfect match or mismatch
			my $perfect_flag = 0;
			for(@MD)
			{
				if($_ == 0){ $perfect_flag = 1;last;}
			}
				$perfect++ if($perfect_flag == 1);
				$mismatch++ if ($perfect_flag == 0);		        
		}
		else
		{ #for PE result,与SE处理方式相同，所read id相同的行放一起处理，只是会成对出现，所以还要区分read1 read2
			$total_mapped++; #read1加一次
			$total_mapped++; #read2也要加一次
			my @reads = split /\n/,$tmpLine;
			my @read1 = ();
			my @read2 = ();
			for(my $i=0;$i<@reads;$i++)
			{ #分别将read1 read2存入不同数组
				push @read1,$reads[$i] if ($i%2!=0);
				push @read2,$reads[$i] if ($i%2==0);
			}

			#对read1进行处理，方式与SE一样
			my @read1_gene = ();			
			my @read1_MD = ();
			for(@read1)
			{ 
 		        my @column = split /\t/,$_;
		        my @MD_Z = grep /MD:Z/,@column;
		        my $MD_Z =($MD_Z[0]=~tr/[ATCG]/[ATCG]/);
		        push @read1_gene, $gene2tr{$column[2]} if (exists $gene2tr{$column[2]});
				push @read1_gene, $column[2] if (!exists $gene2tr{$column[2]});
		        push @read1_MD, $MD_Z;
			}
			#read1:unique or multiple
			if (@read1_gene == 1)
			{
					$unique++;
			}
			else
			{
		        if (!$gene2tr)
				{
	                    $multiple++;
		        }
				else
				{
		                my %count;
		                my @uniq_gene = grep { ++$count{ $_ } < 2; } @read1_gene; #去掉重复基因ID
		                if(@uniq_gene == 1){ $unique++;}else{$multiple++;}
		        }
			}
			
			#read1:perfect match or mismatch
			my $read1_perfect_flag = 0;
			for(@read1_MD)
			{
		        if($_ == 0){ $read1_perfect_flag = 1;last;}
			}
				$perfect++ if($read1_perfect_flag == 1);
				$mismatch++ if ($read1_perfect_flag == 0);

			#对read2进行处理，方式与SE一样
			my @read2_gene = ();
			my @read2_MD = ();
			for(@read2)
			{
		        my @column = split /\t/,$_;
		        my @MD_Z = grep /MD:Z/,@column;
		        my $MD_Z =($MD_Z[0]=~tr/[ATCG]/[ATCG]/);
		        push @read2_gene, $gene2tr{$column[2]} if (exists $gene2tr{$column[2]});
				push @read2_gene, $column[2] if (!exists $gene2tr{$column[2]});
		        push @read2_MD, $MD_Z;
			}
			#read2:unique or multiple
			if (@read2_gene == 1)
			{
		        $unique++;
			}
			else
			{
		        if (!$gene2tr)
				{
		                $multiple++;
		        }
				else
				{
		                my %count;
		                my @uniq_gene = grep { ++$count{ $_ } < 2; } @read2_gene;
		                if(@uniq_gene == 1){ $unique++;}else{$multiple++;}
		        }
			}
			
			#read2:perfect match or mismatch
			my $read2_perfect_flag = 0;
			for(@read2_MD)
			{
		        if($_ == 0){ $read2_perfect_flag = 1;last;}
			}
				$perfect++ if($read2_perfect_flag == 1);
				$mismatch++ if ($read2_perfect_flag == 0);
		}
                $tmpLine = ""; #清空
				$tmpLine .= "$_\n"; #将新起一行的ID存起来
	}
	else
	{
			$oldReadID = $tmp[0];
			$tmpLine .= "$_\n"; #将相同readID的比对结果存放起来
	}
}
close BAM;

if( $tmpLine ne "")
{
	if($seqType eq "SE")
	{ #for SE result
		$total_mapped++;
		my @read1 = split /\n/,$tmpLine;
		my @gene = ();
		my @MD = ();
		for(@read1)
		{ #处理所有相同readID的行
			my @column = split /\t/,$_;
			my @MD_Z = grep /MD:Z/,@column;
			my $MD_Z =($MD_Z[0]=~tr/[ATCG]/[ATCG]/);
			if(defined $gene2tr && exists $gene2tr{$column[2]})
			{
				push @gene,$gene2tr{$column[2]};
			}
			elsif((defined $gene2tr && !exists $gene2tr{$column[2]}) || undef $gene2tr)
			{
				push @gene,$column[2];
			}
				push @MD, $MD_Z;
		}
			
		#unique or multiple	
		if (@gene == 1)
		{
			$unique++;				
		}
		else
		{
			if (!$gene2tr)
			{
				$multiple++;
			}
			else
			{
				my %count;
				my @uniq_gene = grep { ++$count{ $_ } < 2; } @gene; #去掉重复基因ID
				if(@uniq_gene == 1){ $unique++;}else{$multiple++;}
			}
		}
	
		#perfect match or mismatch
		my $perfect_flag = 0;
		for(@MD)
		{
			if($_ == 0){ $perfect_flag = 1;last;}
		}
			$perfect++ if($perfect_flag == 1);
			$mismatch++ if ($perfect_flag == 0);		        
	}
	else
	{ #for PE result,与SE处理方式相同，所read id相同的行放一起处理，只是会成对出现，所以还要区分read1 read2
		$total_mapped++; #read1加一次
		$total_mapped++; #read2也要加一次
		my @reads = split /\n/,$tmpLine;
		my @read1 = ();
		my @read2 = ();
		for(my $i=0;$i<@reads;$i++)
		{ #分别将read1 read2存入不同数组
			push @read1,$reads[$i] if ($i%2!=0);
			push @read2,$reads[$i] if ($i%2==0);
		}

		#对read1进行处理，方式与SE一样
		my @read1_gene = ();			
		my @read1_MD = ();
		for(@read1)
		{ 
			my @column = split /\t/,$_;
			my @MD_Z = grep /MD:Z/,@column;
			my $MD_Z =($MD_Z[0]=~tr/[ATCG]/[ATCG]/);
			push @read1_gene, $gene2tr{$column[2]} if (exists $gene2tr{$column[2]});
			push @read1_gene, $column[2] if (!exists $gene2tr{$column[2]});
			push @read1_MD, $MD_Z;
		}
			#read1:unique or multiple
			if (@read1_gene == 1)
			{
			$unique++;
			}
			else
			{
				if (!$gene2tr)
				{
					$multiple++;
				}
				else
				{
					my %count;
					my @uniq_gene = grep { ++$count{ $_ } < 2; } @read1_gene; #去掉重复基因ID
					if(@uniq_gene == 1){ $unique++;}else{$multiple++;}
				}
			}
			
			#read1:perfect match or mismatch
			my $read1_perfect_flag = 0;
			for(@read1_MD)
			{
				if($_ == 0){ $read1_perfect_flag = 1;last;}
			}
				$perfect++ if($read1_perfect_flag == 1);
				$mismatch++ if ($read1_perfect_flag == 0);

			#对read2进行处理，方式与SE一样
			my @read2_gene = ();
			my @read2_MD = ();
			for(@read2)
			{
				my @column = split /\t/,$_;
				my @MD_Z = grep /MD:Z/,@column;
				my $MD_Z =($MD_Z[0]=~tr/[ATCG]/[ATCG]/);
				push @read2_gene, $gene2tr{$column[2]} if (exists $gene2tr{$column[2]});
				push @read2_gene, $column[2] if (!exists $gene2tr{$column[2]});
				push @read2_MD, $MD_Z;
			}
			#read2:unique or multiple
			if (@read2_gene == 1)
			{
				$unique++;
			}
			else
			{
				if (!$gene2tr)
				{
					$multiple++;
				}
				else
				{
					my %count;
					my @uniq_gene = grep { ++$count{ $_ } < 2; } @read2_gene;
					if(@uniq_gene == 1){ $unique++;}else{$multiple++;}
				}
			}
			
			#read2:perfect match or mismatch
			my $read2_perfect_flag = 0;
			for(@read2_MD)
			{
				if($_ == 0){ $read2_perfect_flag = 1;last;}
			}
				$perfect++ if($read2_perfect_flag == 1);
				$mismatch++ if ($read2_perfect_flag == 0);
	}
}


$total_reads = (scalar keys %reads)*2 if ($seqType eq "PE");
$total_reads = (scalar keys %reads) if ($seqType eq "SE");
$unmap=$total_reads-$total_mapped;
print OUT"Total Reads\t$total_reads\nTotal Bases\t".$total_reads*$length."\n\nTotal Mapped Reads(%)\t$total_mapped\t".&per($total_mapped)."\n\n"."Perfect Match(%)\t$perfect\t".&per($perfect)."\nMismatch(%)\t$mismatch\t".&per($mismatch)."\n\nUnique Match(%)\t$unique\t".&per($unique)."\nMulti-position Match(%)\t$multiple\t".&per($multiple)."\n\nTotal Unmapped Reads(%)\t$unmap\t".&per($unmap)."\n";
close OUT;
`rm -rf $key.merge.bam`;
sub per{
	my $member=shift;
	my $per=sprintf("%.2f",$member/$total_reads*100);
	return $per;
}
sub Gene2Tr{
	my $file=shift;
	my %hash;
	open Gene2Tr,"< $file";
	while(<Gene2Tr>){
		chomp;
		my @tmp=split /\t/,$_;
		$hash{$tmp[1]}=$tmp[0];
	}
	close Gene2Tr;
	return %hash;
}
