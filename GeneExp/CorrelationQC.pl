#!/usr/bin/perl -w
# licong@genomics.cn
# Usage : perl CorrelationQC.pl All.correlation.stat.xls GroupA:A1,A2,A3 GroupB:B1,B2,B3

my ($head,@h,%if);
my $cor=shift;
open IN,"$cor";chomp($head=<IN>);
@h=split /\s+/,$head;
while(<IN>){
	my @a=split;
	for(1..$#h){
		$ha{$a[0]}{$h[$_]}=$a[$_];
		$ha{$h[$_]}{$a[0]}=$a[$_];
	}
}
foreach my $i(@ARGV){
	my @a=split /:/,$i;
	my @b=split /,/,$a[1];
#	LABEL: for my $s(@b){
	for my $s(@b){
		for(@b){
#			if($ha{$s}{$_}**2 < 0.9){
#				print "$i;";last LABEL
#			}
			print "Group:$a[0]\t$s\t$_\t$ha{$s}{$_}\n" if(!$if{"$_\t$s"} && not $_ eq $s);
			$if{"$s\t$_"}=1;
		}
	}
}
			
		
