#!/usr/bin/perl -w
use strict;
use Getopt::Long;

#read option
my ($output,$fqlist,$IDcolumn,$Exprcolumn);
GetOptions(
       "list=s"=>\$fqlist,
       "output=s"=>\$output,
       "IDcolumn=i"=>\$IDcolumn,
       "Exprcolumn=i"=>\$Exprcolumn
);
$IDcolumn ||= 1;
$Exprcolumn ||= 5;
#read directory
my (@xls,@Sname);
open IN,$fqlist or die "Cannot open my  file1 $fqlist:$!\n";
while (<IN>) {
      chomp;
      my @line_1 = split(/\s/);
      push(@Sname,$line_1[0]);
      push(@xls,$line_1[1]);
}
close IN;

my (%hash,%GeneID_all);
my @sample_name;
my $xls_vol = @xls - 1;
my $IDcol = $IDcolumn - 1;
my $Exprcol = $Exprcolumn -1;
my ($i,$sample);
for ($i=0;$i <= $xls_vol;$i++){
   $sample = $Sname[$i]; 
   push(@sample_name,$sample);
   open IN1,$xls[$i] or die "Can't open my file2 $xls[$i]:$!\n";
   while(<IN1>){
         chomp;
         my @line = split(/\t/,$_);
         my $GeneID = $line[$IDcol];
         if($line[$Exprcol] =~ /[a-zA-Z]/ ){
             my $name2 = $line[$Exprcol];
        }else{
             $hash{$GeneID}{$sample} = $line[$Exprcol];
             $GeneID_all{$GeneID} = $line[$Exprcol];
        }
   }
  close IN1;
}

open LIST,">$output";
print LIST "GeneID\t";
my $num;
for($num = 0;$num <= $xls_vol;$num++){
     print LIST "$sample_name[$num]_expression\t";
}
print LIST "\n";
my ($key1,$key2);
foreach $key1(sort keys %GeneID_all){
      print LIST "$key1\t";
       foreach $key2(@sample_name){
        if (exists $hash{$key1}{$key2}){
           print LIST "$hash{$key1}{$key2}\t";
        }else{
            $hash{$key1}{$key2} = 0.01 ;
            print LIST "$hash{$key1}{$key2}\t";
        }
      }
     print LIST "\n";
}
