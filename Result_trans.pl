#!/usr/bin/perl

use strict;
use warnings;

my $sampleid = $ARGV[0];

# Process genes.results file
open(my $genes_in, "<", "${sampleid}.genes.results") or die "Cannot open genes.results file: $!";
open(my $genes_out, ">", "${sampleid}.gene.fpkm.xls") or die "Cannot create gene.fpkm.xls file: $!";
while (my $line = <$genes_in>) {
    chomp($line);
    my @fields = split("\t", $line);
    if ($fields[6] != 0.00) {
        print $genes_out join("\t", @fields[0,1,2,4,6]) . "\n";
    }
}
close($genes_in);
close($genes_out);

# Process isoforms.results file
open(my $isoforms_in, "<", "${sampleid}.isoforms.results") or die "Cannot open isoforms.results file: $!";
open(my $isoforms_out, ">", "${sampleid}.transcript.fpkm.xls") or die "Cannot create transcript.fpkm.xls file: $!";
while (my $line = <$isoforms_in>) {
    chomp($line);
    my @fields = split("\t", $line);
    if ($fields[6] != 0.00) {
        print $isoforms_out join("\t", @fields[0,1,2,4,6]) . "\n";
    }
}
close($isoforms_in);
close($isoforms_out);