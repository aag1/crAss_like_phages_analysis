#!/usr/bin/perl
use strict;
use warnings;


### BASED ON https://github.com/simroux/VirSorter/blob/master/Scripts/Step_1_contigs_cleaning_and_gene_prediction.pl


### INPUT PARAMETERS
my $fastaF = $ARGV[0];
my $outF = $ARGV[1];
my $min_overlap = $ARGV[2];


### READING INPUT FASTA
open FILE1, '<', $fastaF or die "Failed to open $fastaF: $!\n";

my %seq_base;
my $id_seq = "";

while (my $line = <FILE1>) {

	chomp($line);

	if ($line =~ /^>(\S*)/) { $id_seq = $1 }
	else { $seq_base{$id_seq} .= $line }

}

close FILE1 or die "Failed to close $fastaF: $!\n";


### DETECTION OF CIRCULAR SEQUENCES
open FILE2, '>', $outF or die "Failed to open $outF: $!\n";

for my $k (keys %seq_base) {

    my $s = $seq_base{$k};

    my $prefix = substr($s, 0, $min_overlap);

    if ($s =~ /.+($prefix.*?)$/) {

        my $suffix = $1;

        my $l = length($suffix);

        my $test = substr($s, 0, $l);

        if ($suffix eq $test) {

            print FILE2 "$k\t$l\n";

        }
    }
}

close FILE2 or die "Failed to close $outF: $!\n";

# Note that for $s = 'AAbcAAb....AAbcAAb' and $min_overlap = 2, suffix AAb will be detected.
