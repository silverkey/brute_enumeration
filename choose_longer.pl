#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

# This script is used to normalize a sequence fasta file downloaded by using
# BioMart. For example if you ask for UTR3' in the cases of multiple transcripts
# you will get a sequence for each transcript of the same gene.
# This program filters out all the sequences smaller than 10 bp and retain only
# the longer UTR for each gene.

my $USAGE = "\n\tperl $0 [fasta file]\n\n";

die $USAGE unless scalar(@ARGV) == 1;

my $in = Bio::SeqIO->new(-file => $ARGV[0],
												 -format => 'fasta');

$ARGV[0] =~ s/\.fa//;

my $out = Bio::SeqIO->new(-file => ">$ARGV[0]\_longer.fa",
													-format => 'fasta');

my $HASH;

while(my $seq = $in->next_seq) {
	my $id = $seq->id;
	my $string = $seq->seq;
	push(@{$HASH->{$id}},$string);
}

foreach my $id(keys %$HASH) {
	my @seq = @{$HASH->{$id}};
	my @sorted = sort{length($b) <=> length($a)} @seq;
	my $seq = Bio::Seq->new(-id => $id,
													-seq => $sorted[0]);
	$out->write_seq($seq) if length($sorted[0])>=10;
}

