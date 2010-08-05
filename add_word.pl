#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

# This script is just for testing purposes it will add the word ATCGCGT
# to each sequence in the fasta file you will give as param

my $USAGE = "\n\tperl $0 [fasta file]\n\n";

die $USAGE unless scalar(@ARGV) == 1;

my $in = Bio::SeqIO->new(-file => $ARGV[0],
												 -format => 'fasta');

$ARGV[0] =~ s/\.fa//;

my $out = Bio::SeqIO->new(-file => ">$ARGV[0]\_mod.fa",
													-format => 'fasta');

while(my $seq = $in->next_seq) {
	my $new = $seq->seq();
	$new .= 'ATCGCGT';
	$seq->seq($new);
	$out->write_seq($seq);
}
