#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;

my $in = Bio::SeqIO->new(-file => $ARGV[0],
                         -format => 'fasta');

$ARGV[0] =~ s/\.fa//;

my $out = Bio::SeqIO->new(-file => ">$ARGV[0]\_longer.fa",
                          -format => 'fasta');

my $HASH;

while(my $seq = $in->next_seq) {
  my $id = $seq->id;
  my $string = $seq->seq;
  push(@{$HASH->{$id}},$string) unless $string =~ /^Sequence/;
}

foreach my $id(keys %$HASH) {
  my @seq = @{$HASH->{$id}};
  my @sorted = sort{length($b) <=> length($a)} @seq;
  my $seq = Bio::Seq->new(-id => $id,
                          -seq => $sorted[0]);
  $out->write_seq($seq) if length($sorted[0])>=10;
}
