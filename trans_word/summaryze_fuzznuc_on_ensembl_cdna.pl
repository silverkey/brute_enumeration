#!/usr/bin/perl
use strict;
use warnings;
use Bio::SeqIO;
use Data::Dumper;

my $fuzznuc = 'MM_cdna.ATGCCG.gff';
my $cdna = 'Mus_musculus.NCBIM37.62.cdna.all.fa';
my $out = "$fuzznuc";
$out =~ s/.gff/.xls/;

open(OUT,">$out");
print OUT join("\t",'transcript_id','gene_id','length','motif_occurrences','ratio');
print OUT "\n";

my $info = get_ensembl_cdna_info($cdna);
my $res = get_fuzznuc_gff_res($fuzznuc);

foreach my $t(keys %$info) {
  my $g = $info->{$t}->{gene};
  my $l = $info->{$t}->{length};
  my $o = 0;
  my $r = 'NA';
  if(exists $res->{$t}) {
    $o = $res->{$t};
    $r = $o/$l;
  }
  print OUT join("\t",$t,$g,$l,$o,$r);
  print OUT "\n";
}

sub get_fuzznuc_gff_res {
  my $gff = shift;
  my $href = {};
  open(IN,$gff);
  while(my $row = <IN>) {
    next if $row =~ /^\#/;
    my @f = split(/\t/,$row);
    $res->{$f[0]} ++;
  }
  return $res;
}

sub get_ensembl_cdna_info {
  my $file = shift;
  my $href = {};
  my $seqio = Bio::SeqIO->new(-file => $file,
                              -format => 'fasta');
  while(my $seq = $seqio->next_seq) {
    my $desc = $seq->desc;
    $desc =~ s/^.+gene\:(.+)$/$1/;
    $href->{$seq->id}->{gene} = $desc;
    $href->{$seq->id}->{length} = $seq->length;
  }
  return $href;
}
