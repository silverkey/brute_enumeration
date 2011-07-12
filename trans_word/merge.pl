#!/usr/bin/perl
use strict;
use warnings;
use Data::Dumper;

print join("\t",qw(ensgid mgi sym ugid ntr max_occ occ ntissues thyroid brain morula stages),"\n");

my $ugfile = 'unigene_mm_190_profile.txt';
my $ensgugfile = 'ENSGID_UGID';
my $enstugfile = 'ENSTRID_UGID';
my $ensntrfile = 'ENSGID_TRC_BIO';
my $ensmgisymfile = 'ENSGID_MGIID_SYM';
my $countsfile = 'MM_cdna.ATGCCG.xls';

my $unigene = read_unigene($ugfile);
my $ensgug = read_ens_ug($ensgugfile);
my $enstug = read_ens_ug($enstugfile);
my $ensgntr = read_ens_ntr($ensntrfile);
my $ensmgisym = read_ens_mgi_sym($ensmgisymfile);
my $counts = read_counts($countsfile);

foreach my $eid (keys %$counts) {
  print_info($eid);
}

sub print_info {
  my $id = shift;
  my $mgi = $ensmgisym->{$id}->{mgi};
  my $symbol = $ensmgisym->{$id}->{symbol};
  my $ugid = get_ug($id);
  my $ntr = $ensgntr->{$id};
  my $max = $counts->{$id}->{max};
  my $occ = $counts->{$id}->{counts};
  my $ntis = get_ntis($ugid);
  my $thyr = check_tissue($ugid,'thyroid');
  my $brain = check_tissue($ugid,'brain');
  my $morula = check_stage($ugid,'morula');
  my $stages = get_stages($ugid);
  print join("\t",$id,$mgi,$symbol,$ugid,$ntr,$max,$occ,$ntis,$thyr,$brain,$morula,$stages,"\n");
}

sub get_stages {
  my $ids = shift;
  my @ids = split(/\,/,$ids);
  my @tiss;
  foreach my $id(@ids) {
    push(@tiss,keys(%{$unigene->{$id}->{stage}}));
  }
  my %count;
  map { $count{$_}++ } @tiss;
  return join(',',keys(%count));
}

sub check_stage {
  my $ids = shift;
  my $tissue = shift;
  my @ids = split(/\,/,$ids);
  my @tiss;
  foreach my $id(@ids) {
    push(@tiss,keys(%{$unigene->{$id}->{stage}}));
  }
  my %count;
  map { $count{$_}++ } @tiss;
  return 'YES' if exists $count{$tissue};
  return 'NO';
}

sub check_tissue {
  my $ids = shift;
  my $tissue = shift;
  my @ids = split(/\,/,$ids);
  my @tiss;
  foreach my $id(@ids) {
    push(@tiss,keys(%{$unigene->{$id}->{tissue}}));
  }
  my %count;
  map { $count{$_}++ } @tiss;
  return 'YES' if exists $count{$tissue};
  return 'NO';
}

sub get_ntis {
  my $ids = shift;
  my @ids = split(/\,/,$ids);
  my @tiss;
  foreach my $id(@ids) {
    push(@tiss,keys(%{$unigene->{$id}->{tissue}}));
  }
  my %count;
  map { $count{$_}++ } @tiss;
  return scalar(keys(%count));
}

sub get_ug {
  my $id = shift;
  return 'NA' unless exists $ensgug->{$id};
  my $ug = join(',',keys(%{$ensgug->{$id}}));
  return $ug;
}

sub read_unigene {
  my $file = shift;
  my $href = {};
  open(IN,$file);
  my $head = <IN>;
  while(my $row = <IN>) {
    chomp($row);
    my @f = split(/\t/,$row);
    next unless $f[2] >= 1;
    $href->{$f[0]}->{stage}->{$f[1]} = $f[2] if $f[4] eq 'Developmental Stage';
    $href->{$f[0]}->{tissue}->{$f[1]} = $f[2] if $f[4] eq 'Body Sites';
  }
  return $href;
}

sub read_ens_ug {
  my $file = shift;
  my $href = {};
  open(IN,$file);
  my $head = <IN>;
  while(my $row = <IN>) {
    chomp($row);
    my @f = split(/\t/,$row);
    next unless $f[1];
    $href->{$f[0]}->{$f[1]} ++;
  }
  return $href;
}

sub read_ens_ntr {
  my $file = shift;
  my $href = {};
  open(IN,$file);
  my $head = <IN>;
  while(my $row = <IN>) {
    chomp($row);
    my @f = split(/\t/,$row);
    next unless $f[1];
    $href->{$f[0]} = $f[1];
  }
  return $href;
}

sub read_ens_mgi_sym {
  my $file = shift;
  my $href = {};
  open(IN,$file);
  my $head = <IN>;
  while(my $row = <IN>) {
    chomp($row);
    my @f = split(/\t/,$row);
    $href->{$f[0]}->{mgi} = $f[1] || 'NA'; #if $f[1];
    $href->{$f[0]}->{symbol} = $f[2] || 'NA' #if $f[2];
  }
  return $href;
}

sub read_counts {
  my $file = shift;
  my $href = {};
  open(IN,$file);
  my $head = <IN>;
  while(my $row = <IN>) {
    chomp($row);
    my @f = split(/\t/,$row);
    if(exists $href->{$f[1]}->{counts}){
      $href->{$f[1]}->{counts} .= ",$f[3]";
    }
    else {
      $href->{$f[1]}->{counts} = $f[3];
    }
    if(exists $href->{$f[1]}->{max}){
      $href->{$f[1]}->{max} = $f[3] if $f[3] > $href->{$f[1]}->{max};
    }
    else {
      $href->{$f[1]}->{max} = $f[3];
    }
  }
  return $href;
}
