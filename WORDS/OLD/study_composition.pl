#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use DBI;
use DBD::Sqlite;

# DEFAULTS
my $WINDOW = 500;
my $COMPOSITION = 'CTGTG,GTTTG,AAGCA,GCTTT,CTGTT,TCAGA,CTGAA';
my $DIR =  '/Volumes/PACKY/FERG/';
my $DB = 'WORDS';

my $result = GetOptions ('db|db=s' => \$DB,
                         'window|w=i' => \$WINDOW,
                         'composition|c=s' => \$COMPOSITION,
                         'directory|d=s' => \$DIR);

chdir($DIR);

my $DBH = DBI->connect("dbi:SQLite:dbname=$DB","","");

my $INCREMENT = get_last_increment_id();

# POPULATE TABLE COMPOSITION
my $COMPOSITION_ID = populate_table_composition();

# SELECT ALL POSITION FOR ONE WORD
my @WORDS = split(/\,/,$COMPOSITION);

my $POSITIONS = get_all_position($WORDS[0]);

my $OUT_TABLE = "COMPOSITION_ANALYSIS_ID_$COMPOSITION_ID";
open(OUT,">$OUT_TABLE");

foreach my $srn(keys %$POSITIONS) {
  foreach my $position(keys %{$POSITIONS->{$srn}}) {
    my $positive = analyze_position($srn,$position);
    if($positive) {
      $INCREMENT ++;
      my $start = $positive->{start};
      my $end = $positive->{end};
      print OUT "$INCREMENT\t$COMPOSITION_ID\t$srn\t$start\t$end\n";
      print "$INCREMENT\t$COMPOSITION_ID\t$srn\t$start\t$end\n";
    }
  }
}

close(OUT);

populate_window();

sub analyze_position {
  my $srn = shift;
  my $position = shift;
  print "Analyzing $srn\:$position\n";
  my $w = {};

  # BUILD CENTRAL WINDOW
  $w->{c}->{srn} = $srn;
  $w->{c}->{start} = $position - ($WINDOW/2);
  $w->{c}->{end} = $position + ($WINDOW/2);

  # BUILD LEFT WINDOW
  $w->{l}->{srn} = $srn;
  $w->{l}->{start} = $position - $WINDOW;
  $w->{l}->{end} = $position;

  # BUILD RIGHT WINDOW
  $w->{r}->{srn} = $srn;
  $w->{r}->{start} = $position + $WINDOW;
  $w->{r}->{end} = $position;

  # SCAN WINDOWS
  if(scan_window($w->{c})) {
    return $w->{c};
  }
  if(scan_window($w->{l})) {
    return $w->{l};
  }
  if(scan_window($w->{r})) {
    return $w->{r}
  }
  return undef;
}

sub scan_window {
  my $w = shift;
  for(my $i=1;$i<=$#WORDS;$i++) {
    my $word = $WORDS[$i];
    my $counts = look_for_word($word,$w->{srn},$w->{start},$w->{end});
    return undef unless $counts;
  }
  return 1;
}

sub look_for_word {
  my $word = shift;
  my $srn = shift;
  my $start = shift;
  my $end = shift;
  my $sth = $DBH->prepare("SELECT count(*) c FROM match_$word WHERE srn = \'$srn\' AND position >= $start AND position <= $end");
  $sth->execute;
  my $res = $sth->fetchrow_hashref;
  return $res->{c};
}

sub populate_table_composition {
  my $sth = $DBH->prepare_cached("INSERT INTO composition (composition,length) VALUES (?,?)");
  $sth->execute($COMPOSITION,$WINDOW);

  my $id = $DBH->last_insert_id('','','','');
  return $id;
}

sub get_all_position {
  my $word = shift;
  my $href;
  my $sth = $DBH->prepare("SELECT DISTINCT srn, position FROM match_$word WHERE srn=\'25\' ORDER BY position");
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    $href->{$row->{srn}}->{$row->{position}} ++;
  }
  return $href;
}

sub get_last_increment_id {
  my $sth = $DBH->prepare('SELECT max(id) max_id FROM window');
  $sth->execute;
  my $row = $sth->fetchrow_hashref;
  return $row->{max_id};
}

sub populate_window {
  system("sqlite3 -separator \"\t\" $DB \'.import $OUT_TABLE window\'");
}
