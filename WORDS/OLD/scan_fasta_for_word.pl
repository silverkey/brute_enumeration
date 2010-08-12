#!/usr/bin/perl

use strict;
use warnings;

use Cwd;
use Data::Dumper;
use Bio::SeqIO;
use Getopt::Long;
use DBI;
use DBD::Sqlite;

# Actually the database to populate is so constituted:
# CREATE TABLE word(id integer primary key autoincrement, word varchar key);
# CREATE TABLE match(id integer primary key autoincrement, word_id integer key, srn varchar key, position integer key, strand enum(1,-1));

# DEFAULTS
my $FASTA = 'Danio_rerio.Zv8.59.dna.toplevel.fa';
my $DIR = '/Volumes/PACKY/FERG/';
my $DB = 'WORDS';

my $WORD = 'AAAAA';

my $result = GetOptions ('fasta|f=s' => \$FASTA,
                         'directory|d=s' => \$DIR,
                         'db|db=s' => \$DB,
                         'word|w=s' => \$WORD);

chdir($DIR);

my $DBH = DBI->connect("dbi:SQLite:dbname=$DB","","");

my $word_id = populate_table_word();
print "INSERTED ID $word_id\n";

my $TMP = 'TMP.txt';
open(TMP,">$TMP");

my $SEQIO = Bio::SeqIO->new(-file => $FASTA,
                            -format => 'fasta');

while(my $seq = $SEQIO->next_seq) {
  my $srn = $seq->id;
  print "Arrived at $srn\n";
  next unless $seq->length >= 10000;
  print "SCANNING $srn *************************** \n";
#  my $position = get_positions($seq->seq);
}

sub get_positions {
  my $string = shift;
  while(my $match = $string =~ /$WORD/g) {
      print "Match at position ".(length($`)+1)."\n";
  }
}

sub populate_table_word {

  my $sth = $DBH->prepare_cached(qq{
    INSERT INTO word
    (word) VALUES (?)
  });
  $sth->execute($WORD);

  my $id = $DBH->last_insert_id('','','','');
  return $id;
}
