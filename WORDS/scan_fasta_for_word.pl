#!/usr/bin/perl

use strict;
use warnings;

use Data::Dumper;
use Getopt::Long;
use DBI;
use DBD::Sqlite;

# Actually the database to populate is so constituted:
# CREATE TABLE word(id integer primary key autoincrement, word varchar key);
# CREATE TABLE match(id integer primary key autoincrement, word_id integer key, srn varchar key, position integer key, strand varchar key);

# DEFAULTS
my $FASTA = 'Danio_rerio.Zv8.59.dna.toplevel.fa';
my $DIR = '/Volumes/PACKY/FERG/';
my $DB = 'WORDS';
my $WORD = 'CTGAA';

my $result = GetOptions ('fasta|f=s' => \$FASTA,
                         'directory|d=s' => \$DIR,
                         'db|db=s' => \$DB,
                         'word|w=s' => \$WORD);

chdir($DIR);

my $SRN = {};
my $INCREMENT;

my $DBH = DBI->connect("dbi:SQLite:dbname=$DB","","");

my $WORD_ID = 7;
#my $WORD_ID = populate_table_word();

#my $OUT = run_fuzznuc();
my $OUT = "FUZZNUC_$WORD_ID";

parse_fuzznuc_out();

populate_table_match();

populate_table_region();

system('rm *.temp');

sub populate_table_word {
  my $sth = $DBH->prepare("INSERT INTO word (word) VALUES ($WORD)");
  $sth->execute();
  my $id = $DBH->last_insert_id('','','','');
  return $id;
}

sub run_fuzznuc {
  my $out = "FUZZNUC_$WORD_ID";
  my $command = "fuzznuc -sequence $FASTA -pattern $WORD -rformat2 excel -outfile $out -complement";
  system($command);
  return $out;
}

sub parse_fuzznuc_out {
  open(IN,$OUT);
  while(my $line = <IN>) {
    $INCREMENT ++;
    next if $line =~ /^SeqName/;
    chomp($line);
    my($srn,$start,$end,$score,$strand,$pattern,$mismatch) = split(/\t/,$line);
    if($srn =~ /\_/) {
      push(@{$SRN->{scaffold}},"$INCREMENT\t$srn\t$start\t$strand");
    }
    else {
      push(@{$SRN->{$srn}},"$INCREMENT\t$start\t$strand");
    }
  }
  close(OUT);
}

sub populate_table_match {
  foreach my $tabname(keys %$SRN) {
    create_table($tabname);
    create_file($tabname);
    load_data($tabname);
  }
}

sub create_table {
  my $name = shift;
  my $tabname = "match_$WORD\_$name";
  $DBH->do("DROP TABLE IF EXISTS $tabname");
  my $sql;
  if($name eq 'scaffold') {
    $sql = "CREATE TABLE $tabname(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            srn VARCHAR KEY,
            position INTEGER KEY,
            strand VARCHAR KEY)";
  }
  else {
    $sql = "CREATE TABLE $tabname(
            id INTEGER PRIMARY KEY AUTOINCREMENT,
            position INTEGER KEY,
            strand VARCHAR KEY)";
  }
  $DBH->do($sql);
}
    
sub create_file {
  my $name = shift;
  open(OUT,">$name\.temp");
  my $res = $SRN->{$name};
  foreach my $row(@$res) {
    print OUT $row."\n";
  }
  close(OUT);
}

sub load_data {
  my $name = shift;
  my $tabname = "match_$WORD\_$name";
  my $filename = "$name\.temp";
  system("sqlite3 -separator \"\t\" $DB \'.import $filename $tabname\'");
}

sub populate_table_region {
  my $sql = "CREATE TABLE region(
             id INTEGER PRIMARY KEY AUTOINCREMENT,
             name VARCHAR KEY)";
  $DBH->do($sql);
  foreach my $name(keys %$SRN) {
    my $sth = $DBH->prepare("INSERT INTO region (name) VALUES (\'$name\')");
    $sth->execute;
  }
}
