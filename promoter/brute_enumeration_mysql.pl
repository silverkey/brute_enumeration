#!/usr/bin/perl
use strict;
use warnings;
use Cwd;
use File::Copy;
use Bio::Seq;
use DBI;
use Bio::SeqIO;

# This program take in input a fasta file and produce [param] shuffled file of each sequence.
# Than it looks for the word of length [param] in real and shuffled sequences.
# Than it build statistics based on the composition of the real sequences and the shuffled ones
# It produce a directory with a lot of files.... the most important one is the file which
# name starts with 'stat_'. IN THIS FILE THERE ARE THE P-VALUE COMPUTED BY USING THE Z-SCORE
# CALCULATED ON THE DISTRIBUITION OF THE SHUFFLED SEQUENCES FOR EACH MOTIF.

#-----------------------------
# INITIALIZATION
#-----------------------------
my $usage = "\nUsage: perl $0\n
 - [fasta file] 
 - [number of shufflings] 
 - [length of words] 
 - [create db (1 or 0)] \n\n";       
die $usage unless scalar(@ARGV) == 4;
my $fasta = $ARGV[0];
my $nsh = $ARGV[1];
my $wl = $ARGV[2];
my $create = $ARGV[3];
die "You set no shufflings while cting the database de-novo,
therefore the db must exists and contains already shuffling!\n"
if $nsh<1 && $create==1;
# DIRECTORIES WORK
my $startdir = getcwd();
my $workdir = $startdir.'/'.$fasta.'_NSH_'.$nsh.'_WL_'.$wl;
$workdir =~ s/\.fa//;
mkdir($workdir) or die "\nmkdir failed: $!\n";
chdir($workdir) or die "\nchdir failed: $!\n";
copy($startdir.'/'.$fasta,$workdir.'/'.$fasta)
or die "\ncopy failed: $!\n";
# OUTPUT TABLES NAME DEFINITION
my $counts_tab = "$fasta\_$nsh\_$wl\_counts.xls";
$counts_tab =~ s/\.fasta$//;
$counts_tab =~ s/\.fa$//;
my $occurrences_tab = "$fasta\_$nsh\_$wl\_occurrences.txt";
$occurrences_tab =~ s/\.fasta$//;
$occurrences_tab =~ s/\.fa$//;
# DATABASE SPECIFIC SETTINGS
my $db = $fasta.'_NSH_'.$nsh.'_WL_'.$wl;
$db =~ s/\.fa//;
my $usr = 'mysql_dev';
my $pwd = 'riiGbs';
my $host = 'localhost';

#-----------------------------
# ANALYSIS
#-----------------------------
create_db($db,$usr,$pwd,$host) if($create);
my $dbh = connect_to_db($db,$usr,$pwd,$host);
shuffle_fasta($fasta,$nsh);
populate_seq($fasta,$dbh);
populate_shuffled($fasta,$nsh,$dbh);
my($word,$last_inserted) = count_real($dbh,$wl);
count_shuffled($dbh,$wl,$word,$last_inserted,$nsh);
analize_res($dbh,$counts_tab);
write_log();
calculate_p_value($counts_tab);
system('rm *shuffled*');

#-----------------------------
# ROUTINES
#-----------------------------
=head2 shuffle_fasta

Title   : shuffle_fasta

Usage   : shuffle_fasta($fasta, $shuffling)

Function: create the shuffled versions of fasta

Returns : nothing, the file are created in the working dir
          named shuffled_[n].fa where n is the number of the cicle

Args    : 1) The fasta to shuffle
          2) The number of shifflings

Note    : This function use the program "shuffle" from SQUID
          the Sean Eddy's C toolkit. Much faster than perl.
          http://selab.janelia.org/software.html
          
          The seed used to generate the random numbers are coming
          from a counter into the script corresponding to the cicle.
          
          Whatever you need probably Sean Eddy already did....
          and it is working better of course ;-)

=cut

sub shuffle_fasta {
  my $fasta = shift;
  my $nsh = shift;
  my $fastaname = "$fasta";
  $fastaname =~ s/\.fasta$//;
  $fastaname =~ s/\.fa$//;
  my $x;
  for($x=1;$x<=$nsh;$x++) {
    my $sfasta = "$fastaname\_shuffled\_$x\.fa";
    system("shuffle --seed $x $fasta > $sfasta");
  }
}

sub populate_seq {
  my $fasta = shift;
  my $dbh = shift;
  my $tmp = "$fasta\.table\.tmp";
  my $seqio = Bio::SeqIO->new(-file => $fasta,
                              -format => 'fasta');
  open(OUT,">$tmp");
  my $count = 1;
  while(my $seq = $seqio->next_seq) {
    print OUT "$count\t".$seq->id."\t".$seq->seq."\n";
    $count ++;
  }
  close(OUT);
  my $command = "LOAD DATA LOCAL INFILE '$tmp' INTO TABLE seq";
  $dbh->do($command);
  system("rm $tmp");
}

sub populate_shuffled {
  my $source = shift;
  my $nsh = shift;
  my $dbh = shift;
  my $fastaname = "$source";
  $fastaname =~ s/\.fasta$//;
  $fastaname =~ s/\.fa$//;
  my $x;
  for($x=1;$x<=$nsh;$x++) {
    my $sfasta = "$fastaname\_shuffled\_$x\.fa";
    my $tmp = "$sfasta\.table\.tmp";
    my $seqio = Bio::SeqIO->new(-file => $sfasta,
                                -format => 'fasta');
    open(OUT,">$tmp");
    my $count = 1;
    while(my $seq = $seqio->next_seq) {
      print OUT $x."\t$count\t".$seq->seq."\n";
      $count ++;
    }
    close(OUT);
    my $command = "LOAD DATA LOCAL INFILE '$tmp' INTO TABLE shuffled";
    $dbh->do($command);
    system("rm $tmp");
  }
}

=head2 count_real

Title   : count_real

Usage   : count_real($dbh,$wl)

Function: count the occurrences of all the motif in the given fasta

Returns : it populate the tables word and seq_occurrences and returns
          the href with the words and the id of the last word.

Args    : 1) DBH
          2) The length of the word to count

Note    :

=cut

sub count_real {
  my $dbh = shift;
  my $wl = shift;
  my $words = {};
  my $words_counter = 0;

  my $occfile = 'occurrences.tmp';
  my $wordfile = 'words.tmp';
  open(OCCURRENCES,">$occfile");
  open(WORDS,">$wordfile");
    
  my $sth = $dbh->prepare('SELECT * FROM seq');
  $sth->execute;
  
  while(my $row = $sth->fetchrow_hashref) {
    my $id = $row->{id};
    my $string = $row->{seq};
    my $start = 0;
    my $last = length($string) - $wl;
    
    while($start <= $last) {
      my $substring = substr($string,$start,$wl);
      $start ++;
      next if $substring =~ /N/;

      unless (exists($words->{$substring})) {
        $words_counter ++;
        $words->{$substring} = $words_counter;
        print WORDS "$words_counter\t$substring\n";
      }
      print OCCURRENCES "$id\t".$words->{$substring}."\t$start\n";
    }
  }
  close(OCCURRENCES);
  close(WORDS);

  my $command = "LOAD DATA LOCAL INFILE '$occfile' INTO TABLE seq_occurrence";
  $dbh->do($command);
  system("rm $occfile");

  return($words,$words_counter);
}

=head2 count_shuffled

Title   : count_shuffled

Usage   : count_shuffled($dbh,$wl,$words_href,$last_inserted)

Function: count the occurrences of all the motif in the shuffled seq
          and generate a summary used to populate the table
          shuffled_counts

Returns : it populate the tables word and shuffled_counts.

Args    : 1) DBH
          2) The length of the word to count
          3) The words href started in the count_real function
          4) The id of the last word in the words_href
          5) Number of shufflings

Note    :

=cut

sub count_shuffled {
  my $dbh = shift;
  my $wl = shift;
  my $words = shift;
  my $words_counter = shift;
  my $nsh = shift;
  my $shuffling;

  my $wordfile = 'words.tmp';
  open(WORDS,">>$wordfile");
  my $summaryfile = 'summary.tmp';
  open(SUMMARY,">$summaryfile");
    
  for($shuffling=1;$shuffling<=$nsh;$shuffling++) {
    my $occhref = {};
    my $seqhref = {};
    my $sth = $dbh->prepare("SELECT * FROM shuffled WHERE shuffling = $shuffling");
    $sth->execute;  
    while(my $row = $sth->fetchrow_hashref) {
      my $id = $row->{seq_id};
      my $string = $row->{seq};
      my $start = 0;
      my $last = length($string) - $wl;      
      while($start <= $last) {
        my $substring = substr($string,$start,$wl);
        $start ++;
        next if $substring =~ /N/;  
        unless (exists($words->{$substring})) {
          $words_counter ++;
          $words->{$substring} = $words_counter;
          print WORDS "$words_counter\t$substring\n";
        }
        $occhref->{$words->{$substring}} ++;
        $seqhref->{$words->{$substring}}->{$id} ++;
      }
    }
    foreach my $word(keys %$occhref) {
      my $nocc = $occhref->{$word};
      my $nseq = scalar(keys %{$seqhref->{$word}});
      print SUMMARY "$shuffling\t$word\t$nocc\t$nseq\n";
    }
  }
  close(SUMMARY);
  my $command1 = "LOAD DATA LOCAL INFILE '$summaryfile' INTO TABLE shuffled_counts";
  $dbh->do($command1);
  system("rm $summaryfile");
  close(WORDS);
  my $command2 = "LOAD DATA LOCAL INFILE '$wordfile' INTO TABLE word";
  $dbh->do($command2);
  system("rm $wordfile");
}

sub connect_to_db {
  my $db = shift;
  my $usr = shift;
  my $pwd = shift;
  my $host = shift;
  my $dsn = 'dbi:mysql:'.$db;
  $dsn .= ':'.$host if $host; # IN THE CURRENT DBI POD VERSION THERE IS THE '@' IN THE PLACE OF ':'
  my $dbh = DBI->connect($dsn,$usr,$pwd,{PrintWarn=>1,PrintError=>1,RaiseError=>1}) or die $DBI::errstr;
  return $dbh;
}

sub create_db {
  my $db = shift;
  my $usr = shift;
  my $pwd = shift;
  my $host = shift;
  my $dbh = connect_to_db('',$usr,$pwd,$host);
  my $create = "CREATE DATABASE $db";
  $dbh->do($create);
  $dbh->disconnect;
  sleep(1);
  my @data = <DATA>;
  my $schema = join('',@data);
  open(SCHEMA,">schema.sql");
  print SCHEMA $schema;
  close(SCHEMA);
  system("mysql -u$usr -p$pwd -h$host $db < schema.sql");
  sleep(1);
  system('rm schema.sql');
}

=head2 analize_res

Title   : analize_res

Usage   : analize_res()

Function:

Returns :

Args    :

Note    :

=cut

sub analize_res {
  my $dbh = shift;
  my $counts_tab = shift;
  open(OUT,">$counts_tab");  
  my @header = qw(word n.occ n.seq z.occ p.occ adp.occ z.seq p.seq adp.seq);
  print OUT join("\t",@header)."\n";
  
  my $words = get_all_words($dbh);
  my $occreal = get_all_real_occ($dbh);
  my $seqreal = get_all_real_seq($dbh);
  my $countsshuff = get_all_shuff_count($dbh);

  foreach my $wordid(keys %$words) {
    
    my $word = $words->{$wordid};
    
    my $r = ($occreal->{$wordid} || 1);
    my $num_seq = ($seqreal->{$wordid} || 1);

    my $mean_c = ($countsshuff->{$wordid}->{AVGO} || 1);
    my $sd_c = ($countsshuff->{$wordid}->{DEVO} || 1);
    my $diff_c = $r - $mean_c;
    my $z_score_c = $diff_c / $sd_c;
    
    my $mean_g = ($countsshuff->{$wordid}->{AVGS} || 1);
    my $sd_g = ($countsshuff->{$wordid}->{DEVS} || 1);
    my $diff_g = $num_seq - $mean_g;
    my $z_score_g = $diff_g / $sd_g;
    
    print OUT "$word\t$r\t$num_seq\t";
    print OUT "$z_score_c\tNA\tNA\t"; #if $z_score_c;
    print OUT "$z_score_g\tNA\tNA\n"; #if $z_score_g;
  }
  close(OUT);
}

sub calculate_p_value { 
  my $counts_tab = shift; 
  my $script = 'options(echo=FALSE)'."\n";

  $script .=
  '
  f<-read.table(file="'.$counts_tab.'",header=T)
  f$p.occ = 1-pnorm(abs(f$z.occ))
  f$adp.occ <- nrow(f)*(1-pnorm(abs(f$z.occ)))
  f$p.seq <- 1-pnorm(abs(f$z.seq))
  f$adp.seq <- nrow(f)*(1-pnorm(abs(f$z.seq)))
  write.table(f,file="stat_'.$counts_tab.'.csv",sep=",",row.names=F,quote=F)
  ';
  
  my $filename = "SCRIPT.R";
  open(SCRIPT,">$filename");
  print SCRIPT $script;
  system("R CMD BATCH --no-save $filename");
}

sub get_all_words {
  my $dbh = shift;
  my $href = {};
  my $sth = $dbh->prepare('SELECT * FROM word');
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    $href->{$row->{id}} = $row->{word};
  }
  return $href;
}

sub get_all_real_occ {
  my $dbh = shift;
  my $href = {};
  my $sth = $dbh->prepare('SELECT word_id, COUNT(*) C FROM seq_occurrence GROUP BY word_id');
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    $href->{$row->{word_id}} = $row->{C};
  }
  return $href;
}

sub get_all_real_seq {
  my $dbh = shift;
  my $href = {};
  my $sth = $dbh->prepare('SELECT word_id, COUNT(DISTINCT seq_id) S FROM seq_occurrence GROUP BY word_id');
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    $href->{$row->{word_id}} = $row->{S};
  }
  return $href;
}

sub get_all_shuff_count {
  my $dbh = shift;
  my $href = {};
  my $sth = $dbh->prepare('SELECT word_id, AVG(sequences) AVGS, STD(sequences) DEVS, AVG(occurrences) AVGO, STD(occurrences) DEVO FROM shuffled_counts GROUP BY word_id;');
  $sth->execute;
  while(my $row = $sth->fetchrow_hashref) {
    $href->{$row->{word_id}}->{AVGS} = $row->{AVGS};
    $href->{$row->{word_id}}->{DEVS} = $row->{DEVS};
    $href->{$row->{word_id}}->{AVGO} = $row->{AVGO};
    $href->{$row->{word_id}}->{DEVO} = $row->{DEVO};
  }
  return $href;
}

sub write_log {
  my $file_log = "$fasta\_$nsh\_$wl\_enumeration.log";
  $file_log =~ s/\.fa//;
  open(LOG,">$file_log");
  print LOG "\nSequences analyzed with parameters:\n\n";
  print LOG "FASTA FILE:\t$fasta\n";
  print LOG "N SHUFFLING:\t$nsh\n";
  print LOG "WORD LENGTH:\t$wl\n";
}

__DATA__

CREATE TABLE seq(
  id INT(5) UNSIGNED NOT NULL,
  name VARCHAR(30) NOT NULL,
  seq TEXT NOT NULL,
  PRIMARY KEY id(id),
  KEY name(name)
);

CREATE TABLE shuffled(
  shuffling INT(5) UNSIGNED NOT NULL,
  seq_id INT(5) UNSIGNED NOT NULL,
  seq TEXT NOT NULL,
  KEY shuffling(shuffling),
  KEY seq_id(seq_id)
);

CREATE TABLE seq_occurrence(
  seq_id INT(5) UNSIGNED NOT NULL,
  word_id INT(20) UNSIGNED NOT NULL,
  pos INT(10) UNSIGNED NOT NULL,
  KEY seq_id(seq_id),
  KEY word_id(word_id),
  KEY pos(pos)
);

CREATE TABLE shuffled_counts(
  shuffling INT(5) UNSIGNED NOT NULL,
  word_id INT(20) UNSIGNED NOT NULL,
  occurrences INT(20) UNSIGNED NOT NULL,
  sequences INT(10) UNSIGNED NOT NULL,
  KEY shuffling(shuffling),
  KEY word_id(word_id),
  KEY sequences(sequences),
  KEY occurrences(occurrences)
);

CREATE TABLE word(
  id INT(20) UNSIGNED NOT NULL,
  word VARCHAR(100) NOT NULL,
  PRIMARY KEY id(id),
  KEY word(word)
);
