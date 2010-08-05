#!/usr/bin/perl

use strict;
use warnings;

use Cwd;
use Data::Dumper;
use File::Copy;
use Statistics::Descriptive;
use Bio::SeqIO;
use Bio::Seq;

# This program is the wrapper to run brute_enumeration.pl and compare_words.pl together.
# The idea is to compare two list of sequences while testing for both:
# 1) the enrichment of a particular word in the single files
# 2) the enrichment of a particular word in one file in respect to the other
# This is important because with motif finding algorithm we may have the problem of:
# - to be too stringent coming out with no words
# - to be too little stringent coming out with an improbable number of words
# - to come out with two different long list of word from two distinct set of sequences without
#   to know if something is enriched in one and not enriched in the other
# - to choose a background composition of our set of sequences and/or to create it by using
#   some really obscure program.
# ubik.pl overcome all these problems by using a not-so-stringent classical statics based
# on the randomization of the sequences, the z-scores and the proportional test.
# The combination of the statistics on 2 different levels (single files + comparisons) helps
# in mantaining low and higly significant the number of possible results.
# This program produce a lot of folders and files. The final one and the most important,
# probably the only interesting for the final user is the file which name start with
# 'stat_ubik_'.
# The name 'ubik' comes from one of the SCI-FI fixation of mine ;-)

my $USAGE = "\n\tperl $0 [fasta file 1] [fasta file 2] [length of words] [minimum gene number] [number of shuffling]\n\n";
die $USAGE unless scalar(@ARGV) == 5;

my $FASTA1 = $ARGV[0];
my $FASTA2 = $ARGV[1];
my $WL = $ARGV[2];
my $MIN_GEN = $ARGV[3];
my $NUM_SHUFF = $ARGV[4];

my $file_cou_compare = "stat_$FASTA1\_$FASTA2\_$WL\_counts.xls.csv";
$file_cou_compare =~ s/\.fa//g;

my $SWD = getcwd();

my $NWD1 = $SWD.'/'.$FASTA1.'_NSH_'.$NUM_SHUFF.'_WL_'.$WL;
$NWD1 =~ s/\.fa//;

my $NWD2 = $SWD.'/'.$FASTA2.'_NSH_'.$NUM_SHUFF.'_WL_'.$WL;
$NWD2 =~ s/\.fa//;

my $file_cou_brute1 = "stat_$FASTA1\_$NUM_SHUFF\_$WL\_counts.xls.csv";
$file_cou_brute1 =~ s/\.fa//;

my $file_cou_brute2 = "stat_$FASTA2\_$NUM_SHUFF\_$WL\_counts.xls.csv";
$file_cou_brute2 =~ s/\.fa//;


system("perl brute_enumeration.pl $FASTA1 $NUM_SHUFF $WL $MIN_GEN");

system("perl brute_enumeration.pl $FASTA2 $NUM_SHUFF $WL $MIN_GEN");

system("perl compare_words.pl $FASTA1 $FASTA2 $WL $MIN_GEN");


my $F1 = open_file("$NWD1\/$file_cou_brute1");
my $F2 = open_file("$NWD2\/$file_cou_brute2");
my $C = open_file("$file_cou_compare");

my $STAT = "stat_ubik_$ARGV[0]\_$ARGV[1]\_$ARGV[2]\_$ARGV[3]\_$ARGV[4]\.csv";
open(UBIK,">$STAT");

my $STAT_STRING = '"","seq_on_file1","seq_with_word_file1","num_occ_file1","seq_on_file2","seq_with_word_file2","num_occ_file2","p.val_comparison","adj.p.val_comparison",';
$STAT_STRING .= '"p.val_shuff_file1","adj.p.val_shuff_file1","p.val_shuff_file2","adj.p.val_shuff_file2","advise"'."\n";

print UBIK $STAT_STRING;

foreach my $word(keys %$F1) {

	my @string = @{$C->{$word}};

	push(@string,$F1->{$word}->[7]);
	push(@string,$F1->{$word}->[8]);

	if(exists $F2->{$word}) {
		push(@string,$F2->{$word}->[7]);
		push(@string,$F2->{$word}->[8]);
	}
	else {
		push(@string,'"NA"');
		push(@string,'"NA"');
	}

	if($string[8] <= 0.05 && $string[10] <= 0.05) {
		if($string[12] eq '"NA"' || $string[12] >0.05) {
			$string[13] = '"GOOD enriched file 1"';
		}
		else {
			$string[13] = '"GOOD ? enriched in all comparison"';
		}
	}
	elsif($string[12] <= 0.05 && $string[10] <= 0.05) {
		if($string[8] eq '"NA"' || $string[8] >0.05) {
			$string[13] = '"GOOD enriched file 2"';
		}
		else {
			$string[13] = '"GOOD ? enriched in all comparisons"';
		}
	}
	else {
		$string[13] = '"BAD"';
	}

	my $string = join(',',@string);
	print UBIK $string."\n";
}

foreach my $word(keys %$F2) {

	next if exists $F1->{$word};

	my @string = @{$C->{$word}};

	push(@string,'"NA"'); # because F1 doesn't exists
	push(@string,'"NA"'); # because F1 doesn't exists

	push(@string,$F2->{$word}->[7]);
	push(@string,$F2->{$word}->[8]);

	if($string[8] <= 0.05 && $string[12] <= 0.05) {
		if($string[10] eq '"NA"' || $string[10] >0.05) {
			$string[13] = '"GOOD"';
		}
	}
	else {
		$string[13] = '"BAD"';
	}

	my $string = join(',',@string);
	print UBIK $string."\n";
}

sub open_file {
	my $file = shift;
	my $href;
	open(F,$file) or die "\n\n$! $file\n\n";
	my $h = <F>;
	while(my $row = <F>) {
		chomp($row);
		my @field = split(/\,/,$row);
		$href->{$field[0]} = \@field;
	}
	return $href;
}
