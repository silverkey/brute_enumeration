#!/usr/bin/perl

use strict;
use warnings;

use Cwd;
use Data::Dumper;
use File::Copy;
use Statistics::Descriptive;
use Bio::SeqIO;
use Bio::Seq;

# This program take in input 2 fasta containing different sequences and look for all the words
# of size [param] present into the sequences annd calculate statistics for testing the enrichent
# of each word in one group in respect to the other. It produces a lot of files, the most important
# of them is the one which name start with 'stat_'

my $USAGE = "\n\tperl $0 [fasta file 1] [fasta file 2] [length of words] [minimum gene number]\n\n";
die $USAGE unless scalar(@ARGV) == 4;

my $FASTA1 = $ARGV[0];
my $FASTA2 = $ARGV[1];
my $WL = $ARGV[2];
my $MIN_GEN = $ARGV[3];

my $file_cou = "$FASTA1\_$FASTA2\_$WL\_counts.xls";
$file_cou =~ s/\.fa//g;

my($COUNTS1,$GENES1,$NUMSEQ1) = counts($FASTA1);

my($COUNTS2,$GENES2,$NUMSEQ2) = counts($FASTA2);

analize_res();

prop_test();

write_log();



=head2

 Title   : counts

 Usage   : counts()

 Function: count the occurrences of all the motif in the given fasta

 Returns : it populate the $COUNTS hashref

 Args    :

 Note    :

=cut

sub counts {

	my $FASTA = shift;
	my $COUNTS;
	my $GENES;
	my $NUMSEQ;

	my $file_occ = "$FASTA\_$WL\_occurrences.txt";
	$file_occ =~ s/\.fa//g;
	open(OCC,">$file_occ");
	print OCC "ref_id\tword\tposition\n";

	my $seqio = Bio::SeqIO->new(-file => $FASTA,
															-format => 'fasta');

	while(my $seq = $seqio->next_seq) {

		$NUMSEQ ++;

		my $string = $seq->seq;
		my $id = $seq->id;

		my $start = 0;
		my $last = length($string) - $WL;

		while($start <= $last) {

			my $substring = substr($string,$start,$WL);

			$start ++;

			next if $substring =~ /N/;

			$COUNTS->{$substring} ++;
			$GENES->{$substring}->{$id} = 1;

			print OCC "$id\t$substring\t".($start-1)."\n";
		}
	}
	close(OCC);
	return($COUNTS,$GENES,$NUMSEQ);
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

	open(COU,">$file_cou");
	print COU "\ttot_gene_file1\tgene_with_word_file1\tnum_occ_file1\ttot_gene_file2\tgene_with_word_file2\tnum_occ_file2\tp.value\tadj.p.value\n";

	foreach my $word(keys %$COUNTS1) {

		my $counts1 = $COUNTS1->{$word} || '0';
		my $num_seq1 = scalar(values %{$GENES1->{$word}}) || '0';

		my $counts2 = $COUNTS2->{$word} || '0';
		my $num_seq2 = scalar(values %{$GENES2->{$word}}) || '0';

		if($num_seq1 >= $MIN_GEN || $num_seq2 >= $MIN_GEN) {
			print COU "$word\t$NUMSEQ1\t$num_seq1\t$counts1\t$NUMSEQ2\t$num_seq2\t$counts2\tNA\tNA\n";
		}
	}

	foreach my $word(keys %$COUNTS2) {

		next if exists $COUNTS1->{$word};

		my $counts1 = $COUNTS1->{$word} || '0';
		my $num_seq1 = scalar(values %{$GENES1->{$word}}) || '0';

		my $counts2 = $COUNTS2->{$word} || '0';
		my $num_seq2 = scalar(values %{$GENES2->{$word}}) || '0';

		if($num_seq1 >= $MIN_GEN || $num_seq2 >= $MIN_GEN) {
			print COU "$word\t$NUMSEQ1\t$num_seq1\t$counts1\t$NUMSEQ2\t$num_seq2\t$counts2\tNA\tNA\n";
		}
	}

	close(COU);
}

sub prop_test {

	my $script = 'options(echo=FALSE)'."\n";

	$script .= 
'
f<-read.table(file="'.$file_cou.'")
for(i in 1:nrow(f)) {
    test<-prop.test(c(f$gene_with_word_file1[i],f$gene_with_word_file2[i]),c(f$tot_gene_file1[i],f$tot_gene_file2[i]))
    f$p.value[i]<-test$p.value
    f$adj.p.value[i]<-(test$p.value*nrow(f))
}
write.csv(f,file="stat_'.$file_cou.'.csv")
';

  my $filename = "SCRIPT_$ARGV[0]\_$ARGV[1]\_$ARGV[2]\_$ARGV[3]";

	open(SCRIPT,">$filename");
	print SCRIPT $script;

	system("R CMD BATCH --no-save $filename");
	system("rm $filename\*");
}

sub write_log {
	my $file_log = "$FASTA1\_$FASTA2\_$WL\_compare.log";
	$file_log =~ s/\.fa//g;
	open(LOG,">$file_log");
	print LOG "FASTA FILE 1:\t$FASTA1\n";
	print LOG "FASTA FILE 2:\t$FASTA2\n";
	print LOG "WORD LENGTH:\t$WL\n\n\n";
}
