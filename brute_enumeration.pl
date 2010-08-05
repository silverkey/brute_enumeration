#!/usr/bin/perl

use strict;
use warnings;

use Cwd;
use Data::Dumper;
use File::Copy;
use Statistics::Descriptive;
use Bio::SeqIO;
use Bio::Seq;

# This program take in input a fasta file and produce [param] shuffled file of each sequence.
# Than it looks for the word of length [param] in real and shuffled sequences.
# Than it build statistics based on the composition of the real sequences and the shuffled ones
# It produce a directory with a lot of files.... the most important one is the file which
# name starts with 'stat_'. IN THIS FILE THERE ARE THE P-VALUE COMPUTED BY USING THE Z-SCORE
# CALCULATED ON THE DISTRIBUITION OF THE SHUFFLED SEQUENCES FOR EACH MOTIF.

my $USAGE = "\n\tperl $0 [fasta file] [number of shuffling] [length of words] [minimum number of seq]\n\n";
die $USAGE unless scalar(@ARGV) == 4;

my $FASTA = $ARGV[0];
my $NUM_SH = $ARGV[1];
my $WL = $ARGV[2];
my $MIN_SEQ = $ARGV[3];

my $SWD = getcwd();
my $NWD = $SWD.'/'.$FASTA.'_NSH_'.$NUM_SH.'_WL_'.$WL;
$NWD =~ s/\.fa//;

my $NUMSEQ;

my $COUNTS;
my $GENES;
my $SH_GENES;

mkdir($NWD) or die "\nmkdir failed: $!\n";
chdir($NWD) or die "\nchdir failed: $!\n";
copy($SWD.'/'.$FASTA,$NWD.'/'.$FASTA) or die "\ncopy failed: $!\n";

my $file_cou = "$FASTA\_$NUM_SH\_$WL\_counts.xls";
$file_cou =~ s/\.fa//;

my $file_occ = "$FASTA\_$NUM_SH\_$WL\_occurrences.txt";
$file_occ =~ s/\.fa//;


shuffle_fasta();

count_real();

count_shuffled();

analize_res();

write_log();

calculate_p_value();



=head2 shuffle_fasta

Title   : shuffle_fasta

Usage   : shuffle_fasta()

Function: create the shuffled versions of fasta

Returns : nothing, the file are created in the working dir
named shuffled_[n].fa where n is the number of the cicle

Args    :

Note    : 

=cut

sub shuffle_fasta {
	
  my $x;
	
  for($x=1;$x<=$NUM_SH;$x++) {
    
    my $seqio = Bio::SeqIO->new(-file => $FASTA,
		-format => 'fasta');
		
    my $file_shu = "$FASTA\_shuffled\_$x\.fa";
    $file_shu =~ s/\.fa//;
		
    my $seqio_sh = Bio::SeqIO->new(-file => ">$file_shu",
		-format => 'fasta');
		
    while(my $seq = $seqio->next_seq) {
			
      my @seq = split(//,$seq->seq);
			
      my $char = scalar(@seq)-1;
      my @order = (0..$char);
			
      my @new_order1 = _shuffle(\@order);
      my @new_order2 = _shuffle(\@new_order1);
      my @new_order = _shuffle(\@new_order2);
			
      my $string;
      foreach my $inx(@new_order) {
        $string .= $seq[$inx];
      }
			
      $string =~ s/N//g;
			
      my $seq_sh = Bio::Seq->new(-id => $seq->id,
			-seq => $string);
			
      $seqio_sh->write_seq($seq_sh);
			
      $NUMSEQ ++ if $x == 1;
    }
  }
}

=head2 _shuffle

Title   : _shuffle

Usage   : my @new_order = _shuffle(\@order)

Function: take in input an array of elements and returns a new array
in which the order of the elements has been shuffled

Returns : an array

Args    : -1 an array_ref of elements

Note    : called by shuffle()

=cut

sub _shuffle {
	return @_ if !@_ || ref $_ [0] eq 'ARRAY' && !@{$_ [0]};
	my $array = @_ == 1 && ref $_ [0] eq 'ARRAY' ? shift : [@_];
	for (my $i = @$array; -- $i;) {
		my $r = int rand ($i + 1);
		($array -> [$i], $array -> [$r]) = ($array -> [$r], $array -> [$i]);
	}
	wantarray ? @$array : $array;
}

=head2 count_real

Title   : count_real

Usage   : count_real()

Function: count the occurrences of all the motif in the given fasta

Returns : it populate the $COUNTS hashref

Args    :

Note    :

=cut

sub count_real {
	
  open(OCC,">$file_occ");
  print OCC "ref_id\tword\tposition\n";
	
  my $seqio = Bio::SeqIO->new(-file => $FASTA,
	-format => 'fasta');
	
  while(my $seq = $seqio->next_seq) {
		
    my $string = $seq->seq;
    my $id = $seq->id;
		
    my $start = 0;
    my $last = length($string) - $WL;
		
    while($start <= $last) {
			
      my $substring = substr($string,$start,$WL);
			
      $start ++;
			
      next if $substring =~ /N/;
			
      $COUNTS->{$substring}->{'R'} ++;
      $GENES->{$substring}->{$id} = 1;
			
      print OCC "$id\t$substring\t".($start-1)."\n";
    }
  }
}

=head2 count_shuffled

Title   : count_shuffled

Usage   : count_shuffled()

Function: 

Returns : 

Args    :

Note    :

=cut

sub count_shuffled {
	
  my $x;
	
  for($x=1;$x<=$NUM_SH;$x++) {
		
    my $file = "$FASTA\_shuffled_$x\.fa";
    $file =~ s/\.fa//;
		
    my $seqio = Bio::SeqIO->new(-file => $file,
		-format => 'fasta');
		
    while(my $seq = $seqio->next_seq) {
			
      my $string = $seq->seq;
      my $id = $seq->id;
			
      my $start = 0;
      my $last = length($string) - $WL;
			
      while($start <= $last) {
				
        my $substring = substr($string,$start,$WL);
				
        $start ++;
				
        next if $substring =~ /N/;
        next unless exists $COUNTS->{$substring};
				
        $COUNTS->{$substring}->{$x} ++;
        $SH_GENES->{$substring}->{$x}->{$id} = 1;
      }
    }
  }
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
  print COU "\tnum_occ\tnum_seq\tz_score_occ\tp.value_occ\tadj.p.value_occ\tz_score_num_seq\tp.value_num_seq\tadj.p.value_num_seq\n";
	
  foreach my $word(keys %$COUNTS) {
		
    my $stat_c = Statistics::Descriptive::Full->new();
    my $stat_g = Statistics::Descriptive::Full->new();
		
    my $r = $COUNTS->{$word}->{'R'};
    my $num_seq = scalar(keys %{$GENES->{$word}});
		
    next if $num_seq < $MIN_SEQ;
		
    my $x;
    my @s_c;
    my @s_g;
    my $z_score_c;
    my $z_score_g;
		
    for($x=1;$x<=$NUM_SH;$x++) {
      if(exists $COUNTS->{$word}->{$x}) {
				
        my $count_c = $COUNTS->{$word}->{$x};
        push(@s_c,$count_c);
				
        my $count_g = scalar(keys %{$SH_GENES->{$word}->{$x}});
        push(@s_g,$count_g);
      }
    }
    
    $stat_c->add_data(@s_c);
    my $mean_c = $stat_c->mean;
    my $sd_c = $stat_c->standard_deviation;
    my $diff_c = $r - $mean_c;
    $z_score_c = $diff_c / ($sd_c || 1);
		
    $stat_g->add_data(@s_g);
    my $mean_g = $stat_g->mean;
    my $sd_g = $stat_g->standard_deviation;
    my $diff_g = $num_seq - $mean_g;
    $z_score_g = $diff_g / ($sd_g || 1);
		
    print COU "$word\t$r\t$num_seq\t";
    print COU "$z_score_c\tNA\tNA\t"; #if $z_score_c;
    print COU "$z_score_g\tNA\tNA\n"; #if $z_score_g;
  }
  close(COU);
}

sub calculate_p_value {
	
  my $script = 'options(echo=FALSE)'."\n";
	
  $script .=
	'
	f<-read.table(file="'.$file_cou.'")
	
	for(i in 1:nrow(f)) {
		
		f$p.value_occ[i] <- 1-pnorm(abs(f$z_score_occ[i]))
		f$adj.p.value_occ[i] <- nrow(f)*(1-pnorm(abs(f$z_score_occ[i])))
		
		f$p.value_num_seq[i] <- 1-pnorm(abs(f$z_score_num_seq[i]))
		f$adj.p.value_num_seq[i] <- nrow(f)*(1-pnorm(abs(f$z_score_num_seq[i])))
	}
	
	write.csv(f,file="stat_'.$file_cou.'.csv")
	';
  
  my $filename = "SCRIPT_$ARGV[0]\_$ARGV[1]\_$ARGV[2]\_$ARGV[3]";
	
  open(SCRIPT,">$filename");
  print SCRIPT $script;
  
  system("R CMD BATCH --no-save $filename");
	#  system("rm $filename\*");
}

sub write_log {
  my $file_log = "$FASTA\_$NUM_SH\_$WL\_enumeration.log";
  $file_log =~ s/\.fa//;
  open(LOG,">$file_log");
  print LOG "\n$NUMSEQ sequences analyzed with parameters:\n\n";
  print LOG "FASTA FILE:\t$FASTA\n";
  print LOG "N SHUFFLING:\t$NUM_SH\n";
  print LOG "WORD LENGTH:\t$WL\n";
  print LOG "MINIMUM NUMBER OF SEQ:\t$MIN_SEQ\n\n\n";
}
