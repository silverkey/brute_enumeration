#!/usr/bin/perl
use strict;
use warnings;

my $usage = "\n\n\tperl $0 [.profiles file] [version number]\n\n\n";
die $usage unless scalar(@ARGV) == 2;

my $profile = $ARGV[0];
my $version = $ARGV[1];
my $radix = $profile;
$radix =~ s/\.profiles//;

open(PROFILE,$profile);
my $cluster;
my $class;
my $group;
my $count;
my $tot;

my $organism = lc($radix);
open(TAB,">unigene_$organism\_$version\_profile.txt");
print TAB "ug_id\tsource\tcount\tsource_tot\tclass\n";

while(my $row = <PROFILE>) {
	if($row =~ /^\>/) {
		chomp($row);
		my @info = split(/\|/,$row);
		$cluster = $info[0];
		$cluster =~ s/\> /$radix\./;
		$class = $info[1];
		next;
	}
	if($row =~ /^(.+)\t(\d+) \/ (\d+)\s+$/) {
		$group = $1;
		$count = $2;
		$tot = $3;
		print TAB "$cluster\t$group\t$count\t$tot\t$class\n";
	}
	else {
		print TAB "DEAD FOR PROBLEM DURING  PARSING ON LINE: $row!\n";
		die "BAD LINE: $row";
	}
}
