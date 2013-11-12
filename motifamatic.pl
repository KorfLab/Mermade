#!/usr/bin/perl
use strict; use warnings;
use DataBrowser;

die "usage: motifamatic.pl <motif> <seqs> <threshold>\n" unless @ARGV == 3;

my ($mfile, $sfile, $threshold) = @ARGV;

my $motif = read_motif($mfile);

my $len = @$motif;
my $exp = 0.25 ** $len;

open(my $fh, $sfile) or die;
while (<$fh>) {
	chomp;
	my $dna = $_;
	my $score = 1;
	my $max_s = -1e300;
	my $max_i = -1;
	
	for (my $i = 0; $i < length($dna) - @$motif + 1; $i++) {
		my $subseq = substr($dna, $i, $len);
				
		my $score = 1;
		for (my $j = 0; $j < $len; $j++) {
			my $nt = substr($subseq, $j, 1);
			my $s = $motif->[$j]{$nt};
			$s = 0.1 if not defined $s; # N penalty
			$score *= $s;
		}
				
		if ($score > $max_s) {
			$max_s = $score;
			$max_i = $i;
		}
	}
	
	$max_s = log($max_s / $exp) / log(2);
	
	printf "%s %d %f\n", $dna, $max_i, $max_s if $max_s > $threshold;
}

sub read_motif {
	my ($file) = @_;
	
	my @motif;
	open(my $fh, $file) or die;
	my $head = <$fh>;
	for (my $i = 0; $i < 4; $i++) {
		my $line = <$fh>;
		chomp $line;
		my ($nt, @val) = split(/\s+/, $line);
		for (my $j = 0; $j < @val; $j++) {
			$motif[$j]{$nt} = $val[$j] + 1; # pseudocount
		}
	}
	
	for (my $i = 0; $i < @motif; $i ++) {
		my $total = 0;
		foreach my $nt ('A', 'C', 'G', 'T') {
			$total += $motif[$i]{$nt};
		}
		foreach my $nt ('A', 'C', 'G', 'T') {
			$motif[$i]{$nt} /= $total;
		}
	}
	
	return \@motif;
}
