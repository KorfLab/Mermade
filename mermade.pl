#!/usr/bin/perl
use strict; use warnings 'FATAL' => 'all';
use DataBrowser;
use Getopt::Std;
use vars qw($opt_n $opt_m $opt_b $opt_p);
getopts('m:e:p:');

my $MOTIFS = 5;
my $SLOPE  = 0.1;
my $OFFSET = 1;

my @Alph = qw(A C G T);

die "
usage: mermade.pl [options] <kstats file>
options:
  -n <int>   number of motifs to create [$MOTIFS]
  -m <float> edit distance slope [$SLOPE]
  -b <int>   edit distance offset [$OFFSET]
  -p <file>  png file base name (note: requires ImageMagick 'convert')
" unless @ARGV;

$MOTIFS = $opt_n if $opt_n;
$SLOPE  = $opt_m if $opt_m;
$OFFSET = $opt_b if $opt_b;
my $PNG = $opt_p;

print "mermade output\n";

# read k-mer stats file
my %kstats;
my $header = <>;
while (<>) {
	my ($kmer, $fg, $bg, $ratio, $poisson) = split;
	$kstats{$kmer} = {
		fg => $fg,
		bg => $bg,
		ratio => $ratio,
		poisson => $poisson,
	};
}

my $kmers = keys %kstats;
if ($kmers < 1) {exit}

# determine allowable edit distance of neighbors
my ($KMER) = keys %kstats;
my $K = length $KMER;
my $EDIT = int($K * $SLOPE) + $OFFSET;

# build motifs
my @motif;
foreach my $k1 (sort {$kstats{$b}{poisson} <=> $kstats{$a}{poisson}} keys %kstats) {

	next if defined $kstats{$k1}{friended};
	$kstats{$k1}{friended} = 1;

	# gather friends
	my @friend;
	foreach my $k2 (keys %kstats) {
		if (distance($k1, $k2) <= $EDIT) {
			push @friend, $k2 ;
			$kstats{$k2}{friended} = 1;
		}
	}

	# create motif from counts of kmer
	my @count;
	foreach my $kmer (@friend) {
		for (my $i = 0; $i < length($kmer); $i++) {
			$count[$i]{substr($kmer, $i, 1)} += $kstats{$kmer}{fg}; # change?
		}
	}

	for (my $i = 0; $i < length($k1); $i++) {
		foreach my $nt (@Alph) {
			$count[$i]{$nt} = 0 if not defined $count[$i]{$nt};
		}
	}
	
	push @motif, \@count;	
	last if @motif == $MOTIFS;
}

# report motifs
my $mcount = 0;
foreach my $motif (@motif) {
	$mcount++;
	my $mid = "MERMADE.$K.$mcount";
	
	# text output to STDOUT
	print $mid, "\n";
	display($motif);
}

exit(0);

###############################################################################
# subroutines
###############################################################################

# Not edit distance - does not account for indels
sub distance {
	my ($k1, $k2) = @_;
		
	my $mismatch = 0;
	for (my $i = 0; $i < length($k1); $i++) {
		$mismatch++ if substr($k1, $i, 1) ne substr($k2, $i, 1);
	}
	
	return $mismatch;
}

sub display {
	my ($motif) = @_;
	
	foreach my $nt (@Alph) {
		print $nt;
		for (my $i = 0; $i < @$motif; $i++) {
			no warnings;
			printf " %7d", $motif->[$i]{$nt};
			use warnings;
		}
		print "\n";
	}
		
	print "i";
	my $sum = 0;
	for (my $i = 0; $i < @$motif; $i++) {
		my $info = info($motif->[$i]);
		printf "   %.3f", $info;
		$sum += $info;
	}
	printf "\nSum=%.3f bits\n", $sum;
}

sub info {
	my ($c) = @_;
		
	my %p;
	my $total;
	foreach my $k (keys %$c) {$total += $c->{$k}}
	foreach my $k (keys %$c) {
		$p{$k} = $c->{$k} / $total if $c->{$k} > 0;
	}
	
	my $H = 0;
	foreach my $k (keys %p) {
		$H += $p{$k} * log($p{$k});
	}
	$H /= log(2);
	
	return 2 + $H;
}

