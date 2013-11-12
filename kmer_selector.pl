#!/usr/bin/perl
# kmer_selector.pl
use strict;
use warnings;
use FileHandle;
use Getopt::Std;
use DataBrowser;

my %opt = (
	'c' => 0,
	'f' => 0,
	'p' => 10,
	'r' => 2,
	'd' => 0,
	'q' => 0,
);

getopts('c:f:p:r:h:dq', \%opt);

die "
usage: kmer_selector.pl [options] <foreground> <background>
options:
  -c <int>   count threshold      [$opt{c}]
  -f <float> frequency threshold  [$opt{f}]
  -p <float> poisson threshold    [$opt{p}]
  -r <float> ratio threshold      [$opt{r}]
  -q         quiet
  -d         duplicate and remove reverse complement words
  -help
" unless @ARGV == 2;

my ($fcount, $ftotal) = count_kmers($ARGV[0]);
my ($bcount, $btotal) = count_kmers($ARGV[1]);

print "kmer_selector output\n";

exit if $ftotal == 0;

# standardize keys
foreach my $kmer (keys %$fcount) {
	$bcount->{$kmer} = 0 if not defined $bcount->{$kmer};
}
foreach my $kmer (keys %$bcount) {
	$fcount->{$kmer} = 0 if not defined $fcount->{$kmer};
}

# normalize counts (reduce counts of larger file)
normalize($fcount, $ftotal, $bcount, $btotal);

# annotate the k-mers
my %enriched;
my %log;
foreach my $seq (keys %$fcount) {
	my $fore = $fcount->{$seq};
	if ($fore < $opt{c}) {
		$log{'foreground_counts_low'}++;
		next;
	}
	
	my $back = $bcount->{$seq};
	if ($back < $opt{c}) {
		$log{'background_counts_low'}++;
		next;
	}
	
	my $ratio = (1 + $fore ) / (1 + $back);
	if ($ratio < $opt{r}) {
		$log{'ratio_low'}++;
		next;
	}
	
	my $poisson = poisson($fore, $back);
	if ($poisson < $opt{p}) {
		$log{'poisson_low'}++;
		next;
	}
	
	# keep
	$enriched{$seq} = {
		fore => $fore,
		back => $back,
		ratio => $ratio,
		poisson => $poisson,
	};
	
	$log{'kept'}++;
}

foreach my $kmer (sort {$enriched{$b}{poisson} <=> $enriched{$a}{poisson}} keys %enriched) {
	printf "%s\t%d\t%d\t%.3f\t%.3f\n",
		$kmer,
		$enriched{$kmer}{fore},
		$enriched{$kmer}{back},
		$enriched{$kmer}{ratio},
		$enriched{$kmer}{poisson};
}

# logging to STDERR
print STDERR "kmer-selector.pl log\n" unless $opt{'q'};
foreach my $message (keys %log) {
	print STDERR "\t$message $log{$message}\n" unless $opt{'q'};
}

################################################################################

sub anti {
	my ($kmer) = @_;
	$kmer =~ tr/ACGT/TGCA/;
	return scalar reverse $kmer;
}

sub poisson {
	my ($n, $l) = @_;
	
	$n++; $l++; # pseudocount protection
	
	my $lnfac = (0.5 * log(2 * 3.14159265358979))
		+ ($n + 0.5) * log($n)
		- $n + 1 / (12 * $n)
		- 1 / (360 * ($n ** 3)); # stirling approx of ln factorial
	
	my $lnp = $n * log($l) - $l - $lnfac; # ln poisson

	if ($n < $l) {return $lnp}
	else         {return -$lnp}
	
	return -$lnp;
}

sub count_kmers {
	my ($file) = @_;
	
	my %count;
	my $total = 0;
	open(my $fh, $file) or die;
	my $header = <$fh>;
	while (<$fh>) {
		my ($kmer, $count) = split;
		$count{$kmer} += $count;
		if ($opt{d}) {
			$count{anti($kmer)} += $count;
		}
		$total += $count;
	}
	
	if ($opt{d}) {
		my @key = sort keys %count;
		my %keep;
		foreach my $kmer (@key) {
			next unless exists $count{$kmer};
			$keep{$kmer} = $count{$kmer};
			delete $count{$kmer};
			delete $count{anti($kmer)};
		}
		%count = %keep;
	}
		
	return \%count, $total;
}

sub normalize {
	my ($c1, $t1, $c2, $t2) = @_;
	
	$t2 = 1 if ($t2 == 0);
	
	my $r = $t1 / $t2;
	
	if ($r > 1) {
		foreach my $k (keys %$c1) {$c1->{$k} /= $r}
	} else {
		foreach my $k (keys %$c2) {$c2->{$k} *= $r}
	}
		
}

__END__