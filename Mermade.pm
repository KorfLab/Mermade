package Mermade;
use strict;
use warnings 'FATAL' => 'all';

my @Alph = qw(A C G T);


sub distance { # Not edit distance - does not account for indels
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
