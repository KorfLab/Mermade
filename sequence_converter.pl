#!/usr/bin/perl
# sequence_converter.pl

use strict; use warnings;

die "usage: $0 <file>" unless @ARGV;
my $file = $ARGV[0];

if ($file =~/gz$/) {
	open(IN,"gunzip -c $file |") or die
} else {
	open(IN,$file) or die
}

# Scan first 5 lines and store into an array
my @line;
for (my $i = 0; $i < 4; $i++) {
	my $string = <IN>;
	chomp $string;
	push (@line, $string);
}
close(IN);

# Reopen file
if ($file =~/gz$/) {
	open(IN,"gunzip -c $file |") or die
} else {
	open(IN,$file) or die
}

# Assign file type
my $type;
if ($line[0] =~/^@/) {
	# FASTQ
	while (1) {
		my $id   = <IN>;
		last if not defined $id;
		my $read = <IN>;
		my $plus = <IN>;
		my $qual = <IN>;
		print $read;
		
	}
} elsif ($line[0] =~ /^\>/) {
	if ($line[1] =~ /^\>/) {
		# ELAND extended
		while (<IN>) {
			my ($id, $read, $qual) = split;
			print $read, "\n";
		}
	} else {
		# FASTA
		while (1) {
			my $id   = <IN>;
			last if not defined $id;
			my $read = <IN>;
			print $read;
		}
	}
} elsif ($line[0] =~ /^[ACGTN]+$/ and $line[1] =~ /^[ACGTN]+$/) {
	# raw
	while (<IN>) {print}
}
