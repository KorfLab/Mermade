#!/usr/bin/perl
use strict; use warnings;
use FileHandle;
use Getopt::Std;
my %opt = (
	'd' => 'barcodes',
	'b' => 3,
	'm' => 2,
	'n' => 0,
	's' => 21,
	'v' => 'AGATCGGAAG',
	'f' => 'log.txt',
);
getopts('d:b:m:n:s:v:f:h', \%opt);

my $usage = "
usage: de-barcode.pl [options] STDIN
options:
  -d <string> output directory name [$opt{d}]
  -b <int>    barcode length [$opt{b}]
  -m <int>    pre-barcode length [$opt{m}]
  -n <int>    post-barcode length [$opt{n}]
  -s <int>    sequence length [$opt{s}]
  -v <string> verification string [$opt{v}]
  -f <string> name of the log file [$opt{f}]
notes:
  * input is from STDIN (eg. gunzip -c sequences.gz | de-barcode.pl)
  * this program expects RAW format with 1 line for sequence
  * barcodes with N are skipped
";
die $usage if $opt{h} or @ARGV;

my ($DIR, $BLEN, $PRE, $POST, $SLEN, $FLANK, $LOG) = ($opt{d}, $opt{b}, $opt{'m'}, $opt{n}, $opt{'s'}, $opt{v}, $opt{f});

my %fh;
my %bc_count;
my %flank_count;
my %seq_count;
my $flank_skipped = 0;
my $total_sequences = 0;

while (my $seq = <>) {
	chomp $seq;
	my $bc = substr($seq, $PRE, $BLEN);
	my $sq = substr($seq, $PRE + $BLEN + $POST, $SLEN);
	my $flank = substr($seq, $PRE + $BLEN + $POST + $SLEN, length($opt{v}));
	next if $bc =~ /N/;
	$total_sequences++;
	$flank_count{$flank}++;
	if ($flank ne $FLANK) {
		$flank_skipped++;
		next;
	};
	
	if (not defined $fh{$bc}) {
		$fh{$bc} = new FileHandle;
		$fh{$bc}->open(">$DIR/$bc.txt") or die "Could not open $DIR/$bc.txt";
	}
	$bc_count{$bc}++;
	$seq_count{$sq}++;
	$fh{$bc}->print($sq, "\n") if $seq_count{$sq} == 1;

}

# logfile

open(my $logfile, ">$LOG") or die;
print $logfile "\nSequences per Barcode\n";
foreach my $bc (sort keys %bc_count) {
	print $logfile $bc, "\t", $bc_count{$bc}, "\n";
}

if ($total_sequences) {
	printf $logfile "\nReads thrown out due to mismatched flank %.3f\n", $flank_skipped/$total_sequences;
} else {
	printf $logfile "\nNo reads kept :(\n";
}

print $logfile "\nTop 10 Flank Counts (over all barcodes)\n";
my $count = 0;
foreach my $flank (sort {$flank_count{$b} <=> $flank_count{$a}} keys %flank_count) {
	print $logfile "$flank $flank_count{$flank}\n";
	last if $count++ == 10;
}

print $logfile "\nTop 10 Redundant Sequences (over all barcodes)\n";
$count = 0;
foreach my $sq (sort {$seq_count{$b} <=> $seq_count{$a}} keys %seq_count) {
	if ($seq_count{$sq} > 1) {
		print $logfile "$sq $seq_count{$sq}\n";
	}
	last if $count++ == 10;
}