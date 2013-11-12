#!/usr/bin/perl
use strict; use warnings 'FATAL' => 'all';
use sigtrap;
use Getopt::Std;
use vars qw($opt_o $opt_a $opt_b $opt_d $opt_k $opt_K $opt_l $opt_c $opt_p $opt_r $opt_n $opt_e $opt_m $opt_q $opt_v);
getopts('o:d:b:k:K:c:p:r:n:l:m:v:qa');

BEGIN {
	die "MERMADE environment variable must point to a directory\n"
		unless exists $ENV{MERMADE} and -d $ENV{MERMADE};
}

#####################
# usage and options #
#####################

my $OUTFILE      = "bind_n_seq.db";
my $BLEN         = 3;
my $VERIFICATION = 'AGATCGGAAG';
my $MIN_KMER     = 6;
my $MAX_KMER     = 12;
my $MIN_COUNT    = 50;
my $MIN_POISSON  = 2;
my $MIN_RATIO    = 2;
my $KMER_REPORT  = 50;
my $EDIT         = 0.2;
my $MOTIF_REPORT = 5;
my $WDIR         = "wdir";
my $BACKGROUND   = "$ENV{MERMADE}/background.txt";

die "
usage: run_mermade.pl [options] <sequence file> <experiment table>
option:
  -o <path>   output file name [$OUTFILE]
  -d <dir>    working directory [$WDIR]
  -b <file>   background file [$BACKGROUND]
  -l <int>    barcode length [$BLEN]
  -v <string> verification string [$VERIFICATION]
  -k <int>    minimum k-mer [$MIN_KMER]
  -K <int>    maximum k-mer [$MAX_KMER]
  -c <int>    count threshold for k-mers [$MIN_COUNT]
  -p <float>  poission threshold for k-mers [$MIN_POISSON]
  -r <float>  ratio threshold for k-mers [$MIN_RATIO]
  -n <int>    number of k-mers to report [$KMER_REPORT]
  -m <int>    number of motifs to report [$MOTIF_REPORT]
  -q          quiet mode (do not report progress to STDERR)
  -a		  turn on reverse complement remover
" unless @ARGV == 2;
my ($INFILE, $TABLE) = @ARGV;

$OUTFILE      = $opt_o if $opt_o;
$WDIR         = $opt_d if $opt_d;
$BACKGROUND   = $opt_b if $opt_b;
$BLEN         = $opt_l if $opt_l;
$MIN_KMER     = $opt_k if $opt_k;
$MAX_KMER     = $opt_K if $opt_K;
$MIN_COUNT    = $opt_c if $opt_c;
$MIN_POISSON  = $opt_p if $opt_p;
$MIN_RATIO    = $opt_r if $opt_r;
$KMER_REPORT  = $opt_n if $opt_n;
$MOTIF_REPORT = $opt_m if $opt_m;
$VERIFICATION = $opt_v if $opt_v;
my $QUIET     = $opt_q;
my $ANTI	  = $opt_a;

########################################
# executables and other required files #
########################################

my $converter = "$ENV{MERMADE}/sequence_converter.pl";
my $debarcode = "$ENV{MERMADE}/de-barcode.pl";
my $kcounter  = "$ENV{MERMADE}/kmer_counter.pl";
my $kselector = "$ENV{MERMADE}/kmer_selector.pl";
my $mermade   = "$ENV{MERMADE}/mermade.pl";
my $expander  = "$ENV{MERMADE}/motif_expander.pl";
my $creator   = "$ENV{MERMADE}/db_creator.pl";
foreach my $file ($BACKGROUND, $converter, $debarcode, $kselector) {
	die "required file ($file) not found\n" unless -e $file;
}

system("mkdir $WDIR") unless -d $WDIR;
system("mkdir $WDIR/barcodes") unless -d "$WDIR/barcodes";
system("mkdir $WDIR/counts") unless -d "$WDIR/counts";
system("mkdir $WDIR/kstats") unless -d "$WDIR/kstats";
system("mkdir $WDIR/motifs") unless -d "$WDIR/motifs";
system("mkdir $WDIR/mstats") unless -d "$WDIR/mstats";

#############################
# split input into 64 files #
#############################

print STDERR "converting and splitting files... " unless $QUIET;
my @files = `ls $WDIR/barcodes`;
if (@files < 4) {
	run("$converter $INFILE | $debarcode -d $WDIR/barcodes -f $WDIR/log.txt -b $BLEN -v $VERIFICATION");
	print STDERR "done\n" unless $QUIET;
} else {
	print STDERR "previously done\n" unless $QUIET;
}

#####################################
# check file sizes, skip tiny files #
#####################################

print STDERR "checking file sizes... " unless $QUIET;
my @barcode;
my @file = `ls $WDIR/barcodes/*`;
my %size;
chomp @file;
my $total = 0;
foreach my $file (@file) {
	if (not -e $file) {die "can't get size of $file"}
	$size{$file} = -s $file;
	$total += $size{$file};
}
foreach my $file (@file) {
	my $frac = $size{$file} / $total;
	if ($frac > 0.001) { # should be a command line option
		push @barcode, $file;
	} else {
		print STDERR "skipped $file\n" unless $QUIET;
	}
}
print STDERR "done\n" unless $QUIET;

############################
# create background counts #
############################

print STDERR "creating background counts" unless $QUIET;
for (my $k = $MIN_KMER; $k <= $MAX_KMER; $k++) {
	print STDERR "." unless $QUIET;
	unless (-s "$WDIR/counts/bg.$k.counts") {
		run("$kcounter -k $k $BACKGROUND > $WDIR/counts/bg.$k.counts");
	}
}
print STDERR " done\n" unless $QUIET;

#############
# main loop #
#############

foreach my $file (@barcode) {
	my ($bc) = $file =~ /(\w+)\.txt$/;
	
	print STDERR "processing $bc:" unless $QUIET;
	
	for (my $k = $MIN_KMER; $k <= $MAX_KMER; $k++) {
		print STDERR " $k" unless $QUIET;
		
		# count k-mers
		print STDERR "c" unless $QUIET;
		my $fg = "$WDIR/counts/$bc.$k.counts";
		run("$kcounter -k $k $file > $fg") unless -s $fg;
	
		# select seed k-mers
		print STDERR "s" unless $QUIET;
		my $bg = "$WDIR/counts/bg.$k.counts";
		my $ks = "$WDIR/kstats/$bc.$k.kstats";
		run("$kselector -dq $fg $bg > $ks") unless -s $ks;
		next unless -s $ks > 100;
		
		# build seed motifs
		print STDERR "m" unless $QUIET;
		my $mm = "$WDIR/motifs/$bc.$k.mermade";
		run("$mermade -p $WDIR/motifs/$bc $ks > $mm") unless -s $mm;
		next unless -s $mm > 100;
		
		# expand motifs
		print STDERR "e" unless $QUIET;
		my $ex = "$WDIR/mstats/$bc.$k";
		my $ff = "$WDIR/barcodes/$bc.txt";
		my $bf = $BACKGROUND;
		run("$expander -o $ex $mm $ff $bf") unless -d $ex;
		
	}
	print STDERR " done\n" unless $QUIET;
	
}

###################
# create database #
###################

print STDERR "creating database..." unless $QUIET;
run("$creator $WDIR $TABLE $OUTFILE");
print STDERR "done\n" unless $QUIET;


##############################################################################
# subroutines
##############################################################################

sub run {
	my ($cmd) = @_;
#	print STDERR "\n$cmd\n";
	system($cmd) == 0 or die "$cmd failed\n";
}


__END__
