#!/usr/bin/perl
use strict; use warnings 'FATAL' => 'all';
use DataBrowser;
use Getopt::Std;
use vars qw($opt_e $opt_t $opt_o);
getopts('e:t:o:');

my $EXPAND = 2;
my $THRESHOLD = 1.0;
my $ODIR = "html";

die "
usage: motif_expander.pl [options] <mermade output> <fg> <bg>
options:
  -e  <int>    expansion in bp [$EXPAND]
  -t  <float>  threshold in bits per position [$THRESHOLD]
  -o  <name>   output directory name [$ODIR]
" unless @ARGV == 3;

my ($mfile, $ffile, $bfile) = @ARGV;

$EXPAND = $opt_e if $opt_e;
$THRESHOLD = $opt_t if $opt_t;
$ODIR = $opt_o if $opt_o;

###########
# globals #
###########
my $SLEN = `head -1 $ffile`;
chomp $SLEN;
$SLEN = length($SLEN);
my $MLEN = `head -2 $mfile`;
($MLEN) = $MLEN =~ /MERMADE\.(\d+)\./;
die "no motifs found\n" if (not defined $MLEN);
$THRESHOLD *= $MLEN;
my $motifamatic = "$ENV{MERMADE}/motifamatic";
my @alph = qw( A C G T );

###############
# read motifs #
###############
open(my $mfh, $mfile) or die;
my $header = <$mfh>;
my $id;
my $length;
my %motif;
while (<$mfh>) {
	chomp;
	if (/^(MERMADE.(\d+).+)/) {
		$id = $1;
		$length = $2;
	} elsif (/^[ACGT]/) {
		my ($nt, @f) = split;
		for (my $i = 0; $i < @f; $i++) {
			$motif{$id}[$i]{$nt} = $f[$i];
		}
	}
}
close $mfh;

#############
# main loop #
#############
my @motif;
foreach my $mid (sort keys %motif) {

	# create temporary motif file
	my $mfile = "/tmp/motif_expander.$$.motif";
	open(my $out, ">$mfile") or die;
	print $out ">$mid length=$length\n";
	foreach my $nt (@alph) {
		print $out $nt;
		for (my $i = 0; $i < $length; $i++) {
			printf $out "%6d", $motif{$mid}[$i]{$nt};
		}
		print $out "\n";
	}
	close $out;
	
	# run motifamatic on fg and bg
	my $ratio = (-s $ffile) / (-s $bfile);
	my $fseqs = get_matching_seqs($mfile, $ffile);
	my $bseqs = get_matching_seqs($mfile, $bfile);
	my $fcount = @$fseqs;
	my $bcount = int(@$bseqs * $ratio);
	
	# run seqlogo to create its graphical representation
	my $fasta = "/tmp/motif_expander.$$.fasta";
	open(my $fout, ">$fasta") or die;
	print $fout ">motif_expander motifs\n";
	for (my $i = 0; $i < @$fseqs; $i++) {
		print $fout ">seq$i\n", $fseqs->[$i], "\n";
	}

	close $fout;
	my $png = `seqlogo -f $fasta -F PNG -k1 -h4 -w10 -cbnY 2> /dev/null`;
	my $pdf = `seqlogo -f $fasta -F PDF -k1 -h8 -w20 -cbnY 2> /dev/null`;
	
	# note: running seqlogo twice because some PDFs generated by seqlogo do not convert to PNG without errors in ImageMagick:convert
	
	push @motif, {
		name => $mid,
		fcount => $fcount,
		bcount => $bcount,
		ratio  => ($fcount + 1) / ($bcount + 1),
		seed   => `cat $mfile`,
		png    => $png,
		pdf    => $pdf,
	};

	unlink $fasta;
	unlink $mfile;
	
}

##########
# output #
##########
system("mkdir $ODIR") unless -d $ODIR;

foreach my $motif (sort {$b->{ratio} <=> $a->{ratio}} @motif) {
	my $mid = $motif->{name};
	my $bcount = $motif->{bcount};
	my $fcount = $motif->{fcount};
	my $ratio = ($fcount +1) / ($bcount +1);
	
	open(my $out, ">$ODIR/$mid.mstats") or die;
	print $out "$mid\t$fcount\t$bcount\t$ratio\n";
	close $out;
	
	open(OUT, ">$ODIR/$mid.png") or die;
	print OUT $motif->{png};
	close OUT;
	
	open(OUT, ">$ODIR/$mid.pdf") or die;
	print OUT $motif->{pdf};
	close OUT;
}



###############
# subroutines #
###############

sub get_matching_seqs {
	my ($motif_file, $seq_file) = @_;
	
	my $padding = '-' x $EXPAND;

	# run motifamatc to get sequences over threshold
	my @seq;
	open(my $fh, "$motifamatic $motif_file $seq_file $SLEN $THRESHOLD |") or die;
	while (<$fh>) {
		my ($seq, $off, $s) = split;
		$seq = $padding . $seq . $padding;
		$seq = substr($seq, $off, $MLEN + $EXPAND * 2);
		push @seq, $seq;
	}
	
	return \@seq;
}

BEGIN {
	die "MERMADE environment not set" unless defined $ENV{MERMADE} and -d $ENV{MERMADE};
	my $seqlogo = `which seqlogo`;
	die "seqlogo not in path" if not $seqlogo;
}


