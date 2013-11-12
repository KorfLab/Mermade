#!/usr/bin/perl
use strict; use warnings;
use DBI qw(:sql_types);

die "usage: db_creator.pl <working directory> <experiment file> <output name>\n" unless @ARGV == 3;
my ($DIR, $EXP, $DB) = @ARGV;

##############################################################################
# create tables 
##############################################################################

my @table = ("
CREATE TABLE Kmer (
	bc    TEXT,
	k     INTEGER,
	kmer  TEXT,
	fg    INTEGER,
	bg    INTEGER,
	ratio REAL,
	score REAL,
	PRIMARY KEY (bc, kmer)
);
",
"
CREATE TABLE Motif (
	bc    TEXT,
	k     INTEGER,
	kmer  TEXT,
	pwm   TEXT,
	bits  REAL,
	fg    INTEGER,
	bg    INTEGER,
	ratio REAL,
	png   BLOB,
	pdf   BLOB,
	PRIMARY KEY (bc, kmer)
);
",
"
CREATE TABLE Experiment (
	bc    TEXT PRIMARY KEY,
	title TEXT,
	info  TEXT
);
");

my $dbh = DBI->connect("dbi:SQLite:$DB","","");


foreach my $table (@table) {
	my $sth = $dbh->prepare($table);
	$sth->execute();
}

##############################################################################
# load Experiment table 
##############################################################################

open(IN, $EXP) or die "ERROR: problem opening $EXP";
while (<IN>) {
	next unless /^[ACGT]{3,5}/;
	my ($bc, $title, $desc) = split(/,/, $_);
	my $sth = $dbh->prepare("INSERT INTO Experiment VALUES (\"$bc\", \"$title\", \"$desc\");");
	$sth->execute();
}
close IN;

##############################################################################
# load Kmer table - uses a temp file and .import for speed 
##############################################################################

my $tmpfile = "/tmp/db_creator.$$.txt";
open(OUT, ">$tmpfile") or die;
my @kfile = `ls $DIR/kstats/*`;
chomp @kfile;
foreach my $file (@kfile) {
	my ($bc, $k) = $file =~ /(\w+)\.(\d+)\.kstats/;
	open(IN, $file) or die;
	my $header = <IN>;
	while (<IN>) {
		print OUT "$bc\t$k\t$_";
	}
	close IN;
}
close OUT;
open(SQL, "| sqlite3 $DB") or die;
print SQL ".separator \"\t\"\n";
print SQL ".import $tmpfile Kmer\n";
close SQL;
unlink $tmpfile;

##############################################################################
# load Motif table - typical INSERT syntax
##############################################################################

my @file = `ls $DIR/motifs/*`;
chomp @file;
foreach my $file (@file) {
	my ($bc, $k) = $file =~ /(\w+)\.(\d+)\.mermade$/;
	
	open(IN, $file) or die;
	my $header = <IN>;
	while (my $title_line = <IN>) {
		my $aline = <IN>;
		my $cline = <IN>;
		my $gline = <IN>;
		my $tline = <IN>;
		my $iline = <IN>;
		my $sline = <IN>;
		
		my ($n) = $title_line =~ /(\d+)$/;
		my $pwm = $aline . $cline . $gline. $tline . $iline;
		my ($sum) = $sline =~ /Sum=(\S+)/;
		
		my ($kmer, $fg, $bg, $ratio, $score) = read_kstats($bc, $k, $n);		
		my ($pwm_fg, $pwm_bg, $pwm_ratio) = read_mstats($bc, $k, $n);
		my $png = `cat $DIR/mstats/$bc.$k/MERMADE.$k.$n.png`;
		my $pdf = `cat $DIR/mstats/$bc.$k/MERMADE.$k.$n.pdf`;
				
		my $sth = $dbh->prepare("INSERT INTO Motif VALUES (\"$bc\", $k, \"$kmer\", \"$pwm\", $sum, $pwm_fg, $pwm_bg, $pwm_ratio, ?, ?);");
		$sth->bind_param(1, $png, SQL_BLOB);
		$sth->bind_param(2, $pdf, SQL_BLOB);
		$sth->execute();
	}
	close IN;
}


##############################################################################
# subroutines
##############################################################################

sub read_kstats {
	my ($bc, $k, $n) = @_;
	
	open(my $fh, "$DIR/kstats/$bc.$k.kstats");
	my @col;
	for (my $i = 0; $i <= $n; $i++) {
		my $line = <$fh>;
		@col = split(/\s+/, $line);
	}
	close $fh;
	
	return @col;
}

sub read_mstats {
	my ($bc, $k, $n) = @_;
	
	my $line = `cat $DIR/mstats/$bc.$k/MERMADE.$k.$n.mstats`;
	chomp $line;
	my ($id, $fg, $bg, $ratio) = split(/\s+/, $line);
	
	return $fg, $bg, $ratio;
}

__END__
