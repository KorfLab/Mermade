#!/usr/bin/perl
use strict; use warnings;
use DBI qw(:sql_types);

die "usage: $0 <database> <# of motifs> <output dir> <barcodes>" unless @ARGV >= 4;
my ($database, $motifs, $outdir, @barcode) = @ARGV;

if (-d $outdir) {die "will not overwrite $outdir, remove first\n"}
system("mkdir $outdir") == 0 or die "unable to create $outdir\n";
my $pdir = "PNGs";
my $pdf_dir = "PDFs";
my $html = "index.html";
my $html_title = "$outdir MERMADE report";
my $style = "style.css";

my $dbh = DBI->connect("dbi:SQLite:$database","","");
my $sth = $dbh->prepare("");
$sth->execute;

create_style($style);

chdir $outdir;

# Create HTML filehandle
open(HTML, ">$html");
print HTML "<html>
<head>
<title>$html_title</title>
<link rel=\"stylesheet\" type=\"text/css\" href=\"$style\" media=\"screen\" />
</head>
<body>
";

foreach my $bc (@barcode) {
	# Print table headers
	print HTML "<h1>$bc</h1>\n";
	print HTML "<table>
	<tr><td>BC
	<td>K
	<td>Seed
	<td>Foreground
	<td>Background
	<td>Ratio
	<td>PNG (click for PDF)
	<td>Position Weight Matrix
	</tr>\n";
	
	# Create image directories
	system("mkdir -p $pdir/$bc");
	system("mkdir -p $pdf_dir/$bc");
	
	# Prepare and execute database query
	$sth = $dbh->prepare("SELECT bc,k,kmer,pwm,bits,fg,bg,ratio,png,pdf FROM Motif WHERE Motif.bc = \"$bc\" ORDER BY Motif.ratio DESC limit $motifs");
	$sth->execute;
	
	while (my $ref = $sth->fetchrow_hashref) {
		# Retrieve and write PDF
		open(PNG, ">$pdir/$bc/$ref->{kmer}.png");
		print PNG $ref->{png};
		close(PNG);
		
		# Retrieve and write PDF
		open(PDF, ">$pdf_dir/$bc/$ref->{kmer}.pdf");
		print PDF $ref->{pdf};
		close(PDF);
		
		# Populate table
		print HTML "<tr>\n";
		print HTML "<td>$bc\n";
		print HTML "<td>$ref->{k}\n";
		print HTML "<td>$ref->{kmer}\n";
		print HTML "<td>$ref->{fg}
		<td>$ref->{bg}\n";
		printf HTML "<td>%.3f\n",$ref->{ratio};
		print HTML "<td><a href=\"$pdf_dir/$bc/$ref->{kmer}.pdf\"><img src=\"$pdir/$bc/$ref->{kmer}.png\">\n";
		print HTML "<td><pre>$ref->{pwm}</pre>\n";
		print HTML "</tr>\n";
	}
	print HTML "</table>\n";

}

print HTML "</body>
</html>";
close(HTML);

# Creates CSS3 stylesheet (For prettier tables)
sub create_style {
	chdir $outdir;
	my $fh = shift;
	open(STYLE, ">$fh");
	print STYLE '* {
	margin: 0px;
	padding: 0px;
}

h1 {
	font-family: Arial;
	font-size:18pt;
	margin: 8px 8px 0px 8px;
}

table {
	border-collapse: collapse;
	margin: 8px;
	padding: 0px;
	font-family: Arial;
	font-size: 10pt;
	text-align: center;
}

tr:nth-child(2n+1){
	background-color: #EEEEEE;
}

td {
	padding: 5px 10px 5px 10px;
	border-style: solid;
	border-width: 1px;

}';
	close(STYLE);
	return 1;
}

__END__