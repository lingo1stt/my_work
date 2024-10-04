## this script is made to transform fasta file into ideal form for the main scripts
use Getopt::Std;
getopt("fo");

unless ($opt_f && $opt_o){
	print "-f fasta file\n";
	print "-o output file name";
	die();
}
#line variable represent the length
$line = 0;
# seq variable represent the name of the sequence
$seq = "";

open(my $out, ">", $opt_o) or die "Cannot open output file '$opt_o': $!\n";
open FILE, $opt_f;
while (<FILE>){
	chomp;
	$line++;
	if (/>/){
		if ($line>1){print $out "$seq\n";}
		print $out "$_\n";
		$seq = "";
		next;
	}
	$seq .= $_;
}
print $out "$seq\n";
close FILE;
