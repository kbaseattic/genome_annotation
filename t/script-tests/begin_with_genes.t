use strict;
use warnings;

#	Command line tests for genome_annotation
#
#	Test the methods that begin with a genome-typed-object
#	that has gene features defined
#	Begin by downloading MIT9313.genomeTO
#	Call CDSs on the input file and call it MIT9313.genome.genesTO
#	Send output to MIT9313.genome.annotationTO
#
#	Three tests
#	1.	Did the test execute without errors
#	2.	Does the output file exist
#	3.	Is the output file non-empty

use Test::More;
use Data::Dumper;
use Getopt::Long;
use LWP::UserAgent;

my $debug=0;
my $getoptResult=GetOptions(
	'debug'	=>	\$debug,
	
);


# Commands that take a genomeTO with gene features as input
my @annotation_methods = qw(
	genomeTO_to_feature_data
        assign_functions_to_CDSs
	annotate_genome
	genomeTO_to_coding_regions
	genomeTO_to_html
	genomeTO_to_reconstructionTO
	merge_features
);

my $infile  = "MIT9313.genomeTO";
my $outgene = "MIT9313.genome.genesTO";
my $outfile = "MIT9313.genome.annotationTO";
my $ret     = '';
my $cmd     = '';

#  Download test data
unlink $infile  if -e $infile;
unlink $outgene if -e $outgene;
unlink $outfile if -e $outfile;


my $ua = LWP::UserAgent->new();
my $res = $ua->get("http://www.kbase.us/docs/build/MIT9313.genomeTO", 
		   ":content_file" => $infile);

print "Downloaded test data\n";
#
#	Was the data download successful
#
ok($res->is_success, "Downloaded test data");
ok(-e $infile ,"Does the infile exist");

note("Do gene calling for these methods");

# Call genes
	$cmd = "call_CDSs --input $infile --output $outgene ";
	eval { $ret = `$cmd`; };
	is($@, '',     "call_CDSs returns without error");

note("Test the happy cases for annotation methods");
#
#	Loop through the methods (command-line commands)
#	Test each one with the genomeTO as input
#

foreach my $method (@annotation_methods) {
	note ("Testing $method-------------------------");
	$cmd = "$method --input $outgene --output $outfile ";
	eval { $ret = `$cmd`; };
	is($@, '',     "use valid data: $method returns without error");

	ok(-e $outfile ,"Does the outfile exist");
	ok(-s $outfile > 0,"Is the outfile non zero");
	unlink $outfile if (-e $outfile);
}


done_testing();
#unlink $infile if -e $infile;


