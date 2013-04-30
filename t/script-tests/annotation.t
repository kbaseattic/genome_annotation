use strict;
use warnings;

use Test::More;
use Data::Dumper;
use Getopt::Long;
use LWP::UserAgent;

my $debug=0;
my $getoptResult=GetOptions(
	'debug'	=>	\$debug,
	
);


#     methods that take genomeTO and returns genomeTO
my @annotation_methods = qw(
        assign_functions_to_CDSs
	annotate_genome
	genomeTO_strep_pneumo_repeats
	genomeTO_strep_suis_repeats
	genomeTO_to_coding_regions
	genomeTO_to_feature_data
);

my $infile  = "MIT9313.genomeTO";
my $outgene = "MIT9313.genome.genesTO";
my $outfile = "MIT9313.genome.annotationTO";
my $ret     = '';
my $cmd     = '';


#  Download test data
unlink $infile if -e $infile;

my $ua = LWP::UserAgent->new();
my $res = $ua->get("http://www.kbase.us/docs/build/MIT9313.genomeTO", 
		   ":content_file" => "MIT9313.genomeTO");

print "Downloaded test data\n";

note("Test the happy cases for annotation methods");

# Call genes
	$cmd = "call_CDSs --input $infile --output $outgene ";
	eval { $ret = `$cmd`; };
	is($@, '',     "call_CDSs returns without error");

note("Test the happy cases for annotation methods");

foreach my $method (@annotation_methods) {
	$cmd = "$method --input $infile --output $outfile ";
	eval { $ret = `$cmd`; };
	is($@, '',     "use valid data: $method returns without error");

	ok(-e $outfile ,"Does the outfile exist");
	ok(-s $outfile > 0,"Is the outfile non zero");
	unlink $outfile if (-e $outfile);
}


done_testing();
#unlink $infile if -e $infile;


