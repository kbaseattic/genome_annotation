use strict;
use warnings;
#	Command line tests for genome_annotation
#
#	Test the methods that begin with a genome-typed-object
#	that has dna/contigs defined
#	Begin by downloading MIT9313.genomeTO
#	Send output to MIT9313.genome.genesTO
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
my $localServer=0;
my $getoptResult=GetOptions(
        'debug' =>      \$debug,
);

# Commands (including Gene calling) that take a genomeTO as input
my @genecall_methods = qw(
	genomeTO_strep_pneumo_repeats
	genomeTO_strep_suis_repeats
	call_selenoproteins
	call_pyrrolysoproteins
	call_RNAs
	call_CDSs_by_glimmer
	call_CDSs
);

# the call_CDSs_by_projection needs data that isn't yet computable by a service - it
# comes from the genomeTT_to_coding_regions script. See
# http://kbase.science.energy.gov/developer-zone/tutorials/iris/constructing-rast2-in-the-iris-environment/

my $infile  = "MIT9313.genomeTO";
my $outfile = "MIT9313.genome.genesTO";
my $ret     = '';
my $cmd     = '';


#  Download test data
unlink $infile  if -e $infile;
unlink $outfile if -e $outfile;

my $ua = LWP::UserAgent->new();
my $res = $ua->get("http://www.kbase.us/docs/build/MIT9313.genomeTO",
                   ":content_file" => $infile);
#
#	Was the data download successful
#
ok($res->is_success, "Downloaded test data");
ok(-e $infile ,"Does the infile exist");

note("Test the happy cases for gene calling methods");
#
#	Loop through the methods (command-line commands)
#	Test each one with the genomeTO as input
#
foreach my $method (@genecall_methods) {
	note ("Testing $method-------------------------");
	$cmd = "$method --input $infile --output $outfile ";
	eval { $ret = `$cmd`; };
	is($@, '',     "use valid data: $method returns without error");

	ok(-e $outfile ,"Does the outfile exist");
	ok(-s $outfile > 0,"Is the outfile non zero");
	unlink $outfile if (-e $outfile);
}


done_testing();
#unlink $infile if -e $infile;

