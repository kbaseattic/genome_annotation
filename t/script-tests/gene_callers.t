use strict;
use warnings;

use Test::More;
use Data::Dumper;
use Getopt::Long;
use LWP::UserAgent;


my $debug=0;
my $localServer=0;
my $getoptResult=GetOptions(
        'debug' =>      \$debug,
);

# Gene calling methods that takes a genomeTO as input
my @genecall_methods = qw(
	call_selenoproteins
	call_pyrrolysoproteins
	call_RNAs
	call_CDSs_by_glimmer
	call_CDSs_by_projection
	call_CDSs
);

# the call_CDSs_by_projection needs data that isn't yet computable by a service - it
# comes from the genomeOT_to_coding_regions script. See
# http://kbase.science.energy.gov/developer-zone/tutorials/iris/constructing-rast2-in-the-iris-environment/

my $infile  = "MIT9313.genomeTO";
my $outfile = "MIT9313.genome.genesTO";
my $ret     = '';
my $cmd     = '';


#  Test 6 - Download test data
unlink "MIT9313.genomeTO" if -e "MIT9313.genomeTO";
my $ua = LWP::UserAgent->new();
my $res = $ua->get("http://www.kbase.us/docs/build/MIT9313.genomeTO",
                   ":content_file" => $infile);

ok($res->is_success, "Downloaded test data");
ok(-e $infile ,"Does the infile exist");

note("Test the happy cases for gene calling methods");

foreach my $method (@genecall_methods) {
	$cmd = "$method --input $infile --output $outfile ";
	eval { $ret = `$cmd`; };
	is($@, '',     "use valid data: $method returns without error");

	ok(-e $outfile ,"Does the outfile exist");
	ok(-s $outfile > 0,"Is the outfile non zero");
	unlink $outfile if (-e $outfile);
}


done_testing();
#unlink $infile if -e $infile;

