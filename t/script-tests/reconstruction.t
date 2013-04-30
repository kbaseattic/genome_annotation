use strict;
use warnings;

use Test::More;
use Data::Dumper;
use Getopt::Long;
use LWP::UserAgent;

use Bio::KBase::GenomeAnnotation::Client;

my $debug=0;
my $getoptResult=GetOptions(
        'debug' =>      \$debug,
);

my $infile = "MIT9313.genome.annotated.reconstructionTO";
my $outfile = "MIT9313.genome.annotated.reconstructionTO.out";
my $cmd;
my $ret;

#  Test 3 - Can the object do all of the methods that take
#           reconstruction typed objects as input.
my @reconstruction_methods = qw(
        reconstructionTO_to_roles
        reconstructionTO_to_subsystems
);

#  Test 6 - Download test data
unlink $infile if -e $infile;

my $ua = LWP::UserAgent->new();
my $res = $ua->get("http://www.kbase.us/docs/build/MIT9313.genome.annotated.reconstructionTO",
                   ":content_file" => "MIT9313.genome.annotated.reconstructionTO");

ok($res->is_success, "Downloaded test data");

# Create a genome typed object

note("Test the happy cases for reconstruction methods");

foreach my $method (@reconstruction_methods) {
	$cmd = "$method --input $infile --output $outfile ";
	eval { $ret = `$cmd`; };
	is($@, '',     "use valid data: $method returns without error");

	ok(-e $outfile ,"Does the outfile exist");
	ok(-s $outfile > 0,"Is the outfile non zero");
	unlink $outfile if (-e $outfile);
}

done_testing();
unlink "MIT9313.genome.annotated.reconstructionTO";


