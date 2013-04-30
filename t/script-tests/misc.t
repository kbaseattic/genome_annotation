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



#  Test 3 - Can the object do all of the misc_methods
# misc_methods that take a genomeTO as input and return something else
my @misc_methods = qw(
        genomeTO_to_reconstructionTO
        genomeTO_to_feature_data
);

my $infile  = "MIT9313.genomeTO.annotated";
my $outfile = "MIT9313.genomeTO.annotated.out";
my $ret     = '';
my $cmd     = '';

#  Test 6 - Download test data
unlink $infile if -e $infile;

my $ua = LWP::UserAgent->new();
my $res = $ua->get("http://www.kbase.us/docs/build/MIT9313.genomeTO.annotated",
                   ":content_file" => $infile);

ok($res->is_success, "Downloaded test data");
ok(-e $infile ,"Does the infile exist");

note("Test the happy cases for misc misc_methods");

foreach my $method (@misc_methods) {
	$cmd = "$method --input $infile --output $outfile ";
	eval { $ret = `$cmd`; };
	is($@, '',     "use valid data: $method returns without error");

	ok(-e $outfile ,"Does the outfile exist");
	ok(-s $outfile > 0,"Is the outfile non zero");
#	unlink $outfile if (-e $outfile);
}

done_testing();
#unlink $infile;

