use strict;
use warnings;
#	Command line tests for genome_annotation
#
#	Test the methods that begin with a reconstructionTO
#	Begin by downloading MIT9313.genome.annotated.reconstructionTO
#	Send output to MIT9313.genome.annotated.reconstructionTO.out 
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
        'debug' =>      \$debug,
);

my $infile = "MIT9313.genome.annotated.reconstructionTO";
my $outfile = "MIT9313.genome.annotated.reconstructionTO.out";
my $cmd;
my $ret;

# Commands (including Gene calling) that take a
#           reconstruction typed objects as input.
my @reconstruction_methods = qw(
        reconstructionTO_to_roles
        reconstructionTO_to_subsystems
);

#  Download test data
unlink $infile  if -e $infile;
unlink $outfile if -e $outfile;

my $ua = LWP::UserAgent->new();
my $res = $ua->get("http://www.kbase.us/docs/build/MIT9313.genome.annotated.reconstructionTO",
                   ":content_file" => $infile);

#
#	Was the data download successful
#
ok($res->is_success, "Downloaded test data");
ok(-e $infile ,"Does the infile exist");

note("Test the happy cases for reconstruction methods");
#
#	Loop through the methods (command-line commands)
#	Test each one with the reconstructionTO as input
#

foreach my $method (@reconstruction_methods) {
	$cmd = "$method --input $infile --output $outfile ";
	eval { $ret = `$cmd`; };
	is($@, '',     "use valid data: $method returns without error");

	ok(-e $outfile ,"Does the outfile exist");
	ok(-s $outfile > 0,"Is the outfile non zero");
	unlink $outfile if (-e $outfile);
}

done_testing();
unlink $infile;


