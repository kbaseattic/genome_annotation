use strict;
use warnings;
#       Command line tests for genome_annotation
#
#       Test the misc commands

use Test::More;
use Data::Dumper;
use Getopt::Long;
use LWP::UserAgent;

my $debug=0;
my $getoptResult=GetOptions(
        'debug' =>      \$debug,
);

my $tmpfile = "testing.tmp.txt";
my $outfile = "testing.tmp";
my $ret     = '';
my $cmd     = '';

#  Delete before running the test
unlink $outfile if -e $outfile;

open (OUT,">$tmpfile") || die "Did not create tmp output file"; 
while (<DATA>)
{
	print OUT;
}
close OUT;

#
#-----------------------------------------------------------------
#
my $method = "cluster_objects";
note ("Testing $method---------------------------------------");

$cmd = "$method < $tmpfile > $outfile";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e $outfile ,"Does the outfile exist");
ok(-s $outfile > 0,"Is the outfile non zero");
unlink $outfile if (-e $outfile);
#
#-----------------------------------------------------------------
#
$method = "file_to_spreadsheet";
note ("Testing $method---------------------------------------");

$cmd = "$method -f $outfile < $tmpfile ";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e $outfile ,"Does the outfile exist");
ok(-s $outfile > 0,"Is the outfile non zero");
unlink $outfile if (-e $outfile);

#
#-----------------------------------------------------------------
#
$method = "cs_to_genome";
note ("Testing $method---------------------------------------");

$cmd = "$method 'kb|g.2323' > $outfile ";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e $outfile ,"Does the outfile exist");
ok(-s $outfile > 0,"Is the outfile non zero");
unlink $outfile if (-e $outfile);

#
#-----------------------------------------------------------------
#
$method = "sort_by_loc";
note ("Testing $method---------------------------------------");

$cmd = "$method < $tmpfile > tmp1 ";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e 'tmp1' ,"Does the outfile exist");
ok(-s 'tmp1' > 0,"Is the outfile non zero");

$cmd = "diff $tmpfile tmp1 > $outfile";
eval { $ret = `$cmd`; };
is($@, '',     "DIFF returns without error");

ok(-e $outfile ,"Does the diff outfile exist");
ok(-s $outfile > 0,"Is the diff outfile non zero");

$method = "sort_by_id";
note ("Testing $method---------------------------------------");

$cmd = "$method < $tmpfile > tmp2 ";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e 'tmp1' ,"Does the outfile exist");
ok(-s 'tmp1' > 0,"Is the outfile non zero");

$cmd = "diff $tmpfile tmp2 > $outfile";
eval { $ret = `$cmd`; };
is($@, '',     "DIFF returns without error");

ok(-e $outfile ,"Does the diff outfile exist");
ok(-s $outfile > 0,"Is the diff outfile non zero");

#
#	Output from the two sorts needs to be different
#
$cmd = "diff tmp1 tmp2 > $outfile";
eval { $ret = `$cmd`; };
is($@, '',     "DIFF returns without error");

ok(-e $outfile ,"Does the diff outfile exist");
ok(-s $outfile > 0,"Is the diff outfile non zero");
unlink $outfile if (-e $outfile);
unlink 'tmp1' if (-e 'tmp1');
unlink 'tmp2' if (-e 'tmp2');

#
#-----------------------------------------------------------------
#	NOT TESTED
#	regions_around
#	rel2tabs
#	tabs2rel
#-----------------------------------------------------------------
#

open (OUT,">$tmpfile") || die "Did not create tmp output file";        
print OUT ">test\ntacgtacagctagctagctagctagctagtcttgagctgactatcgatcgtagctagctagctagctgatcgtattatatattatgatcgtagctagctagctgactactgatcgtagctagctactagtg\n";
close OUT;

$method = "fasta_to_genome";
note ("Testing $method---------------------------------------");

$cmd = "$method 'Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)' Bacteria 11 < $tmpfile > $outfile ";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e $outfile ,"Does the outfile exist");
ok(-s $outfile > 0,"Is the outfile non zero");
unlink $outfile if (-e $outfile);

#
#-----------------------------------------------------------------
#


done_testing();
unlink $tmpfile;
unlink $outfile;

__DATA__
Col1	Col2
92	kb|g.2620
89	kb|g.416
89	kb|g.415
76	kb|g.873
76	kb|g.3149
92	kb|g.220
89	kb|g.46
89	kb|g.45
76	kb|g.83
76	kb|g.349
