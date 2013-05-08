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

$cmd = "$method -c 2 < $tmpfile > tmp1 ";
#print "DEBUG: $cmd\n";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e 'tmp1' ,"Does the outfile exist");
ok(-s 'tmp1' > 0,"Is the outfile non zero");

$cmd = "diff $tmpfile tmp1 > $outfile";
#print "DEBUG: $cmd\n";
eval { $ret = `$cmd`; };
is($@, '',     "DIFF returns without error");

ok(-e $outfile ,"Does the diff outfile exist");
ok(-s $outfile > 0,"Is the diff outfile non zero");

$method = "sort_by_id";
note ("Testing $method---------------------------------------");

$cmd = "$method -c 1 < $tmpfile > tmp2 ";
#print "DEBUG: $cmd\n";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e 'tmp1' ,"Does the outfile exist");
ok(-s 'tmp1' > 0,"Is the outfile non zero");

$cmd = "diff $tmpfile tmp2 > $outfile";
#print "DEBUG: $cmd\n";
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
#
$method = "tabs2rel";
note ("Testing $method---------------------------------------");

$cmd = "$method < $tmpfile > $outfile ";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e $outfile ,"Does the outfile exist");
ok(-s $outfile > 0,"Is the outfile non zero");

#
#-----------------------------------------------------------------
#
$method = "rel2tabs";
note ("Testing $method---------------------------------------");

$cmd = "$method < $outfile > $tmpfile ";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e $tmpfile ,"Does the outfile exist");
ok(-s $tmpfile > 0,"Is the outfile non zero");

#
#-----------------------------------------------------------------
#	NOT TESTED
#	regions_around
#-----------------------------------------------------------------
#

open (OUT,">$tmpfile") || die "Did not create tmp output file";        
print OUT ">test\ntacgtacagctagctagctagctagctagtcttgagctgactatcgatcgtagctagctagctagctgatcgtattatatattatgatcgtagctagctagctgactactgatcgtagctagctactagtg\n";
close OUT;

$method = "fasta_to_genome";
note ("Testing $method---------------------------------------");

$cmd = "$method 'Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)' Bacteria 11 < $tmpfile > $outfile ";
print "DEBUG: $cmd\n";
eval { $ret = `$cmd`; };
is($@, '',     "use valid data: $method returns without error");

ok(-e $outfile ,"Does the outfile exist");
ok(-s $outfile > 0,"Is the outfile non zero");
#unlink $outfile if (-e $outfile);

#
#-----------------------------------------------------------------
#


done_testing();
#unlink $tmpfile;
unlink $outfile;

__DATA__
kb|g.22250.CDS.2733286	kb|g.19976.c.0_24099+588	62
kb|g.22250.CDS.2733262	kb|g.19976.c.0_20799+1359	60
kb|g.22250.CDS.2733278	kb|g.19976.c.0_21959+615	91
kb|g.22250.CDS.2733279	kb|g.19976.c.0_23968+249	60
kb|g.22250.CDS.2620	kb|g.19976.c.0_25099+588	92
kb|g.22250.CDS.416	kb|g.19976.c.0_26099+588	89
kb|g.22250.CDS.415	kb|g.19976.c.0_27099+588	89
kb|g.22250.CDS.873	kb|g.19976.c.0_28099+588	76
kb|g.22250.CDS.3149	kb|g.19976.c.0_29099+588	76
kb|g.22250.CDS.220	kb|g.19976.c.0_30099+588	92
kb|g.22250.CDS.46	kb|g.19976.c.0_31099+588	89
kb|g.22250.CDS.45	kb|g.19976.c.0_32099+588	89
kb|g.22250.CDS.83	kb|g.19976.c.0_33099+588	76
kb|g.22250.CDS.349	kb|g.19976.c.0_33099+588	76
