#
# Given a set of json dumps of the spgene_ref core, make a dbfile lookup table
#

use Data::Dumper;
use File::Slurp;
use JSON::XS;
use DB_File;
use strict;
use File::Temp;

if (@ARGV != 1)
{
    die "Usage: $0 data-dir\n";
}

my $dir = shift;

my $idx_file = "$dir/spgenes.idx";

my $t1 = File::Temp->new;
my $t2 = File::Temp->new;
my $t3 = File::Temp->new;
my $t4 = File::Temp->new;
$t1->close();
$t2->close();
$t3->close();
$t4->close();
my $hdrs = "-H 'Accept: application/solr+json'  -H 'Content-type: application/solrquery+x-www-form-urlencoded'";

system("curl -o $t1 $hdrs 'https://www.patricbrc.org/api/sp_gene_ref?q=*%3A*&rows=25000&wt=json&indent=true'");
system("curl -o $t2 $hdrs 'https://www.patricbrc.org/api/sp_gene_ref?q=*%3A*&rows=25000&start=25000&wt=json&indent=true'");
system("curl -o $t3 $hdrs 'https://www.patricbrc.org/api/sp_gene_ref?q=*%3A*&rows=25000&start=50000&wt=json&indent=true'");
system("curl -o $t4 $hdrs 'https://www.patricbrc.org/api/sp_gene_ref?q=*%3A*&rows=25000&start=75000&wt=json&indent=true'");

my @in = ($t1, $t2, $t3, $t4);

my %idx;
tie %idx, 'DB_File', $idx_file, O_RDWR | O_CREAT, 0644, $DB_HASH or die "Cannot tie $idx_file\n";

for my $in (@in)
{
    my $txt = read_file("$in");
    $txt or die "Cannot read $in: $!";
    my $j = eval { decode_json($txt); };

    if (!$j || $@)
    {
	die "Cannot decode json from $in: $@\n";
    }
    my $docs = $j->{response}->{docs};
    ref($docs) eq 'ARRAY' or die "Invalid json in $in\n"; 
    for my $doc (@$docs)
    {
	$doc->{classification} = join(",", @{$doc->{classification}}) if ref($doc->{classification});
	$doc->{classification} = lc($doc->{classification});
	
	$idx{$doc->{source}, $doc->{source_id}} = encode_json($doc);
    }
}
