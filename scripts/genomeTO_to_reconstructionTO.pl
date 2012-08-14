use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;

use Getopt::Long;
use strict;
use Data::Dumper;

my $input_file;
my $output_file;
my $url = "http://bio-data-1.mcs.anl.gov/services/genome_annotation";

my $rc = GetOptions('url=s'     => \$url,
		    'input=s' 	=> \$input_file,
		    'output=s'  => \$output_file,
		    );

my $usage = "genomeTO_to_reconstructionTO [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> reconstruction-file]";

@ARGV == 0 or die "Usage: $usage\n";

my $anno_server = Bio::KBase::GenomeAnnotation::Client->new($url);
my $in_fh;
if ($input_file)
{
    open($in_fh, "<", $input_file) or die "Cannot open $input_file: $!";
}
else
{
    $in_fh = \*STDIN;
}

my $out_fh;
if ($output_file)
{
    open($out_fh, ">", $output_file) or die "Cannot open $output_file: $!";
}
else
{
    $out_fh = \*STDOUT;
}
my $json = JSON::XS->new;

my $input_genome;
{
    local $/;
    undef $/;
    my $input_genome_txt = <$in_fh>;
    $input_genome = $json->decode($input_genome_txt);
}

my $output_reconstruction = $anno_server->genomeTO_to_reconstructionTO($input_genome);

$json->pretty(1);
print $out_fh $json->encode($output_reconstruction);
close($out_fh);
