
use strict;
use gjoseqlib;
use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;

use Getopt::Long;

my $input_file;
my $output_file;
my $url = "http://bio-data-1.mcs.anl.gov/services/genome_annotation";

my $rc = GetOptions('url=s'     => \$url,
		    'input=s' 	=> \$input_file,
		    'output=s'  => \$output_file,
		    );

my $usage = "assign_functions_to_CDSs [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> genome-file]";

@ARGV == 0 or die "Usage: $usage\n";

my $kbase_server = Bio::KBase::GenomeAnnotation::Client->new($url);

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

my $output_genome = $kbase_server->assign_functions_to_CDSs($input_genome);

$json->pretty(1);
print $out_fh $json->encode($output_genome);
close($out_fh);
