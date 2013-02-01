
=head1 genomeTO_to_strep_suis_repeats


This routine takes as input a file containing a JSON-encoded
genomeTO.  It invokes a tool to locate a set of control sites (often called "repeats")
in Streptococcus suisnia genomes.  This tool locates the sites and adds the features
representing them to the genomeTO.

Example:

    genomeTO_to_strep_suis_repeats < input > output

=cut

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

my $usage = "genomeTO_to_strep_suis_repeats [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> extended-genome-file]";

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

my $genomeTO;
{
    local $/;
    undef $/;
    my $input_genome_txt = <$in_fh>;
    $genomeTO = $json->decode($input_genome_txt);
}
use JSON::XS;

#...$resultTO appears to differ from $genomeTO, so apparently =NOT= "In Place"!
my $resultTO = $kbase_server->get_strep_pneumo_repeats($genomeTO);  #...Modified "in place"

$json->pretty(1);
print $out_fh $json->encode($resultTO);
close($out_fh);
