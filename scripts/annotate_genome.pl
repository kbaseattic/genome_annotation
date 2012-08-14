
=head1 annotate_genome


This routine takes as input a file containing a JSON-encoded
genomeTO.  It attempts a basic annotation task.  This involves
calling protein-encoding genes using Glimmer, assign functions using a Kmer-based strategy,
calling tRNAs using a tool built by Sean Eddy, calling rRNAs by a tool written by Niels Larsen,
as well as a set of enhancements that were supplied by Gary Olsen.  In other words, it attempts
to use a set of common tools to rapidly produce a fairly high-quality initial annotation.

Example:

    annotate_genome < input.genomeTO > output.genomeTO

=cut


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

my $usage = "annotate_genome [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> genome-file]";

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

my $output_genome = $anno_server->annotate_genome($input_genome);

$json->pretty(1);
print $out_fh $json->encode($output_genome);
close($out_fh);
