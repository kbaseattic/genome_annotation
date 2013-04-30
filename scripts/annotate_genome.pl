use strict;

=head1 NAME

annotate_genome

=head1 SYNOPSIS

annotate_genome [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> genome-file]

=head1 DESCRIPTION

This routine takes as input a file containing a JSON-encoded genomeTO.
It attempts to produce a basic automated genome annotation.
This involves calling protein-encoding genes using Glimmer,
assign functions using a Kmer-based strategy,
calling tRNAs using a tool built by Sean Eddy,
calling rRNAs by a tool written by Niels Larsen,
as well as a set of enhancements that were supplied by Gary Olsen.
In other words, it attempts to use a set of common tools to rapidly produce
a fairly high-quality initial annotation.

Example:

    annotate_genome < input.genomeTO > output.genomeTO

where 'input.genomeTO' is an initial KBase object of 'genome' type,
and 'output.genomeTO' is 'input.genomeTO' augmented by the feature calls
and feature annotations.

=head1 COMMAND-LINE OPTIONS

Usage: annotate_genome --url service-url < genome-file > genome-file
Usage: annotate_genome --url service-url --input genome-file --output genome-file

    --url    --- Optional URL for alternate KBase server (D: http://bio-data-1.mcs.anl.gov/services/genome_annotation)

    --input  --- Option to read genome-typed-object from input file instead of from STDIN

    --output --- Option to write enhanced genome-typed-object to output file instead of STDOUT

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut


use gjoseqlib;
use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;

use Getopt::Long;

my $help;
my $input_file;
my $output_file;
# my $url = "http://bio-data-1.mcs.anl.gov/services/genome_annotation";
my $url = "https://kbase.us/services/genome_annotation";

my $rc = GetOptions('help'      => \$help,
		    'url=s'     => \$url,
		    'input=s' 	=> \$input_file,
		    'output=s'  => \$output_file,
		    );

my $usage = "annotate_genome [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> genome-file]";

if (!$rc || $help || @ARGV != 0) {
    seek(DATA, 0, 0);
    while (<DATA>) {
	last if /^=head1 COMMAND-LINE /;
    }
    while (<DATA>) {
	last if (/^=/);
	print $_;
    }
    exit($help ? 0 : 1);
}

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

__DATA__
