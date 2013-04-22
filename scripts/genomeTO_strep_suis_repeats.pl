
=head1 NAME

genomeTO_to_strep_suis_repeats

=head1 SYNOPSIS

genomeTO_to_strep_suis_repeats [-url] < input > output

=head1 DESCRIPTION

This routine takes as input a file containing a JSON-encoded
genomeTO.  It invokes a tool to locate a set of control sites (often called "repeats")
in Streptococcus suis genomes.  This tool locates the sites and adds the features
representing them to the genomeTO.

Example:

    genomeTO_to_strep_suis_repeats < input > output

=head1 COMMAND-LINE OPTIONS

Usage: genomeTO_to_strep_suis_repeats [--url service-url]  < genome-file  > extended-genome-file
Usage: genomeTO_to_strep_suis_repeats --input genome-file  --output genome-file [--url service-url]

    --url    --- Optional URL for alternate KBase server (D: http://bio-data-1.mcs.anl.gov/services/genome_annotation)

    --input  --- Option to read genome-typed-object from input file instead of from STDIN

    --output --- Option to write enhanced genome-typed-object to output file instead of STDOUT

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut

use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;

use Getopt::Long;
use strict;
use Data::Dumper;

my $input_file;
my $output_file;
# my $url = "http://bio-data-1.mcs.anl.gov/services/genome_annotation";
my $url = "https://kbase.us/services/genome_annotation";

my $help;
my $rc = GetOptions('url:s'     => \$url,
		    'help'      => \$help,
		    'input:s' 	=> \$input_file,
		    'output:s'  => \$output_file,
		    );

my $usage = "genomeTO_to_strep_suis_repeats [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> extended-genome-file]";

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

#...Result appears to differ from $genomeTO, so apparently =NOT= "In Place"!
my $resultTO = $kbase_server->get_strep_suis_repeats($genomeTO);

$json->pretty(1);
print $out_fh $json->encode($resultTO);
close($out_fh);

__DATA__
