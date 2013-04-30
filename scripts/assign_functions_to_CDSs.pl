
=head1 NAME

assign_functions_to_CDSs.pl

=head1 SYNOPSIS

assign_functions_to_CDSs [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Reads a genome-typed-object file, and annotates each coding DNA sequaence (CDS)
with the predicted function based on FIGfams and Kmers.

Example:

    assign_functions_to_CDSs  < genome.file  > genome.with.CDS-functions

=head1 COMMAND-LINE OPTIONS

Usage: assign_functions_to_CDSs [--url service-url]  < genome-file > genome-file
Usage: assign_functions_to_CDSs --input genome-file  --output genome-file  --url service-url

    --url    --- Optional URL for alternate KBase server (D: http://bio-data-1.mcs.anl.gov/services/genome_annotation)

    --input  --- Option to read genome-typed-object from input file instead of from STDIN

    --output --- Option to write enhanced genome-typed-object to output file instead of STDOUT

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut

use strict;
use gjoseqlib;
use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;

use Getopt::Long;

my $input_file;
my $output_file;
# my $url = "http://bio-data-1.mcs.anl.gov/services/genome_annotation";
my $url = "https://kbase.us/services/genome_annotation";

my $help;
my $rc = GetOptions('help'      => \$help,
		    'url=s'     => \$url,
		    'input=s' 	=> \$input_file,
		    'output=s'  => \$output_file,
		    );

my $usage = "assign_functions_to_CDSs [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> genome-file]";

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

__DATA__
