use strict;

=head1 NAME

annotate_plant_transcripts

=head1 SYNOPSIS

annotate_plant_transcripts [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> genome-file]

=head1 DESCRIPTION

This routine takes as input a file containing a JSON-encoded KBaseGenomes.Genome object.
The genome object should contain a set of protein sequences, and this
routine attempts to annotate the protein sequences directly with metabolic
functions.

Example:

    annotate_plant_transcripts < input.genome > output.genome

where 'input.genome' is an initial KBase object of 'genome' type,
and 'output.genome' is 'input.genome' augmented by the feature annotations.

=head1 COMMAND-LINE OPTIONS

Usage: annotate_plant_transcripts --url service-url < genome-file > genome-file
Usage: annotate_plant_transcripts --url service-url --input genome-file --output genome-file

    --url    --- Optional URL for alternate KBase server (D: http://bio-data-1.mcs.anl.gov/services/genome_annotation)

    --input  --- Option to read genome-typed-object from input file instead of from STDIN

    --output --- Option to write enhanced genome-typed-object to output file instead of STDOUT

=head1 AUTHORS

seaver

=cut

use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;

use Getopt::Long;

my $help;
my $input_file;
my $output_file;
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

#Check to make sure genome is from supported domain
my $trouble = 0;
if ($input_genome->{domain} !~ m/^Plant/o) {
    $trouble = 1;
    die "Application \'annotate_plant_transcripts\' can only annotate plant transcripts\n";
}

my $output_genome = $anno_server->annotate_proteins_kmer_v1($input_genome,{dataset_name=>"Release70",kmer_size=>8});

$json->pretty(1);
print $out_fh $json->encode($output_genome);
close($out_fh);
