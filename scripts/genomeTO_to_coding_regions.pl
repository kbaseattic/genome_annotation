
=head1 NAME

genomeTO_to_coding_regions

=head1 SYNOPSIS

genomeTO_to_coding_regions [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> extended-genome-file]

=head1 DESCRIPTION

This routine takes as input a file containing a JSON-encoded
genomeTO.  It invokes a kmer-based service which attempts to
locate portions of coding regions based on the notion
of signature kmers.  In the process it proposes functional
assignments of the encoded proteins.  It adds the estimates 
to the hash that makes up the genomeTO as a field with the
key 'DNA_kmer_data'.  The field is a pointer to a list of 5-tuples, each
containing [Contig,Beg,End,KmerHits,Function].

Example:

    genomeTO_to_coding_regions < input > output

=head1 COMMAND-LINE OPTIONS

Usage: genomeTO_to_coding_regions [--url service-url]  < genome-file  > extended-genome-file]
Usage: genomeTO_to_coding_regions  --input genome-file  --output genome-file  [--url service-url]

    --url    --- Optional URL for alternate KBase server (D: http://bio-data-1.mcs.anl.gov/services/genome_annotation)

    --input  --- optional name of genome-typed-object input file

    --output --- optional name of genome-typed-object output file

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
my $rc = GetOptions('help'      => \$help,
		    'url=s'     => \$url,
		    'input=s' 	=> \$input_file,
		    'output=s'  => \$output_file,
		    );

my $usage = "genomeTO_to_coding_regions [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> extended-genome-file]";

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

use CloseGenomes;

use JSON::XS;
my $tmp = $input_genome->{contigs};
my @raw_contigs = map { [$_->{id},'',$_->{dna}] }  @$tmp;
my $contigs = \@raw_contigs;
my $close = &CloseGenomes::get_coding_regions($contigs, {});

$input_genome->{DNA_kmer_data} = $close;
$json->pretty(1);
print $out_fh $json->encode($input_genome);
close($out_fh);

__DATA__
