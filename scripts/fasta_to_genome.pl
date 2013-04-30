
=head1 NAME

fasta_to_genome

=head1 SYNOPSIS

fasta_to_genome [--input fasta-file] [--output output-file] [--source source] [--source-id source-id] scientific-name domain genetic-code [< contig-fasta] [> genome-data]

=head1 DESCRIPTION

This command takes as input a file containing contigs for a prokaryotic genome
and it creates as output a genomeTO

Example:

    fasta_to_genome  scientific-name domain genetic-code < contigs > genomeTO

=head1 COMMAND-LINE OPTIONS

Usage: fasta_to_genome [--input fasta-file] [--output output-file] [--source source] [--source-id source-id] scientific-name domain genetic-code [< contig-fasta] [> genome-data]


    scientific-name --- This is the full name you want associated with the genome
                        (e.g., 'Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)')

    domain          --- This should be 'Bacteria' or 'Archaea' for prokaryotes.

    genetic-code    --- This is normally 11 for prokaryotes, but is 4 for a group including the <i>Mycoplasmas</i>.
                        (Check with NCBI if you are not sure which code to use.)
    
=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut


use strict;
use gjoseqlib;
use Bio::KBase::IDServer::Client;
use JSON::XS;

use Getopt::Long;

my $help;
my $scientific_name;
my $domain;
my $genetic_code;
my $source;
my $source_id;
my $input_file;
my $output_file;

my $rc = GetOptions('help'        => \$help,
		    'source=s' 	  => \$source,
		    'source-id=s' => \$source_id,
		    'input=s' 	  => \$input_file,
		    'output=s'    => \$output_file,
		    );

my $usage = "fasta_to_genome [--input fasta-file] [--output output-file] [--source source] [--source-id source-id] scientific-name domain genetic-code [< contig-fasta] [> genome-data]";

if (!$rc || $help || @ARGV != 3) {
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


my $scientific_name = shift;
my $domain = shift;
my $genetic_code = shift;

# my $id_server = Bio::KBase::IDServer::Client->new('http://bio-data-1.mcs.anl.gov/services/idserver');
my $id_server = Bio::KBase::IDServer::Client->new('https://kbase.us/services/idserver');

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

my $genome_idx = $id_server->allocate_id_range('kb|g', 1);
my $genome_id = "kb|g.$genome_idx";

my $contigs = [];
my $genome = {
    id 		    => $genome_id,
    scientific_name => $scientific_name,
    domain 	    => $domain,
    genetic_code    => $genetic_code,
    contigs         => $contigs,
    features        => [],
    defined($source) ? (source => $source) : (),
    defined($source_id) ? (source_id => $source_id) : (),
};
    
while (my($id, $def, $seq) = read_next_fasta_seq($in_fh))
{
    push(@$contigs, { id => $id, dna => $seq });
}

close($in_fh);

my $jdumper = JSON::XS->new;
$jdumper->pretty(1);
print $out_fh $jdumper->encode($genome);
close($out_fh);

__DATA__
