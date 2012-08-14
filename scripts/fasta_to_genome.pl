
=head1 fasta_to_genome

This command takes as input a file containing contigs for a prokaryotic genome
and it creates as output a genomeTO

Example:

    fasta_to_genome scientific-name domain genetic-code < contigs > genomeTO

=over 4

=item scientific-name

This is the full name you want associated with the genome (e.g., 'Buchnera aphidicola str. Tuc7 (Acyrthosiphon pisum)')

=item domain

This should be 'Bacteria' or 'Archaea' for prokaryotes.

=item genetic-code

This is normally 11 for prokaryotes, but is 4 for a group including the <i>Mycoplasmas</i>.  Check NCBI
if you are not sure which to use.

=back
=cut


use strict;
use gjoseqlib;
use Bio::KBase::IDServer::Client;
use JSON::XS;

use Getopt::Long;

my $scientific_name;
my $domain;
my $genetic_code;
my $source;
my $source_id;
my $input_file;
my $output_file;

my $rc = GetOptions('source=s' 	  => \$source,
		    'source-id=s' => \$source_id,
		    'input=s' 	  => \$input_file,
		    'output=s'    => \$output_file,
		    );

my $usage = "fasta_to_genome [--input fasta-file] [--output output-file] [--source source] [--source-id source-id] scientific-name domain genetic-code [< contig-fasta] [> genome-data]";

@ARGV == 3 or die "Usage: $usage\n";
my $scientific_name = shift;
my $domain = shift;
my $genetic_code = shift;

my $id_server = Bio::KBase::IDServer::Client->new('http://bio-data-1.mcs.anl.gov/services/idserver');

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
