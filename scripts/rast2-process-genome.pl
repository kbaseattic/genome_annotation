use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);

=head1 NAME

rast2-process-genome

=head1 SYNOPSIS

  rast2-process-genome --scientific-name "Genus species"
                       --domain Bacteria
		       --genetic-code 11
		       --output-format genbank
		       < genome.contigs.fa > genome.gbk

  rast2-process-genome --batch-input-directory my.genomes.dir

  rast2-process-genome --batch-input-file my.genomes.file

=head1 DESCRIPTION

Annotate bacterial genomes using the RAST2 pipeline. Eukaryotic genomes are not supported.

This program may be used in one of two modes.

In immediate mode a single genome is processed in real time. The rast2-process-genome
script will not return until processing has completed.

In batch mode one or more genomes are submitted for processing by the backend
annotation services. The rast2-process-genome script returns after the data has been
uploaded, and emits a job identifier that may be used to query the status of the
computation and to retrieve the results when they have been completed.

=head2 Input Formats and Metadata

The RAST2 pipeline takes as input genomic contigs in L<FASTA|http://en.wikipedia.org/wiki/FASTA_format> format. Each contig
must have a unique identifier on each definition line.

In immediate mode, the contigs are provided either from the standard input (by default) or
via a file specfied using the C<--input> parameter.

In batch mode, the location of the contig is defined by the batch-mode description files. See
L<BATCH FORMAT> for more details.

The RAST2 pipeline also requires a minimal amount of metadata in order to complete its analysis.

=over 4

=item *

The scientific name of the genome. This takes the form minimally of "Genus species" and
may include a strain specifier. For the most part the pipeline does not interpret the
name, but some components of the annotation pipeline may produce better results if given
an accurate scientific name.

=item *

The genetic code of the genome. This is used to use the proper DNA to protein translation table.
Valid values for this option are 11 for Archaea, most Bacteria, most Virii, and some Mitochondria;
and 4 for Mycoplasmaea, Spiroplasmaea, Ureoplasmaea, and Fungal Mitochondria.

=item *

The phylogenetic domain of the genome. Again, there are components of the annotation pipeline that
produce more accurate results if the domain is specified accurately.  Valid values for this option
are Bacteria and Archaea.

=back

In immediate mode the genome metadata must be specified using the C<--scientific-name>, C<--domain>,
and C<--genetic-code> parameters.

In batch mode the metadata is defined by the metadata files in the input batch formats.

=head1 BATCH FORMAT

The RAST2 pipeline supports the processing of numerous genomes in a single run. Since this
may be time consuming, the processing is done in the background on the RAST2 service. Once the
genomes are submitted, the rast2-process-genome script will return a job identifier which may be used by
the L<rast2-status> script.

To use the batch format, one must specify the set of input data files to be processed
as well as the metadata for each. We support two different mechanisms for providing this information.

First is the directory-of-genomes mechanism. Here, the user provides a directory with the following structure:

 dir/genome-1/meta
             /contigs.fa
    /genome-2/meta
             /contigs.fa

In other words, a directory contains a set of directories, one per genome to be processed. In
each of the genome directories we have two files; contigs.fa contains the DNA contigs to
be processed and meta contains a set of key-value pairs:

 scientific-name Genome name
 genetic-code value
 domain value

The values in these pairs correspond to the values as described in L<Input Formats and Metadata>. The use
of this mechanism is triggered by the C<--batch-input-directory> parameter.

The second mechanism is triggere by the C<--batch-input-file> parameter. Here, the user specifies
a file containing the following tab-delimited values:

=over 4

=item *

B<Identifier>. This is an identifier for the genome, unique to the input file, used
to distinguish the processed output genomes.

=item *

B<Contigs-file>. This is the pathname on the local machine to the contigs to be processed.

=item *

B<Scientific-name>. The scientific name for this genome.

=item *

B<Genetic-code>. The genetic code for this genome.

=item *

B<Domain>. The domain of this genome.

=back

=cut
 
my @options = (['scientific-name=s', 'Scientific name for the genome to be processed.'],
	       ['domain=s', 'Domain of the genome to be processed (Bacteria or Archaea)'],
	       ['genetic-code=i', 'Genetic code of the genome to be processed (4 or 11)'],
	       ['output-format=s@', 'Output format. This option may be repeated to generate multiple output formats in batch mode. Defaults to genbank.',
		{ default => ['genbank'] }],
	       ['input|i=s', 'Process a single genome from the given file. Input defaults to standard input.'],
	       ['output|o=s', 'Output file to which a single processed genome will be written.'],
	       ['batch-input-directory=s', 'Process a batch of genomes from the given directory.'],
	       ['batch-input-file=s', 'Process a batch of genomes defined by the given file.'],
	       ['workflow=s', 'Workflow definition for this genome or batch.'],
	       ['url=s', 'URL to genome annotation service.'],
	       ['help|h', 'Show this help message'],
	       );

my($opt, $usage) = describe_options("%c %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $client = get_annotation_client($opt);

my $workflow;

if ($opt->workflow)
{
    die "User specification of workflow is not supported yet.\n";
}
else
{
    $workflow = $client->default_workflow();
}

#print Dumper($workflow);
#print Dumper($opt);

my $genome = load_input($opt);

my $res = $client->run_pipeline($genome, $workflow);

for my $format (@{$opt->output_format})
{
    my $txt = $client->export_genome($res, $format, []);
    my $fh = get_output_fh($opt);
    print $fh $txt;
    last;
}
    

   
  
