
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-create-genome

=head1 SYNOPSIS

rast-create-genome [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Create a new empty genome object.

If a GenBank file is specified using the --from-genbank parameter, initialize the
genome object data and metadata from the specified GenBank file.

The RAST2 pipeline requires a minimal amount of metadata in order to complete its analysis:

=over 4

=item *

The scientific name of the genome. This takes the form minimally of "Genus species" and
may include a strain specifier. For the most part the pipeline does not interpret the
name, but some components of the annotation pipeline may produce better results if given
an accurate scientific name.

=item *

The phylogenetic domain of the genome. Again, there are components of the annotation pipeline that
produce more accurate results if the domain is specified accurately.  Valid values for this option
are Bacteria and Archaea.

=item *

The genetic code of the genome. This is used to use the proper DNA to protein translation table.
Valid values for this option are 11 for Archaea, most Bacteria, most Virii, and some Mitochondria;
and 4 for Mycoplasmaea, Spiroplasmaea, Ureoplasmaea, and Fungal Mitochondria.

=back



=head1 COMMAND-LINE OPTIONS

rast-create-genome [-ho] [long options...] > output
	-o --output            file to which the output is to be written
	-h --help              print usage message and exit
        --from-genbank gb-file
	--url                  URL for the genome annotation service
	--genome-id            Genome identifier
	--scientific-name      Scientific name (Genus species strain) for the
	                       genome
	--domain               Domain (Bacteria/Archaea/Virus/Eukaryota) for
	                       the genome
	--genetic-code         Genetic code for the genome (probably 11 or 4
	                       for bacterial genomes)

	--source               Source (external database) name for this genome
	--source-id            Identifier for this genome in the source
	                       (external database)
	--contigs              Fasta file containing DNA contig data


=cut

my @options = (options_genome_out(), options_genome_metadata(), options_contigs(), options_help());
my($opt, $usage) = describe_options("rast-create-genome %o > output",
				    ["from-genbank=s", "Create from this genbank file"],
				    @options);

print($usage->text), exit if $opt->help;

my $client = get_annotation_client($opt);

my $genome_out;

if ($opt->from_genbank)
{
    my $txt;
    if (open(F, "<", $opt->from_genbank))
    {
	undef $/;
	$txt = <F>;
	close(F);
    }
    else
    {
	die "Cannot open genbank file " . $opt->from_genbank . ": $!\n";
    }
    $genome_out = $client->create_genome_from_genbank($txt);
}
else
{
    $genome_out = $client->create_genome(get_params_for_genome_metadata($opt));
    
    if ($opt->{contigs})
    {
	$genome_out = $client->add_contigs($genome_out, get_params_for_contigs($opt));
    }
}

write_output($genome_out, $opt);
