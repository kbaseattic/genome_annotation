
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-set-metadata

=head1 SYNOPSIS

rast-set-metadata [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Set the given metadata fields in the genome object.

=head1 COMMAND-LINE OPTIONS

rast-set-metadata [-hio] [long options...] < input > output
	-i --input             file from which the input is to be read
	-o --output            file to which the output is to be written
	-h --help              print usage message and exit
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
=cut

my @options = (options_common(), options_genome_metadata());

my($opt, $usage) = describe_options("rast-set-metadata %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->set_metadata($genome_in, get_params_for_genome_metadata($opt));

write_output($genome_out, $opt);
