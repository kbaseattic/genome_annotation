
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
use GenomeTypeObject;

=head1 NAME

rast-export-SEED

=head1 SYNOPSIS

rast-export-SEED [--input genome-file] directory [< stdin]

=head1 DESCRIPTION

Export the given genome as a SEED directory.

=head1 COMMAND-LINE OPTIONS

rast-export-genome [-hio] [long options...] format < input > output
	-i --input           file from which the input is to be read
	-o --output          file to which the output is to be written
	-h --help            print usage message and exit
	--url                URL for the genome annotation service
	--feature-type       Include this feature type in output. If no
	                     feature-types specified, include all feature
	                     types
	Supported formats: 
	  genbank         Genbank format
	  genbank_merged  Genbank format as single merged locus, suitable for Artemis
	  feature_data    Tabular form of feature data
	  protein_fasta   Protein translations in fasta format
	  contig_fasta    Contig DNA in fasta format
	  feature_dna     Feature DNA sequences in fasta format
	  gff             GFF format
	  embl            EMBL format

=cut

my @options = (options_genome_in(), options_help() );

my($opt, $usage) = describe_options("%c %o output-directory < input > output",
				    @options);

print($usage->text), exit  if $opt->help;
print($usage->text), exit 1  unless @ARGV == 1;

my $dir = shift;

-d $dir or mkdir($dir) or die "Cannot mkdir $dir: $!";


my $gobj = GenomeTypeObject->initialize(load_input($opt));
$gobj->write_seed_dir($dir, { map_CDS_to_peg => 1 });
