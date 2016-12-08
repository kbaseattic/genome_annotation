
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);

=head1 NAME

rast-export-genome

=head1 SYNOPSIS

rast-export-genome [--input genome-file] [--output exported-data] format [< genome-file] [> exportdata]

=head1 DESCRIPTION

Export the given genome using the specified format. Format is one of

=over 4

=item genbank

A genbank file.

=item embl

An EMBL file.

=item gff

A GFF3 file.

=back

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
          spreadsheet_txt RAST-style spreadsheet (tab-separated text format)
	  spreadsheet_xls RAST-style spreadsheet (Excel XLS format)
          seed_dir        SEED directory (in tar.gz format)
	  feature_data    Tabular form of feature data
	  protein_fasta   Protein translations in fasta format
	  contig_fasta    Contig DNA in fasta format
	  feature_dna     Feature DNA sequences in fasta format
	  patric_features PATRIC features.tab format
          patric_specialty_genes PATRIC spgenes.tab format
          patric_genome_metadata PATRIC genome metadata format
	  gff             GFF format
	  embl            EMBL format

=cut

my @valid_format = map { $_->[0] } @Bio::KBase::GenomeAnnotation::CmdHelper::export_formats;
my %valid_format = map { $_->[0] => 1 } @Bio::KBase::GenomeAnnotation::CmdHelper::export_formats;

my @options = (options_common(), options_export(), [], options_export_formats());

my($opt, $usage) = describe_options("rast-export-genome %o format < input > output",
				    @options);

print($usage->text), exit  if $opt->help;
die($usage->text) unless @ARGV == 1;

my $format = lc(shift);
$format =~ s/-/_/g;

if (!$valid_format{$format})
{
    die "Invalid format $format. Valid formats are @valid_format\n";
}
	       
my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $formatted = $client->export_genome($genome_in, $format, $opt->feature_type);

write_text_output($formatted, $opt);
