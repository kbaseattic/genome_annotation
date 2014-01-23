
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast2-export-genome

=head1 SYNOPSIS

rast2-export-genome [--input genome-file] [--output exported-data] format [< genome-file] [> exportdata]

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

rast2-export-genome [-io] [long options...] format < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit
	--url           URL for the genome annotation service

=cut

my @valid_format = qw(gff genbank embl protein_fasta contig_fasta);
my %valid_format = map { $_ => 1 } @valid_format;

my @options = (options_common());

my($opt, $usage) = describe_options("rast2-export-genome %o format < input > output",
				    @options);

print($usage->text), exit  if $opt->help;
print($usage->text), exit 1  unless @ARGV == 1;

my $format = lc(shift);
$format =~ s/-/_/g;

if (!$valid_format{$format})
{
    die "Invalid format $format. Valid formats are @valid_format\n";
}
	       
my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $formatted = $client->export_genome($genome_in, $format);

write_text_output($formatted, $opt);
