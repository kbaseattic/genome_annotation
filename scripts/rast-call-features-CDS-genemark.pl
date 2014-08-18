
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast2-call-features-CDS-genemark

=head1 SYNOPSIS

rast2-call-features-CDS-genemark [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Call features from the contigs in the given genome using the Genemark gene caller.

=head1 COMMAND-LINE OPTIONS

rast2-annotate-proteins-kmer-v2 [-io] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit

=cut

my @options = options_common();


my($opt, $usage) = describe_options("%c %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->call_features_CDS_genemark($genome_in);

write_output($genome_out, $opt);
