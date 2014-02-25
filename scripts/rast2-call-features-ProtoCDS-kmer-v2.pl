
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast2-call-features-ProtoCDS-kmer-v2

=head1 SYNOPSIS

rast2-call-features-ProtoCDS-kmer-v2 [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Call features from the contigs in the given genome using the Prodigal gene caller.

=head1 COMMAND-LINE OPTIONS

rast2-annotate-proteins-kmer-v2 [-io] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit
	--min-hits      minimum number of Kmer hits required for a call to be
	                made
	--max-gap       maximum size of a gap allowed for a call to be made

=cut

my @options = (options_common(), options_kmer_v2());


my($opt, $usage) = describe_options("rast2-call-features-ProtoCDS-kmer-v2 %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->call_features_ProtoCDS_kmer_v2($genome_in, get_params_for_kmer_v2($opt));

write_output($genome_out, $opt);
