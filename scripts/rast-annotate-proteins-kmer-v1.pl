
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-annotate-proteins-kmer-v1 

=head1 SYNOPSIS

rast-annotate-proteins-kmer-v1 [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Annotate proteins using the SEED version-2 kmers.

=head1 COMMAND-LINE OPTIONS

rast-annotate-proteins-kmer-v1 [-io] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit
	--min-hits      minimum number of Kmer hits required for a call to be
	                made
	--max-gap       maximum size of a gap allowed for a call to be made

=cut

my @options = (options_common(), options_kmer_v1());


my($opt, $usage) = describe_options("rast-annotate-proteins-kmer-v1 %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->annotate_proteins_kmer_v1($genome_in, get_params_for_kmer_v1($opt));

write_output($genome_out, $opt);
