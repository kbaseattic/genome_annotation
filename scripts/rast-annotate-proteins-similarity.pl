
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-annotate-proteins-similarity 

=head1 SYNOPSIS

rast-annotate-proteins-similarity [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Annotate proteins using the SEED version-2 kmers.

=head1 COMMAND-LINE OPTIONS

rast-annotate-proteins-similarity [-io] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit

=cut

my @options = (options_common(), options_similarity());


my($opt, $usage) = describe_options("%c %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->annotate_proteins_similarity($genome_in, get_params_for_similarity($opt));

write_output($genome_out, $opt);
