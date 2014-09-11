
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-resolve-overlapping-features

=head1 SYNOPSIS

rast-resolve-overlapping-features [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Resolve overlapping features in this genome.

=head1 COMMAND-LINE OPTIONS

rast-resolve-overlapping-features [-hio] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	-h --help       print usage message and exit
	--url           URL for the genome annotation service

=cut

my @options = options_common();

my($opt, $usage) = describe_options("rast-resolve-overlapping-features %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $params = {};
my $genome_out = $client->resolve_overlapping_features($genome_in, $params);

write_output($genome_out, $opt);
