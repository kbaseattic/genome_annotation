
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-call-features-insertion-sequences

=head1 SYNOPSIS

rast-call-features-insertion-sequences [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Call insertion sequences.

=head1 COMMAND-LINE OPTIONS

rast-call-features-insertion-sequences.pl [-hio] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	-h --help       print usage message and exit
	--url           URL for the genome annotation service

=cut

my @options = options_common();

my($opt, $usage) = describe_options("%c %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->call_features_insertion_sequences($genome_in);

write_output($genome_out, $opt);
