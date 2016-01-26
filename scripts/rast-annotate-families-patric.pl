
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-annotate-families-patric 

=head1 SYNOPSIS

rast-annotate-families-patric [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Determin PATRIC family membership for the given genome.

=head1 COMMAND-LINE OPTIONS

rast-annotate-families-patric [-io] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit

=cut

my @options = (options_common());


my($opt, $usage) = describe_options("rast-annotate-families-patric %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->annotate_families_patric($genome_in);

write_output($genome_out, $opt);
