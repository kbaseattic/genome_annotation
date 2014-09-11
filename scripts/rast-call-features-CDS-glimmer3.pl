use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);

=head1 NAME

rast-call-features-CDS-glimmer3

=head1 SYNOPSIS

rast-call-features-CDS-glimmer3 [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Call features from the contigs in the given genome using the Glimmer3 gene caller.

=head1 COMMAND-LINE OPTIONS

rast-call-features-CDS-glimmer3 [-io] [long options...] < input > output
	-i --input              file from which the input is to be read
	-o --output             file to which the output is to be written
	--help                  print usage message and exit
	--url                   URL for the genome annotation service
	--min-training-len      Minimum size of a contig to be used for
	                        training glimmer3
=cut

my @options = (options_common(), options_glimmer3());


my($opt, $usage) = describe_options("rast-call-features-CDS-glimmer3 %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->call_features_CDS_glimmer3($genome_in, get_params_for_glimmer3($opt));

write_output($genome_out, $opt);
