
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-call-features-crispr

=head1 SYNOPSIS

rast-call-features-crispr [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Call CRISPR features in the given genome.

=head1 COMMAND-LINE OPTIONS

rast-call-features-crispr [-hio] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	-h --help       print usage message and exit
	--url           URL for the genome annotation service
    
=cut

my @options = (options_common());

my($opt, $usage) = describe_options("rast-call-features-crispr %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->call_features_crispr($genome_in);

write_output($genome_out, $opt);
