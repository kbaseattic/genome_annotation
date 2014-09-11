
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-create-genome-from-RAST

=head1 SYNOPSIS

rast-create-genome-from-RAST [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Create a genome object based on a RAST job.

=head1 COMMAND-LINE OPTIONS

rast-annotate-proteins-kmer-v2 [-io] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit
	--min-hits      minimum number of Kmer hits required for a call to be
	                made
	--max-gap       maximum size of a gap allowed for a call to be made

=cut

my @options = (options_genome_out(), options_help());

my($opt, $usage) = describe_options("rast-create-genome-from-RAST %o job-number-or-genome-id < input > output",
				    @options);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 1;

my $spec = shift;

my $client = get_annotation_client($opt);

my $genome_out = $client->create_genome_from_RAST($spec);

write_output($genome_out, $opt);
