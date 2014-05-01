
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
use gjoseqlib;

=head1 NAME

rast2-add-contigs

=head1 SYNOPSIS

rast2-add-contigs [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Add the given contig data to the genome object.

=head1 COMMAND-LINE OPTIONS

rast2-annotate-proteins-kmer-v2 [-io] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit
	--min-hits      minimum number of Kmer hits required for a call to be
	                made
	--max-gap       maximum size of a gap allowed for a call to be made

=cut

my @options = (options_common(), options_contigs(1));

my($opt, $usage) = describe_options("rast2-add-contigs %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

if (!$opt->contigs)
{
    print "A value for the --contigs flag must be provided.\n";
    print($usage->text);
    exit 1;
}

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->add_contigs($genome_in, get_params_for_contigs($opt));

write_output($genome_out, $opt);
