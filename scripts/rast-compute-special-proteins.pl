
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-compute-special-proteins

=head1 SYNOPSIS

rast-compute-special-proteins [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Compute the instances of the specialty proteins in the given genome.

=head1 COMMAND-LINE OPTIONS

rast-annotate-proteins-kmer-v2 [-io] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit
	--min-hits      minimum number of Kmer hits required for a call to be
	                made
	--max-gap       maximum size of a gap allowed for a call to be made

=cut

my @options = (options_genome_in(), options_help());

my($opt, $usage) = describe_options("%c %o < input > output",
				    ['db=s@', "Database name to search (option may be repeated). Defaults to all available databases"],
				    ['output|o=s', 'Output file'],
				    @options);

print($usage->text), exit if $opt->help;

my $db_list = $opt->db || [];

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $out = $client->compute_special_proteins($genome_in, $db_list);

my $out_fh;
if ($opt->output)
{
    open($out_fh, ">", $opt->output) or die "Cannot write output file ". $opt->output . ": $!\n";
}
else
{
    $out_fh = \*STDOUT;
}

print $out_fh join("\t", @$_), "\n" foreach @$out;
