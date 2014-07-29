
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast2-add-features

=head1 SYNOPSIS

rast2-add-features [--input genome-file] [--output genome-file] features-file [< genome-file] [> genome-file]

=head1 DESCRIPTION

Add a set of features to the genome.

=head1 COMMAND-LINE OPTIONS

rast2-annotate-proteins-kmer-v2 [-io] [long options...] < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	--help          print usage message and exit
	--min-hits      minimum number of Kmer hits required for a call to be
	                made
	--max-gap       maximum size of a gap allowed for a call to be made

=cut

my @options = options_common();

push(@options,
     [],
     ["Each line in the features-file contains the following fields:"],
     [" id           ID of the feature. A new feature ID will be assigned for this feature"],
     [" location     Location of the feature on the contig, in the format ContigID_<start-pos>[+-]<length>"],
     [" feature-type Type of this feature (CDS, rna, etc.)"],
     [" function     Function assigned to this feature."],
     [" aliases      Comma-separated list of aliases for this feature"],
    );
     

my($opt, $usage) = describe_options("rast2-add-features %o features-file < input > output",
				    @options);


print($usage->text), exit if ($opt->help || @ARGV != 1);

my $features_file = shift;

open(F, "<", $features_file) or die "Cannot open $features_file: $!\n";
my @features;
while (<F>)
{
    chomp;
    my($id, $loc, $type, $func, $aliases) = split(/\t/);
    push(@features, [$id, $loc, $type, $func, $aliases]);
}
close(F);

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->add_features($genome_in, \@features);

write_output($genome_out, $opt);
