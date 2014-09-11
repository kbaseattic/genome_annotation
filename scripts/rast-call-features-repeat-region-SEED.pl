
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-call-features-repeat-regions-SEED

=head1 SYNOPSIS

rast-call-features-repeat-regions-SEED [--min-identity iden] [--min-length len] [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Call features from the contigs in the given genome using the Prodigal gene caller.

=head1 COMMAND-LINE OPTIONS

rast-call-features-repeat-regions-SEED [-io] [long options...] < input > output

            -i --input      file from which the input is to be read
            -o --output     file to which the output is to be written
            --help          print usage message and exit
            --url           URL for the genome annotation service
            --call-5S       Call 5S RNA features
            --call-SSU      Call SSU RNA features
            --call-LSU      Call LSU RNA features
    
=cut

my @options = (options_common(), options_repeat_regions_seed());

my($opt, $usage) = describe_options("rast-call-features-repeat-regions-SEED %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my %params;
$params{$_} = $opt->{$_} foreach grep { exists($opt->{$_}) } qw(min_length min_identity);

my $genome_out = $client->call_features_repeat_region_SEED($genome_in, \%params);

write_output($genome_out, $opt);
