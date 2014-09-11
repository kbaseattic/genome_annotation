
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast2-call-features-rRNA-SEED

=head1 SYNOPSIS

rast2-call-features-rRNA-SEED [--call-5S] [--call-LSU] [--call-SSU] [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Call features from the contigs in the given genome using the Prodigal gene caller.

=head1 COMMAND-LINE OPTIONS

rast2-call-features-rRNA-SEED [-io] [long options...] < input > output

            -i --input      file from which the input is to be read
            -o --output     file to which the output is to be written
            --help          print usage message and exit
            --url           URL for the genome annotation service
            --call-5S       Call 5S RNA features
            --call-SSU      Call SSU RNA features
            --call-LSU      Call LSU RNA features
    
=cut

my @options = (options_common(), options_rrna_seed());

my($opt, $usage) = describe_options("rast2-call-features-rRNA-SEED %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my @ftypes = qw(5S LSU SSU);
my @req_ftypes = grep { $opt->{"call_" . lc($_)} } @ftypes;
if (@req_ftypes == 0)
{
    push(@req_ftypes, 'ALL');
}

my $genome_out = $client->call_features_rRNA_SEED($genome_in, \@req_ftypes);

write_output($genome_out, $opt);
