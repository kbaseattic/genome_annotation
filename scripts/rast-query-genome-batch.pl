
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
use Bio::KBase::HandleService;
use JSON::XS;

=head1 NAME

rast-query-genome-batch

=head1 SYNOPSIS

rast-query-genome-batch batch-id
    
=head1 DESCRIPTION
    
=head1 COMMAND-LINE OPTIONS

rast-process-genome-batch.pl [-h] [long options...] directory-of-genome-objects
    --workflow      File containing genome processing workflow
                    specification
    -h --help       print usage message and exit
    --url           URL for the genome annotation service

=cut

my @options = (options_help());

my($opt, $usage) = describe_options("%c %o batch-id",
				    @options);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 1;

my $batch_id = shift;

my $client = get_annotation_client($opt);

my $status = $client->pipeline_batch_query($batch_id);

