
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
use JSON::XS;

=head1 NAME

rast-get-default-workflow

=head1 SYNOPSIS

rast-get-default-workflow [--input genome-file] [--output genome-file] [< genome-file] [> genome-file]

=head1 DESCRIPTION

Retrieve default workflow.

=head1 COMMAND-LINE OPTIONS

rast-get-default-workflow.pl [-ho] [long options...] > output
	-o --output     file to which the output is to be written
	-h --help       print usage message and exit
	--url           URL for the genome annotation service

=cut

my @options = (options_genome_out(), options_help());

my($opt, $usage) = describe_options("%c %o  > output",
				    @options);

print($usage->text), exit if $opt->help;

my $client = get_annotation_client($opt);

my $wf = $client->default_workflow();

my $fh = get_output_fh($opt);

my $json = JSON::XS->new->pretty(1);
print $fh $json->encode($wf);
