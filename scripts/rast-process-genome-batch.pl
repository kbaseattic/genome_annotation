
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
use Bio::KBase::HandleService;
use JSON::XS;

=head1 NAME

rast-submit-genome-batch

=head1 SYNOPSIS

rast-submit-genome-batch directory-of-genome-objects
    
=head1 DESCRIPTION

Submit a batch of genomes for processing with RASTtk.

The specified directory will contain the genomes to be processed. Each file in the
directory contains a genome typed object.

By default the standard RASTtk pipeline will be run on each genome. If the
C<--workflow> flag is provided, the specified workflow file will be used instead.
    
=head1 COMMAND-LINE OPTIONS

rast-process-genome-batch.pl [-h] [long options...] directory-of-genome-objects
    --workflow      File containing genome processing workflow
                    specification
    -h --help       print usage message and exit
    --url           URL for the genome annotation service

=cut

my @options = (options_workflow_specification(), options_help());

my($opt, $usage) = describe_options("%c %o directory-of-genome-objects",
				    @options);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 1;

my $dir = shift;

my $client = get_annotation_client($opt);
my $hservice = Bio::KBase::HandleService->new();

my $workflow;

if ($opt->workflow)
{
    open(F, "<", $opt->workflow) or die "Cannot open workflow file " . $opt->workflow . ": $!\n";
    local $/;
    undef $/;
    $workflow = decode_json(scalar <F>);
    close(F);
}
else
{
    $workflow = $client->default_workflow();
}

die Dumper($workflow);

#
# For each file in the directory, we push to the handle service and retain
# the handle. The handles are then sent to the genome annotation service
# to begin processing.
#

opendir(D, $dir) or die "Cannot open genome directory $dir: $!\n";

my @files = sort { $a <=> $b } grep { !/^\./ && -f "$dir/$_" } readdir(D);

my @handles;

for my $file (@files)
{
    my $path = "$dir/$file";

    my $handle = $hservice->upload($path);
    print "Uploaded $path: " . Dumper($handle);
    push(@handles, $handle);
}

my $batch_id = $client->pipeline_batch_start(\@handles, $workflow);

