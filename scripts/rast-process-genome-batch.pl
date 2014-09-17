
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
				    @options,
				    ['offset=i' => 'Submit genomes starting with this one. Zero-based, order based on alphabetical sort of contents of directory', { default => 0 }],
				    ['count=i' => 'Submit this many genomes, starting at the given offset.',
				    {  default => -1 }],
				    ['verbose|v' => 'Show verbose details during submission'],
				   );

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 1;

my $dir = shift;

my $client = get_annotation_client($opt);
my $hservice = Bio::KBase::HandleService->new();

my $workflow;

if ($opt->workflow)
{
    $workflow = parse_json_file($opt->workflow);
}
else
{
    $workflow = $client->default_workflow();
}

#
# For each file in the directory, we push to the handle service and retain
# the handle. The handles are then sent to the genome annotation service
# to begin processing.
#

opendir(D, $dir) or die "Cannot open genome directory $dir: $!\n";

my @files = sort { $a <=> $b } grep { !/^\./ && -f "$dir/$_" } readdir(D);

my @genomes;

if ($opt->offset > 0)
{
    splice(@files, 0, $opt->offset);
}

if ($opt->count != -1)
{
    splice(@files, $opt->count);
}

if ($opt->verbose)
{
    my $json = JSON::XS->new->pretty(1);
    print "Submitting workflow: \n";
    print $json->encode($workflow);
}

for my $file (@files)
{
    my $path = "$dir/$file";

    my $gobj = parse_json_file($path);
    my $handle = $hservice->upload($path);
    print "Uploaded $path: " . Dumper($handle, $gobj->{id}) if $opt->verbose;
    push(@genomes, { genome_id => $gobj->{id}, data => $handle, filename => $file });
}

my $batch_id = $client->pipeline_batch_start(\@genomes, $workflow);
print "$batch_id\n";

sub parse_json_file
{
    my($file) = @_;
    local $/;
    undef $/;
    my $fh;
    open($fh, "<", $file) or die "Cannot open $file: $!";
    return decode_json(scalar <$fh>);
}
