
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
use Bio::KBase::HandleService;
use JSON::XS;

=head1 NAME

rast-download-genome-batch

=head1 SYNOPSIS

rast-download-genome-batch batch-id
    
=head1 DESCRIPTION
    
=head1 COMMAND-LINE OPTIONS

rast-process-genome-batch.pl [-h] [long options...] directory-of-genome-objects
    --workflow      File containing genome processing workflow
                    specification
    -h --help       print usage message and exit
    --url           URL for the genome annotation service
    
=cut

my @options = (options_help());

our($opt, $usage) = describe_options("%c %o batch-id output-directory",
				    ["verbose|v", "Verbose output"],
				    @options);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 2;

my $batch_id = shift;
my $dir = shift;

-d $dir or die "Output directory $dir does not exist\n";

my $client = get_annotation_client($opt);
my $hservice = Bio::KBase::HandleService->new();

my $batch_status = $client->pipeline_batch_status($batch_id);
my $status = $batch_status->{details};

for my $ent (@$status)
{
    my $gid = $ent->{genome_id};
    my $filename = $ent->{filename};
    
    if ($ent->{status} ne 'completed' && 0)
    {
	print STDERR "Skipping $gid: not complete (status is $ent->{status})\n" if $opt->verbose;
	next;
    }

    my $out_file = $gid;
    $out_file = $filename if $filename;
    
    save($ent->{stdout}, "$dir/$out_file.stdout");
    save($ent->{stderr}, "$dir/$out_file.stderr");
    save($ent->{output}, "$dir/$out_file.gto");
}

sub save
{
    my($h, $file) = @_;
    if (!ref($h))
    {
	return;
    }
    my $url = get_url($h);

    print "Save $url to $file\n" if $opt->verbose;
    unlink($file);
    $hservice->download($h, $file);
}

sub get_url
{
    my($h) = @_;
    return "$h->{url}/node/$h->{id}?download";
}
     


