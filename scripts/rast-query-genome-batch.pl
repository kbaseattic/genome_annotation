
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

rast-query-genome-batch.pl [-hr] [long options...] batch-id
    -r --readable   Show output in human-readable form
    --summary       Show summary of job status only
    --raw           Print the raw AWE status
    -h --help       print usage message and exit
    --url           URL for the genome annotation service
=cut

my @options = (options_help());

my($opt, $usage) = describe_options("%c %o batch-id",
				    ["readable|r", 'Show output in human-readable form'],
				    ["summary", 'Show summary of job status only'],
				    ["raw", "Print the raw AWE status"],
				    @options);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 1;

my $batch_id = shift;

my $client = get_annotation_client($opt);

my $batch_status = $client->pipeline_batch_status($batch_id);

if ($opt->raw)
{
    my $json = JSON::XS->new->pretty(1);
    print $json->encode($batch_status);
    exit 0;
}

my $status = $batch_status->{details};

my %stats;

$stats{$_->{status}}++ foreach @$status;

if ($opt->readable || $opt->summary)
{
    print "Status counts:\n";
    for my $s (sort { $a cmp $b } keys %stats)
    {
	printf "    %-14s %4d\n", "$s:", $stats{$s};
    }
    print "\n";

    print "Batch submit time:     $batch_status->{submit_date}\n";
    print "Batch start time:      $batch_status->{start_date}\n";
    print "Batch completion time: $batch_status->{completion_date}\n";
    print "\n";
}

exit if $opt->summary;

for my $ent (@$status)
{
    if ($opt->readable)
    {
	print "$ent->{genome_id}:\n";
	print "    status:          $ent->{status}\n";
	print "    creation date:   $ent->{creation_date}\n";
	print "    start date:      $ent->{start_date}\n";
	print "    completion date: $ent->{completion_date}\n";
	print "    stdout:          " . get_url($ent->{stdout}) . "\n" if $ent->{stdout};
	print "    stderr:          " . get_url($ent->{stderr}) . "\n" if $ent->{stderr};
	print "    output:          " . get_url($ent->{output}) . "\n" if $ent->{output};
	print "\n";
    }
    else
    {
	print join("\t", @$ent{qw(genome_id status creation_date completion_date)},
		   map { get_url($_) } @$ent{qw(stdout stderr output)},
		   ), "\n";
				   
    }
}


sub get_url
{
    my($h) = @_;
    
    return ref($h) ? "$h->{url}/node/$h->{id}?download" : "";
}
     


