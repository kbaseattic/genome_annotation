
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
				    ["readable|r", 'Show output in human-readable form'],
				    @options);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if @ARGV != 1;

my $batch_id = shift;

my $client = get_annotation_client($opt);

my $status = $client->pipeline_batch_status($batch_id);

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
     


