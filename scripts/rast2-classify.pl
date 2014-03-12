
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use gjoseqlib;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast2-classify 

=head1 SYNOPSIS

rast2-classify > classifiers

=head1 DESCRIPTION

Enumerate the available kmer classifiers.

=head1 COMMAND-LINE OPTIONS

rast2-classify [-h] [long options...] classifier < input > output
	-h --help   print usage message and exit
	--url       URL for the genome annotation service

=cut

my @options = (options_common(), options_classifier());

my($opt, $usage) = describe_options("%c %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

@ARGV == 1 or die $usage->text;

my $classifier = shift;

my $client = get_annotation_client($opt);

my $bin_only = 1;

my $det_fh;
if ($opt->detailed_output_file)
{
    open($det_fh, ">", $opt->detailed_output_file) or die "Cannot write $opt->{detailed_output_file}: $!";
    $bin_only = 0;
}

my $unc_fh;
if ($opt->unclassified_output_file)
{
    open($unc_fh, ">", $opt->unclassified_output_file) or die "Cannot write $opt->{unclassified_output_file}: $!";
    $bin_only = 0;
}

my $in_fh = get_input_fh($opt);
my $out_fh = get_output_fh($opt);

my $dna = [];
while (my($id, $def, $seq) = gjoseqlib::read_next_fasta_seq($in_fh))
{
    push(@$dna, [$id, $seq]);
}


my $bins;
if ($bin_only)
{
    $bins = $client->classify_into_bins($classifier,$dna);
}
else
{
    my($raw, $unassigned);
    ($bins, $raw, $unassigned) = $client->classify_full($classifier,$dna);
    print $det_fh $raw;
    close($det_fh);
    print $unc_fh "$_\n" foreach @$unassigned;
    close($unc_fh);
}
    
for my $g (sort { $bins->{$b} <=> $bins->{$a} } keys %$bins)
{
    print $out_fh "$g\t$bins->{$g}\n";
}
