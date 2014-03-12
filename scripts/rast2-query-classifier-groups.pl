
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast2-query-classifier-groups 

=head1 SYNOPSIS

rast2-query-classifier-groups classifier-name > classifiers

=head1 DESCRIPTION

Enumerate the available kmer classifiers.

=head1 COMMAND-LINE OPTIONS

rast2-query-classifier-groups [-io] [long options...] < input > output
	--help          print usage message and exit

=cut

my @options = (options_help());

my($opt, $usage) = describe_options("%c %o classifier-name",
				    @options);

print($usage->text), exit if $opt->help;

@ARGV == 1 or die $usage->text;

my $classifier = shift;

my $client = get_annotation_client($opt);

my $groups = $client->query_classifier_groups($classifier);

for my $gname (sort { my($aa) = $a =~ /group(\d+)/; my($bb) = $b =~ /group(\d+)/;
		      defined($aa) && defined($bb) ? ($aa <=> $bb) : ($a cmp $b) } keys %$groups)
{
    print join("\t", $gname, @{$groups->{$gname}}), "\n";
}

