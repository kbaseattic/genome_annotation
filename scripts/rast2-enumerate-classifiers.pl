
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast2-enumerate-classifiers 

=head1 SYNOPSIS

rast2-enumerate-classifiers > classifiers

=head1 DESCRIPTION

Enumerate the available kmer classifiers.

=head1 COMMAND-LINE OPTIONS

rast2-enumerate-classifiers [-h] [long options...] < input > output
	-h --help   print usage message and exit
	--url       URL for the genome annotation service

=cut

my @options = (options_help());

my($opt, $usage) = describe_options("%c %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $client = get_annotation_client($opt);

my $classifiers = $client->enumerate_classifiers();

print "$_\n" foreach @$classifiers;
