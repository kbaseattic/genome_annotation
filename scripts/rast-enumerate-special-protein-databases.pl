
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-enumerate-special-protein-databases

=head1 SYNOPSIS

rast-enumerate-special-protein-databases 

=head1 DESCRIPTION

Enumerate the available special protein databases.

=head1 COMMAND-LINE OPTIONS

rast-enumerate-special-protein-databases [-h] [long options...] < input > output
	-h --help   print usage message and exit
	--url       URL for the genome annotation service

=cut

my @options = options_help();

my($opt, $usage) = describe_options("%c %o < input > output",
				    @options);

print($usage->text), exit if $opt->help;

my $client = get_annotation_client($opt);

my $out = $client->enumerate_special_protein_databases();

print "$_\n" foreach @$out;
