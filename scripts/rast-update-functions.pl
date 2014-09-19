
use strict;
use Getopt::Long::Descriptive;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::CmdHelper qw(:all);
=head1 NAME

rast-update-functions

=head1 SYNOPSIS

rast-update-functions [--input genome-file] [--output genome-file] functions-file [< genome-file] [> genome-file]

=head1 DESCRIPTION

Update feature functions.

=head1 COMMAND-LINE OPTIONS

rast-update-functions.pl [-hio] [long options...] functions-file < input > output
	-i --input      file from which the input is to be read
	-o --output     file to which the output is to be written
	-h --help       print usage message and exit
	--url           URL for the genome annotation service

=cut

my @options = options_common();

my($opt, $usage) = describe_options("%c %o functions-file < input > output",
				    @options);

print($usage->text), exit if $opt->help;
print($usage->text), exit 1 if (@ARGV != 1);

my $functions_file = shift;

open(F, "<", $functions_file) or die "Cannot open $functions_file: $!\n";
my @functions;
while (<F>)
{
    chomp;
    my($id, $func) = split(/\t/);
    push(@functions, [$id, $func]);
}
close(F);


my $genome_in = load_input($opt);

my $client = get_annotation_client($opt);

my $genome_out = $client->update_functions($genome_in, \@functions, {});

write_output($genome_out, $opt);
