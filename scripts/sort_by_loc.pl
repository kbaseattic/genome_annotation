use strict;
use Data::Dumper;
use Carp;

#
# This is a SAS Component
#

=head1 NAME

sort_by_loc

=head1 SYNOPSIS

sort_by_loc  [-c N]  < input  > sorted.input

=head1 DESCRIPTION

This routine takes as input a file in which one column is composed
of locations, and it outputs a sorted version (sorted by locations).

Example:

    sort_by_loc [-c N] < input > output

The standard input should be a tab-separated table (i.e., each line
is a tab-separated set of fields).  Normally, the last field in each
line would contain the identifer. If another column contains the location,
use

    -c N

where N is the column (from 1) that contains the location.

This is a pipe command. The input is normally taken from the standard input,
and the output is written to standard output.

=head1 COMMAND-LINE OPTIONS

Usage:  sort_by_loc  [-c N]  < input  > sorted.input
Usage   sort_by_loc  [-c N]  --input input-file  --output sorted-input

    -c N     --- The number of the column (from 1) that contains the location

    --input  --- Option to read from input file instead of from STDIN

    --output --- Option to write to output file instead of STDOUT

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut

my $usage = "usage: sort_by_loc [-c N] < input  > sorted.input ";

my $column;
my $input_file;
my $output_file;
use Bio::KBase::CDMI::CDMIClient;
use Bio::KBase::Utilities::ScriptThing;
my $kbO = Bio::KBase::CDMI::CDMIClient->new_for_script('c=i' => \$column,
						       'input=s' 	=> \$input_file,
						       'output=s'  => \$output_file);

if (! $kbO) { print STDERR $usage; exit }

my $in_fh;
if ($input_file)
{
    open($in_fh, "<", $input_file) or die "Cannot open $input_file: $!";
}
else
{
    $in_fh = \*STDIN;
}

my $out_fh;
if ($output_file)
{
    open($out_fh, ">", $output_file) or die "Cannot open $output_file: $!";
}
else
{
    $out_fh = \*STDOUT;
}

my @tuples = Bio::KBase::Utilities::ScriptThing::GetBatch($in_fh, 1000000, $column);
@tuples = sort { &by_loc($a->[0],$b->[0]) } @tuples;
foreach $_ (@tuples)
{
    print $out_fh $_->[1],"\n";
}

sub by_loc {
    my($x,$y) = @_;

    if ($x =~ /^(\S+)_(\d+)[\+\-]\d+$/)
    {
	my($c1,$s1) = ($1,$2);
	if ($y =~ /^(\S+)_(\d+)[\+\-]\d+$/)
	{
	    my($c2,$s2) = ($1,$2);
	    return (($c1 cmp $c2) or ($s1 <=> $s2));
	}
    }
    return 0;
}

__DATA__
