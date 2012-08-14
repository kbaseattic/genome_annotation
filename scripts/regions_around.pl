use strict;
use Data::Dumper;
use Carp;

#
# This is a SAS Component
#

=head1 regions_around


This routine takes as input a file in which one column is composed
of fids, and it outputs a sorted version (sorted by fid locations).


Example:

    regions_around [-c N] -base Type -upstream Upstream -downstream Downstream < input > output

The standard input should be a tab-separated table (i.e., each line
is a tab-separated set of fields).  Normally, the last field in each
line would contain the identifer. If another column contains the identifier
use

    -c N

where N is the column (from 1) that contains the subsystem.

This is a pipe command. The input is taken from the standard input, and the
output is to the standard output.

=cut

my $usage = "usage: regions_around [-c N] -base Type -upstream Upstream -downstream Downstream < input > output";

my $column;
my $input_file;
my $output_file;
my $base;
my $upstream;
my $downstream;

use Bio::KBase::CDMI::CDMIClient;
use Bio::KBase::Utilities::ScriptThing;
my $kbO = Bio::KBase::CDMI::CDMIClient->new_for_script('c=i'          => \$column,
						       'base=s'       => \$base,
						       'upstream=i'   => \$upstream,
						       'downstream=i' => \$downstream,
						       'input=s'      => \$input_file,
						       'output=s'     => \$output_file);

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
@tuples = map { [&extract_region($_->[0],$base,$upstream,$downstream),$_->[1]] } @tuples;
foreach $_ (@tuples)
{
    print $out_fh $_->[1],"\t",$_->[0],"\n";
}

sub extract_region {
    my($loc,$base,$upstream,$downstream) = @_;

    my($contig,$beg,$strand,$len);
    if ($loc =~ /^(\S+)_(\d+)([\-\+])(\d+)$/)
    {
	($contig,$beg,$strand,$len) = ($1,$2,$3,$4);
    }
    elsif ($loc =~ /^(\S+)_(\d+)$/)
    {
	($contig,$beg,$strand,$len) = ($1,$2,'+',1);
    }
    else
    {
	print STDERR "badly formed input location: $loc\n";
	return $loc;
    }

    my $ln;
    if ($base =~/^s/i)
    {
	$ln  = $upstream + $downstream + 1;
	if ($strand eq '+')
	{
	    $b   = $beg - $upstream;
	}
	else
	{
	    $b   = $beg + $upstream;
	}
    }	    
    elsif ($base =~ /^e/i)
    {
	$ln  = $upstream + $downstream + 1;
	if ($strand eq '+')
	{
	    $b   = ($beg + $len - 1) - $upstream;
	}
	else
	{
	    $b   = ($beg - ($len - 1)) + $upstream;
	}
    }	    
    elsif ($base =~ /^i/i)
    {
	$ln  = $upstream + $downstream + $len + 1;
	if ($strand eq '+')
	{
	    $b   = $beg - $upstream;
	}
	else
	{
	    $b   = $beg + $upstream;
	}
    }	    
    return "$contig\_$b$strand$ln";
}

