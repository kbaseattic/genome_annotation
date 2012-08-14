use strict;
use Data::Dumper;
use Carp;

#
# This is a SAS Component
#

=head1 sort_by_ids


This routine takes as input a file in which one column is composed
of fids, and it outputs a sorted version (sorted by fid ID).


Example:

    sort_by_id [-c N] < input > output

The standard input should be a tab-separated table (i.e., each line
is a tab-separated set of fields).  Normally, the last field in each
line would contain the identifer. If another column contains the identifier
use

    -c N

where N is the column (from 1) that contains the subsystem.

This is a pipe command. The input is taken from the standard input, and the
output is to the standard output.

=cut

my $usage = "usage: sort_by_id [-c N] < input  > sorted.input ";

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
@tuples = sort { &by_fid($a->[0],$b->[0]) } @tuples;
foreach $_ (@tuples)
{
    print $out_fh $_->[1],"\n";
}

sub by_fid {
    my($x,$y) = @_;

    if ($x =~ /^kb\|g\.(\d+)\.([^\.]+)\.(\d+)$/)
    {
	my($g1,$t1,$n1) = ($1,$2,$3);
	if ($y =~ /^kb\|g\.(\d+)\.([^\.]+)\.(\d+)$/)
	{
	    my($g2,$t2,$n2) = ($1,$2,$3);
	    return ($g1 <=> $g2) or ($t1 cmp $t2) or ($n1 <=> $n2);
	}
    }
    return 0;
}

