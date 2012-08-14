use strict;

my($f1,$f2);
my $usage = "usage: a_and_b File1 File2 > intersection";

@ARGV == 2 or die $usage;

(
 ($f1 = shift @ARGV) && (-s $f1) &&
 ($f2 = shift @ARGV) && (-s $f2)
)
    || die $usage;

my %f2H = map { $_ ? ($_ => 1) : () } `cat $f2`;
my %keep;
foreach my $x (grep { $f2H{$_} } `cat $f1`)
{
    $keep{$x} = 1;
}

foreach $_ (sort keys(%keep))
{
    print $_;
}
