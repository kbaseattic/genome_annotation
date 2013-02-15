use strict;

#
# This is a SAS Component
#

=head1 NAME

a_not_b

=head1 SYNOPSIS

a_not_b File1 File2 > set-theoretic-difference

=head1 DESCRIPTION

Finds the lines in File1 that are B<NOT> present in File2.

Example:

    a_not_b File1 File2 > set-theoretic-difference

where File1 and File2 are lists of assigned functions.

=head1 COMMAND-LINE OPTIONS

Usage: a_not_b File1 File2 > set-theoretic-difference

File1 --- Name of a text file

File2 --- Name of text-file to be compared.

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut

my($f1,$f2);
my $usage = "usage: a_and_b File1 File2 > intersection";

my $help;
use Getopt::Long;
my $rc = GetOptions('help' => \$help);

if (!$rc || $help || @ARGV < 2) {
    seek(DATA, 0, 0);
    while (<DATA>) {
	last if /^=head1 COMMAND-LINE /;
    }
    while (<DATA>) {
	last if (/^=/);
	print $_;
    }
    exit($help ? 0 : 1);
}

(
 ($f1 = shift @ARGV) && (-s $f1) &&
 ($f2 = shift @ARGV) && (-s $f2)
)
    || die $usage;

my %f2H = map { $_ ? ($_ => 1) : () } `cat $f2`;
my %keep;
foreach my $x (grep { ! defined($f2H{$_}) } `cat $f1`)
{
    $keep{$x} = 1;
}

foreach $_ (sort keys(%keep))
{
    print $_;
}

__DATA__
