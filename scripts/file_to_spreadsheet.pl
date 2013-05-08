#!/usr/bin/perl -w

use Getopt::Std;
use Spreadsheet::WriteExcel;
use strict;
use Data::Dumper;
use Carp;


#
# This is a SAS Component
#


=head1 NAME

file_to_spreadsheet

=head1 SYNOPSIS

file_to_spreadsheet -f filename [-u hyperlink_template]

=head1 DESCRIPTION

Writes the contents of the tab separated file on STDIN to a spreadsheet

The output is an xls file

Example: file_to_spreadsheet -f test.xls -u < test.txt 

=head1 COMMAND-LINE OPTIONS

Usage: file_to_spreadsheet -f test.xls [-u http://pubseed.theseed.org/?page=Annotation&feature=PEG] < test.txt

    -f --- The filename for the output xls format spreadsheet.

    -u --- An optional hyperlink template; PEG will be replaced by any FIG-IDs in the input file

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut


my $help;
our ($opt_f, $opt_u);
use Getopt::Long;
my $rc = GetOptions('help' => \$help,
		    'f=s'  => \$opt_f,
		    'u:s'  => \$opt_u,
    );

if (!$rc || $help || @ARGV != 0) {
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



# Create a new Excel workbook
my $workbook = Spreadsheet::WriteExcel->new($opt_f);
if (!$workbook) {die "$opt_f workbook failure"};

# Add a worksheet
my $worksheet = $workbook->add_worksheet();
if (!$worksheet) {die "$opt_f worksheet failure"};

my $row = 0;
my $col = 0;
my @ctot;
while (<STDIN>) {
	chomp;
	my @line = split("\t", $_);

	$col = 0;
	foreach my $element (@line) {
	    $ctot[$col] += length($element);
	    if ($opt_u && $element =~ /fig\|/)
	    {
		my $link = $opt_u;
		$link =~ s/PEG/$element/;
		
		$worksheet->write_url($row, $col, $link, $element);
	    }
	    else
	    {
		$worksheet->write($row, $col, $element);
	    }
	    $col++;
	}
	$row++;
}
for my $col (0..$#ctot)
{
    my $avg = int($ctot[$col] / $row);
    if ($avg > 5)
    {
	$worksheet->set_column($col, $col, $avg);
    }
}

__DATA__
