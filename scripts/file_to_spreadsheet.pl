#!/usr/bin/perl -w
use Getopt::Std;
use Spreadsheet::WriteExcel;
use strict;
use Data::Dumper;
use Carp;


#
# This is a SAS Component
#


=head1 svr_file_to_spreadsheet -f filename 

Writes the contents of the tab separated file on STDIN to a spreadsheet

The output is an xls file

Example: svr_file_to_spreadsheet -f test.xls  < test.txt 

=head2 Command-Line Options

=over 4

=item -f 

The file name given to the output xls format spreadsheet.

=back

=head2 Output Format

The output is a file in xls format

=cut


our ($opt_f, $opt_u);
getopt("f:u:");
if (!$opt_f) {
	die "Usage: No File Specified (-f missing)\n";
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
