#
# This is a SAS Component
#
#
# Copyright (c) 2003-2006 University of Chicago and Fellowship
# for Interpretations of Genomes. All Rights Reserved.
#
# This file is part of the SEED Toolkit.
# 
# The SEED Toolkit is free software. You can redistribute
# it and/or modify it under the terms of the SEED Toolkit
# Public License. 
#
# You should have received a copy of the SEED Toolkit Public License
# along with this program; if not write to the University of Chicago
# at info@ci.uchicago.edu or the Fellowship for Interpretation of
# Genomes at veronika@thefig.info or download a copy from
# http://www.theseed.org/LICENSE.TXT.
#


# usage: rel2tabs < 2-column-table > 1-column.version";

=head1 NAME                                                                                                                                

rel2tabs                                                                                                                                 

=head1 SYNOPSIS                                                                                                                         

rel2tabs < 2-column-table > 1-column.version                                                                                                                           
=head1 DESCRIPTION                                                                                                                   

detailed_description_of_purpose_of_script                                                                                    

Example:                                                                                                                                   

   example_of_use                                                                                                                         

example_description                                                                                                                       

=head1 COMMAND-LINE OPTIONS                                                                                              

Usage: rel2tabs < 2-column-table > 1-column.version                                                                                                               
=head1 AUTHORS                                                                                                                         

L<The SEED Project|http://www.theseed.org>

=cut

use Carp;
use strict;
use Getopt::Long;

$| = 1;

my $help;
my $rc = GetOptions('help' => \$help);

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



$_ = <STDIN>;
while (defined($_) && ($_ =~ /^(\S+)/))
{
    my @set = ();
    my $curr = $1;
    while ($_ && ($_ =~ /^(\S+)\t(\S+)/) && ($1 eq $curr))
    {
	push(@set,$2);
	$_ = <STDIN>;
    }
    print join("\t", @set),"\n";
}

__DATA__
