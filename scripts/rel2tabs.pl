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


# usage: rel2tabs < 2-column-table > multi-column.version";

=head1 NAME                                                                                                                                

rel2tabs                                                                                                                                 

=head1 SYNOPSIS                                                                                                                         

rel2tabs < 2-column-table > multi-column.version                                                                                                                           
=head1 DESCRIPTION                                                                                                                   

detailed_description_of_purpose_of_script                                                                                    
Given:

1   kb|g.22250.peg.2733286
1   kb|g.19976.c.0_2409901+588
1   CDS
1   hypothetical protein
2   kb|g.22250.peg.2733262
2   kb|g.19976.c.0_2379935+1359
2   CDS
2   Chloride channel protein
3   kb|g.22250.peg.27332783   kb|g.19976.c.0_2395935+615
3   CDS
3   hypothetical protein4   kb|g.22250.peg.2733279
4   kb|g.19976.c.0_2396827+2494   CDS
                                       
rel2tabs will produce:

kb|g.22250.peg.2733286  kb|g.19976.c.0_2409901+588  CDS hypothetical
kb|g.22250.peg.2733262  kb|g.19976.c.0_2379935+1359 CDS Chloride
kb|g.22250.peg.2733278  kb|g.19976.c.0_2395935+615  CDS hypothetical
kb|g.22250.peg.2733279  kb|g.19976.c.0_2396827+249  CDS hypothetical


Example:                                                                                                                                   

   example_of_use                                                                                                                         

example_description                                                                                                                       

=head1 COMMAND-LINE OPTIONS                                                                                              

Usage: rel2tabs < 2-column-table > multi-column.version                                                                                                               
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
