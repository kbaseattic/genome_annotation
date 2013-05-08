# -*- perl -*-
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

use Getopt::Long;

=head1 NAME                                                                                                                                
cluster_objects                                                                                                                                 
=head1 SYNOPSIS                                                                                                                         
cluster_objects < related > sets

=head1 DESCRIPTION                                                                                                                   
cluster related objects into sets
input is a two column table of related objects (id's)
Using transitive closure, cluster_objects produces all sets of related objects.

Example:                                                                                                                                   
Computing protein sets for the genomes in an OTU

   use strict;
   use Data::Dumper;

   my $usage = "usage: pan_genome OTU > protein.families";
   my $otu = shift @ARGV;
   if (! defined($otu)) { die $usage };

   my @genomes = map { chop; $_ } `echo $otu | get_relationship_IsCollectionOf -to id | cut -f2`;

   open(CLUSTER,"| cluster_objects") || die "could not open clustering";
   foreach my $x (@genomes)
   {
       foreach my $y (@genomes) {
           if ($x lt $y) {
               my @output = `corresponds "$x" "$y" -a 0.6 2> /dev/null`;
               foreach $_ (@output) { 
                   chomp;
                   my($peg1,$sc,$peg2) = split(/\t/,$_);
                   print CLUSTER "$peg1\t$peg2\n";
               }
           }
       }
   }
   close(CLUSTER);

For otu 512, this produces

kb|g.1997.peg.115   kb|g.3378.peg.1544  kb|g.3460.peg.1206  kb|g.927.peg.1035
kb|g.1997.peg.1288  kb|g.3378.peg.1435  kb|g.3460.peg.1811  kb|g.927.peg.436
kb|g.1997.peg.1047  kb|g.3378.peg.361   kb|g.3460.peg.935   kb|g.927.peg.316
kb|g.3378.peg.1151  kb|g.3460.peg.1464  kb|g.927.peg.136
kb|g.1997.peg.1662  kb|g.3378.peg.563   kb|g.3460.peg.93    kb|g.927.peg.1468
kb|g.1997.peg.952   kb|g.3378.peg.1558  kb|g.3460.peg.1788  kb|g.927.peg.724
kb|g.1997.peg.466   kb|g.3378.peg.1417  kb|g.3460.peg.1355  kb|g.927.peg.1073
kb|g.1997.peg.1649  kb|g.3378.peg.243   kb|g.3460.peg.773   kb|g.927.peg.1157
kb|g.1997.peg.849   kb|g.3378.peg.1530  kb|g.3460.peg.1726  kb|g.927.peg.19
kb|g.1997.peg.1449  kb|g.3378.peg.222   kb|g.3460.peg.1424  kb|g.927.peg.1506
kb|g.1997.peg.1417  kb|g.3378.peg.28    kb|g.3460.peg.1264  kb|g.927.peg.311
kb|g.1997.peg.1031  kb|g.927.peg.425
...


=head1 COMMAND-LINE OPTIONS                                                                                              
Usage: cluster_objects < related > sets

=head1 AUTHORS                                                                                                                         
L<The SEED Project|http://www.theseed.org>

=cut

use Carp;
use strict;

# usage: cluster_objects < related > sets

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


my %to_cluster;
my %in_cluster;

my $nxt = 1;



while (defined($_ = <STDIN>))
{
    if ($_ =~ /^(\S[^\t]*\S)\t(\S[^\t]*\S)/)
    {
	my $obj1 = $1;
	my $obj2 = $2;
	my $in1 = $to_cluster{$obj1};
	my $in2 = $to_cluster{$obj2};

	if (defined($in1) && defined($in2) && ($in1 != $in2))
	{
	    push(@{$in_cluster{$in1}},@{$in_cluster{$in2}});
	    foreach $_ (@{$in_cluster{$in2}})
	    {
		$to_cluster{$_} = $in1;
	    }
	    delete $in_cluster{$in2};
	}
	elsif ((! defined($in1)) && defined($in2))
	{
	    push(@{$in_cluster{$in2}},$obj1);
	    $to_cluster{$obj1} = $in2;
	}
	elsif ((! defined($in2)) && defined($in1))
	{
	    push(@{$in_cluster{$in1}},$obj2);
	    $to_cluster{$obj2} = $in1;
	}
	elsif ((! defined($in1)) && (! defined($in2)))
	{
	    $to_cluster{$obj1} = $to_cluster{$obj2} = $nxt;
	    $in_cluster{$nxt} = [$obj1,$obj2];
	    $nxt++;
	}
    }
}

foreach my $cluster (keys(%in_cluster))
{
    my $set = $in_cluster{$cluster};
    print join("\t",sort @$set),"\n";
}

__DATA__
