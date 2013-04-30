
=head1 NAME

merge_features

=head1 SYNOPSIS 

merge_features [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> extended-genome-file]";

=head1 DESCRIPTION

This routines take a genomeTO as input.  It is presumed that
all of the called features are in the genomeTO, as well as many
that might be redundant (called by alternative methods).  This program
is supposed to select that set of features that are to be kept.

Example:

    merge_features < input > output

=head1 COMMAND-LINE OPTIONS

Usage:  merge_features  [--url service-url] < input > output
Usage:  merge_features  --input genome-file  --output genome-file  [--url service-url]

    --url    --- Optional URL for alternate KBase server (D: http://bio-data-1.mcs.anl.gov/services/genome_annotation)

    --input  --- Option to read genome-typed-object from input file instead of from STDIN

    --output --- Option to write enhanced genome-typed-object to output file instead of STDOUT

=head1 AUTHORS

L<The SEED Project|http://www.theseed.org>

=cut

use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;

use Getopt::Long;
use strict;
use Data::Dumper;

my $input_file;
my $output_file;
# my $url = "http://bio-data-1.mcs.anl.gov/services/genome_annotation";
my $url = "https://kbase.us/services/genome_annotation";

my $help;
my $rc = GetOptions('help'      => \$help,
		    'url=s'     => \$url,
		    'input=s' 	=> \$input_file,
		    'output=s'  => \$output_file,
		    );

my $usage = "merge_features [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> extended-genome-file]";

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

my $anno_server = Bio::KBase::GenomeAnnotation::Client->new($url);

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
my $json = JSON::XS->new;

my $input_genome;
{
    local $/;
    undef $/;
    my $input_genome_txt = <$in_fh>;
    $input_genome = $json->decode($input_genome_txt);
}

my $all_features = $input_genome->{features};
my $keep         = &create_merged_features($all_features);
$input_genome->{features} = $keep;

use JSON::XS;
$json->pretty(1);
print $out_fh $json->encode($input_genome);
close($out_fh);

sub create_merged_features {
    my($all_features) = @_;

    my @overlap_sets = &build_overlap_sets($all_features);
    my @multiple = grep { @$_ > 1 } @overlap_sets;
    my @single  = grep { @$_ == 1 } @overlap_sets;
    my @to_keep = map { $_->[0] } @single;

# This has a number of problems.  Basically, it is not picking the "best" from a set of CDSs -
# just the first (which may be the longest -- I am not even sure)

    foreach my $overlap_set (@multiple)
    {
	my $grouped_overlap_set = &group($overlap_set);
	my $supported = &supported_sets($grouped_overlap_set);
	my $j;
	for ($j=0; ($j < @$grouped_overlap_set); $j++)
	{
	    my $group = $grouped_overlap_set->[$j];
	    my $type  = $group->[0]->{type};
	    if ($type ne 'CDS')
	    {
		push(@to_keep,@$group);
	    }
	    elsif ($supported->{$j})
	    {
		push(@to_keep,$group->[$supported->{$j}]);  # push the first supported version
		                                                  # You really want "the best" supported version
	    }
	    elsif ((($j > 0) || (&overlap_sz($group->[$j-1],$group->[$j]) < 100)) &&
		   (($j < (@$group - 1)) || (&overlap_sz($group->[$j],$group->[$j+1]) < 100)))
	    {
		push(@to_keep,$group->[0]);   # again, this takes the first, when you want the best
	    }
	}
    }
    return \@to_keep;
}

sub supported_sets {
    my($grouped) = @_;

    my $supported = {};
    my $groupI;
    for ($groupI = 0; ($groupI < @$grouped); $groupI++)
    {
	my $group = $grouped->[$groupI];
	my $i;
	for ($i=0; ($i < @$group) && (! &supported_call($group->[$i])); $i++) {}
	if ($i < @$group)
	{
	    $supported->{$groupI} = $;
	}
    }
    return $supported;
}

sub supported_call {
    my($f) = @_;
    if ($f->{type} =~ /rna/i) { return 1 }
    if ($f->{type} eq 'CDS')
    {
	my $ann = $f->{annotations};
	my $i;
	for ($i=0; ($i < @$ann) && (! &record_of_solid_kmers($ann->[$i])); $i++) {}
	return ($i < @$ann);
    }
    return 0;
}

sub record_of_solid_kmers {
    my($annot) = @_;
    
    return (($annot->[0] =~ /^Set function.* hits=(\d+)/s) && ($1 > 6)) ||
	   (($annot->[0] =~ /^called by projection.* kmers=(\d+)/s) && ($1 > 6));
}

sub group {
    my($overlap_set) = @_;

    my $grouped = [];
    my $i=0;
    while ($i < @$overlap_set)
    {
	my $group = [$overlap_set->[$i]];
	$i++;
	while (($i < @$overlap_set) && &groupable($group,$overlap_set->[$i]))
	{
	    push(@$group,$overlap_set->[$i]);
	    $i++;
	}
	push(@$grouped,$group);
    }
    return $grouped;
}

sub groupable {
    my($group,$x) = @_;

    if ($group->[0]->{type} ne $x->{type}) { return 0 }
    if ($group->[0]->{type} eq 'CDS')
    {
	if (&stop($group->[0]->{location}) eq &stop($x->{location}))
	{
	    return 1;
	}
    }
    return 0;
}

sub stop {
    my($loc) = @_;

    my $last = $loc->[-1];
    return ($last->[2] eq '+') ? ($last->[1] + ($last->[3] - 1)) : 
	                         ($last->[1] - ($last->[3] - 1));
}

sub build_overlap_sets {
    my($all_features) = @_;

    my @sorted = sort { &compare_locs($a->{location},$b->{location}) } @$all_features;
    my @sets;
    while (my $f = shift @sorted)
    {
	my $set = [$f];
	while ((@sorted > 0) && &overlaps($set,$sorted[0]))
	{
	    my $x = shift @sorted;
	    push(@$set,$x);
	}
	push(@sets,$set);
    }
    return @sets;
}

sub overlaps {
    my($set,$x) = @_;

    my $i;
    for ($i=0; ($i < @$set) && (! &overlaps_entry($set->[$i],$x)); $i++) {}
    return ($i < @$set);
}

use SeedUtils;

sub overlaps_entry {
    my($x,$y) = @_;

    my $loc1 = $x->{location};
    my $loc2 = $y->{location};
    if ($loc1->[0]->[0] ne $loc2->[0]->[0]) { return 0 }
    my($b1,$e1) = &bounds($loc1);
    my($b2,$e2) = &bounds($loc2);
    return &SeedUtils::between($b1,$b2,$e1) || &SeedUtils::between($b2,$b1,$e2);
}
						   
sub bounds {  ### assumes regions are on same contig and strand
    my($loc) = @_;

    my $b = $loc->[0]->[1];
    my $e = ($loc->[0]->[2] eq '+') ? ($loc->[-1]->[1] + ($loc->[-1]->[3] - 1)) : 
	                              ($loc->[-1]->[1] - ($loc->[-1]->[3] - 1));
    return ($b < $e) ? ($b,$e) : ($e,$b);
}

sub compare_locs {
    my($loc1,$loc2) = @_;
    return (($loc1->[0]->[0] cmp $loc2->[0]->[0]) or ($loc1->[0]->[1] <=> $loc2->[0]->[1]));
}

sub overlap_sz {
    my($x,$y) = @_;
    my $c1 = $x->{location}->[0]->[0];
    my $c2 = $y->{location}->[0]->[0];
    if ($c1 ne $c2) { return 0 }
    my($b1,$e1) = &bounds($x->{location});
    my($b2,$e2) = &bounds($y->{location});

    if ($e1 < $b1) { ($b1,$e1) = ($e1,$b1) }
    if ($e2 > $b2) { ($b2,$e2) = ($e2,$b2) }

    my $minE = ($e1 > $e2) ? $e2 : $e1;
    my $maxB = ($b1 > $b2) ? $b1 : $b2;

    if ($maxB > $minE) { return 0 }
    return ($minE - $maxB) + 1;
}

__DATA__
