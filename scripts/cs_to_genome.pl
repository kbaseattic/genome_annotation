
use strict;
use gjoseqlib;
use Bio::KBase::IDServer::Client;
use Bio::KBase::CDMI::Client;
use JSON::XS;
use Data::Dumper;

use Getopt::Long;

my $scientific_name;
my $domain;
my $genetic_code;
my $source;
my $source_id;
my $input_file;
my $output_file;


use Getopt::Long;

=head1 NAME                                                                                                                                
cs_to_genome                                                                                                                                 
=head1 SYNOPSIS                                                                                                                         
cs_to_genome [--output output-file] KBase-genome-id > genome-object                                                                                                                           
=head1 DESCRIPTION                                                                                                                   
detailed_description_of_purpose_of_script                                                                                    
Example:                                                                                                                                   
   example_of_use                                                                                                                         
example_description                                                                                                                       
=head1 COMMAND-LINE OPTIONS                                                                                              
Usage: cs_to_genome [--output output-file] KBase-genome-id > genome-object                                                                                                               
=head1 AUTHORS                                                                                                                         
L<The SEED Project|http://www.theseed.org>

=cut

use Carp;
use strict;

$| = 1;

my $help;

my $rc = GetOptions('output=s'    => \$output_file,
                    'help' => \$help,
		    );

if (!$rc || $help || @ARGV != 1) {
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



my $usage = "cs_to_genome [--output output-file] KBase-genome-id > genome-object";

@ARGV == 1 or die "Usage: $usage\n";

my $src_genome_id = shift;

# my $id_server = Bio::KBase::IDServer::Client->new('http://bio-data-1.mcs.anl.gov/services/idserver');
# my $cs = Bio::KBase::CDMI::Client->new('http://bio-data-1.mcs.anl.gov/services/cdmi_api');
my $id_server = Bio::KBase::IDServer::Client->new('https://kbase.us/services/idserver');
my $cs = Bio::KBase::CDMI::Client->new('https://kbase.us/services/cdmi_api');


my $out_fh;
if ($output_file)
{
    open($out_fh, ">", $output_file) or die "Cannot open $output_file: $!";
}
else
{
    $out_fh = \*STDOUT;
}

#my $genome_idx = $id_server->allocate_id_range('kb|g', 1);
#my $genome_id = "kb|g.$genome_idx";
my $genome_id = "kb|g.22186";

my $ginfo = $cs->genomes_to_genome_data([$src_genome_id]);
$ginfo = $ginfo->{$src_genome_id};

my $domain = $ginfo->{taxonomy};
$domain =~ s/;.*$//;

my $contigs = [];
my $features = [];
my $genome = {
    id 		    => $genome_id,
    scientific_name => $ginfo->{scientific_name},
    domain 	    => $domain,
    genetic_code    => $ginfo->{genetic_code},
    contigs         => $contigs,
    features        => $features,
};

my $fres = $cs->get_relationship_IsOwnerOf([$src_genome_id], [], [], ['id', 'feature_type', 'source_id', 'function']);
my $fitems = [ map { $_->[2] } @$fres ];

my $fids = [ map { $_->{id} } @$fitems ];
my $trans = $cs->get_relationship_Produces($fids, [], [], ['id', 'sequence']);
my $f2trans = { map { $_->[1]->{from_link} => $_->[2] } @$trans };

my $locs = $cs->get_relationship_IsLocatedIn($fids, [], ['begin', 'dir', 'len', 'ordinal', 'to_link'], []);
my $f2loc = {};
for my $l (@$locs)
{
    my $id = $l->[1]->{from_link};
    my $item = $l->[1];
    
    $f2loc->{$id}->[$item->{ordinal}] = $item;
}

my $c = $cs->genomes_to_contigs([$src_genome_id]);
my $src_contigs = $c->{$src_genome_id};
my $cseqs = $cs->contigs_to_sequences($src_contigs);
my $cprefix = "$genome_id.c";
my $cstart = $id_server->allocate_id_range($cprefix, scalar @$src_contigs) + 0;

my %cmap;
for my $ctg (@$src_contigs)
{
    my $nctg = "$cprefix.$cstart";
    $cstart++;
    push(@$contigs, { id => $nctg, dna => $cseqs->{$ctg} });
    $cmap{$ctg} = $nctg;
}

my %items;
for my $item (@$fitems)
{
    push(@{$items{$item->{feature_type}}}, $item);
}

for my $type (keys %items)
{
    my $fids = $items{$type};
    my $n = @$fids;

    my $prefix = "$genome_id.$type";
    my $start = $id_server->allocate_id_range($prefix, $n) + 0;

    for my $item (sort { my $afid = $a->{id};
			 my $bfid = $b->{id};
			 my $aloc = $f2loc->{$afid}->[0];
			 my $bloc = $f2loc->{$bfid}->[0];
			 $aloc->{to_link} cmp $bloc->{to_link} or
			     $aloc->{begin} <=> $bloc->{begin} }
		  @$fids)
    {
	my $feature = {};
	my $fid = $item->{id};
	my $nfid = "$prefix.$start";
	$start++;
	my $locs = $f2loc->{$fid};

	my @nloc;
	for my $loc (@$locs)
	{
	    my $nctg = $cmap{$loc->{to_link}};

	    push(@nloc, [$nctg, $loc->{begin}, $loc->{dir}, $loc->{len}]);
	}

	if (exists($f2trans->{$fid}))
	{
	    my $trans = $f2trans->{$fid};
	    $feature->{protein_translation} = $trans->{sequence};
	}
	$feature->{location} = [@nloc];
	$feature->{function} = $item->{function};
	$feature->{id} = $nfid;
	$feature->{type} = $type;
	$feature->{annotations} = [];
	$feature->{aliases} = [];
	push(@$features, $feature);
    }
}
    
my $jdumper = JSON::XS->new;
$jdumper->pretty(1);
print $out_fh $jdumper->encode($genome);
close($out_fh);


__DATA__
