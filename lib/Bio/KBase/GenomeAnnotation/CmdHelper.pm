package Bio::KBase::GenomeAnnotation::CmdHelper;

use strict;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;
use File::Slurp qw(read_file write_file);
use base 'Exporter';
use gjoseqlib;
use List::Util 'max';

our @EXPORT_OK = qw(load_input write_output get_annotation_client write_text_output
		    get_params_for_kmer_v1 get_params_for_kmer_v2 get_params_for_glimmer3 get_params_for_genome_metadata
		    get_params_for_contigs
		    options_help options_common options_kmer_v1 options_kmer_v2 options_rrna_seed options_glimmer3
		    options_repeat_regions_seed options_export options_classifier options_genome_metadata
		    options_genome_in options_genome_out options_contigs options_export_formats
		    get_input_fh get_output_fh 
		   );
our %EXPORT_TAGS = (all => \@EXPORT_OK);

sub options_genome_in
{
    return (['input|i=s', 'file from which the input is to be read'])
	    
}

sub options_genome_out
{
    return (['output|o=s', 'file to which the output is to be written']);
	    
}

sub options_common
{
    return (options_genome_in(),
	    options_genome_out(),
	    &options_help());
	    
}

sub options_contigs
{
    return(['contigs=s', 'Fasta file containing DNA contig data']);
}
	

sub options_help
{
    return (['help|h', 'print usage message and exit'],
	    ['url=s', 'URL for the genome annotation service'],
	    );
	    
}

sub options_kmer_v1
{
    return(['kmer-size=i' => "kmer size", { default => 8 } ],
	   ['dataset-name=s' => "kmer dataset name"],
	   ['score-threshold' => 'score threshold'],
	   ['hit-threshold' => 'hit threshold'],
	   ['sequential-hit-threshold' => 'sequential-hit threshold'],
	   ['min-hits=i', 'minimum number of Kmer hits required for a call to be made'],
	   ['max-gap=i',  'maximum size of a gap allowed for a call to be made'],
	   ['min-size=i' => 'minimum size of DNA feature to call', { default => 48 }],
	  );
}

sub options_kmer_v2
{
    return (['min-hits=i', 'minimum number of Kmer hits required for a call to be made'],
	    ['max-gap=i',  'maximum size of a gap allowed for a call to be made'],
	    );
	    
}

sub options_glimmer3
{
    return (['min-training-len=i',  'Minimum size of a contig to be used for training glimmer3', { default => 2000 }],
	    );
	    
}

sub options_rrna_seed
{
    return(['call-5S', 'Call 5S RNA features'],
	   ['call-SSU', 'Call SSU RNA features'],
	   ['call-LSU', 'Call LSU RNA features'],
	   );
}

sub options_repeat_regions_seed
{
    return(['min-identity=f', 'minimum BLAST idendity'],
	   ['min-length=i', 'minimum length'],
	   );
}

sub options_export
{
    return (['feature-type=s@', 'Include this feature type in output. If no feature-types specified, include all feature types', { default => [] }]);
}


sub options_classifier
{
    return (["detailed-output-file|d=s" => "File to write detailed output (reads and hit information)"],
	    ["unclassified-output-file|u=s" => "File to write unclassified read IDs to"]);

}

sub options_genome_metadata
{
    return (['genome-id=s' => "Genome identifier"],
	    ['scientific-name=s' => "Scientific name (Genus species strain) for the genome"],
	    ['domain=s' => "Domain (Bacteria/Archaea/Virus/Eukaryota) for the genome"],
	    ['genetic-code=i' => "Genetic code for the genome (probably 11 or 4 for bacterial genomes)"],
	    ['source=s' => "Source (external database) name for this genome"],
	    ['source-id=s' => "Identifier for this genome in the source (external database)"],
	   );
}

#
# This format list is copied from the backend rast_export.pl code. Changes in the
# export format list happen slowly, so we can manually keep this list up to date. Much
# simpler than requiring a roundtrip to the service to find the formats.
#
# Changes also need to be propagated to the documentation in the API spec.
#   
our @export_formats = ([genbank => "Genbank format"],
		       [genbank_merged => "Genbank format as single merged locus, suitable for Artemis"],
		       [feature_data => "Tabular form of feature data"],
		       [protein_fasta => "Protein translations in fasta format"],
		       [contig_fasta => "Contig DNA in fasta format"],
		       [feature_dna => "Feature DNA sequences in fasta format"],
		       [gff => "GFF format"],
		       [embl => "EMBL format"]);

sub options_export_formats
{
    my @options;
    push(@options, ["Supported formats:"]);

    my $len = max(map { length($_->[0]) } @export_formats);
    $len++;
    for my $fmt (@export_formats)
    {
	my($name, $desc) = @$fmt;
	push(@options, [sprintf("  %-${len}s $desc", $name)]);
    }
    return @options;
}

sub get_annotation_client
{
    my($opts) = @_;
    my $client = Bio::KBase::GenomeAnnotation::Client->new($opts->{url});
    return $client;
}

sub get_input_fh
{
    my($opts) = @_;

    my $fh;
    if ($opts->{input})
    {
	open($fh, "<", $opts->{input}) or die "Cannot open input file $opts->{input} :$!";
    }
    else
    {
	$fh = \*STDIN;
    }
    return $fh;
}	

sub get_output_fh
{
    my($opts) = @_;

    my $fh;
    if ($opts->{output})
    {
	open($fh, ">", $opts->{output}) or die "Cannot open input file $opts->{output} :$!";
    }
    else
    {
	$fh = \*STDOUT;
    }
    return $fh;
}	

sub load_input
{
    my($opts) = @_;

    my $fh;
    if ($opts->{input})
    {
	open($fh, "<", $opts->{input}) or die "Cannot open input file $opts->{input} :$!";
    }
    else
    {
	$fh = \*STDIN;
    }

    my $text = read_file($fh);
    undef $fh;
    my $obj = decode_json($text);
    return $obj;
}

sub write_output
{
    my($genome, $opts) = @_;

    my $coder = JSON::XS->new->pretty;
    my $text = $coder->encode($genome);
    if ($opts->{output})
    {
	write_file($opts->{output}, \$text);
    }
    else
    {
	write_file(\*STDOUT, \$text);
    }
}

sub write_text_output
{
    my($text, $opts) = @_;

    if ($opts->{output})
    {
	write_file($opts->{output}, \$text);
    }
    else
    {
	write_file(\*STDOUT, \$text);
    }
}

sub get_params_for_kmer_v2
{
    my($opt) = @_;
    my $params = {};
    $params->{min_hits} = $opt->{min_hits} if defined($opt->{min_hits});
    $params->{max_gap} = $opt->{max_gap} if defined($opt->{max_gap});
    return $params;
}

sub get_params_for_kmer_v1
{
    my($opt) = @_;
    my $params = {};
    for my $p (qw(kmer_size dataset_name score_threshold hit_threshold sequential_hit_threshold max_gaps min_hits min_size))
    {
	$params->{$p} = $opt->{$p} if defined($opt->{$p});
    }
    return $params;
}

sub get_params_for_glimmer3
{
    my($opt) = @_;
    my $params = {};
    for my $p (qw(min_training_len))
    {
	$params->{$p} = $opt->{$p} if defined($opt->{$p});
    }
    return $params;
}

sub get_params_for_genome_metadata
{
    my($opt) = @_;
    my $params = {};
    for my $p (qw(genome_id scientific_name domain genetic_code source source_id))
    {
	my $to = ($p eq 'genome_id') ? 'id' : $p;
	$params->{$to} = $opt->{$p} if defined($opt->{$p});
    }
    return $params;
}

sub get_params_for_contigs
{
    my($opt) = @_;

    if ($opt->{contigs})
    {
	my $fh;
	open($fh, "<", $opt->{contigs}) or die "Cannot open contigs data file $opt->{contigs}: $!";
	my @ctgs;
	while (my($id, $def, $seq) = read_next_fasta_seq($fh))
	{
	    push(@ctgs, { id => $id, dna => $seq });
	}
	close($fh);
	return \@ctgs;
    }
    return undef;
}
