package Bio::KBase::GenomeAnnotation::CmdHelper;

use strict;
use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;
use File::Slurp qw(read_file write_file);
use base 'Exporter';
our @EXPORT_OK = qw(load_input write_output get_annotation_client write_text_output
		    get_params_for_kmer_v1 get_params_for_kmer_v2 get_params_for_glimmer3
		    options_help options_common options_kmer_v1 options_kmer_v2 options_rrna_seed options_glimmer3
		    options_repeat_regions_seed options_export options_classifier
		    get_input_fh get_output_fh
		   );
our %EXPORT_TAGS = (all => \@EXPORT_OK);

sub options_common
{
    return (['input|i=s', 'file from which the input is to be read'],
	    ['output|o=s', 'file to which the output is to be written'],
	    &options_help());
	    
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
