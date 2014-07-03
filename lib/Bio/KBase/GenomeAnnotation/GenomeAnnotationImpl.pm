package Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;
use strict;
use Bio::KBase::Exceptions;
# Use Semantic Versioning (2.0.0-rc.1)
# http://semver.org 
our $VERSION = "0.1.0";

=head1 NAME

GenomeAnnotation

=head1 DESCRIPTION

API Access to the Genome Annotation Service.

Provides support for gene calling, functional annotation, re-annotation. Use to extract annotation in
formation about an existing genome, or to create new annotations.

=cut

#BEGIN_HEADER

use Safe;
use File::Temp;
use File::Slurp;
use Data::Dumper;
use Digest::MD5 'md5_hex';
use Time::HiRes 'gettimeofday';

use Bio::KBase::IDServer::Client;
use Bio::KBase::KmerAnnotationByFigfam::Client;
use KmerClassifier;
#use Bio::KBase::KIDL::Helpers qw(json_to_tempfile tempfile_to_json);
use IPC::Run qw(run);
use JSON::XS;
use File::Slurp;
use Bio::KBase::GenomeAnnotation::Awe;
use Bio::KBase::GenomeAnnotation::Shock;

use Bio::KBase::GenomeAnnotation::Glimmer;
use GenomeTypeObject;
use ANNOserver;
use SeedUtils;
use gjoseqlib;
use StrepRepeats;
use overlap_resolution;

use Bio::KBase::DeploymentConfig;

our $idserver_url = 'https://kbase.us/services/idserver';

sub _get_coder
{
    return JSON::XS->new->ascii->pretty->allow_nonref;
}

sub _call_using_strep_repeats
{
    my($self, $genome_in, $tool) = @_;
    
    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $parsed = StrepRepeats::get_raw_repeats($genome_in, $tool);

    my $event = {
	tool_name => $tool,
	execution_time => scalar gettimeofday,
	parameters => [],
	hostname => $self->{hostname},
    };

    my $idc = Bio::KBase::IDServer::Client->new($idserver_url);
    my $event_id = $genome_in->add_analysis_event($event);

    my $type = 'repeat_unit';

    for my $r (@$parsed)
    {
	my $loc = $r->{location};
	$genome_in->add_feature({
	    -id_client 	     => $idc,
	    -id_prefix 	     => $genome_in->{id},
	    -type 	     => $type,
	    -location 	     => $loc,
	    -function 	     => $r->{function},
	    -annotation      => $r->{annotation},
	    -analysis_event_id 	     => $event_id,
	});
    }

    $genome_in = $genome_in->prepare_for_return();

    return $genome_in;
}


#END_HEADER

sub new
{
    my($class, @args) = @_;
    my $self = {
    };
    bless $self, $class;
    #BEGIN_CONSTRUCTOR

    my $cfg = Bio::KBase::DeploymentConfig->new($ENV{KB_SERVICE_NAME} || "GenomeAnnotation");

    my $dir = $cfg->setting("kmer_v2_data_directory");
    $dir or die "Configuration parameter for kmer_v2_data_directory not set";
    -d $dir or die "Directory $dir for kmer_v2_data_directory does not exist";
	
    $self->{kmer_v2_data_directory} = $dir;

    my $dir = $cfg->setting("kmer_classifier_data_directory");
    #
    # Make these soft errors for now.
    #
    if (!$dir)
    {
	warn "Configuration parameter for kmer_classifier_data_directory not set";
    }
    elsif (! -d $dir)
    {
	warn "Directory $dir for kmer_classifier_data_directory does not exist";
    }
	
    $self->{kmer_classifier_data_directory} = $dir;

    my $i = $cfg->setting("idserver_url");
    $idserver_url = $i if $i;

    $self->{kmer_service_url} = $cfg->setting("kmer_service_url");
    $self->{awe_server} = $cfg->setting("awe-server");
    $self->{shock_server} = $cfg->setting("shock-server");

    print STDERR "kmer_v2_data_directory = $self->{kmer_v2_data_directory}\n";
    print STDERR "idserver = $idserver_url\n";

    my $h = `hostname`;
    chomp $h;
    $self->{hostname} = $h;
    
    #END_CONSTRUCTOR

    if ($self->can('_init_instance'))
    {
	$self->_init_instance();
    }
    return $self;
}

=head1 METHODS



=head2 genome_ids_to_genomes

  $genomes = $obj->genome_ids_to_genomes($ids)

=over 4

=item Parameter and return types

=begin html

<pre>
$ids is a reference to a list where each element is a genome_id
$genomes is a reference to a list where each element is a genomeTO
genome_id is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$ids is a reference to a list where each element is a genome_id
$genomes is a reference to a list where each element is a genomeTO
genome_id is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Given one or more Central Store genome IDs, convert them into genome objects.

=back

=cut

sub genome_ids_to_genomes
{
    my $self = shift;
    my($ids) = @_;

    my @_bad_arguments;
    (ref($ids) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"ids\" (value was \"$ids\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to genome_ids_to_genomes:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'genome_ids_to_genomes');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genomes);
    #BEGIN genome_ids_to_genomes

    

    #END genome_ids_to_genomes
    my @_bad_returns;
    (ref($genomes) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"genomes\" (value was \"$genomes\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to genome_ids_to_genomes:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'genome_ids_to_genomes');
    }
    return($genomes);
}




=head2 create_genome

  $genome = $obj->create_genome($metadata)

=over 4

=item Parameter and return types

=begin html

<pre>
$metadata is a genome_metadata
$genome is a genomeTO
genome_metadata is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
genome_id is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$metadata is a genome_metadata
$genome is a genomeTO
genome_metadata is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
genome_id is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Create a new genome object and assign metadata.

=back

=cut

sub create_genome
{
    my $self = shift;
    my($metadata) = @_;

    my @_bad_arguments;
    (ref($metadata) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"metadata\" (value was \"$metadata\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to create_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create_genome');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome);
    #BEGIN create_genome

    $genome = GenomeTypeObject->new();
    $genome->set_metadata($metadata);

    #
    # If no ID was created, allocate a genome ID from our ID server.
    #
    if (!$genome->{id})
    {
	my $id_prefix = "kb";
	
	my $idc = Bio::KBase::IDServer::Client->new($idserver_url);
	$genome->{id} = "$id_prefix|g." . $idc->allocate_id_range("$id_prefix|g", 1);
    }
    $genome = $genome->prepare_for_return();

    #END create_genome
    my @_bad_returns;
    (ref($genome) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome\" (value was \"$genome\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to create_genome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create_genome');
    }
    return($genome);
}




=head2 create_genome_from_SEED

  $genome = $obj->create_genome_from_SEED($genome_id)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_id is a string
$genome is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_id is a string
$genome is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Create a new genome object based on data from the SEED project.

=back

=cut

sub create_genome_from_SEED
{
    my $self = shift;
    my($genome_id) = @_;

    my @_bad_arguments;
    (!ref($genome_id)) or push(@_bad_arguments, "Invalid type for argument \"genome_id\" (value was \"$genome_id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to create_genome_from_SEED:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create_genome_from_SEED');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome);
    #BEGIN create_genome_from_SEED
    #END create_genome_from_SEED
    my @_bad_returns;
    (ref($genome) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome\" (value was \"$genome\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to create_genome_from_SEED:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create_genome_from_SEED');
    }
    return($genome);
}




=head2 create_genome_from_RAST

  $genome = $obj->create_genome_from_RAST($genome_or_job_id)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_or_job_id is a string
$genome is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_or_job_id is a string
$genome is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Create a new genome object based on a RAST genome.

=back

=cut

sub create_genome_from_RAST
{
    my $self = shift;
    my($genome_or_job_id) = @_;

    my @_bad_arguments;
    (!ref($genome_or_job_id)) or push(@_bad_arguments, "Invalid type for argument \"genome_or_job_id\" (value was \"$genome_or_job_id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to create_genome_from_RAST:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create_genome_from_RAST');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome);
    #BEGIN create_genome_from_RAST

    print STDERR "get RAST : ctx=" . Dumper($ctx);
    $genome = {};
    #END create_genome_from_RAST
    my @_bad_returns;
    (ref($genome) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome\" (value was \"$genome\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to create_genome_from_RAST:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'create_genome_from_RAST');
    }
    return($genome);
}




=head2 set_metadata

  $genome_out = $obj->set_metadata($genome_in, $metadata)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$metadata is a genome_metadata
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
genome_metadata is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$metadata is a genome_metadata
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
genome_metadata is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string


=end text



=item Description

Modify genome metadata.

=back

=cut

sub set_metadata
{
    my $self = shift;
    my($genome_in, $metadata) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($metadata) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"metadata\" (value was \"$metadata\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to set_metadata:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'set_metadata');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN set_metadata

    $genome_out = GenomeTypeObject->initialize_without_indexes($genome_in);
    $genome_out->set_metadata($metadata);
    $genome_out = $genome_out->prepare_for_return();
    
    #end set_metadata
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to set_metadata:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'set_metadata');
    }
    return($genome_out);
}




=head2 add_contigs

  $genome_out = $obj->add_contigs($genome_in, $contigs)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$contigs is a reference to a list where each element is a contig
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$contigs is a reference to a list where each element is a contig
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Add a set of contigs to the genome object.

=back

=cut

sub add_contigs
{
    my $self = shift;
    my($genome_in, $contigs) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($contigs) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"contigs\" (value was \"$contigs\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to add_contigs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'add_contigs');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #END set_metadata
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to set_metadata:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'set_metadata');
    }
    return($genome_out);
}




=head2 add_contigs

  $genome_out = $obj->add_contigs($genome_in, $contigs)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$contigs is a reference to a list where each element is a contig
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$contigs is a reference to a list where each element is a contig
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Add a set of contigs to the genome object.

=back

=cut

sub add_contigs
{
    my $self = shift;
    my($genome_in, $contigs) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($contigs) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"contigs\" (value was \"$contigs\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to add_contigs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'add_contigs');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN add_contigs

    $genome_out = GenomeTypeObject->initialize_without_indexes($genome_in);
    $genome_out->add_contigs($contigs);
    $genome_out = $genome_out->prepare_for_return();

    #END add_contigs
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to add_contigs:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'add_contigs');
    }
    return($genome_out);
}




=head2 add_contigs_from_handle

  $genome_out = $obj->add_contigs_from_handle($genome_in, $contigs)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$contigs is a reference to a list where each element is a contig
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$contigs is a reference to a list where each element is a contig
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Add a set of contigs to the genome object, loading the contigs
from the given handle service handle.

=back

=cut

sub add_contigs_from_handle
{
    my $self = shift;
    my($genome_in, $contigs) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($contigs) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"contigs\" (value was \"$contigs\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to add_contigs_from_handle:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'add_contigs_from_handle');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN add_contigs_from_handle
    #END add_contigs_from_handle
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to add_contigs_from_handle:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'add_contigs_from_handle');
    }
    return($genome_out);
}




=head2 add_features

  $genome_out = $obj->add_features($genome_in, $features)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$features is a reference to a list where each element is a compact_feature
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
compact_feature is a reference to a list containing 5 items:
	0: (id) a string
	1: (location) a string
	2: (feature_type) a string
	3: (function) a string
	4: (aliases) a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$features is a reference to a list where each element is a compact_feature
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
compact_feature is a reference to a list containing 5 items:
	0: (id) a string
	1: (location) a string
	2: (feature_type) a string
	3: (function) a string
	4: (aliases) a string


=end text



=item Description

Add a set of features in tabular form.

=back

=cut

sub add_features
{
    my $self = shift;
    my($genome_in, $features) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($features) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"features\" (value was \"$features\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to add_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'add_features');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN add_features

    $genome_out = GenomeTypeObject->initialize($genome_in);
    $genome_out->add_features_from_list($features);
    $genome_out = $genome_out->prepare_for_return();
    
    #END add_features
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to add_features:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'add_features');
    }
    return($genome_out);
}




=head2 genomeTO_to_reconstructionTO

  $return = $obj->genomeTO_to_reconstructionTO($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a reconstructionTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
reconstructionTO is a reference to a hash where the following keys are defined:
	subsystems has a value which is a variant_subsystem_pairs
	bindings has a value which is a fid_role_pairs
	assignments has a value which is a fid_function_pairs
variant_subsystem_pairs is a reference to a list where each element is a variant_of_subsystem
variant_of_subsystem is a reference to a list containing 2 items:
	0: a subsystem
	1: a variant
subsystem is a string
variant is a string
fid_role_pairs is a reference to a list where each element is a fid_role_pair
fid_role_pair is a reference to a list containing 2 items:
	0: a fid
	1: a role
fid is a string
role is a string
fid_function_pairs is a reference to a list where each element is a fid_function_pair
fid_function_pair is a reference to a list containing 2 items:
	0: a fid
	1: a function
function is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a reconstructionTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
reconstructionTO is a reference to a hash where the following keys are defined:
	subsystems has a value which is a variant_subsystem_pairs
	bindings has a value which is a fid_role_pairs
	assignments has a value which is a fid_function_pairs
variant_subsystem_pairs is a reference to a list where each element is a variant_of_subsystem
variant_of_subsystem is a reference to a list containing 2 items:
	0: a subsystem
	1: a variant
subsystem is a string
variant is a string
fid_role_pairs is a reference to a list where each element is a fid_role_pair
fid_role_pair is a reference to a list containing 2 items:
	0: a fid
	1: a role
fid is a string
role is a string
fid_function_pairs is a reference to a list where each element is a fid_function_pair
fid_function_pair is a reference to a list containing 2 items:
	0: a fid
	1: a function
function is a string


=end text



=item Description



=back

=cut

sub genomeTO_to_reconstructionTO
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to genomeTO_to_reconstructionTO:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'genomeTO_to_reconstructionTO');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN genomeTO_to_reconstructionTO

    use ANNOserver;
    my $annoO = ANNOserver->new;
    my %in_models = map { chop; ($_ => 1) } `all_roles_used_in_models`;
    my %bindings_to_roles;

    my $features = $genomeTO->{features};
    my @role_fid_tuples;
    my $assignments = [];
    foreach my $fidH (@$features)
    {
	my $fid = $fidH->{id};
	my $f = $fidH->{function};
	if ($f)
	{
	    
	    push(@$assignments,[$fid,$f]);
	    foreach my $role (&SeedUtils::roles_of_function($f))
	    {
		push(@role_fid_tuples,[$role,$fid]);
		if ($in_models{$role}) { $bindings_to_roles{$role}->{$fid} = 1 }
	    }
	}
    }
    my $mr = $annoO->metabolic_reconstruction({-roles => \@role_fid_tuples});
    my $sub_vars = [];
    my $bindings = [];
    my %subsys;
    foreach my $tuple (@$mr)
    {
	my($sub_var,$role,$fid) = @$tuple;
	my($sub,$var) = split(/:/,$sub_var);
	if ($var !~ /\*?(0|-1)\b/)
	{
	    $subsys{$sub} = $var;
	    $bindings_to_roles{$role}->{$fid} = 1;
	}
    }
    foreach my $role (keys(%bindings_to_roles))
    {
	my $roles = $bindings_to_roles{$role};
	my @fids = keys(%$roles);
	foreach my $fid (@fids)
	{
	    push(@$bindings,[$fid,$role]);
	}
    }
    my @sv = map { [$_,$subsys{$_}] } keys(%subsys);
    $return = {
	subsystems => \@sv,
	assignments => $assignments,
	bindings   => $bindings,
    };
    
    #END genomeTO_to_reconstructionTO
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to genomeTO_to_reconstructionTO:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'genomeTO_to_reconstructionTO');
    }
    return($return);
}




=head2 genomeTO_to_feature_data

  $return = $obj->genomeTO_to_feature_data($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a fid_data_tuples
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
fid_data_tuples is a reference to a list where each element is a fid_data_tuple
fid_data_tuple is a reference to a list containing 4 items:
	0: a fid
	1: a md5
	2: a location
	3: a function
fid is a string
md5 is a string
function is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a fid_data_tuples
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
fid_data_tuples is a reference to a list where each element is a fid_data_tuple
fid_data_tuple is a reference to a list containing 4 items:
	0: a fid
	1: a md5
	2: a location
	3: a function
fid is a string
md5 is a string
function is a string


=end text



=item Description



=back

=cut

sub genomeTO_to_feature_data
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to genomeTO_to_feature_data:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'genomeTO_to_feature_data');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN genomeTO_to_feature_data

    my $feature_data = [];
    my $features = $genomeTO->{features};
    foreach my $feature (@$features)
    {
	my $fid = $feature->{id};
	my $loc = join(",",map { my($contig,$beg,$strand,$len) = @$_; 
				 "$contig\_$beg$strand$len" 
			       } @{$feature->{location}});
	my $type = $feature->{type};
	my $func = $feature->{function};
	my $md5 = "";
	$md5 = md5_hex(uc($feature->{protein_translation})) if $feature->{protein_translation};
	my $aliases = join(",",@{$feature->{aliases}});
	push(@$feature_data,[$fid,$loc,$type,$func,$aliases,$md5]);
    }
    $return = $feature_data;
    #END genomeTO_to_feature_data
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to genomeTO_to_feature_data:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'genomeTO_to_feature_data');
    }
    return($return);
}




=head2 reconstructionTO_to_roles

  $return = $obj->reconstructionTO_to_roles($reconstructionTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$reconstructionTO is a reconstructionTO
$return is a reference to a list where each element is a role
reconstructionTO is a reference to a hash where the following keys are defined:
	subsystems has a value which is a variant_subsystem_pairs
	bindings has a value which is a fid_role_pairs
	assignments has a value which is a fid_function_pairs
variant_subsystem_pairs is a reference to a list where each element is a variant_of_subsystem
variant_of_subsystem is a reference to a list containing 2 items:
	0: a subsystem
	1: a variant
subsystem is a string
variant is a string
fid_role_pairs is a reference to a list where each element is a fid_role_pair
fid_role_pair is a reference to a list containing 2 items:
	0: a fid
	1: a role
fid is a string
role is a string
fid_function_pairs is a reference to a list where each element is a fid_function_pair
fid_function_pair is a reference to a list containing 2 items:
	0: a fid
	1: a function
function is a string

</pre>

=end html

=begin text

$reconstructionTO is a reconstructionTO
$return is a reference to a list where each element is a role
reconstructionTO is a reference to a hash where the following keys are defined:
	subsystems has a value which is a variant_subsystem_pairs
	bindings has a value which is a fid_role_pairs
	assignments has a value which is a fid_function_pairs
variant_subsystem_pairs is a reference to a list where each element is a variant_of_subsystem
variant_of_subsystem is a reference to a list containing 2 items:
	0: a subsystem
	1: a variant
subsystem is a string
variant is a string
fid_role_pairs is a reference to a list where each element is a fid_role_pair
fid_role_pair is a reference to a list containing 2 items:
	0: a fid
	1: a role
fid is a string
role is a string
fid_function_pairs is a reference to a list where each element is a fid_function_pair
fid_function_pair is a reference to a list containing 2 items:
	0: a fid
	1: a function
function is a string


=end text



=item Description



=back

=cut

sub reconstructionTO_to_roles
{
    my $self = shift;
    my($reconstructionTO) = @_;

    my @_bad_arguments;
    (ref($reconstructionTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"reconstructionTO\" (value was \"$reconstructionTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to reconstructionTO_to_roles:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'reconstructionTO_to_roles');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN reconstructionTO_to_roles

    my $bindings = $reconstructionTO->{bindings};
    my %roles = map { ($_->[1] => 1) } @$bindings;
    $return = [sort keys(%roles)];

    #END reconstructionTO_to_roles
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to reconstructionTO_to_roles:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'reconstructionTO_to_roles');
    }
    return($return);
}




=head2 reconstructionTO_to_subsystems

  $return = $obj->reconstructionTO_to_subsystems($reconstructionTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$reconstructionTO is a reconstructionTO
$return is a variant_subsystem_pairs
reconstructionTO is a reference to a hash where the following keys are defined:
	subsystems has a value which is a variant_subsystem_pairs
	bindings has a value which is a fid_role_pairs
	assignments has a value which is a fid_function_pairs
variant_subsystem_pairs is a reference to a list where each element is a variant_of_subsystem
variant_of_subsystem is a reference to a list containing 2 items:
	0: a subsystem
	1: a variant
subsystem is a string
variant is a string
fid_role_pairs is a reference to a list where each element is a fid_role_pair
fid_role_pair is a reference to a list containing 2 items:
	0: a fid
	1: a role
fid is a string
role is a string
fid_function_pairs is a reference to a list where each element is a fid_function_pair
fid_function_pair is a reference to a list containing 2 items:
	0: a fid
	1: a function
function is a string

</pre>

=end html

=begin text

$reconstructionTO is a reconstructionTO
$return is a variant_subsystem_pairs
reconstructionTO is a reference to a hash where the following keys are defined:
	subsystems has a value which is a variant_subsystem_pairs
	bindings has a value which is a fid_role_pairs
	assignments has a value which is a fid_function_pairs
variant_subsystem_pairs is a reference to a list where each element is a variant_of_subsystem
variant_of_subsystem is a reference to a list containing 2 items:
	0: a subsystem
	1: a variant
subsystem is a string
variant is a string
fid_role_pairs is a reference to a list where each element is a fid_role_pair
fid_role_pair is a reference to a list containing 2 items:
	0: a fid
	1: a role
fid is a string
role is a string
fid_function_pairs is a reference to a list where each element is a fid_function_pair
fid_function_pair is a reference to a list containing 2 items:
	0: a fid
	1: a function
function is a string


=end text



=item Description



=back

=cut

sub reconstructionTO_to_subsystems
{
    my $self = shift;
    my($reconstructionTO) = @_;

    my @_bad_arguments;
    (ref($reconstructionTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"reconstructionTO\" (value was \"$reconstructionTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to reconstructionTO_to_subsystems:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'reconstructionTO_to_subsystems');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN reconstructionTO_to_subsystems

    my $subsys_pairs = $reconstructionTO->{subsystems};
    $return = $subsys_pairs;
    
    #END reconstructionTO_to_subsystems
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to reconstructionTO_to_subsystems:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'reconstructionTO_to_subsystems');
    }
    return($return);
}




=head2 assign_functions_to_CDSs

  $return = $obj->assign_functions_to_CDSs($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Given a genome object populated with contig data, perform gene calling
and functional annotation and return the annotated genome.

=back

=cut

sub assign_functions_to_CDSs
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to assign_functions_to_CDSs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'assign_functions_to_CDSs');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN assign_functions_to_CDSs
    my $features = $genomeTO->{features};
    my %to;
    my $i;
    my @prots;
    for ($i=0; ($i < @$features); $i++)
    {
	$to{$features->[$i]->{id}} = $i;
	my $fid = $features->[$i];
	my $translation;
	if (defined($translation = $fid->{protein_translation}))
	{
	    my $id = $fid->{id};
	    push(@prots,[$id,'',$translation]);
	}
    }
    my $anno = ANNOserver->new();
    my $handle = $anno->assign_function_to_prot(-input => \@prots,
						-kmer => 8,
						-scoreThreshold => 3,
						-seqHitThreshold => 3);
    while (my $res = $handle->get_next())
    {
	my($id, $function, $otu, $score, $nonoverlap_hits, $overlap_hits, $details, $fam) = @$res;
	$features->[$to{$id}]->{function} = $function;
	push(@{$features->[$to{$id}]->{annotations}},
	     ["Set function to\n$function\nby assign_function_to_CDSs with otu=$otu score=$score nonoverlap=$nonoverlap_hits hits=$overlap_hits figfam=$fam",
	      'genome annotation service',
	      time
	     ]);
    }
    $return = $genomeTO;
    #END assign_functions_to_CDSs
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to assign_functions_to_CDSs:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'assign_functions_to_CDSs');
    }
    return($return);
}




=head2 annotate_genome

  $return = $obj->annotate_genome($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub annotate_genome
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genome');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN annotate_genome
    my $genome = $genomeTO;
    my $anno = ANNOserver->new();

    if ($genome->{genetic_code} !~ /^\d+$/)
    {
	die "Genome has invalid genetic code $genome->{genetic_code}";
    }
    if ($genome->{scientific_name} eq '')
    {
	die "Genome does not have a scientific name defined";
    }
    if ($genome->{domain} eq '')
    {
	die "Genome does not have a domain defined";
    }

    #
    # Reformat the contigs for use with the ANNOserver.
    #
    my @contigs;
    foreach my $gctg (@{$genome->{contigs}})
    {
	push(@contigs, [$gctg->{id}, undef, $gctg->{dna}]);
    }

    #
    # Call genes.
    #
    print STDERR "Call genes...\n";
    my $peg_calls = $anno->call_genes(-input => \@contigs, -geneticCode => $genome->{genetic_code});
    print STDERR "Call genes...done\n";


    #
    # Call RNAs
    #
    my($genus, $species, $strain) = split(/\s+/, $genome->{scientific_name}, 3);
    $species = "sp" if $species eq '';
    print STDERR "Call rnas '$genus' '$species' '$strain' '$genome->{domain}'...\n";
    my $rna_calls = $anno->find_rnas(-input => \@contigs, -genus => $genus, -species => $species,
				     -domain => $genome->{domain});
    print STDERR "Call rnas...done\n";

    my($fasta_rna, $rna_locations) = @$rna_calls;

    my %feature_loc;
    my %feature_func;
    my %feature_anno;
    
    for my $ent (@$rna_locations)
    {
	my($loc_id, $contig, $start, $stop, $func) = @$ent;
	my $len = abs($stop - $start) + 1;
	my $strand = ($stop > $start) ? '+' : '-';
	$feature_loc{$loc_id} = [$contig, $start, $strand, $len];
	$feature_func{$loc_id} = $func if $func;
    }

    my($fasta_proteins, $protein_locations) = @$peg_calls;

    my $features = $genome->{features};
    if (!$features)
    {
	$features = [];
	$genome->{features} = $features;
    }

    #
    # Assign functions for proteins.
    #

    my $prot_fh;
    open($prot_fh, "<", \$fasta_proteins) or die "Cannot open the fasta string as a filehandle: $!";
    my $handle = $anno->assign_function_to_prot(-input => $prot_fh,
						-kmer => 8,
						-scoreThreshold => 3,
						-seqHitThreshold => 3);
    while (my $res = $handle->get_next())
    {
	my($id, $function, $otu, $score, $nonoverlap_hits, $overlap_hits, $details, $fam) = @$res;
	$feature_func{$id} = $function;
	$feature_anno{$id} = "Set function to\n$function\nby assign_function_to_prot with otu=$otu score=$score nonoverlap=$nonoverlap_hits hits=$overlap_hits figfam=$fam";
    }
    close($prot_fh);
    
    for my $ent (@$protein_locations)
    {
	my($loc_id, $contig, $start, $stop) = @$ent;
	my $len = abs($stop - $start) + 1;
	my $strand = ($stop > $start) ? '+' : '-';
	$feature_loc{$loc_id} = [$contig, $start, $strand, $len];
    }

    my $id_server = Bio::KBase::IDServer::Client->new($idserver_url);

    #
    # Create features for PEGs
    #
    my $n_pegs = @$protein_locations;
    my $protein_prefix = "$genome->{id}.CDS";
    my $peg_id_start = $id_server->allocate_id_range($protein_prefix, $n_pegs) + 0;
    print STDERR "allocated CDS id start $peg_id_start for $n_pegs CDSs\n";

    open($prot_fh, "<", \$fasta_proteins) or die "Cannot open the fasta string as a filehandle: $!";
    my $next_id = $peg_id_start;
    while (my($id, $def, $seq) = read_next_fasta_seq($prot_fh))
    {
	my $loc = $feature_loc{$id};
	my $kb_id = "$protein_prefix.$next_id";
	$next_id++;
	my $annos = [];
	push(@$annos, ['Initial gene call performed by call_genes',
		       'genome annotation service',
		       time
		       ]);
	if ($feature_anno{$id})
	{
	    push(@$annos, [$feature_anno{$id}, 'genome annotation service', time]);
	}
	my $feature = {
	    id => $kb_id,
	    location => [$loc],
	    type => 'CDS',
	    protein_translation => $seq,
	    aliases => [],
	    $feature_func{$id} ? (function => $feature_func{$id}) : (),
	    annotations => $annos,
	};
	push(@$features, $feature);
    }
    close($prot_fh);

    #
    # Create features for RNAs
    #
    my $n_rnas = @$rna_locations;
    my $rna_prefix = "$genome->{id}.rna";
    my $rna_id_start = $id_server->allocate_id_range($rna_prefix, $n_rnas) + 0;
    print STDERR "allocated id start $rna_id_start for $n_rnas nras\n";

    my $rna_fh;
    open($rna_fh, "<", \$fasta_rna) or die "Cannot open the fasta string as a filehandle: $!";
    $next_id = $rna_id_start;
    while (my($id, $def, $seq) = read_next_fasta_seq($rna_fh))
    {
	my $loc = $feature_loc{$id};
	my $kb_id = "$rna_prefix.$next_id";
	$next_id++;
	my $feature = {
	    id => $kb_id,
	    location => [$loc],
	    type => 'rna',
	    $feature_func{$id} ? (function => $feature_func{$id}) : (),
	    aliases => [],
	    annotations => [ ['Initial RNA call performed by find_rnas', 'genome annotation service', time] ],
	};
	push(@$features, $feature);
    }

    $return = $genome;
    
    #END annotate_genome
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_genome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_genome');
    }
    return($return);
}




=head2 call_selenoproteins

  $return = $obj->call_selenoproteins($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub call_selenoproteins
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_selenoproteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_selenoproteins');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_selenoproteins

    my $idc = Bio::KBase::IDServer::Client->new($idserver_url);

    my $coder = _get_coder();
    
    my $genomeTO_json = $coder->encode($genomeTO);

    my $genomeOut_json;
    my $stderr;

    my $tmp = File::Temp->new();
    print $tmp $genomeTO_json;
    close($tmp);

    my $ok = run(['rast_call_special_proteins',
		  '--seleno',
		  '--input', $tmp,
		  '--id-server', $idc->{url}],
		 '>', \$genomeOut_json,
		 '2>', \$stderr);

    undef $tmp;
    undef $genomeTO;

    if ($ok)
    {
	$return = $coder->decode($genomeOut_json);
    }
    else
    {
	die "rast_call_special_proteins failed: $stderr";
    }
    
    #END call_selenoproteins
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_selenoproteins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_selenoproteins');
    }
    return($return);
}




=head2 call_pyrrolysoproteins

  $return = $obj->call_pyrrolysoproteins($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub call_pyrrolysoproteins
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_pyrrolysoproteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_pyrrolysoproteins');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_pyrrolysoproteins
    my $idc = Bio::KBase::IDServer::Client->new($idserver_url);

    my $coder = _get_coder();
    
    my $genomeTO_json = $coder->encode($genomeTO);

    my $genomeOut_json;
    my $stderr;

    my $tmp = File::Temp->new();
    print $tmp $genomeTO_json;
    close($tmp);

    my $ok = run(['rast_call_special_proteins',
		  '--pyrro',
		  '--input', $tmp,
		  '--id-server', $idc->{url}],
		 '>', \$genomeOut_json,
		 '2>', \$stderr);

    undef $tmp;
    undef $genomeTO;

    if ($ok)
    {
	$return = $coder->decode($genomeOut_json);
    }
    else
    {
	die "rast_call_special_proteins failed: $stderr";
    }

    #END call_pyrrolysoproteins
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_pyrrolysoproteins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_pyrrolysoproteins');
    }
    return($return);
}




=head2 call_features_selenoprotein

  $return = $obj->call_features_selenoprotein($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Given a genome typed object, call selenoprotein features.

=back

=cut

sub call_features_selenoprotein
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_selenoprotein:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_selenoprotein');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_selenoprotein

    $return = $self->call_selenoproteins($genomeTO);
    
    #END call_features_selenoprotein
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_selenoprotein:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_selenoprotein');
    }
    return($return);
}




=head2 call_features_pyrrolysoprotein

  $return = $obj->call_features_pyrrolysoprotein($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Given a genome typed object, call pyrrolysoprotein features.

=back

=cut

sub call_features_pyrrolysoprotein
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_pyrrolysoprotein:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_pyrrolysoprotein');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_pyrrolysoprotein

    $return = $self->call_pyrrolysoproteins($genomeTO);
    
    #END call_features_pyrrolysoprotein
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_pyrrolysoprotein:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_pyrrolysoprotein');
    }
    return($return);
}




=head2 call_features_rRNA_SEED

  $genome_out = $obj->call_features_rRNA_SEED($genome_in, $types)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$types is a reference to a list where each element is a rna_type
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
rna_type is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$types is a reference to a list where each element is a rna_type
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
rna_type is a string


=end text



=item Description

Given a genome typed object, find instances of ribosomal RNAs in
the genome.
The types parameter is used to select the types of RNAs to
call. It is a list of strings where each value is one of
   "5S"
   "SSU"
   "LSU"
or "ALL" to choose all available rRNA types.

=back

=cut

sub call_features_rRNA_SEED
{
    my $self = shift;
    my($genome_in, $types) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($types) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"types\" (value was \"$types\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_rRNA_SEED:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_rRNA_SEED');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_rRNA_SEED
    my $coder = JSON::XS->new;
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my %types = map { lc($_) => 1 } @$types;

    my @type_args;

    if ($types{all} || @$types == 0)
    {
	# Don't need to set type arg since tool defaults to calling all.
    }
    else
    {
	for my $type (qw(5S SSU LSU))
	{
	    if ($types{lc($type)})
	    {
		push(@type_args, "-$type");
	    }
	}
    }
print STDERR Dumper(\@type_args, \%types);
    my @cmd = ("rast_call_rRNAs", "--input", $tmp_in, "--output", $tmp_out,
	       "--id-prefix", $genome_in->{id}, "--id-server", $idserver_url, @type_args);
    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "error calling rRNAs: $rc\non command @cmd";
    }

    $genome_out = $coder->decode(scalar read_file("" . $tmp_out));
    #END call_features_rRNA_SEED
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_rRNA_SEED:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_rRNA_SEED');
    }
    return($genome_out);
}




=head2 call_features_tRNA_trnascan

  $genome_out = $obj->call_features_tRNA_trnascan($genome_in)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Given a genome typed object, find instances of tRNAs in
the genome.

=back

=cut

sub call_features_tRNA_trnascan
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_tRNA_trnascan:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_tRNA_trnascan');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_tRNA_trnascan
    my $coder = JSON::XS->new;
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @cmd = ("rast_call_tRNAs", "--input", $tmp_in, "--output", $tmp_out,
	       "--id-prefix", $genome_in->{id}, "--id-server", $idserver_url);
    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "error calling tRNAs: $rc\non command @cmd";
    }

    $genome_out = $coder->decode(scalar read_file("" . $tmp_out));
    #END call_features_tRNA_trnascan
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_tRNA_trnascan:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_tRNA_trnascan');
    }
    return($genome_out);
}




=head2 call_RNAs

  $genome_out = $obj->call_RNAs($genome_in)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Given a genome typed object, find instances of all RNAs we currently
have support for detecting.

=back

=cut

sub call_RNAs
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_RNAs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_RNAs');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_RNAs
    ##################################################
    use ANNOserver;
    my $anno = ANNOserver->new();
    #
    # Reformat the contigs for use with the ANNOserver.
    #
    my @contigs;
    foreach my $gctg (@{$genome_in->{contigs}})
    {
	push(@contigs, [$gctg->{id}, undef, $gctg->{dna}]);
    }

    if ($genome_in->{scientific_name} eq '')
    {
	die "Genome does not have a scientific name defined";
    }
    if ($genome_in->{domain} eq '')
    {
	die "Genome does not have a domain defined";
    }

    #
    # Call RNAs
    #
    my($genus, $species, $strain) = split(/\s+/, $genome_in->{scientific_name}, 3);
    $species = "sp" if $species eq '';
    print STDERR "Call rnas '$genus' '$species' '$strain' '$genome_in->{domain}'...\n";
    my $rna_calls = $anno->find_rnas(-input => \@contigs, -genus => $genus, -species => $species,
				     -domain => $genome_in->{domain});
    print STDERR "Call rnas...done\n";
    my($fasta_rna, $rna_locations) = @$rna_calls;

    my %feature_loc;
    my %feature_func;
    my %feature_anno;
    
    for my $ent (@$rna_locations)
    {
	my($loc_id, $contig, $start, $stop, $func) = @$ent;
	my $len = abs($stop - $start) + 1;
	my $strand = ($stop > $start) ? '+' : '-';
	$feature_loc{$loc_id} = [$contig, $start, $strand, $len];
	$feature_func{$loc_id} = $func if $func;
    }
    my $features = $genome_in->{features};
    if (!$features)
    {
	$features = [];
	$genome_in->{features} = $features;
    }

    my $id_server = Bio::KBase::IDServer::Client->new($idserver_url);

    #
    # Create features for RNAs
    #
    my $n_rnas = @$rna_locations;
    my $rna_prefix = "$genome_in->{id}.rna";
    my $rna_id_start = $id_server->allocate_id_range($rna_prefix, $n_rnas) + 0;
    print STDERR "allocated id start $rna_id_start for $n_rnas nras\n";

    my $rna_fh;
    open($rna_fh, "<", \$fasta_rna) or die "Cannot open the fasta string as a filehandle: $!";
    my $next_id = $rna_id_start;
    while (my($id, $def, $seq) = read_next_fasta_seq($rna_fh))
    {
	my $loc = $feature_loc{$id};
	my $kb_id = "$rna_prefix.$next_id";
	$next_id++;
	my $feature = {
	    id => $kb_id,
	    location => [$loc],
	    type => 'rna',
	    $feature_func{$id} ? (function => $feature_func{$id}) : (),
	    aliases => [],
	    annotations => [ ['Initial RNA call performed by find_rnas', 'genome annotation service', time] ],
	};
	push(@$features, $feature);
    }
    $genome_out = $genome_in;

    #END call_RNAs
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_RNAs:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_RNAs');
    }
    return($genome_out);
}




=head2 call_features_CDS_glimmer3

  $return = $obj->call_features_CDS_glimmer3($genomeTO, $params)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$params is a glimmer3_parameters
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
glimmer3_parameters is a reference to a hash where the following keys are defined:
	min_training_len has a value which is an int

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$params is a glimmer3_parameters
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
glimmer3_parameters is a reference to a hash where the following keys are defined:
	min_training_len has a value which is an int


=end text



=item Description



=back

=cut

sub call_features_CDS_glimmer3
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_glimmer3:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_CDS_glimmer3');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_glimmer3

    my $genome_in = GenomeTypeObject->initialize($genomeTO);
    my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();

    my $gparams = { %$params };
    $gparams->{genetic_code} = $genome_in->{genetic_code};
    my $calls = Bio::KBase::GenomeAnnotation::Glimmer::call_genes_with_glimmer($sequences_file, $gparams);

    my $trans_table = SeedUtils::genetic_code($genome_in->{genetic_code});
    
    my $event = {
	tool_name => "glimmer3",
	execution_time => scalar gettimeofday,
	parameters => [ map { join("=", $_, $gparams->{$_}) } sort keys %$gparams ],
	hostname => $self->{hostname},
    };

    my $idc = Bio::KBase::IDServer::Client->new($idserver_url);
    my $event_id = $genome_in->add_analysis_event($event);
    my $type = 'CDS';

    my $id_prefix = $genome_in->{id};
    my $typed_prefix = join(".", $id_prefix, $type);

    my $count = @$calls;
    my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, $count);

    for my $call (@$calls)
    {
	my($fid, $contig, $begin, $end, $dna) = @$call;

	my($strand, $len);
	my $fix_start = 1;
	if ($begin < $end)
	{
	    $fix_start = 0 if $begin <= 3;

	    $strand = '+';
	    $len = $end - $begin + 1;
	}
	else
	{
	    my $cobj = $genome_in->find_contig($contig);
	    if (ref $cobj)
	    {
		my $clen = length($cobj->{dna});
		$fix_start = 0 if $begin > ($clen - 3);
	    }

	    $strand = '-';
	    $len = $begin - $end + 1;
	}

	my $loc = [[$contig, $begin, $strand, $len]];

	my $trans = SeedUtils::translate($dna, $trans_table, $fix_start);
	$trans =~ s/\*$//;

	my $id = join(".", $typed_prefix, $cur_id_suffix);
	$cur_id_suffix++;
	
	$genome_in->add_feature({
	    -id		     => $id,
	    -type 	     => $type,
	    -location 	     => $loc,
	    -analyis_event_id 	     => $event_id,
	    -annotator => 'glimmer3',
	    -protein_translation => $trans,
	});
    }
				
    $return = $genome_in;
    $return = $return->prepare_for_return();
    
    #END call_features_CDS_glimmer3
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_glimmer3:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_CDS_glimmer3');
    }
    return($return);
}




=head2 call_features_CDS_prodigal

  $return = $obj->call_features_CDS_prodigal($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub call_features_CDS_prodigal
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_prodigal:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_CDS_prodigal');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_prodigal

    my $coder = JSON::XS->new;
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genomeTO));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @cmd = ("rast_call_CDSs_using_prodigal", "--input", $tmp_in, "--output", $tmp_out,
	       "--id-prefix", $genomeTO->{id}, "--id-server", $idserver_url);
    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "error calling CDSs: $rc\non command @cmd";
    }

    $return = $coder->decode(scalar read_file("" . $tmp_out));
    
    #END call_features_CDS_prodigal
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_prodigal:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_CDS_prodigal');
    }
    return($return);
}




=head2 call_features_CDS_SEED_projection

  $return = $obj->call_features_CDS_SEED_projection($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub call_features_CDS_SEED_projection
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_SEED_projection:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_CDS_SEED_projection');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_SEED_projection
    die "Not implemented";
    #END call_features_CDS_SEED_projection
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_SEED_projection:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_CDS_SEED_projection');
    }
    return($return);
}




=head2 call_features_CDS_FragGeneScan

  $return = $obj->call_features_CDS_FragGeneScan($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub call_features_CDS_FragGeneScan
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_FragGeneScan:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_CDS_FragGeneScan');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_FragGeneScan
    die "Not implemented";
    #END call_features_CDS_FragGeneScan
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_CDS_FragGeneScan:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_CDS_FragGeneScan');
    }
    return($return);
}




=head2 call_features_repeat_region_SEED

  $genome_out = $obj->call_features_repeat_region_SEED($genome_in, $params)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$params is a repeat_region_SEED_parameters
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
repeat_region_SEED_parameters is a reference to a hash where the following keys are defined:
	min_identity has a value which is a float
	min_length has a value which is an int

</pre>

=end html

=begin text

$genome_in is a genomeTO
$params is a repeat_region_SEED_parameters
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
repeat_region_SEED_parameters is a reference to a hash where the following keys are defined:
	min_identity has a value which is a float
	min_length has a value which is an int


=end text



=item Description



=back

=cut

sub call_features_repeat_region_SEED
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_repeat_region_SEED:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_repeat_region_SEED');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_repeat_region_SEED

    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();
    my $output_file = File::Temp->new();

    my @opts;
    push(@opts, "-i", $params->{min_identity}) if defined($params->{min_identity});
    push(@opts, "-l", $params->{min_length}) if defined($params->{min_length});

    my @cmd = ("svr_big_repeats", @opts);

    my $tmpdir = File::Temp->newdir(undef, CLEANUP => 1);

    print STDERR "seq file $sequences_file\n";
    my $ok = run(\@cmd,
		 init => sub { print STDERR "init chdir to $tmpdir\n"; chdir($tmpdir) or die $!; },
		 "<", $sequences_file,
		 "|",
		 ["svr_condense_repeats"],
		 ">", $output_file,
		);
    print STDERR "ok=$ok\n" . `cat $output_file`;

    unlink($sequences_file);

    if (!$ok)
    {
	die "Error running svr_big_repeats: @cmd\n";
    }

    close($output_file);
    my($res_fh);
    open($res_fh, "<", $output_file) or die "Cannot open svr_big_repeats file $output_file: $!";

    my $event = {
	tool_name => "svr_big_repeats | svr_condense_repeats",
	execution_time => scalar gettimeofday,
	parameters => \@opts,
	hostname => $self->{hostname},
    };

    my $idc = Bio::KBase::IDServer::Client->new($idserver_url);
    my $event_id = $genome_in->add_analysis_event($event);

    # olson@bio-data-1:~/FIGdisk/dist/releases/dev2$ svr_condense_repeats < r
    #    NC_000913       15377   16741   99.71
    #    NC_000913       19796   20564   98.83

    my $type = 'repeat';
    my $function = 'repeat region';

    while(<$res_fh>)
    {
	chomp;
	my($contig, $left, $right, $iden) = split(/\t/);

	next unless $left =~ /^\d+$/ && $right =~ /^\d+$/;

	my $confidence = $iden / 100.0;

	my $quality = {
	    existence_confidence => $confidence,
	};

	my $len = 1 + $right - $left;
	my $loc = [[$contig, $left, '+', $len]];
	$genome_in->add_feature({
	    -id_client 	     => $idc,
	    -id_prefix 	     => $genome_in->{id},
	    -type 	     => $type,
	    -location 	     => $loc,
	    -function 	     => $function,
	    -analysis_event_id 	     => $event_id,
	    -quality_measure => $quality,
	});
    }

    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();

    #END call_features_repeat_region_SEED
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_repeat_region_SEED:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_repeat_region_SEED');
    }
    return($genome_out);
}




=head2 call_features_prophage_phispy

  $genome_out = $obj->call_features_prophage_phispy($genome_in)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub call_features_prophage_phispy
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_prophage_phispy:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_prophage_phispy');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_prophage_phispy
    my $coder = JSON::XS->new;
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @cmd = ("rast_call_prophage_using_phispy", "--input", $tmp_in, "--output", $tmp_out,
	       "--id-prefix", $genome_in->{id}, "--id-server", $idserver_url);
    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "error calling prophages: $rc\non command @cmd";
    }

    $genome_out = $coder->decode(scalar read_file("" . $tmp_out));
    #END call_features_prophage_phispy
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_prophage_phispy:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_prophage_phispy');
    }
    return($genome_out);
}




=head2 call_features_scan_for_matches

  $genome_out = $obj->call_features_scan_for_matches($genome_in, $pattern, $feature_type)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$pattern is a string
$feature_type is a string
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$pattern is a string
$feature_type is a string
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub call_features_scan_for_matches
{
    my $self = shift;
    my($genome_in, $pattern, $feature_type) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (!ref($pattern)) or push(@_bad_arguments, "Invalid type for argument \"pattern\" (value was \"$pattern\")");
    (!ref($feature_type)) or push(@_bad_arguments, "Invalid type for argument \"feature_type\" (value was \"$feature_type\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_scan_for_matches:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_scan_for_matches');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_scan_for_matches
    die "Not implemented";
    #END call_features_scan_for_matches
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_scan_for_matches:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_scan_for_matches');
    }
    return($genome_out);
}




=head2 annotate_proteins_kmer_v1

  $return = $obj->annotate_proteins_kmer_v1($genomeTO, $params)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$params is a kmer_v1_parameters
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$params is a kmer_v1_parameters
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int


=end text



=item Description



=back

=cut

sub annotate_proteins_kmer_v1
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_proteins_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_proteins_kmer_v1');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN annotate_proteins_kmer_v1

    my $n_proteins_per_call = 100;

    my $genome_in = GenomeTypeObject->initialize($genomeTO);

    my $kmer_service = Bio::KBase::KmerAnnotationByFigfam::Client->new($self->{kmer_service_url});


    if (!defined($params->{dataset_name}))
    {
	$params->{dataset_name} = $kmer_service->get_default_dataset_name();
    }

    my $event = {
	tool_name => "KmerAnnotationByFigfam",
	execution_time => scalar gettimeofday,
	parameters => [ map { join("=", $_, $params->{$_}) } sort keys %$params ],
	hostname => $self->{hostname},
    };

    my $event_id = $genome_in->add_analysis_event($event);

    my @proteins;

    my $do_anno = sub {
	my($proteins) = @_;
	my $res = $kmer_service->annotate_proteins($proteins, $params);
	for my $hit (@$res)
	{
	    my($id, $func, $otu, $score, $nonover, $over, $details) = @$hit;
	    if ($func)
	    {
		$genome_in->update_function("GenomeAnnotationImpl", $id, $func, $event_id);
	    }
	}
    };

    for my $feature ($genome_in->features)
    {
	my $trans = $feature->{protein_translation};
	next unless $trans;

	push(@proteins, [$feature->{id}, $trans]);
	if (@proteins >= $n_proteins_per_call)
	{
	    $do_anno->(\@proteins);
	    @proteins = ();
	}
    }
    $do_anno->(\@proteins) if @proteins;
    
    $return = $genome_in->prepare_for_return();

    #END annotate_proteins_kmer_v1
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_proteins_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_proteins_kmer_v1');
    }
    return($return);
}




=head2 annotate_proteins_kmer_v2

  $genome_out = $obj->annotate_proteins_kmer_v2($genome_in, $params)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$params is a kmer_v2_parameters
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int

</pre>

=end html

=begin text

$genome_in is a genomeTO
$params is a kmer_v2_parameters
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int


=end text



=item Description



=back

=cut

sub annotate_proteins_kmer_v2
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_proteins_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_proteins_kmer_v2');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN annotate_proteins_kmer_v2

    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $sequences_file = $genome_in->extract_protein_sequences_to_temp_file();
    my $output_file = File::Temp->new();

    my $min_hits = 5;
    my $max_gap = 200;

    if (defined($params->{min_hits}))
    {
	$min_hits = $params->{min_hits};
    }
	
    if (defined($params->{max_gap}))
    {
	$max_gap = $params->{max_gap};
    }

    my @params = ("-a", "-g", $max_gap, "-m", $min_hits, "-d", $self->{kmer_v2_data_directory});

    my @cmd = ("kmer_search", @params);
    my $ok = run(\@cmd,
		 "<", $sequences_file,
		 ">", $output_file);

    unlink($sequences_file);

    if (!$ok)
    {
	die "Error running kmer_search: @cmd\n";
    }

    close($output_file);
    my($res_fh);
    open($res_fh, "<", $output_file) or die "Cannot open kmer_search output file $output_file: $!";

    my $event = {
	tool_name => "kmer_search",
	execution_time => scalar gettimeofday,
	parameters => \@params,
	hostname => $self->{hostname},
    };

    my $event_id = $genome_in->add_analysis_event($event);
    
    while(<$res_fh>)
    {
	chomp;
	my($fid, $function) = split(/\t/);
	$genome_in->update_function("GenomeAnnotationImpl", $fid, $function, $event_id);
    }

    $genome_out = $genome_in->prepare_for_return();

    #END annotate_proteins_kmer_v2
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_proteins_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_proteins_kmer_v2');
    }
    return($genome_out);
}




=head2 resolve_overlapping_features

  $genome_out = $obj->resolve_overlapping_features($genome_in, $params)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$params is a resolve_overlapping_features_parameters
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
resolve_overlapping_features_parameters is a reference to a hash where the following keys are defined:
	placeholder has a value which is an int

</pre>

=end html

=begin text

$genome_in is a genomeTO
$params is a resolve_overlapping_features_parameters
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
resolve_overlapping_features_parameters is a reference to a hash where the following keys are defined:
	placeholder has a value which is an int


=end text



=item Description



=back

=cut

sub resolve_overlapping_features
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to resolve_overlapping_features:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'resolve_overlapping_features');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN resolve_overlapping_features

    $genome_in = GenomeTypeObject->initialize($genome_in);

    $genome_out = overlap_resolution::resolve_overlapping_features($genome_in, $params);

    $genome_out = $genome_out->prepare_for_return();

    #END resolve_overlapping_features
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to resolve_overlapping_features:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'resolve_overlapping_features');
    }
    return($genome_out);
}




=head2 call_features_ProtoCDS_kmer_v1

  $return = $obj->call_features_ProtoCDS_kmer_v1($genomeTO, $params)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$params is a kmer_v1_parameters
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$params is a kmer_v1_parameters
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int


=end text



=item Description



=back

=cut

sub call_features_ProtoCDS_kmer_v1
{
    my $self = shift;
    my($genomeTO, $params) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_ProtoCDS_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_ProtoCDS_kmer_v1');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_ProtoCDS_kmer_v1

    my $genome_in = GenomeTypeObject->initialize($genomeTO);
    
    my $min_hits = 5;
    my $max_gap = 200;

    if (defined($params->{min_hits}))
    {
	$min_hits = $params->{min_hits};
    }
	
    if (defined($params->{max_gap}))
    {
	$max_gap = $params->{max_gap};
    }

    my $event = {
	tool_name => "KmerAnnotationByFigfam",
	execution_time => scalar gettimeofday,
	parameters => [ map { join("=", $_, $params->{$_}) } sort keys %$params ],
	hostname => $self->{hostname},
    };

    my $idc = Bio::KBase::IDServer::Client->new($idserver_url);
    my $event_id = $genome_in->add_analysis_event($event);

    my $type = 'protoCDS';
    my $id_prefix = $genome_in->{id};
    my $typed_prefix = join(".", $id_prefix, $type);

    my $kmer_service = Bio::KBase::KmerAnnotationByFigfam::Client->new($self->{kmer_service_url});
    if (!defined($params->{dataset_name}))
    {
	$params->{dataset_name} = $kmer_service->get_default_dataset_name();
    }

    for my $ctg ($genome_in->contigs)
    {
	my $hits = $kmer_service->call_genes_in_dna([[$ctg->{id}, $ctg->{dna}]], $params);
	# print STDERR Dumper($hits);

	my $count = @$hits;
	my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, $count);

	for my $hit (@$hits)
	{
	    my($nhits, $id, $begin, $end, $function, $otu) = @$hit;
	    my $quality = {
		existence_confidence => 0.5,
		hit_count => 0 + $nhits,
	    };
	    my($strand, $len);
	    if ($begin < $end)
	    {
		$strand = '+';
		$len = $end - $begin + 1;
	    }
	    else
	    {
		$strand = '-';
		$len = $begin - $end + 1;
	    }
	    
	    my $fid = join(".", $typed_prefix, $cur_id_suffix);
	    $cur_id_suffix++;

	    my $loc = [[$id, $begin, $strand, $len]];
	    $genome_in->add_feature({
		-id		     => $fid,
		-type 	     => $type,
		-location 	     => $loc,
		-function 	     => $function,
		-analysis_event_id 	     => $event_id,
		-quality_measure => $quality,
	    });
	}
    }

    $return = $genome_in;
    $return = $return->prepare_for_return();

    #END call_features_ProtoCDS_kmer_v1
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_ProtoCDS_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_ProtoCDS_kmer_v1');
    }
    return($return);
}




=head2 call_features_ProtoCDS_kmer_v2

  $genome_out = $obj->call_features_ProtoCDS_kmer_v2($genome_in, $params)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$params is a kmer_v2_parameters
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int

</pre>

=end html

=begin text

$genome_in is a genomeTO
$params is a kmer_v2_parameters
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int


=end text



=item Description

RAST-style kmers

=back

=cut

sub call_features_ProtoCDS_kmer_v2
{
    my $self = shift;
    my($genome_in, $params) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($params) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"params\" (value was \"$params\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_ProtoCDS_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_ProtoCDS_kmer_v2');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_ProtoCDS_kmer_v2

    $genome_in = GenomeTypeObject->initialize($genome_in);
    
    my $sequences_file = $genome_in->extract_contig_sequences_to_temp_file();
    my $output_file = File::Temp->new();

    my $min_hits = 5;
    my $max_gap = 200;

    if (defined($params->{min_hits}))
    {
	$min_hits = $params->{min_hits};
    }
	
    if (defined($params->{max_gap}))
    {
	$max_gap = $params->{max_gap};
    }

    my @params = ("-g", $max_gap, "-m", $min_hits, "-d", $self->{kmer_v2_data_directory});

    my @cmd = ("kmer_search", @params);
    my $ok = run(\@cmd,
		 "<", $sequences_file,
		 ">", $output_file);

    unlink($sequences_file);

    if (!$ok)
    {
	die "Error running kmer_search: @cmd\n";
    }

    close($output_file);
    my($res_fh);
    open($res_fh, "<", $output_file) or die "Cannot open kmer_search output file $output_file: $!";

    my $event = {
	tool_name => "kmer_search",
	execution_time => scalar gettimeofday,
	parameters => \@params,
	hostname => $self->{hostname},
    };

    my $idc = Bio::KBase::IDServer::Client->new($idserver_url);
    my $event_id = $genome_in->add_analysis_event($event);

#	NC_002952       2902278 2902407 -       1       37      LSU ribosomal protein L34p      116.145798
#	Contig
#	Left
#	Right
#	Strand
#	Frame
#	Score-1
#	Function
#	Score-2
#    

    my $type = 'protoCDS';
    my $id_prefix = $genome_in->{id};
    my $typed_prefix = join(".", $id_prefix, $type);

    my @hits;

    while(<$res_fh>)
    {
	chomp;
	my($contig, $left, $right, $strand, $frame, $hit_count, $function, $weighted_hit_count) = split(/\t/, $_);

	next unless $left =~ /^\d+$/ && $right =~ /^\d+$/;

	push(@hits, $_);
    }
    close($res_fh);

    my $count = @hits;
    my $cur_id_suffix = $idc->allocate_id_range($typed_prefix, $count);

    for my $hit (@hits)
    {
	my($contig, $left, $right, $strand, $frame, $hit_count, $function, $weighted_hit_count) = split(/\t/, $hit);

	my $confidence = 1 - 0.5 ** ($weighted_hit_count / 3);

	my $quality = {
	    existence_confidence => $confidence,
	    hit_count => 0 + $hit_count,
	    weighted_hit_count => 0 + $weighted_hit_count,
	};

	my $id = join(".", $typed_prefix, $cur_id_suffix);
	$cur_id_suffix++;

	my $begin = 0 + (($strand eq '+') ? $left : $right);
	my $len = 1 + $right - $left;
	my $loc = [[$contig, $begin, $strand, $len]];
	$genome_in->add_feature({
	    -id              => $id,
	    -type 	     => $type,
	    -location 	     => $loc,
	    -function 	     => $function,
	    -analysis_event_id 	     => $event_id,
	    -quality_measure => $quality,
	});
    }

    $genome_out = $genome_in;
    $genome_out = $genome_out->prepare_for_return();

    #END call_features_ProtoCDS_kmer_v2
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_ProtoCDS_kmer_v2:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_ProtoCDS_kmer_v2');
    }
    return($genome_out);
}




=head2 annotate_proteins

  $return = $obj->annotate_proteins($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Ross's new kmers

=back

=cut

sub annotate_proteins
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to annotate_proteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_proteins');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN annotate_proteins
    $return = $self->assign_functions_to_CDSs($genomeTO);
    #END annotate_proteins
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to annotate_proteins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'annotate_proteins');
    }
    return($return);
}




=head2 estimate_crude_phylogenetic_position_kmer

  $position_estimate = $obj->estimate_crude_phylogenetic_position_kmer($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$position_estimate is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$position_estimate is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Determine close genomes.

=back

=cut

sub estimate_crude_phylogenetic_position_kmer
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to estimate_crude_phylogenetic_position_kmer:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'estimate_crude_phylogenetic_position_kmer');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($position_estimate);
    #BEGIN estimate_crude_phylogenetic_position_kmer
    die "Not implemented";
    #END estimate_crude_phylogenetic_position_kmer
    my @_bad_returns;
    (!ref($position_estimate)) or push(@_bad_returns, "Invalid type for return variable \"position_estimate\" (value was \"$position_estimate\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to estimate_crude_phylogenetic_position_kmer:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'estimate_crude_phylogenetic_position_kmer');
    }
    return($position_estimate);
}




=head2 find_close_neighbors

  $return = $obj->find_close_neighbors($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub find_close_neighbors
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to find_close_neighbors:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'find_close_neighbors');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN find_close_neighbors
    die "Not implemented";
    #END find_close_neighbors
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to find_close_neighbors:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'find_close_neighbors');
    }
    return($return);
}




=head2 call_features_strep_suis_repeat

  $return = $obj->call_features_strep_suis_repeat($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Interface to Strep repeats and "boxes" tools

=back

=cut

sub call_features_strep_suis_repeat
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_strep_suis_repeat:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_strep_suis_repeat');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_strep_suis_repeat

    $return = $self->_call_using_strep_repeats($genomeTO, "suis_repeat_annotation");

    #END call_features_strep_suis_repeat
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_strep_suis_repeat:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_strep_suis_repeat');
    }
    return($return);
}




=head2 call_features_strep_pneumo_repeat

  $return = $obj->call_features_strep_pneumo_repeat($genomeTO)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genomeTO is a genomeTO
$return is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub call_features_strep_pneumo_repeat
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_strep_pneumo_repeat:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_strep_pneumo_repeat');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_strep_pneumo_repeat

    $return = $self->_call_using_strep_repeats($genomeTO, "pneumococcal_repeat_annotation");

    #END call_features_strep_pneumo_repeat
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_strep_pneumo_repeat:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_strep_pneumo_repeat');
    }
    return($return);
}




=head2 call_features_crispr

  $genome_out = $obj->call_features_crispr($genome_in)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description



=back

=cut

sub call_features_crispr
{
    my $self = shift;
    my($genome_in) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_crispr:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_crispr');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN call_features_crispr

    my $coder = JSON::XS->new;
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @cmd = ("rast_call_crisprs", "--input", $tmp_in, "--output", $tmp_out,
	       "--id-prefix", $genome_in->{id}, "--id-server", $idserver_url);
    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "error calling rRNAs: $rc\non command @cmd";
    }

    $genome_out = $coder->decode(scalar read_file("" . $tmp_out));

    #END call_features_crispr
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_crispr:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_crispr');
    }
    return($genome_out);
}




=head2 export_genome

  $exported_data = $obj->export_genome($genome_in, $format, $feature_types)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$format is a string
$feature_types is a reference to a list where each element is a string
$exported_data is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string

</pre>

=end html

=begin text

$genome_in is a genomeTO
$format is a string
$feature_types is a reference to a list where each element is a string
$exported_data is a string
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string


=end text



=item Description

Export genome typed object to one of the supported output formats:
genbank, embl, or gff.
If feature_types is a non-empty list, limit the output to the given
feature types.

=back

=cut

sub export_genome
{
    my $self = shift;
    my($genome_in, $format, $feature_types) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (!ref($format)) or push(@_bad_arguments, "Invalid type for argument \"format\" (value was \"$format\")");
    (ref($feature_types) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"feature_types\" (value was \"$feature_types\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to export_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'export_genome');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($exported_data);
    #BEGIN export_genome

    my $coder = JSON::XS->new;
    my $tmp_in = File::Temp->new();
    write_file($tmp_in, $coder->encode($genome_in));

    close($tmp_in);
    my $tmp_out = File::Temp->new();
    close($tmp_out);

    my @type_flag = map { ("--feature-type", $_) } @$feature_types;

    my @cmd = ("rast_export_genome", @type_flag, "--input", $tmp_in, "--output", $tmp_out, $format);

    my $rc = system(@cmd);
    if ($rc != 0)
    {
	die "error exporting genome: $rc\non command @cmd";
    }

    $exported_data = scalar read_file("" . $tmp_out);

    #END export_genome
    my @_bad_returns;
    (!ref($exported_data)) or push(@_bad_returns, "Invalid type for return variable \"exported_data\" (value was \"$exported_data\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to export_genome:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'export_genome');
    }
    return($exported_data);
}




=head2 enumerate_classifiers

  $return = $obj->enumerate_classifiers()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a reference to a list where each element is a string

</pre>

=end html

=begin text

$return is a reference to a list where each element is a string


=end text



=item Description

Enumerate the available classifiers. Returns the list of identifiers for
the classifiers.

=back

=cut

sub enumerate_classifiers
{
    my $self = shift;

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN enumerate_classifiers

    $return = [];
    my $dir = $self->{kmer_classifier_data_directory};
    for my $ent (<$dir/*/groups>)
    {
	print STDERR "Try $ent\n";
	my($name) = $ent =~ m,([^/]+)/groups$,;
	push(@$return, $name);
    }
    #END enumerate_classifiers
    my @_bad_returns;
    (ref($return) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to enumerate_classifiers:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'enumerate_classifiers');
    }
    return($return);
}




=head2 query_classifier_groups

  $return = $obj->query_classifier_groups($classifier)

=over 4

=item Parameter and return types

=begin html

<pre>
$classifier is a string
$return is a reference to a hash where the key is a string and the value is a reference to a list where each element is a genome_id
genome_id is a string

</pre>

=end html

=begin text

$classifier is a string
$return is a reference to a hash where the key is a string and the value is a reference to a list where each element is a genome_id
genome_id is a string


=end text



=item Description

Query the groups included in the given classifier. This is a
mapping from the group name to the list of genome IDs included
in the group. Note that these are genome IDs native to the
system that created the classifier; currently these are
SEED genome IDs that may be translated using the
source IDs on the Genome entity.

=back

=cut

sub query_classifier_groups
{
    my $self = shift;
    my($classifier) = @_;

    my @_bad_arguments;
    (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument \"classifier\" (value was \"$classifier\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to query_classifier_groups:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'query_classifier_groups');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN query_classifier_groups

    my $cobj = KmerClassifier->new("$self->{kmer_classifier_data_directory}/$classifier");
    $return = $cobj->group_membership_hash();

    #END query_classifier_groups
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to query_classifier_groups:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'query_classifier_groups');
    }
    return($return);
}




=head2 query_classifier_taxonomies

  $return = $obj->query_classifier_taxonomies($classifier)

=over 4

=item Parameter and return types

=begin html

<pre>
$classifier is a string
$return is a reference to a hash where the key is a string and the value is a string

</pre>

=end html

=begin text

$classifier is a string
$return is a reference to a hash where the key is a string and the value is a string


=end text



=item Description

Query the taxonomy strings that this classifier maps.

=back

=cut

sub query_classifier_taxonomies
{
    my $self = shift;
    my($classifier) = @_;

    my @_bad_arguments;
    (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument \"classifier\" (value was \"$classifier\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to query_classifier_taxonomies:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'query_classifier_taxonomies');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN query_classifier_taxonomies
    #END query_classifier_taxonomies
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to query_classifier_taxonomies:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'query_classifier_taxonomies');
    }
    return($return);
}




=head2 classify_into_bins

  $return = $obj->classify_into_bins($classifier, $dna_input)

=over 4

=item Parameter and return types

=begin html

<pre>
$classifier is a string
$dna_input is a reference to a list where each element is a reference to a list containing 2 items:
	0: (id) a string
	1: (dna_data) a string
$return is a reference to a hash where the key is a string and the value is an int

</pre>

=end html

=begin text

$classifier is a string
$dna_input is a reference to a list where each element is a reference to a list containing 2 items:
	0: (id) a string
	1: (dna_data) a string
$return is a reference to a hash where the key is a string and the value is an int


=end text



=item Description

Classify a dataset, returning only the binned output.

=back

=cut

sub classify_into_bins
{
    my $self = shift;
    my($classifier, $dna_input) = @_;

    my @_bad_arguments;
    (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument \"classifier\" (value was \"$classifier\")");
    (ref($dna_input) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"dna_input\" (value was \"$dna_input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to classify_into_bins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'classify_into_bins');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN classify_into_bins

    my $cobj = KmerClassifier->new("$self->{kmer_classifier_data_directory}/$classifier");

    my $tmp = File::Temp->new;
    close($tmp);
    open(TMP, ">", $tmp) or die "cannot write tempfile $tmp: $!";
    for my $ent (@$dna_input)
    {
	gjoseqlib::print_alignment_as_fasta(\*TMP, [$ent->[0], undef, $ent->[1]]);
    }
    close(TMP);
    my($bins, $missed) = $cobj->classify("" . $tmp);

    $return = $bins;

    #END classify_into_bins
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to classify_into_bins:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'classify_into_bins');
    }
    return($return);
}




=head2 classify_full

  $return_1, $raw_output, $unassigned = $obj->classify_full($classifier, $dna_input)

=over 4

=item Parameter and return types

=begin html

<pre>
$classifier is a string
$dna_input is a reference to a list where each element is a reference to a list containing 2 items:
	0: (id) a string
	1: (dna_data) a string
$return_1 is a reference to a hash where the key is a string and the value is an int
$raw_output is a string
$unassigned is a reference to a list where each element is a string

</pre>

=end html

=begin text

$classifier is a string
$dna_input is a reference to a list where each element is a reference to a list containing 2 items:
	0: (id) a string
	1: (dna_data) a string
$return_1 is a reference to a hash where the key is a string and the value is an int
$raw_output is a string
$unassigned is a reference to a list where each element is a string


=end text



=item Description

Classify a dataset, returning the binned output along with the raw assignments and the list of
sequences that were not assigned.

=back

=cut

sub classify_full
{
    my $self = shift;
    my($classifier, $dna_input) = @_;

    my @_bad_arguments;
    (!ref($classifier)) or push(@_bad_arguments, "Invalid type for argument \"classifier\" (value was \"$classifier\")");
    (ref($dna_input) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"dna_input\" (value was \"$dna_input\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to classify_full:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'classify_full');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return_1, $raw_output, $unassigned);
    #BEGIN classify_full

    my $cobj = KmerClassifier->new("$self->{kmer_classifier_data_directory}/$classifier");

    my $tmp = File::Temp->new;
    close($tmp);
    open(TMP, ">", $tmp) or die "cannot write tempfile $tmp: $!";

    my $raw_tmp = File::Temp->new();
    for my $ent (@$dna_input)
    {
	gjoseqlib::print_alignment_as_fasta(\*TMP, [$ent->[0], undef, $ent->[1]]);
    }
    close(TMP);
    my($bins, $missed) = $cobj->classify("" . $tmp, $raw_tmp);
    close($raw_tmp);

    $return_1 = $bins;
    $unassigned = $missed;
    $raw_output = read_file("" . $raw_tmp);

    #END classify_full
    my @_bad_returns;
    (ref($return_1) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return_1\" (value was \"$return_1\")");
    (!ref($raw_output)) or push(@_bad_returns, "Invalid type for return variable \"raw_output\" (value was \"$raw_output\")");
    (ref($unassigned) eq 'ARRAY') or push(@_bad_returns, "Invalid type for return variable \"unassigned\" (value was \"$unassigned\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to classify_full:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'classify_full');
    }
    return($return_1, $raw_output, $unassigned);
}




=head2 default_workflow

  $return = $obj->default_workflow()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a workflow
workflow is a reference to a hash where the following keys are defined:
	stages has a value which is a reference to a list where each element is a pipeline_stage
pipeline_stage is a reference to a hash where the following keys are defined:
	name has a value which is a string
	condition has a value which is a string
	repeat_region_SEED_parameters has a value which is a repeat_region_SEED_parameters
	glimmer3_parameters has a value which is a glimmer3_parameters
	kmer_v1_parameters has a value which is a kmer_v1_parameters
	kmer_v2_parameters has a value which is a kmer_v2_parameters
repeat_region_SEED_parameters is a reference to a hash where the following keys are defined:
	min_identity has a value which is a float
	min_length has a value which is an int
glimmer3_parameters is a reference to a hash where the following keys are defined:
	min_training_len has a value which is an int
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int

</pre>

=end html

=begin text

$return is a workflow
workflow is a reference to a hash where the following keys are defined:
	stages has a value which is a reference to a list where each element is a pipeline_stage
pipeline_stage is a reference to a hash where the following keys are defined:
	name has a value which is a string
	condition has a value which is a string
	repeat_region_SEED_parameters has a value which is a repeat_region_SEED_parameters
	glimmer3_parameters has a value which is a glimmer3_parameters
	kmer_v1_parameters has a value which is a kmer_v1_parameters
	kmer_v2_parameters has a value which is a kmer_v2_parameters
repeat_region_SEED_parameters is a reference to a hash where the following keys are defined:
	min_identity has a value which is a float
	min_length has a value which is an int
glimmer3_parameters is a reference to a hash where the following keys are defined:
	min_training_len has a value which is an int
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int


=end text



=item Description



=back

=cut

sub default_workflow
{
    my $self = shift;

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN default_workflow

    my @stages = (
	      { name => 'call_features_rRNA_SEED' },
	      { name => 'call_features_tRNA_trnascan' },
	      { name => 'call_features_repeat_region_SEED',
		    repeat_region_SEED_parameters => { } },
	      { name => 'call_selenoproteins' },
	      { name => 'call_pyrrolysoproteins' },
	      { name => 'call_features_strep_suis_repeat',
		condition => '$genome->{scientific_name} =~ /^Streptococcus\s/' },
	      { name => 'call_features_strep_pneumo_repeat',
		condition => '$genome->{scientific_name} =~ /^Streptococcus\s/' },
	      { name => 'call_features_crispr' },
	      { name => 'call_features_CDS_prodigal' },
	      { name => 'annotate_proteins_kmer_v2', kmer_v2_parameters => {} },
	      # { name => 'call_features_prophage_phispy' },
		 );
    $return = { stages => \@stages };
    
    #END default_workflow
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to default_workflow:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'default_workflow');
    }
    return($return);
}




=head2 run_pipeline

  $genome_out = $obj->run_pipeline($genome_in, $workflow)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$workflow is a workflow
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
workflow is a reference to a hash where the following keys are defined:
	stages has a value which is a reference to a list where each element is a pipeline_stage
pipeline_stage is a reference to a hash where the following keys are defined:
	name has a value which is a string
	condition has a value which is a string
	repeat_region_SEED_parameters has a value which is a repeat_region_SEED_parameters
	glimmer3_parameters has a value which is a glimmer3_parameters
	kmer_v1_parameters has a value which is a kmer_v1_parameters
	kmer_v2_parameters has a value which is a kmer_v2_parameters
repeat_region_SEED_parameters is a reference to a hash where the following keys are defined:
	min_identity has a value which is a float
	min_length has a value which is an int
glimmer3_parameters is a reference to a hash where the following keys are defined:
	min_training_len has a value which is an int
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int

</pre>

=end html

=begin text

$genome_in is a genomeTO
$workflow is a workflow
$genome_out is a genomeTO
genomeTO is a reference to a hash where the following keys are defined:
	id has a value which is a genome_id
	scientific_name has a value which is a string
	domain has a value which is a string
	genetic_code has a value which is an int
	source has a value which is a string
	source_id has a value which is a string
	quality has a value which is a genome_quality_measure
	contigs has a value which is a reference to a list where each element is a contig
	contigs_handle has a value which is a Handle
	features has a value which is a reference to a list where each element is a feature
	close_genomes has a value which is a reference to a list where each element is a close_genome
	analysis_events has a value which is a reference to a list where each element is an analysis_event
genome_id is a string
genome_quality_measure is a reference to a hash where the following keys are defined:
	frameshift_error_rate has a value which is a float
	sequence_error_rate has a value which is a float
contig is a reference to a hash where the following keys are defined:
	id has a value which is a contig_id
	dna has a value which is a string
	genetic_code has a value which is an int
	cell_compartment has a value which is a string
	replicon_type has a value which is a string
	replicon_geometry has a value which is a string
	complete has a value which is a bool
contig_id is a string
bool is an int
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
feature is a reference to a hash where the following keys are defined:
	id has a value which is a feature_id
	location has a value which is a location
	type has a value which is a feature_type
	function has a value which is a string
	protein_translation has a value which is a string
	aliases has a value which is a reference to a list where each element is a string
	alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
	0: (source) a string
	1: (alias) a string

	annotations has a value which is a reference to a list where each element is an annotation
	quality has a value which is a feature_quality_measure
	feature_creation_event has a value which is an analysis_event_id
feature_id is a string
location is a reference to a list where each element is a region_of_dna
region_of_dna is a reference to a list containing 4 items:
	0: a contig_id
	1: (begin) an int
	2: (strand) a string
	3: (length) an int
feature_type is a string
annotation is a reference to a list containing 4 items:
	0: (comment) a string
	1: (annotator) a string
	2: (annotation_time) an int
	3: an analysis_event_id
analysis_event_id is a string
feature_quality_measure is a reference to a hash where the following keys are defined:
	truncated_begin has a value which is a bool
	truncated_end has a value which is a bool
	existence_confidence has a value which is a float
	frameshifted has a value which is a bool
	selenoprotein has a value which is a bool
	pyrrolysylprotein has a value which is a bool
	overlap_rules has a value which is a reference to a list where each element is a string
	existence_priority has a value which is a float
	hit_count has a value which is a float
	weighted_hit_count has a value which is a float
close_genome is a reference to a hash where the following keys are defined:
	genome has a value which is a genome_id
	closeness_measure has a value which is a float
analysis_event is a reference to a hash where the following keys are defined:
	id has a value which is an analysis_event_id
	tool_name has a value which is a string
	execution_time has a value which is a float
	parameters has a value which is a reference to a list where each element is a string
	hostname has a value which is a string
workflow is a reference to a hash where the following keys are defined:
	stages has a value which is a reference to a list where each element is a pipeline_stage
pipeline_stage is a reference to a hash where the following keys are defined:
	name has a value which is a string
	condition has a value which is a string
	repeat_region_SEED_parameters has a value which is a repeat_region_SEED_parameters
	glimmer3_parameters has a value which is a glimmer3_parameters
	kmer_v1_parameters has a value which is a kmer_v1_parameters
	kmer_v2_parameters has a value which is a kmer_v2_parameters
repeat_region_SEED_parameters is a reference to a hash where the following keys are defined:
	min_identity has a value which is a float
	min_length has a value which is an int
glimmer3_parameters is a reference to a hash where the following keys are defined:
	min_training_len has a value which is an int
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int


=end text



=item Description



=back

=cut

sub run_pipeline
{
    my $self = shift;
    my($genome_in, $workflow) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (ref($workflow) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"workflow\" (value was \"$workflow\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to run_pipeline:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'run_pipeline');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_out);
    #BEGIN run_pipeline

    my %param_defs = (annotate_proteins_kmer_v1 => 'kmer_v1_parameters',
		      annotate_proteins_kmer_v2 => 'kmer_v2_parameters',
		      call_features_repeat_region_SEED => 'repeat_regions_SEED_parameters',
		      call_features_CDS_glimmer3 => 'glimmer3_parameters',
		      call_features_repeat_region_SEED => 'repeat_region_SEED_parameters',
		      );

    my $cur = $genome_in;
    for my $stage (@{$workflow->{stages}})
    {
	my $method = $stage->{name};
	my $condition = $stage->{condition};
	if ($condition)
	{
	    my $safe = Safe->new();
	    my $g = $safe->varglob('genome');
	    $$g = $cur;
	    my $ok = $safe->reval($condition);
	    print STDERR "Condition eval of '$condition' returns $ok\n";
	    if (!$ok)
	    {
		print STDERR "Skipping $method due to condition $condition\n";
		next;
	    }
	}
	my @params;
	if (my $param_def = $param_defs{$method})
	{
	    push(@params, $stage->{$param_def});
	}
	elsif ($method eq 'call_features_rRNA_SEED')
	{
	    # Special case.
	    push(@params, []);
	}
	if ($self->can($method))
	{
	    print STDERR "Call $method with @params\n";
	    print STDERR Dumper($stage);
	    my $out;
	    eval {
		$out = $self->$method($cur, @params);
	    };

	    if ($@)
	    {
		die "Error invoking method $method: $@";
	    }
	    else
	    {
		print STDERR "Finished\n";
		$cur = $out;
	    }
	}
	else
	{
	    die "Trying to call invalid method $method";
	}
    }

    $genome_out = $cur;

    #END run_pipeline
    my @_bad_returns;
    (ref($genome_out) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_out\" (value was \"$genome_out\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to run_pipeline:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'run_pipeline');
    }
    return($genome_out);
}




=head2 pipeline_batch_start

  $batch_id = $obj->pipeline_batch_start($genomes, $workflow)

=over 4

=item Parameter and return types

=begin html

<pre>
$genomes is a reference to a list where each element is a pipeline_batch_input
$workflow is a workflow
$batch_id is a string
pipeline_batch_input is a reference to a hash where the following keys are defined:
	genome_id has a value which is a string
	data has a value which is a Handle
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
workflow is a reference to a hash where the following keys are defined:
	stages has a value which is a reference to a list where each element is a pipeline_stage
pipeline_stage is a reference to a hash where the following keys are defined:
	name has a value which is a string
	condition has a value which is a string
	repeat_region_SEED_parameters has a value which is a repeat_region_SEED_parameters
	glimmer3_parameters has a value which is a glimmer3_parameters
	kmer_v1_parameters has a value which is a kmer_v1_parameters
	kmer_v2_parameters has a value which is a kmer_v2_parameters
repeat_region_SEED_parameters is a reference to a hash where the following keys are defined:
	min_identity has a value which is a float
	min_length has a value which is an int
glimmer3_parameters is a reference to a hash where the following keys are defined:
	min_training_len has a value which is an int
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int

</pre>

=end html

=begin text

$genomes is a reference to a list where each element is a pipeline_batch_input
$workflow is a workflow
$batch_id is a string
pipeline_batch_input is a reference to a hash where the following keys are defined:
	genome_id has a value which is a string
	data has a value which is a Handle
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string
workflow is a reference to a hash where the following keys are defined:
	stages has a value which is a reference to a list where each element is a pipeline_stage
pipeline_stage is a reference to a hash where the following keys are defined:
	name has a value which is a string
	condition has a value which is a string
	repeat_region_SEED_parameters has a value which is a repeat_region_SEED_parameters
	glimmer3_parameters has a value which is a glimmer3_parameters
	kmer_v1_parameters has a value which is a kmer_v1_parameters
	kmer_v2_parameters has a value which is a kmer_v2_parameters
repeat_region_SEED_parameters is a reference to a hash where the following keys are defined:
	min_identity has a value which is a float
	min_length has a value which is an int
glimmer3_parameters is a reference to a hash where the following keys are defined:
	min_training_len has a value which is an int
kmer_v1_parameters is a reference to a hash where the following keys are defined:
	kmer_size has a value which is an int
	dataset_name has a value which is a string
	return_scores_for_all_proteins has a value which is an int
	score_threshold has a value which is an int
	hit_threshold has a value which is an int
	sequential_hit_threshold has a value which is an int
	detailed has a value which is an int
	min_hits has a value which is an int
	min_size has a value which is an int
	max_gap has a value which is an int
kmer_v2_parameters is a reference to a hash where the following keys are defined:
	min_hits has a value which is an int
	max_gap has a value which is an int


=end text



=item Description



=back

=cut

sub pipeline_batch_start
{
    my $self = shift;
    my($genomes, $workflow) = @_;

    my @_bad_arguments;
    (ref($genomes) eq 'ARRAY') or push(@_bad_arguments, "Invalid type for argument \"genomes\" (value was \"$genomes\")");
    (ref($workflow) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"workflow\" (value was \"$workflow\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to pipeline_batch_start:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'pipeline_batch_start');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($batch_id);
    #BEGIN pipeline_batch_start

    #
    # Construct an AWE workflow and submit it to our AWE server.
    #
    # The genomes have already been uploaded - the genomes list here is a
    # list of handle objects.
    #
    # The workflow we create will have one task per genome submitted. The data
    # passed to each task are the serialized form of the handle and the
    # workflow description.
    #

    my $json = JSON::XS->new->pretty(1);
    my $shock = Bio::KBase::GenomeAnnotation::Shock->new($self->{shock_server}, $ctx->token);
    my $awe = Bio::KBase::GenomeAnnotation::Awe->new($self->{awe_server}, $ctx->token);
    
    my $job = Bio::KBase::GenomeAnnotation::Awe::JobDescription->new(pipeline => 'rasttk',
								     name => 'rasttk',
								     project => 'rasttk',
								     user => $ctx->user_id,
								     clientgroups => '');
    
    my $i = 0;
    for my $gspec (@$genomes)
    {
	my($genome_id, $handle) = @$gspec{'genome_id', 'data'};
	
	my $txt = $json->encode([$handle, $workflow]);
	my $node = $shock->put_file_data($txt);
	
	my $in_file = Bio::KBase::GenomeAnnotation::Awe::JobFile->new("pipeinput_$i.json", $shock->server, $node);
	my $out_file = Bio::KBase::GenomeAnnotation::Awe::JobFile->new("pipeoutput_$i.json", $shock->server);
	my $stdout_file = Bio::KBase::GenomeAnnotation::Awe::JobFile->new("stdout_$i.txt", $shock->server);
	my $stderr_file = Bio::KBase::GenomeAnnotation::Awe::JobFile->new("stderr_$i.txt", $shock->server);
	$i++;
	
	my $id = $job->add_task("rast_run_pipeline_local",
				"rast_run_pipeline_local",
				join(" ",
				     $in_file->in_name,
				     $out_file->name,
				     $stdout_file->name,
				     $stderr_file->name),
				[],
				[$in_file],
				[$out_file, $stdout_file, $stderr_file],
				undef, undef, $awe,
			        { genome_id => $genome_id });
    }
    
    $batch_id = $awe->submit($job);
    print STDERR "submitted $batch_id\n";
    
    #END pipeline_batch_start
    my @_bad_returns;
    (!ref($batch_id)) or push(@_bad_returns, "Invalid type for return variable \"batch_id\" (value was \"$batch_id\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to pipeline_batch_start:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'pipeline_batch_start');
    }
    return($batch_id);
}




=head2 pipeline_batch_status

  $genome_status = $obj->pipeline_batch_status($batch_id)

=over 4

=item Parameter and return types

=begin html

<pre>
$batch_id is a string
$genome_status is a pipeline_batch_status_entry
pipeline_batch_status_entry is a reference to a hash where the following keys are defined:
	genome_id has a value which is a string
	status has a value which is a string
	stdout has a value which is a Handle
	stderr has a value which is a Handle
	output has a value which is a Handle
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string

</pre>

=end html

=begin text

$batch_id is a string
$genome_status is a pipeline_batch_status_entry
pipeline_batch_status_entry is a reference to a hash where the following keys are defined:
	genome_id has a value which is a string
	status has a value which is a string
	stdout has a value which is a Handle
	stderr has a value which is a Handle
	output has a value which is a Handle
Handle is a reference to a hash where the following keys are defined:
	file_name has a value which is a string
	id has a value which is a string
	type has a value which is a string
	url has a value which is a string
	remote_md5 has a value which is a string
	remote_sha1 has a value which is a string


=end text



=item Description



=back

=cut

sub pipeline_batch_status
{
    my $self = shift;
    my($batch_id) = @_;

    my @_bad_arguments;
    (!ref($batch_id)) or push(@_bad_arguments, "Invalid type for argument \"batch_id\" (value was \"$batch_id\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to pipeline_batch_status:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'pipeline_batch_status');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($genome_status);
    #BEGIN pipeline_batch_status

    my $awe = Bio::KBase::GenomeAnnotation::Awe->new($self->{awe_server}, $ctx->token);
    
    my $job = $awe->job($batch_id);
    print STDERR Dumper($job);
    
    #END pipeline_batch_status
    my @_bad_returns;
    (ref($genome_status) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"genome_status\" (value was \"$genome_status\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to pipeline_batch_status:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'pipeline_batch_status');
    }
    return($genome_status);
}




=head2 version 

  $return = $obj->version()

=over 4

=item Parameter and return types

=begin html

<pre>
$return is a string
</pre>

=end html

=begin text

$return is a string

=end text

=item Description

Return the module version. This is a Semantic Versioning number.

=back

=cut

sub version {
    return $VERSION;
}

=head1 TYPES



=head2 Handle

=over 4



=item Description

* This is a handle service handle object, used for by-reference
* passing of data files.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
file_name has a value which is a string
id has a value which is a string
type has a value which is a string
url has a value which is a string
remote_md5 has a value which is a string
remote_sha1 has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
file_name has a value which is a string
id has a value which is a string
type has a value which is a string
url has a value which is a string
remote_md5 has a value which is a string
remote_sha1 has a value which is a string


=end text

=back



=head2 bool

=over 4



=item Definition

=begin html

<pre>
an int
</pre>

=end html

=begin text

an int

=end text

=back



=head2 md5

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 md5s

=over 4



=item Definition

=begin html

<pre>
a reference to a list where each element is a md5
</pre>

=end html

=begin text

a reference to a list where each element is a md5

=end text

=back



=head2 genome_id

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 feature_id

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 contig_id

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 feature_type

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 region_of_dna

=over 4



=item Description

A region of DNA is maintained as a tuple of four components:

                the contig
                the beginning position (from 1)
                the strand
                the length

           We often speak of "a region".  By "location", we mean a sequence
           of regions from the same genome (perhaps from distinct contigs).

           Strand is either '+' or '-'.


=item Definition

=begin html

<pre>
a reference to a list containing 4 items:
0: a contig_id
1: (begin) an int
2: (strand) a string
3: (length) an int

</pre>

=end html

=begin text

a reference to a list containing 4 items:
0: a contig_id
1: (begin) an int
2: (strand) a string
3: (length) an int


=end text

=back



=head2 location

=over 4



=item Description

a "location" refers to a sequence of regions


=item Definition

=begin html

<pre>
a reference to a list where each element is a region_of_dna
</pre>

=end html

=begin text

a reference to a list where each element is a region_of_dna

=end text

=back



=head2 analysis_event_id

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 analysis_event

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is an analysis_event_id
tool_name has a value which is a string
execution_time has a value which is a float
parameters has a value which is a reference to a list where each element is a string
hostname has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is an analysis_event_id
tool_name has a value which is a string
execution_time has a value which is a float
parameters has a value which is a reference to a list where each element is a string
hostname has a value which is a string


=end text

=back



=head2 annotation

=over 4



=item Definition

=begin html

<pre>
a reference to a list containing 4 items:
0: (comment) a string
1: (annotator) a string
2: (annotation_time) an int
3: an analysis_event_id

</pre>

=end html

=begin text

a reference to a list containing 4 items:
0: (comment) a string
1: (annotator) a string
2: (annotation_time) an int
3: an analysis_event_id


=end text

=back



=head2 feature_quality_measure

=over 4



=item Description

* The numeric priority of this feature's right to exist. Specialty
* tools will give the features they create a high priority; more generic
* tools will give their features a lower priority. The overlap removal procedure
* will use this priority to determine which of a set of overlapping features
* should be removed.
*
* The intent is that a change of 1 in the priority value represents a factor of 2 in
* preference.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
truncated_begin has a value which is a bool
truncated_end has a value which is a bool
existence_confidence has a value which is a float
frameshifted has a value which is a bool
selenoprotein has a value which is a bool
pyrrolysylprotein has a value which is a bool
overlap_rules has a value which is a reference to a list where each element is a string
existence_priority has a value which is a float
hit_count has a value which is a float
weighted_hit_count has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
truncated_begin has a value which is a bool
truncated_end has a value which is a bool
existence_confidence has a value which is a float
frameshifted has a value which is a bool
selenoprotein has a value which is a bool
pyrrolysylprotein has a value which is a bool
overlap_rules has a value which is a reference to a list where each element is a string
existence_priority has a value which is a float
hit_count has a value which is a float
weighted_hit_count has a value which is a float


=end text

=back



=head2 feature

=over 4



=item Description

A feature object represents a feature on the genome. It contains 
the location on the contig with a type, the translation if it
represents a protein, associated aliases, etc. It also contains
information gathered during the annotation process that is involved
in stages that perform overlap removal, quality testing, etc.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a feature_id
location has a value which is a location
type has a value which is a feature_type
function has a value which is a string
protein_translation has a value which is a string
aliases has a value which is a reference to a list where each element is a string
alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
0: (source) a string
1: (alias) a string

annotations has a value which is a reference to a list where each element is an annotation
quality has a value which is a feature_quality_measure
feature_creation_event has a value which is an analysis_event_id

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a feature_id
location has a value which is a location
type has a value which is a feature_type
function has a value which is a string
protein_translation has a value which is a string
aliases has a value which is a reference to a list where each element is a string
alias_pairs has a value which is a reference to a list where each element is a reference to a list containing 2 items:
0: (source) a string
1: (alias) a string

annotations has a value which is a reference to a list where each element is an annotation
quality has a value which is a feature_quality_measure
feature_creation_event has a value which is an analysis_event_id


=end text

=back



=head2 contig

=over 4



=item Description

circular / linear


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a contig_id
dna has a value which is a string
genetic_code has a value which is an int
cell_compartment has a value which is a string
replicon_type has a value which is a string
replicon_geometry has a value which is a string
complete has a value which is a bool

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a contig_id
dna has a value which is a string
genetic_code has a value which is an int
cell_compartment has a value which is a string
replicon_type has a value which is a string
replicon_geometry has a value which is a string
complete has a value which is a bool


=end text

=back



=head2 close_genome

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome has a value which is a genome_id
closeness_measure has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome has a value which is a genome_id
closeness_measure has a value which is a float


=end text

=back



=head2 genome_quality_measure

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
frameshift_error_rate has a value which is a float
sequence_error_rate has a value which is a float

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
frameshift_error_rate has a value which is a float
sequence_error_rate has a value which is a float


=end text

=back



=head2 genomeTO

=over 4



=item Description

All of the information about particular genome


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a genome_id
scientific_name has a value which is a string
domain has a value which is a string
genetic_code has a value which is an int
source has a value which is a string
source_id has a value which is a string
quality has a value which is a genome_quality_measure
contigs has a value which is a reference to a list where each element is a contig
contigs_handle has a value which is a Handle
features has a value which is a reference to a list where each element is a feature
close_genomes has a value which is a reference to a list where each element is a close_genome
analysis_events has a value which is a reference to a list where each element is an analysis_event

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a genome_id
scientific_name has a value which is a string
domain has a value which is a string
genetic_code has a value which is an int
source has a value which is a string
source_id has a value which is a string
quality has a value which is a genome_quality_measure
contigs has a value which is a reference to a list where each element is a contig
contigs_handle has a value which is a Handle
features has a value which is a reference to a list where each element is a feature
close_genomes has a value which is a reference to a list where each element is a close_genome
analysis_events has a value which is a reference to a list where each element is an analysis_event


=end text

=back



=head2 genome_metadata

=over 4



=item Description

* Genome metadata. We use this structure to define common metadata
* settings used in the API calls below. It's possible this data should
* have been separated in this way in the genome object itself, but there
* is an extant body of code that assumes the current structure of the genome
* object.


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
id has a value which is a genome_id
scientific_name has a value which is a string
domain has a value which is a string
genetic_code has a value which is an int
source has a value which is a string
source_id has a value which is a string

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
id has a value which is a genome_id
scientific_name has a value which is a string
domain has a value which is a string
genetic_code has a value which is an int
source has a value which is a string
source_id has a value which is a string


=end text

=back



=head2 subsystem

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 variant

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 variant_of_subsystem

=over 4



=item Definition

=begin html

<pre>
a reference to a list containing 2 items:
0: a subsystem
1: a variant

</pre>

=end html

=begin text

a reference to a list containing 2 items:
0: a subsystem
1: a variant


=end text

=back



=head2 variant_subsystem_pairs

=over 4



=item Definition

=begin html

<pre>
a reference to a list where each element is a variant_of_subsystem
</pre>

=end html

=begin text

a reference to a list where each element is a variant_of_subsystem

=end text

=back



=head2 fid

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 role

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 function

=over 4



=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 fid_role_pair

=over 4



=item Definition

=begin html

<pre>
a reference to a list containing 2 items:
0: a fid
1: a role

</pre>

=end html

=begin text

a reference to a list containing 2 items:
0: a fid
1: a role


=end text

=back



=head2 fid_role_pairs

=over 4



=item Definition

=begin html

<pre>
a reference to a list where each element is a fid_role_pair
</pre>

=end html

=begin text

a reference to a list where each element is a fid_role_pair

=end text

=back



=head2 fid_function_pair

=over 4



=item Definition

=begin html

<pre>
a reference to a list containing 2 items:
0: a fid
1: a function

</pre>

=end html

=begin text

a reference to a list containing 2 items:
0: a fid
1: a function


=end text

=back



=head2 fid_function_pairs

=over 4



=item Definition

=begin html

<pre>
a reference to a list where each element is a fid_function_pair
</pre>

=end html

=begin text

a reference to a list where each element is a fid_function_pair

=end text

=back



=head2 reconstructionTO

=over 4



=item Description

Metabolic reconstruction
represents the set of subsystems that we infer are present in this genome


=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
subsystems has a value which is a variant_subsystem_pairs
bindings has a value which is a fid_role_pairs
assignments has a value which is a fid_function_pairs

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
subsystems has a value which is a variant_subsystem_pairs
bindings has a value which is a fid_role_pairs
assignments has a value which is a fid_function_pairs


=end text

=back



=head2 fid_data_tuple

=over 4



=item Definition

=begin html

<pre>
a reference to a list containing 4 items:
0: a fid
1: a md5
2: a location
3: a function

</pre>

=end html

=begin text

a reference to a list containing 4 items:
0: a fid
1: a md5
2: a location
3: a function


=end text

=back



=head2 fid_data_tuples

=over 4



=item Definition

=begin html

<pre>
a reference to a list where each element is a fid_data_tuple
</pre>

=end html

=begin text

a reference to a list where each element is a fid_data_tuple

=end text

=back



=head2 compact_feature

=over 4



=item Description

* This tuple defines a compact form for defining features to be batch-loaded
* into a genome object.


=item Definition

=begin html

<pre>
a reference to a list containing 5 items:
0: (id) a string
1: (location) a string
2: (feature_type) a string
3: (function) a string
4: (aliases) a string

</pre>

=end html

=begin text

a reference to a list containing 5 items:
0: (id) a string
1: (location) a string
2: (feature_type) a string
3: (function) a string
4: (aliases) a string


=end text

=back



=head2 rna_type

=over 4



=item Description

[ validate.enum("5S", "SSU", "LSU", "ALL") ]


=item Definition

=begin html

<pre>
a string
</pre>

=end html

=begin text

a string

=end text

=back



=head2 glimmer3_parameters

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
min_training_len has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
min_training_len has a value which is an int


=end text

=back



=head2 repeat_region_SEED_parameters

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
min_identity has a value which is a float
min_length has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
min_identity has a value which is a float
min_length has a value which is an int


=end text

=back



=head2 kmer_v1_parameters

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
kmer_size has a value which is an int
dataset_name has a value which is a string
return_scores_for_all_proteins has a value which is an int
score_threshold has a value which is an int
hit_threshold has a value which is an int
sequential_hit_threshold has a value which is an int
detailed has a value which is an int
min_hits has a value which is an int
min_size has a value which is an int
max_gap has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
kmer_size has a value which is an int
dataset_name has a value which is a string
return_scores_for_all_proteins has a value which is an int
score_threshold has a value which is an int
hit_threshold has a value which is an int
sequential_hit_threshold has a value which is an int
detailed has a value which is an int
min_hits has a value which is an int
min_size has a value which is an int
max_gap has a value which is an int


=end text

=back



=head2 kmer_v2_parameters

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
min_hits has a value which is an int
max_gap has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
min_hits has a value which is an int
max_gap has a value which is an int


=end text

=back



=head2 resolve_overlapping_features_parameters

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
placeholder has a value which is an int

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
placeholder has a value which is an int


=end text

=back



=head2 pipeline_stage

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
name has a value which is a string
condition has a value which is a string
repeat_region_SEED_parameters has a value which is a repeat_region_SEED_parameters
glimmer3_parameters has a value which is a glimmer3_parameters
kmer_v1_parameters has a value which is a kmer_v1_parameters
kmer_v2_parameters has a value which is a kmer_v2_parameters

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
name has a value which is a string
condition has a value which is a string
repeat_region_SEED_parameters has a value which is a repeat_region_SEED_parameters
glimmer3_parameters has a value which is a glimmer3_parameters
kmer_v1_parameters has a value which is a kmer_v1_parameters
kmer_v2_parameters has a value which is a kmer_v2_parameters


=end text

=back



=head2 workflow

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
stages has a value which is a reference to a list where each element is a pipeline_stage

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
stages has a value which is a reference to a list where each element is a pipeline_stage


=end text

=back



=head2 pipeline_batch_input

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome_id has a value which is a string
data has a value which is a Handle

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome_id has a value which is a string
data has a value which is a Handle


=end text

=back



=head2 pipeline_batch_status_entry

=over 4



=item Definition

=begin html

<pre>
a reference to a hash where the following keys are defined:
genome_id has a value which is a string
status has a value which is a string
stdout has a value which is a Handle
stderr has a value which is a Handle
output has a value which is a Handle

</pre>

=end html

=begin text

a reference to a hash where the following keys are defined:
genome_id has a value which is a string
status has a value which is a string
stdout has a value which is a Handle
stderr has a value which is a Handle
output has a value which is a Handle


=end text

=back



=cut

1;
