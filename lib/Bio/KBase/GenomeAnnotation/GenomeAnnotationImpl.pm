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

use File::Temp;
use File::Slurp;
use Data::Dumper;
use Digest::MD5 'md5_hex';
use Time::HiRes 'gettimeofday';
use Data::Structure::Util qw(unbless);

use Bio::KBase::IDServer::Client;
#use Bio::KBase::KIDL::Helpers qw(json_to_tempfile tempfile_to_json);
use IPC::Run qw(run);
use JSON::XS;

use GenomeTypeObject;
use ANNOserver;
use SeedUtils;
use gjoseqlib;
use StrepRepeats;

use Bio::KBase::DeploymentConfig;

our $idserver_url = 'https://kbase.us/services/idserver';

sub _get_coder
{
    return JSON::XS->new->ascii->pretty->allow_nonref;
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

    my $i = $cfg->setting("idserver_url");
    $idserver_url = $i if $i;

    print STDERR "kmer_v2_data_directory = $dir\n";
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
 NOTE: Many of these "transformations" modify the input hash and
       copy the pointer.  Be warned.

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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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

    my $ok = run(['rast_call_special_proteins',
		  '--seleno',
		  '--id-server', $idc->{url}],
		 '<', $genomeTO_json,
		 '>', \$genomeOut_json,
		 '2>', \$stderr);

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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
    
    #
    #...Find the pyrrolysoproteins...
    #
    use find_special_proteins;
    
    # Reformat the contigs into "Gary-tuples"
    my @contigs;
    foreach my $gctg (@{$genomeTO->{contigs}}) {
	push(@contigs, [$gctg->{id}, undef, $gctg->{dna}]);
    }
    
    #...Only difference from 'call_selenoproteins' is 'pyrrolysine' flag, and annotations written
    my $parms   = { contigs => \@contigs, pyrrolysine => 1 };
    my @results = &find_special_proteins::find_selenoproteins( $parms );
    
    #
    # Allocate IDs for PEGs
    #
    my $n_pegs = @results;
    my $protein_prefix = "$genomeTO->{id}.CDS";
    my $id_server = Bio::KBase::IDServer::Client->new($idserver_url);
    my $peg_id_start = $id_server->allocate_id_range($protein_prefix, $n_pegs) + 0;
    my $next_id = $peg_id_start;
    print STDERR "allocated CDS id start $peg_id_start for $n_pegs CDSs\n";
    
    #
    # Create features for PEGs
    #
    my $features = $genomeTO->{features};
    if (!$features)
    {
	$features = [];
	$genomeTO->{features} = $features;
    }
    
    # Reformat result from &find_special_proteins::find_selenoproteins({pyrrolysine => 1}).
    foreach my $feature (@results) {
	my $loc  = $feature->{location};
	my $seq  = $feature->{sequence};
	my $func = $feature->{reference_def};
	
	my ($contig, $start, $stop, $strand) = &SeedUtils::parse_location( $feature->{location} );
	my $len = abs($stop - $start) + 1;
	my $strand = ($stop > $start) ? '+' : '-';
	
	my $kb_id = "$protein_prefix.$next_id";
	++$next_id;
	
	my $annos = [];
	push(@$annos, ["Set function to\n$func\nfor initial gene call performed by call_pyrrolysoproteins",
		       'genome annotation service',
		       time
		       ]);
	
	my $feature = {
	    id => $kb_id,
	    location => [[ $contig, $start, $strand, $len ]],
	    type => 'CDS',
	    protein_translation => $seq,
	    aliases => [],
	    $func ? (function => $func) : (),
	    annotations => $annos,
	};
	push(@$features, $feature);
    }
    $return = $genomeTO;

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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	for my $type (qw('5S' SSU LSU))
	{
	    if ($types{lc($type)})
	    {
		push(@type_args, "-$type");
	    }
	}
    }

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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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

  $return = $obj->call_features_CDS_glimmer3($genomeTO)

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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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

sub call_features_CDS_glimmer3
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_CDS_glimmer3:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_CDS_glimmer3');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_CDS_glimmer3
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	die "Error running kmer_search: @cmd\n";
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
	    -event_id 	     => $event_id,
	    -quality_measure => $quality,
	});
    }

    $genome_out = $genome_in;
    $genome_out->prepare_for_return();
    unbless $genome_out;

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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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

    $genome_out = $genome_in;
    unbless $genome_out;

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




=head2 call_features_ProtoCDS_kmer_v1

  $return = $obj->call_features_ProtoCDS_kmer_v1($genomeTO)

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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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

sub call_features_ProtoCDS_kmer_v1
{
    my $self = shift;
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_ProtoCDS_kmer_v1:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_ProtoCDS_kmer_v1');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_ProtoCDS_kmer_v1
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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

    while(<$res_fh>)
    {
	chomp;
	my($contig, $left, $right, $strand, $frame, $hit_count, $function, $weighted_hit_count) = split(/\t/);

	next unless $left =~ /^\d+$/ && $right =~ /^\d+$/;

	my $confidence = 1 - 0.5 ** ($weighted_hit_count / 3);

	my $quality = {
	    existence_confidence => $confidence,
	    hit_count => 0 + $hit_count,
	    weighted_hit_count => 0 + $weighted_hit_count,
	};


	my $begin = 0 + (($strand eq '+') ? $left : $right);
	my $len = 1 + $right - $left;
	my $loc = [[$contig, $begin, $strand, $len]];
	$genome_in->add_feature({
	    -id_client 	     => $idc,
	    -id_prefix 	     => $genome_in->{id},
	    -type 	     => $type,
	    -location 	     => $loc,
	    -function 	     => $function,
	    -event_id 	     => $event_id,
	    -quality_measure => $quality,
	});
    }

    $genome_out = $genome_in;
    $genome_out->prepare_for_return();
    unbless $genome_out;

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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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

  $return = $obj->call_features_crispr($genomeTO)

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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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
    my($genomeTO) = @_;

    my @_bad_arguments;
    (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genomeTO\" (value was \"$genomeTO\")");
    if (@_bad_arguments) {
	my $msg = "Invalid arguments passed to call_features_crispr:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_crispr');
    }

    my $ctx = $Bio::KBase::GenomeAnnotation::Service::CallContext;
    my($return);
    #BEGIN call_features_crispr
    #END call_features_crispr
    my @_bad_returns;
    (ref($return) eq 'HASH') or push(@_bad_returns, "Invalid type for return variable \"return\" (value was \"$return\")");
    if (@_bad_returns) {
	my $msg = "Invalid returns passed to call_features_crispr:\n" . join("", map { "\t$_\n" } @_bad_returns);
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
							       method_name => 'call_features_crispr');
    }
    return($return);
}




=head2 export_genome

  $exported_data = $obj->export_genome($genome_in, $format)

=over 4

=item Parameter and return types

=begin html

<pre>
$genome_in is a genomeTO
$format is a string
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
	overlap_allowed has a value which is a bool
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
	overlap_allowed has a value which is a bool
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

=back

=cut

sub export_genome
{
    my $self = shift;
    my($genome_in, $format) = @_;

    my @_bad_arguments;
    (ref($genome_in) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument \"genome_in\" (value was \"$genome_in\")");
    (!ref($format)) or push(@_bad_arguments, "Invalid type for argument \"format\" (value was \"$format\")");
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

    my @cmd = ("rast_export_genome", "--input", $tmp_in, "--output", $tmp_out, $format);

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

Is this a real feature?


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
overlap_allowed has a value which is a bool
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
overlap_allowed has a value which is a bool
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
features has a value which is a reference to a list where each element is a feature
close_genomes has a value which is a reference to a list where each element is a close_genome
analysis_events has a value which is a reference to a list where each element is an analysis_event


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



=cut

1;
