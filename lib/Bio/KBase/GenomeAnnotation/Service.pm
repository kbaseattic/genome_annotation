package Bio::KBase::GenomeAnnotation::Service;


use Data::Dumper;
use Moose;
use POSIX;
use JSON;
use Bio::KBase::Log;
use Class::Load qw();
use Config::Simple;
my $get_time = sub { time, 0 };
eval {
    require Time::HiRes;
    $get_time = sub { Time::HiRes::gettimeofday };
};

use Bio::KBase::AuthToken;

extends 'RPC::Any::Server::JSONRPC::PSGI';

has 'instance_dispatch' => (is => 'ro', isa => 'HashRef');
has 'user_auth' => (is => 'ro', isa => 'UserAuth');
has 'valid_methods' => (is => 'ro', isa => 'HashRef', lazy => 1,
			builder => '_build_valid_methods');
has 'loggers' => (is => 'ro', required => 1, builder => '_build_loggers');
has 'config' => (is => 'ro', required => 1, builder => '_build_config');

our $CallContext;

our %return_counts = (
        'genome_ids_to_genomes' => 1,
        'create_genome' => 1,
        'create_genome_from_genbank' => 1,
        'create_genome_from_SEED' => 1,
        'create_genome_from_RAST' => 1,
        'set_metadata' => 1,
        'add_contigs' => 1,
        'add_contigs_from_handle' => 1,
        'add_features' => 1,
        'genomeTO_to_reconstructionTO' => 1,
        'genomeTO_to_feature_data' => 1,
        'reconstructionTO_to_roles' => 1,
        'reconstructionTO_to_subsystems' => 1,
        'assign_functions_to_CDSs' => 1,
        'annotate_genome' => 1,
        'call_selenoproteins' => 1,
        'call_pyrrolysoproteins' => 1,
        'call_features_selenoprotein' => 1,
        'call_features_pyrrolysoprotein' => 1,
        'call_features_insertion_sequences' => 1,
        'call_features_rRNA_SEED' => 1,
        'call_features_tRNA_trnascan' => 1,
        'call_RNAs' => 1,
        'call_features_CDS_glimmer3' => 1,
        'call_features_CDS_prodigal' => 1,
        'call_features_CDS_genemark' => 1,
        'call_features_CDS_SEED_projection' => 1,
        'call_features_CDS_FragGeneScan' => 1,
        'call_features_repeat_region_SEED' => 1,
        'call_features_prophage_phispy' => 1,
        'call_features_scan_for_matches' => 1,
        'annotate_proteins_similarity' => 1,
        'annotate_proteins_kmer_v1' => 1,
        'annotate_proteins_kmer_v2' => 1,
        'resolve_overlapping_features' => 1,
        'propagate_genbank_feature_metadata' => 1,
        'call_features_ProtoCDS_kmer_v1' => 1,
        'call_features_ProtoCDS_kmer_v2' => 1,
        'enumerate_special_protein_databases' => 1,
        'compute_special_proteins' => 1,
        'annotate_special_proteins' => 1,
        'annotate_families_figfam_v1' => 1,
        'annotate_families_patric' => 1,
        'annotate_null_to_hypothetical' => 1,
        'annotate_strain_type_MLST' => 1,
        'compute_cdd' => 1,
        'annotate_proteins' => 1,
        'estimate_crude_phylogenetic_position_kmer' => 1,
        'find_close_neighbors' => 1,
        'call_features_strep_suis_repeat' => 1,
        'call_features_strep_pneumo_repeat' => 1,
        'call_features_crispr' => 1,
        'update_functions' => 1,
        'renumber_features' => 1,
        'export_genome' => 1,
        'enumerate_classifiers' => 1,
        'query_classifier_groups' => 1,
        'query_classifier_taxonomies' => 1,
        'classify_into_bins' => 1,
        'classify_full' => 3,
        'default_workflow' => 1,
        'complete_workflow_template' => 1,
        'run_pipeline' => 1,
        'pipeline_batch_start' => 1,
        'pipeline_batch_status' => 1,
        'pipeline_batch_enumerate_batches' => 1,
        'version' => 1,
);

our %method_authentication = (
        'genome_ids_to_genomes' => 'none',
        'create_genome' => 'none',
        'create_genome_from_genbank' => 'none',
        'create_genome_from_SEED' => 'none',
        'create_genome_from_RAST' => 'none',
        'set_metadata' => 'none',
        'add_contigs' => 'none',
        'add_contigs_from_handle' => 'none',
        'add_features' => 'none',
        'genomeTO_to_reconstructionTO' => 'none',
        'genomeTO_to_feature_data' => 'none',
        'reconstructionTO_to_roles' => 'none',
        'reconstructionTO_to_subsystems' => 'none',
        'assign_functions_to_CDSs' => 'none',
        'annotate_genome' => 'none',
        'call_selenoproteins' => 'none',
        'call_pyrrolysoproteins' => 'none',
        'call_features_selenoprotein' => 'none',
        'call_features_pyrrolysoprotein' => 'none',
        'call_features_insertion_sequences' => 'none',
        'call_features_rRNA_SEED' => 'none',
        'call_features_tRNA_trnascan' => 'none',
        'call_RNAs' => 'none',
        'call_features_CDS_glimmer3' => 'none',
        'call_features_CDS_prodigal' => 'none',
        'call_features_CDS_genemark' => 'none',
        'call_features_CDS_SEED_projection' => 'none',
        'call_features_CDS_FragGeneScan' => 'none',
        'call_features_repeat_region_SEED' => 'none',
        'call_features_prophage_phispy' => 'none',
        'call_features_scan_for_matches' => 'none',
        'annotate_proteins_similarity' => 'none',
        'annotate_proteins_kmer_v1' => 'none',
        'annotate_proteins_kmer_v2' => 'none',
        'resolve_overlapping_features' => 'none',
        'propagate_genbank_feature_metadata' => 'none',
        'call_features_ProtoCDS_kmer_v1' => 'none',
        'call_features_ProtoCDS_kmer_v2' => 'none',
        'enumerate_special_protein_databases' => 'none',
        'compute_special_proteins' => 'none',
        'annotate_special_proteins' => 'none',
        'annotate_families_figfam_v1' => 'none',
        'annotate_families_patric' => 'none',
        'annotate_null_to_hypothetical' => 'none',
        'annotate_strain_type_MLST' => 'none',
        'compute_cdd' => 'none',
        'annotate_proteins' => 'none',
        'estimate_crude_phylogenetic_position_kmer' => 'none',
        'find_close_neighbors' => 'none',
        'call_features_strep_suis_repeat' => 'none',
        'call_features_strep_pneumo_repeat' => 'none',
        'call_features_crispr' => 'none',
        'update_functions' => 'none',
        'renumber_features' => 'none',
        'export_genome' => 'none',
        'enumerate_classifiers' => 'none',
        'query_classifier_groups' => 'none',
        'query_classifier_taxonomies' => 'none',
        'classify_into_bins' => 'none',
        'classify_full' => 'none',
        'default_workflow' => 'none',
        'complete_workflow_template' => 'none',
        'run_pipeline' => 'none',
        'pipeline_batch_start' => 'required',
        'pipeline_batch_status' => 'required',
        'pipeline_batch_enumerate_batches' => 'required',
);


sub _build_valid_methods
{
    my($self) = @_;
    my $methods = {
        'genome_ids_to_genomes' => 1,
        'create_genome' => 1,
        'create_genome_from_genbank' => 1,
        'create_genome_from_SEED' => 1,
        'create_genome_from_RAST' => 1,
        'set_metadata' => 1,
        'add_contigs' => 1,
        'add_contigs_from_handle' => 1,
        'add_features' => 1,
        'genomeTO_to_reconstructionTO' => 1,
        'genomeTO_to_feature_data' => 1,
        'reconstructionTO_to_roles' => 1,
        'reconstructionTO_to_subsystems' => 1,
        'assign_functions_to_CDSs' => 1,
        'annotate_genome' => 1,
        'call_selenoproteins' => 1,
        'call_pyrrolysoproteins' => 1,
        'call_features_selenoprotein' => 1,
        'call_features_pyrrolysoprotein' => 1,
        'call_features_insertion_sequences' => 1,
        'call_features_rRNA_SEED' => 1,
        'call_features_tRNA_trnascan' => 1,
        'call_RNAs' => 1,
        'call_features_CDS_glimmer3' => 1,
        'call_features_CDS_prodigal' => 1,
        'call_features_CDS_genemark' => 1,
        'call_features_CDS_SEED_projection' => 1,
        'call_features_CDS_FragGeneScan' => 1,
        'call_features_repeat_region_SEED' => 1,
        'call_features_prophage_phispy' => 1,
        'call_features_scan_for_matches' => 1,
        'annotate_proteins_similarity' => 1,
        'annotate_proteins_kmer_v1' => 1,
        'annotate_proteins_kmer_v2' => 1,
        'resolve_overlapping_features' => 1,
        'propagate_genbank_feature_metadata' => 1,
        'call_features_ProtoCDS_kmer_v1' => 1,
        'call_features_ProtoCDS_kmer_v2' => 1,
        'enumerate_special_protein_databases' => 1,
        'compute_special_proteins' => 1,
        'annotate_special_proteins' => 1,
        'annotate_families_figfam_v1' => 1,
        'annotate_families_patric' => 1,
        'annotate_null_to_hypothetical' => 1,
        'annotate_strain_type_MLST' => 1,
        'compute_cdd' => 1,
        'annotate_proteins' => 1,
        'estimate_crude_phylogenetic_position_kmer' => 1,
        'find_close_neighbors' => 1,
        'call_features_strep_suis_repeat' => 1,
        'call_features_strep_pneumo_repeat' => 1,
        'call_features_crispr' => 1,
        'update_functions' => 1,
        'renumber_features' => 1,
        'export_genome' => 1,
        'enumerate_classifiers' => 1,
        'query_classifier_groups' => 1,
        'query_classifier_taxonomies' => 1,
        'classify_into_bins' => 1,
        'classify_full' => 1,
        'default_workflow' => 1,
        'complete_workflow_template' => 1,
        'run_pipeline' => 1,
        'pipeline_batch_start' => 1,
        'pipeline_batch_status' => 1,
        'pipeline_batch_enumerate_batches' => 1,
        'version' => 1,
    };
    return $methods;
}

my $DEPLOY = 'KB_DEPLOYMENT_CONFIG';
my $SERVICE = 'KB_SERVICE_NAME';

sub get_config_file
{
    my ($self) = @_;
    if(!defined $ENV{$DEPLOY}) {
        return undef;
    }
    return $ENV{$DEPLOY};
}

sub get_service_name
{
    my ($self) = @_;
    if(!defined $ENV{$SERVICE}) {
        return 'GenomeAnnotation';
    }
    return $ENV{$SERVICE};
}

sub _build_config
{
    my ($self) = @_;
    my $sn = $self->get_service_name();
    my $cf = $self->get_config_file();
    if (!($cf)) {
        return {};
    }
    my $cfg = new Config::Simple($cf);
    my $cfgdict = $cfg->get_block($sn);
    if (!($cfgdict)) {
        return {};
    }
    return $cfgdict;
}

sub logcallback
{
    my ($self) = @_;
    $self->loggers()->{serverlog}->set_log_file(
        $self->{loggers}->{userlog}->get_log_file());
}

sub log
{
    my ($self, $level, $context, $message, $tag) = @_;
    my $user = defined($context->user_id()) ? $context->user_id(): undef; 
    $self->loggers()->{serverlog}->log_message($level, $message, $user, 
        $context->module(), $context->method(), $context->call_id(),
        $context->client_ip(), $tag);
}

sub _build_loggers
{
    my ($self) = @_;
    my $submod = $self->get_service_name();
    my $loggers = {};
    my $callback = sub {$self->logcallback();};
    $loggers->{userlog} = Bio::KBase::Log->new(
            $submod, {}, {ip_address => 1, authuser => 1, module => 1,
            method => 1, call_id => 1, changecallback => $callback,
	    tag => 1,
            config => $self->get_config_file()});
    $loggers->{serverlog} = Bio::KBase::Log->new(
            $submod, {}, {ip_address => 1, authuser => 1, module => 1,
            method => 1, call_id => 1,
	    tag => 1,
            logfile => $loggers->{userlog}->get_log_file()});
    $loggers->{serverlog}->set_log_level(6);
    return $loggers;
}

#
# Override method from RPC::Any::Server::JSONRPC 
# to eliminate the deprecation warning for Class::MOP::load_class.
#
sub _default_error {
    my ($self, %params) = @_;
    my $version = $self->default_version;
    $version =~ s/\./_/g;
    my $error_class = "JSON::RPC::Common::Procedure::Return::Version_${version}::Error";
    Class::Load::load_class($error_class);
    my $error = $error_class->new(%params);
    my $return_class = "JSON::RPC::Common::Procedure::Return::Version_$version";
    Class::Load::load_class($return_class);
    return $return_class->new(error => $error);
}


#override of RPC::Any::Server
sub handle_error {
    my ($self, $error) = @_;
    
    unless (ref($error) eq 'HASH' ||
           (blessed $error and $error->isa('RPC::Any::Exception'))) {
        $error = RPC::Any::Exception::PerlError->new(message => $error);
    }
    my $output;
    eval {
        my $encoded_error = $self->encode_output_from_exception($error);
        $output = $self->produce_output($encoded_error);
    };
    
    return $output if $output;
    
    die "$error\n\nAlso, an error was encountered while trying to send"
        . " this error: $@\n";
}

#override of RPC::Any::JSONRPC
sub encode_output_from_exception {
    my ($self, $exception) = @_;
    my %error_params;
    if (ref($exception) eq 'HASH') {
        %error_params = %{$exception};
        if(defined($error_params{context})) {
            my @errlines;
            $errlines[0] = $error_params{message};
            push @errlines, split("\n", $error_params{data});
            $self->log($Bio::KBase::Log::ERR, $error_params{context}, \@errlines);
            delete $error_params{context};
        }
    } else {
        %error_params = (
            message => $exception->message,
            code    => $exception->code,
        );
    }
    my $json_error;
    if ($self->_last_call) {
        $json_error = $self->_last_call->return_error(%error_params);
    }
    # Default to default_version. This happens when we throw an exception
    # before inbound parsing is complete.
    else {
        $json_error = $self->_default_error(%error_params);
    }
    return $self->encode_output_from_object($json_error);
}

sub trim {
    my ($str) = @_;
    if (!(defined $str)) {
        return $str;
    }
    $str =~ s/^\s+|\s+$//g;
    return $str;
}

sub getIPAddress {
    my ($self) = @_;
    my $xFF = trim($self->_plack_req->header("X-Forwarded-For"));
    my $realIP = trim($self->_plack_req->header("X-Real-IP"));
    my $nh = $self->config->{"dont_trust_x_ip_headers"};
    my $trustXHeaders = !(defined $nh) || $nh ne "true";

    if ($trustXHeaders) {
        if ($xFF) {
            my @tmp = split(",", $xFF);
            return trim($tmp[0]);
        }
        if ($realIP) {
            return $realIP;
        }
    }
    return $self->_plack_req->address;
}

#
# Ping method reflected from /ping on the service.
#
sub ping
{
    my($self, $env) = @_;
    return [ 200, ["Content-type" => "text/plain"], [ "OK\n" ] ];
}


#
# Authenticated ping method reflected from /auth_ping on the service.
#
sub auth_ping
{
    my($self, $env) = @_;

    my $req = Plack::Request->new($env);
    my $token = $req->header("Authorization");

    if (!$token)
    {
	return [401, [], ["Authentication required"]];
    }

    my $auth_token = Bio::KBase::AuthToken->new(token => $token, ignore_authrc => 1);
    my $valid = $auth_token->validate();

    if ($valid)
    {
	return [200, ["Content-type" => "text/plain"], ["OK " . $auth_token->user_id . "\n"]];
    }
    else
    {
	return [403, [], "Authentication failed"];
    }
}

sub call_method {
    my ($self, $data, $method_info) = @_;

    my ($module, $method, $modname) = @$method_info{qw(module method modname)};
    
    my $ctx = Bio::KBase::GenomeAnnotation::ServiceContext->new($self->{loggers}->{userlog},
                           client_ip => $self->getIPAddress());
    $ctx->module($modname);
    $ctx->method($method);
    $ctx->call_id($self->{_last_call}->{id});
    
    my $args = $data->{arguments};

{
    # Service GenomeAnnotation requires authentication.

    my $method_auth = $method_authentication{$method};
    $ctx->authenticated(0);
    if ($method_auth eq 'none')
    {
	# No authentication required here. Move along.
    }
    else
    {
	my $token = $self->_plack_req->header("Authorization");

	if (!$token && $method_auth eq 'required')
	{
	    $self->exception('PerlError', "Authentication required for GenomeAnnotation but no authentication header was passed");
	}

	my $auth_token = Bio::KBase::AuthToken->new(token => $token, ignore_authrc => 1);
	my $valid = $auth_token->validate();
	# Only throw an exception if authentication was required and it fails
	if ($method_auth eq 'required' && !$valid)
	{
	    $self->exception('PerlError', "Token validation failed: " . $auth_token->error_message);
	} elsif ($valid) {
	    $ctx->authenticated(1);
	    $ctx->user_id($auth_token->user_id);
	    $ctx->token( $token);
	}
    }
}
    my $new_isa = $self->get_package_isa($module);
    no strict 'refs';
    local @{"${module}::ISA"} = @$new_isa;
    local $CallContext = $ctx;
    my @result;
    {
	# 
	# Process tag and metadata information if present.
	#
	my $tag = $self->_plack_req->header("Kbrpc-Tag");
	if (!$tag)
	{
	    my ($t, $us) = &$get_time();
	    $us = sprintf("%06d", $us);
	    my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
	    $tag = "S:$self->{hostname}:$$:$ts";
	}
	local $ENV{KBRPC_TAG} = $tag;
	my $kb_metadata = $self->_plack_req->header("Kbrpc-Metadata");
	my $kb_errordest = $self->_plack_req->header("Kbrpc-Errordest");
	local $ENV{KBRPC_METADATA} = $kb_metadata if $kb_metadata;
	local $ENV{KBRPC_ERROR_DEST} = $kb_errordest if $kb_errordest;

	my $stderr = Bio::KBase::GenomeAnnotation::ServiceStderrWrapper->new($ctx, $get_time);
	$ctx->stderr($stderr);

        my $xFF = $self->_plack_req->header("X-Forwarded-For");
        if ($xFF) {
            $self->log($Bio::KBase::Log::INFO, $ctx,
                "X-Forwarded-For: " . $xFF, $tag);
        }
	
        my $err;
        eval {
            $self->log($Bio::KBase::Log::INFO, $ctx, "start method", $tag);
	    local $SIG{__WARN__} = sub {
		my($msg) = @_;
		$stderr->log($msg);
		print STDERR $msg;
	    };

            @result = $module->$method(@{ $data->{arguments} });
            $self->log($Bio::KBase::Log::INFO, $ctx, "end method", $tag);
        };
	
        if ($@)
        {
            my $err = $@;
	    $stderr->log($err);
	    $ctx->stderr(undef);
	    undef $stderr;
            $self->log($Bio::KBase::Log::INFO, $ctx, "fail method", $tag);
            my $nicerr;
            if(ref($err) eq "Bio::KBase::Exceptions::KBaseException") {
                $nicerr = {code => -32603, # perl error from RPC::Any::Exception
                           message => $err->error,
                           data => $err->trace->as_string,
                           context => $ctx
                           };
            } else {
                my $str = "$err";
                $str =~ s/Bio::KBase::CDMI::Service::call_method.*//s; # is this still necessary? not sure
                my $msg = $str;
                $msg =~ s/ at [^\s]+.pm line \d+.\n$//;
                $nicerr =  {code => -32603, # perl error from RPC::Any::Exception
                            message => $msg,
                            data => $str,
                            context => $ctx
                            };
            }
            die $nicerr;
        }
	$ctx->stderr(undef);
	undef $stderr;
    }
    my $result;
    if ($return_counts{$method} == 1)
    {
        $result = [[$result[0]]];
    }
    else
    {
        $result = \@result;
    }
    return $result;
}


sub get_method
{
    my ($self, $data) = @_;
    
    my $full_name = $data->{method};
    
    $full_name =~ /^(\S+)\.([^\.]+)$/;
    my ($package, $method) = ($1, $2);
    
    if (!$package || !$method) {
	$self->exception('NoSuchMethod',
			 "'$full_name' is not a valid method. It must"
			 . " contain a package name, followed by a period,"
			 . " followed by a method name.");
    }

    if (!$self->valid_methods->{$method})
    {
	$self->exception('NoSuchMethod',
			 "'$method' is not a valid method in service GenomeAnnotation.");
    }
	
    my $inst = $self->instance_dispatch->{$package};
    my $module;
    if ($inst)
    {
	$module = $inst;
    }
    else
    {
	$module = $self->get_module($package);
	if (!$module) {
	    $self->exception('NoSuchMethod',
			     "There is no method package named '$package'.");
	}
	
	Class::Load::load_class($module);
    }
    
    if (!$module->can($method)) {
	$self->exception('NoSuchMethod',
			 "There is no method named '$method' in the"
			 . " '$package' package.");
    }
    
    return { module => $module, method => $method, modname => $package };
}

package Bio::KBase::GenomeAnnotation::ServiceContext;

use strict;

=head1 NAME

Bio::KBase::GenomeAnnotation::ServiceContext

head1 DESCRIPTION

A KB RPC context contains information about the invoker of this
service. If it is an authenticated service the authenticated user
record is available via $context->user. The client IP address
is available via $context->client_ip.

=cut

use base 'Class::Accessor';

__PACKAGE__->mk_accessors(qw(user_id client_ip authenticated token
                             module method call_id hostname stderr));

sub new
{
    my($class, $logger, %opts) = @_;
    
    my $self = {
        %opts,
    };
    chomp($self->{hostname} = `hostname`);
    $self->{hostname} ||= 'unknown-host';
    $self->{_logger} = $logger;
    $self->{_debug_levels} = {7 => 1, 8 => 1, 9 => 1,
                              'DEBUG' => 1, 'DEBUG2' => 1, 'DEBUG3' => 1};
    return bless $self, $class;
}

sub _get_user
{
    my ($self) = @_;
    return defined($self->user_id()) ? $self->user_id(): undef; 
}

sub _log
{
    my ($self, $level, $message) = @_;
    $self->{_logger}->log_message($level, $message, $self->_get_user(),
        $self->module(), $self->method(), $self->call_id(),
        $self->client_ip());
}

sub log_err
{
    my ($self, $message) = @_;
    $self->_log($Bio::KBase::Log::ERR, $message);
}

sub log_info
{
    my ($self, $message) = @_;
    $self->_log($Bio::KBase::Log::INFO, $message);
}

sub log_debug
{
    my ($self, $message, $level) = @_;
    if(!defined($level)) {
        $level = 1;
    }
    if($self->{_debug_levels}->{$level}) {
    } else {
        if ($level =~ /\D/ || $level < 1 || $level > 3) {
            die "Invalid log level: $level";
        }
        $level += 6;
    }
    $self->_log($level, $message);
}

sub set_log_level
{
    my ($self, $level) = @_;
    $self->{_logger}->set_log_level($level);
}

sub get_log_level
{
    my ($self) = @_;
    return $self->{_logger}->get_log_level();
}

sub clear_log_level
{
    my ($self) = @_;
    $self->{_logger}->clear_user_log_level();
}

package Bio::KBase::GenomeAnnotation::ServiceStderrWrapper;

use strict;
use POSIX;
use Time::HiRes 'gettimeofday';

sub new
{
    my($class, $ctx, $get_time) = @_;
    my $self = {
	get_time => $get_time,
    };
    my $dest = $ENV{KBRPC_ERROR_DEST} if exists $ENV{KBRPC_ERROR_DEST};
    my $tag = $ENV{KBRPC_TAG} if exists $ENV{KBRPC_TAG};
    my ($t, $us) = gettimeofday();
    $us = sprintf("%06d", $us);
    my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);

    my $name = join(".", $ctx->module, $ctx->method, $ctx->hostname, $ts);

    if ($dest && $dest =~ m,^/,)
    {
	#
	# File destination
	#
	my $fh;

	if ($tag)
	{
	    $tag =~ s,/,_,g;
	    $dest = "$dest/$tag";
	    if (! -d $dest)
	    {
		mkdir($dest);
	    }
	}
	if (open($fh, ">", "$dest/$name"))
	{
	    $self->{file} = "$dest/$name";
	    $self->{dest} = $fh;
	}
	else
	{
	    warn "Cannot open log file $dest/$name: $!";
	}
    }
    else
    {
	#
	# Log to string.
	#
	my $stderr;
	$self->{dest} = \$stderr;
    }
    
    bless $self, $class;

    for my $e (sort { $a cmp $b } keys %ENV)
    {
	$self->log_cmd($e, $ENV{$e});
    }
    return $self;
}

sub redirect
{
    my($self) = @_;
    if ($self->{dest})
    {
	return("2>", $self->{dest});
    }
    else
    {
	return ();
    }
}

sub redirect_both
{
    my($self) = @_;
    if ($self->{dest})
    {
	return(">&", $self->{dest});
    }
    else
    {
	return ();
    }
}

sub timestamp
{
    my($self) = @_;
    my ($t, $us) = $self->{get_time}->();
    $us = sprintf("%06d", $us);
    my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);
    return $ts;
}

sub log
{
    my($self, $str) = @_;
    my $d = $self->{dest};
    my $ts = $self->timestamp();
    if (ref($d) eq 'SCALAR')
    {
	$$d .= "[$ts] " . $str . "\n";
	return 1;
    }
    elsif ($d)
    {
	print $d "[$ts] " . $str . "\n";
	return 1;
    }
    return 0;
}

sub log_cmd
{
    my($self, @cmd) = @_;
    my $d = $self->{dest};
    my $str;
    my $ts = $self->timestamp();
    if (ref($cmd[0]))
    {
	$str = join(" ", @{$cmd[0]});
    }
    else
    {
	$str = join(" ", @cmd);
    }
    if (ref($d) eq 'SCALAR')
    {
	$$d .= "[$ts] " . $str . "\n";
    }
    elsif ($d)
    {
	print $d "[$ts] " . $str . "\n";
    }
	 
}

sub dest
{
    my($self) = @_;
    return $self->{dest};
}

sub text_value
{
    my($self) = @_;
    if (ref($self->{dest}) eq 'SCALAR')
    {
	my $r = $self->{dest};
	return $$r;
    }
    else
    {
	return $self->{file};
    }
}


1;
