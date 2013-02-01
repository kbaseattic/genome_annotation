package Bio::KBase::GenomeAnnotation::Client;

use JSON::RPC::Client;
use strict;
use Data::Dumper;
use URI;
use Bio::KBase::Exceptions;

# Client version should match Impl version
# This is a Semantic Version number,
# http://semver.org
our $VERSION = "0.1.0";

=head1 NAME

Bio::KBase::GenomeAnnotation::Client

=head1 DESCRIPTION



=cut

sub new
{
    my($class, $url) = @_;

    my $self = {
	client => Bio::KBase::GenomeAnnotation::Client::RpcClient->new,
	url => $url,
    };
    my $ua = $self->{client}->ua;	 
    my $timeout = $ENV{CDMI_TIMEOUT} || (30 * 60);	 
    $ua->timeout($timeout);
    bless $self, $class;
    #    $self->_validate_version();
    return $self;
}




=head2 $result = genomeTO_to_reconstructionTO(genomeTO)



=cut

sub genomeTO_to_reconstructionTO
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function genomeTO_to_reconstructionTO (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to genomeTO_to_reconstructionTO:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'genomeTO_to_reconstructionTO');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.genomeTO_to_reconstructionTO",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'genomeTO_to_reconstructionTO',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method genomeTO_to_reconstructionTO",
					    status_line => $self->{client}->status_line,
					    method_name => 'genomeTO_to_reconstructionTO',
				       );
    }
}



=head2 $result = genomeTO_to_feature_data(genomeTO)



=cut

sub genomeTO_to_feature_data
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function genomeTO_to_feature_data (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to genomeTO_to_feature_data:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'genomeTO_to_feature_data');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.genomeTO_to_feature_data",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'genomeTO_to_feature_data',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method genomeTO_to_feature_data",
					    status_line => $self->{client}->status_line,
					    method_name => 'genomeTO_to_feature_data',
				       );
    }
}



=head2 $result = reconstructionTO_to_roles(reconstructionTO)



=cut

sub reconstructionTO_to_roles
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function reconstructionTO_to_roles (received $n, expecting 1)");
    }
    {
	my($reconstructionTO) = @args;

	my @_bad_arguments;
        (ref($reconstructionTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"reconstructionTO\" (value was \"$reconstructionTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to reconstructionTO_to_roles:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'reconstructionTO_to_roles');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.reconstructionTO_to_roles",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'reconstructionTO_to_roles',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method reconstructionTO_to_roles",
					    status_line => $self->{client}->status_line,
					    method_name => 'reconstructionTO_to_roles',
				       );
    }
}



=head2 $result = reconstructionTO_to_subsystems(reconstructionTO)



=cut

sub reconstructionTO_to_subsystems
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function reconstructionTO_to_subsystems (received $n, expecting 1)");
    }
    {
	my($reconstructionTO) = @args;

	my @_bad_arguments;
        (ref($reconstructionTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"reconstructionTO\" (value was \"$reconstructionTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to reconstructionTO_to_subsystems:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'reconstructionTO_to_subsystems');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.reconstructionTO_to_subsystems",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'reconstructionTO_to_subsystems',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method reconstructionTO_to_subsystems",
					    status_line => $self->{client}->status_line,
					    method_name => 'reconstructionTO_to_subsystems',
				       );
    }
}



=head2 $result = annotate_genome(genomeTO)

Given a genome object populated with contig data, perform gene calling
and functional annotation and return the annotated genome.
 NOTE: Many of these "transformations" modify the input hash and
       copy the pointer.  Be warned.

=cut

sub annotate_genome
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function annotate_genome (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_genome:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'annotate_genome');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.annotate_genome",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'annotate_genome',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method annotate_genome",
					    status_line => $self->{client}->status_line,
					    method_name => 'annotate_genome',
				       );
    }
}



=head2 $result = call_selenoproteins(genomeTO)



=cut

sub call_selenoproteins
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function call_selenoproteins (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_selenoproteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'call_selenoproteins');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.call_selenoproteins",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'call_selenoproteins',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method call_selenoproteins",
					    status_line => $self->{client}->status_line,
					    method_name => 'call_selenoproteins',
				       );
    }
}



=head2 $result = call_pyrrolysoproteins(genomeTO)



=cut

sub call_pyrrolysoproteins
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function call_pyrrolysoproteins (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_pyrrolysoproteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'call_pyrrolysoproteins');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.call_pyrrolysoproteins",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'call_pyrrolysoproteins',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method call_pyrrolysoproteins",
					    status_line => $self->{client}->status_line,
					    method_name => 'call_pyrrolysoproteins',
				       );
    }
}



=head2 $result = call_RNAs(genomeTO)



=cut

sub call_RNAs
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function call_RNAs (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_RNAs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'call_RNAs');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.call_RNAs",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'call_RNAs',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method call_RNAs",
					    status_line => $self->{client}->status_line,
					    method_name => 'call_RNAs',
				       );
    }
}



=head2 $result = call_CDSs(genomeTO)



=cut

sub call_CDSs
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function call_CDSs (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_CDSs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'call_CDSs');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.call_CDSs",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'call_CDSs',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method call_CDSs",
					    status_line => $self->{client}->status_line,
					    method_name => 'call_CDSs',
				       );
    }
}



=head2 $result = find_close_neighbors(genomeTO)



=cut

sub find_close_neighbors
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function find_close_neighbors (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to find_close_neighbors:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'find_close_neighbors');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.find_close_neighbors",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'find_close_neighbors',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method find_close_neighbors",
					    status_line => $self->{client}->status_line,
					    method_name => 'find_close_neighbors',
				       );
    }
}



=head2 $result = assign_functions_to_CDSs(genomeTO)



=cut

sub assign_functions_to_CDSs
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function assign_functions_to_CDSs (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to assign_functions_to_CDSs:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'assign_functions_to_CDSs');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.assign_functions_to_CDSs",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'assign_functions_to_CDSs',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method assign_functions_to_CDSs",
					    status_line => $self->{client}->status_line,
					    method_name => 'assign_functions_to_CDSs',
				       );
    }
}



=head2 $result = annotate_proteins(genomeTO)

Given a genome object populated with feature data, reannotate
the features that have protein translations. Return the updated
genome object.

=cut

sub annotate_proteins
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function annotate_proteins (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to annotate_proteins:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'annotate_proteins');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.annotate_proteins",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'annotate_proteins',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method annotate_proteins",
					    status_line => $self->{client}->status_line,
					    method_name => 'annotate_proteins',
				       );
    }
}



=head2 $result = call_CDSs_by_projection(genomeTO)



=cut

sub call_CDSs_by_projection
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function call_CDSs_by_projection (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to call_CDSs_by_projection:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'call_CDSs_by_projection');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.call_CDSs_by_projection",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'call_CDSs_by_projection',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method call_CDSs_by_projection",
					    status_line => $self->{client}->status_line,
					    method_name => 'call_CDSs_by_projection',
				       );
    }
}



=head2 $result = get_strep_suis_repeats(genomeTO)

Interface to Strep repeats and "boxes" tools

=cut

sub get_strep_suis_repeats
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function get_strep_suis_repeats (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to get_strep_suis_repeats:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'get_strep_suis_repeats');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.get_strep_suis_repeats",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'get_strep_suis_repeats',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method get_strep_suis_repeats",
					    status_line => $self->{client}->status_line,
					    method_name => 'get_strep_suis_repeats',
				       );
    }
}



=head2 $result = get_strep_pneumo_repeats(genomeTO)



=cut

sub get_strep_pneumo_repeats
{
    my($self, @args) = @_;

    if ((my $n = @args) != 1)
    {
	Bio::KBase::Exceptions::ArgumentValidationError->throw(error =>
							       "Invalid argument count for function get_strep_pneumo_repeats (received $n, expecting 1)");
    }
    {
	my($genomeTO) = @args;

	my @_bad_arguments;
        (ref($genomeTO) eq 'HASH') or push(@_bad_arguments, "Invalid type for argument 1 \"genomeTO\" (value was \"$genomeTO\")");
        if (@_bad_arguments) {
	    my $msg = "Invalid arguments passed to get_strep_pneumo_repeats:\n" . join("", map { "\t$_\n" } @_bad_arguments);
	    Bio::KBase::Exceptions::ArgumentValidationError->throw(error => $msg,
								   method_name => 'get_strep_pneumo_repeats');
	}
    }

    my $result = $self->{client}->call($self->{url}, {
	method => "GenomeAnnotation.get_strep_pneumo_repeats",
	params => \@args,
    });
    if ($result) {
	if ($result->is_error) {
	    Bio::KBase::Exceptions::JSONRPC->throw(error => $result->error_message,
					       code => $result->content->{code},
					       method_name => 'get_strep_pneumo_repeats',
					      );
	} else {
	    return wantarray ? @{$result->result} : $result->result->[0];
	}
    } else {
        Bio::KBase::Exceptions::HTTP->throw(error => "Error invoking method get_strep_pneumo_repeats",
					    status_line => $self->{client}->status_line,
					    method_name => 'get_strep_pneumo_repeats',
				       );
    }
}



sub version {
    my ($self) = @_;
    my $result = $self->{client}->call($self->{url}, {
        method => "GenomeAnnotation.version",
        params => [],
    });
    if ($result) {
        if ($result->is_error) {
            Bio::KBase::Exceptions::JSONRPC->throw(
                error => $result->error_message,
                code => $result->content->{code},
                method_name => 'get_strep_pneumo_repeats',
            );
        } else {
            return wantarray ? @{$result->result} : $result->result->[0];
        }
    } else {
        Bio::KBase::Exceptions::HTTP->throw(
            error => "Error invoking method get_strep_pneumo_repeats",
            status_line => $self->{client}->status_line,
            method_name => 'get_strep_pneumo_repeats',
        );
    }
}

sub _validate_version {
    my ($self) = @_;
    my $svr_version = $self->version();
    my $client_version = $VERSION;
    my ($cMajor, $cMinor) = split(/\./, $client_version);
    my ($sMajor, $sMinor) = split(/\./, $svr_version);
    if ($sMajor != $cMajor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Major version numbers differ.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor < $cMinor) {
        Bio::KBase::Exceptions::ClientServerIncompatible->throw(
            error => "Client minor version greater than Server minor version.",
            server_version => $svr_version,
            client_version => $client_version
        );
    }
    if ($sMinor > $cMinor) {
        warn "New client version available for Bio::KBase::GenomeAnnotation::Client\n";
    }
    if ($sMajor == 0) {
        warn "Bio::KBase::GenomeAnnotation::Client version is $svr_version. API subject to change.\n";
    }
}

package Bio::KBase::GenomeAnnotation::Client::RpcClient;
use base 'JSON::RPC::Client';

#
# Override JSON::RPC::Client::call because it doesn't handle error returns properly.
#

sub call {
    my ($self, $uri, $obj) = @_;
    my $result;

    if ($uri =~ /\?/) {
       $result = $self->_get($uri);
    }
    else {
        Carp::croak "not hashref." unless (ref $obj eq 'HASH');
        $result = $self->_post($uri, $obj);
    }

    my $service = $obj->{method} =~ /^system\./ if ( $obj );

    $self->status_line($result->status_line);

    if ($result->is_success) {

        return unless($result->content); # notification?

        if ($service) {
            return JSON::RPC::ServiceObject->new($result, $self->json);
        }

        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    elsif ($result->content_type eq 'application/json')
    {
        return JSON::RPC::ReturnObject->new($result, $self->json);
    }
    else {
        return;
    }
}

1;
