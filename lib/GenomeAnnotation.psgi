use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;

use Bio::KBase::GenomeAnnotation::Service;



my @dispatch;

{
    my $obj = Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl->new;
    push(@dispatch, 'GenomeAnnotation' => $obj);
}


my $server = Bio::KBase::GenomeAnnotation::Service->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub {
    my($env) = @_;
    
    my $resp = $server->handle_input($env);

    if ($env->{HTTP_ORIGIN})
    {
	my($code, $hdrs, $body) = @$resp;
	push(@$hdrs, 'Access-Control-Allow-Origin', $env->{HTTP_ORIGIN});
    }
    return $resp;
};

$handler;
