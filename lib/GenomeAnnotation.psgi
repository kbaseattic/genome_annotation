use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;

use Bio::KBase::GenomeAnnotation::Service;
use Plack::Middleware::CrossOrigin;



my @dispatch;

{
    my $obj = Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl->new;
    push(@dispatch, 'GenomeAnnotation' => $obj);
}


my $server = Bio::KBase::GenomeAnnotation::Service->new(instance_dispatch => { @dispatch },
				allow_get => 0,
			       );

my $handler = sub { $server->handle_input(@_) };

$handler = Plack::Middleware::CrossOrigin->wrap( $handler, origins => "*", headers => "*");
