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

my $handler = sub { $server->handle_input(@_) };

$handler;
