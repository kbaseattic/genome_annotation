package Bio::KBase::GenomeAnnotation::Shock;

use strict;
use REST::Client;
use LWP::UserAgent;
use JSON::XS;
use Data::Dumper;

use base 'Class::Accessor';

__PACKAGE__->mk_accessors(qw(ua server rest json auth_token));

sub new
{
    my($class, $server, $auth_token) = @_;

    my $ua = LWP::UserAgent->new();
    my $rest = REST::Client->new(useragent => $ua);
    my @auth_header;
    if ($auth_token)
    {
	@auth_header = (Authorization => "OAuth $auth_token");
	$rest->addHeader(@auth_header);
    }

    my $self = {
	server => $server,
	ua => $ua,
	rest => $rest,
	json => JSON::XS->new,
	(defined($auth_token) ? (auth_token => $auth_token) : ()),
	auth_header => \@auth_header,
    };
    return bless $self, $class;
}

sub get_file
{
    my($self, $node) = @_;

    my $res = $self->rest->GET($self->server . "/node/$node?download");
    if ($self->rest->responseCode != 200)
    {
	die "get_file failed: " . $self->rest->responseContent();
    }
    return $self->rest->responseContent();
}

sub put_file
{
    my($self, $file) = @_;
    my $url = $self->server . "/node";
    my $res = $self->ua->post($url,
			      Content_Type => 'multipart/form-data',
			      Content => [upload => [$file]]);
    if ($res->is_success)
    {
	my $ret = $self->json->decode($res->content);
	my $id = $ret->{data}->{id};
	return $id;
    }
    else
    {
	die "put_file failed: " . $res->content;
    }
}    
	    
sub put_file_data
{
    my($self, $file_data) = @_;
    my $url = $self->server . "/node";
    my $res = $self->ua->post($url,
			      @{$self->{auth_header}},
			      Content_Type => 'multipart/form-data',
			      Content_Length => length($file_data),
			      Content => [upload => [undef, 'user_data', Content => $file_data]]);
    if ($res->is_success)
    {
	my $ret = $self->json->decode($res->content);
	my $id = $ret->{data}->{id};
	# print Dumper($ret);
	return $id;
    }
    else
    {
	die "put_file failed: " . $res->content;
    }
}    
	    

1;
