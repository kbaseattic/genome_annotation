use strict;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;
use Bio::KBase::HandleService;
use JSON::XS;
use File::Slurp qw(read_file write_file);
use File::Temp ':POSIX';
use IO::File;
use Capture::Tiny 'capture';

@ARGV == 2 or @ARGV == 4 or die "Usage: $0 genome-workflow-file output-file [stdout stderr]\n";

our $gwfile = shift;
our $out_file = shift;

my $stdout_file = shift;
my $stderr_file = shift;

if ($stdout_file)
{
    my $stdout_fh = IO::File->new($stdout_file, "w+");
    my $stderr_fh = IO::File->new($stderr_file, "w+");

    capture(\&run, stdout => $stdout_fh, stderr => $stderr_fh);
}
else
{
    run();
}

sub run
{
    my $json = JSON::XS->new->pretty(1);

    my $impl = Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl->new();
    my $hservice = Bio::KBase::HandleService->new();

    my $ctx = ContextObj->new;
    $Bio::KBase::GenomeAnnotation::Service::CallContext = $ctx;

    open(OF, ">", $out_file) or die "Cannot open $out_file: $!";
    
    my($hobj, $wobj);
    {
	my $gtext = read_file($gwfile);
	$gtext or die "Error reading $gwfile: $!";
	my $obj = $json->decode($gtext);
	($hobj, $wobj) = @$obj;
    }
    
    print STDERR Dumper($wobj);
    
    #
    # Our genome object is a handle. We need to pull it down and parse.
    #
    
    my $tmp = tmpnam();
    print STDERR "tmp is $tmp\n";
    $hservice->download($hobj, "" . $tmp);
    
    my $gtext = read_file($tmp);
    my $gobj = $json->decode($gtext);

    my $out;
    eval {
	$out = $impl->run_pipeline($gobj, $wobj);
    };
    if ($@)
    {
	print STDERR "FAILURE running pipeline:\n$@\n";
	print OF $json->encode({failure => $@});
    }
    else
    {
	print OF $json->encode($out);
    }
    close(OF);
}

#
# Do a fairly minor emulation of the call context.
# We will need to at some point properly configure the auth stuff so that
# incoming authentication tokens (via the AWE environment) are propagated
# properly.
# 
package ContextObj;
use strict;

use base 'Class::Accessor';

__PACKAGE__->mk_accessors(qw(user_id client_ip authenticated token
                             module method call_id hostname stderr));

sub new
{
    my($class) = @_;
    my $h = `hostname`;
    chomp $h;
    my $self = { hostname => $h };
    return bless $class, $self;
}
