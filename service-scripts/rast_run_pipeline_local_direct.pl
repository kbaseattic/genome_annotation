use strict;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;
use Bio::KBase::HandleService;
use JSON::XS;
use File::Slurp qw(read_file write_file);
use File::Temp ':POSIX';
use IO::File;
use Capture::Tiny 'capture';
use Getopt::Long::Descriptive;

my($opt, $usage) = describe_options("%c %o [< input] [> output]",
				    ["input|i=s" => "Input file"],
				    ["output|o=s" => "Output file"],
				    ["workflow|w=s" => "Optional workflow definition"],
				    ["help|h" => "Show this help message"]);

print($usage->text), exit if $opt->help;
die($usage->text)  if (@ARGV != 0);

my $json = JSON::XS->new->pretty(1);

my $in_fh;
if ($opt->input)
{
    open($in_fh, "<", $opt->input) or die "Cannot open input file " . $opt->input . ": $!\n";
}
else
{
    $in_fh = \*STDIN;
}
my $out_fh;
if ($opt->output)
{
    open($out_fh, ">", $opt->output) or die "Cannot open output file " . $opt->output . ": $!\n";
}
else
{
    $out_fh = \*STDOUT;
}

my $ctx = ContextObj->new;
my $impl = Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl->new();
$Bio::KBase::GenomeAnnotation::Service::CallContext = $ctx;

my $wf_obj;
if ($opt->workflow)
{
    my $gtext = read_file($opt->workflow);
    $gtext or die "Error reading workflow file " . $opt->workflow . ": $!\n";
    $wf_obj = $json->decode($gtext);
}
else
{
    $wf_obj = $impl->default_workflow();
}
     
print STDERR Dumper($wobj);
    
my $gtext = read_file($in_fh);
my $gobj = $json->decode($gtext);
$gobj or die "Cannot parse input genome object\n";

my $out;
eval {
    $out = $impl->run_pipeline($gobj, $wobj);
};
if ($@)
{
    print STDERR "FAILURE running pipeline:\n$@\n";
    print $out_fh $json->encode({failure => $@});
}
else
{
    print $out_fh $json->encode($out);
}
close($out_fh);

#
# Do a fairly minor emulation of the call context.
# We will need to at some point properly configure the auth stuff so that
# incoming authentication tokens (via the AWE environment) are propagated
# properly.
# 
package ContextObj;
use strict;
use Data::Dumper;
use base 'Class::Accessor';

BEGIN {
    ContextObj->mk_accessors(qw(user_id client_ip authenticated token
				module method call_id hostname stderr));
};

sub new
{
    my($class) = @_;
    my $h = `hostname`;
    chomp $h;
    my $self = { hostname => $h };

    bless $self, $class;

    $self->module("run_pipeline");
    $self->method("unknown");

    my $stderr = ServiceStderrWrapper->new($self);
    $self->stderr($stderr);

    return $self;
}


package ServiceStderrWrapper;

use strict;
use POSIX;
use Time::HiRes 'gettimeofday';

sub new
{
    my($class, $ctx) = @_;
    my $self = {};
    my $dest = $ENV{KBRPC_ERROR_DEST};
    my $tag = $ENV{KBRPC_TAG};
    my ($t, $us) = gettimeofday();
    $us = sprintf("%06d", $us);
    my $ts = strftime("%Y-%m-%dT%H:%M:%S.${us}Z", gmtime $t);

    my $name = join(".", $ctx->module, $ctx->method, $ctx->hostname, $ts);

    if ($dest =~ m,^/,)
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

sub log
{
    my($self, $str) = @_;
    my $d = $self->{dest};
    if (ref($d) eq 'SCALAR')
    {
	$$d .= $str . "\n";
	return 1;
    }
    elsif ($d)
    {
	print $d $str . "\n";
	return 1;
    }
    return 0;
}

sub log_cmd
{
    my($self, @cmd) = @_;
    my $d = $self->{dest};
    my $str;
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
	$$d .= $str . "\n";
    }
    elsif ($d)
    {
	print $d $str . "\n";
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


