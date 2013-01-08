use strict;
use warnings;

use Test::More;
use Data::Dumper;
use Getopt::Long;
use LWP::UserAgent;

use Bio::KBase::GenomeAnnotation::Client;

my $debug=0;
my $localServer=0;
my $getoptResult=GetOptions(
        'debug' =>      \$debug,
        'localServer'   =>      \$localServer,
);

my ($url,$pid);
$url='http://localhost:7050' unless ($localServer);
# Start a server on localhost if desired
($pid, $url) = Server::start('GenomeAnnotation') unless ($url);
my $obj;

#  Test 1 - Can a new object be created without parameters? 
$obj = Bio::KBase::GenomeAnnotation::Client->new();
ok( defined $obj, "Did an object get defined" );               

#  Test 2 - Is the object in the right class?
isa_ok( $obj, 'Bio::KBase::GenomeAnnotation::Client', "Is it in the right class" );   

#  Test 3 - Can the object do all of the methods that take
#           reconstruction typed objects as input.
my @reconstruction_methods = qw(
        reconstructionTO_to_roles
        reconstructionTO_to_subsystems
);
can_ok($obj, @reconstruction_methods);    

#  Test 4 - Can a new object be created with valid parameter? 
my $annotation_server = Bio::KBase::GenomeAnnotation::Client->new($url);
ok( defined $annotation_server, "Did an object get defined" );               

#  Test 5 - Is the object in the right class?
isa_ok( $annotation_server, 'Bio::KBase::GenomeAnnotation::Client', "Is it in the right class" );   

#  Test 6 - Download test data
unlink "MIT9313.genome.annotated.reconstructionTO" if -e "MIT9313.genome.annotated.reconstructionTO";

my $ua = LWP::UserAgent->new();
my $res = $ua->get("http://www.kbase.us/docs/build/MIT9313.genome.annotated.reconstructionTO",
                   ":content_file" => "MIT9313.genome.annotated.reconstructionTO");

ok($res->is_success, "Downloaded test data");

# Create a genome typed object
my $genome_to = ReconstructionTO->new("MIT9313.genome.annotated.reconstructionTO");

my $genomeTO;
my $reconstructionTO;
my %results;

note("Test the happy cases for reconstruction methods");

foreach my $method (@reconstruction_methods) {
	if ($debug)
	{
		ok($results{$method} = $annotation_server->$method($genome_to->{decode}),
			"Test $method");
	} else {
		eval {$results{$method} = $annotation_server->$method($genome_to->{decode}); };
		ok(!$@, "Test $method");
	}
}

done_testing();
Server::stop($pid) if ($pid);
unlink "MIT9313.genome.annotated.reconstructionTO";




# Helper packages Server and GenomeTO
package Server;
use Plack::Runner;
use IO::Socket::INET;

sub start {
  my $service = shift;

  # Find an unused port.
  my $port;
  {
     my $sock = IO::Socket::INET->new('LocalPort' => 0);
     $sock->listen;
     $port = $sock->sockport;
     $sock->close();
  }

  # Fork and create service.
  my $child_pid = fork;
  if ($child_pid == 0)
  {
     die "can not find lib/$service.psgi" unless -e "lib/$service.psgi";
     open STDOUT, "/dev/null";
     open STDERR, "/dev/null";
     my $runner = Plack::Runner->new();
     $runner->parse_options("--listen", "0:$port");
     $runner->run("lib/$service.psgi");
     exit;
  }

  # Wait for server to start.
  sleep 5;
  return ($child_pid, "http://localhost:$port");
}

sub stop {
  my($child_pid, $url) = shift;
  kill 1, $child_pid;
}

package ReconstructionTO;
use JSON -support_by_pp;
use Data::Dumper;

sub new {
        my $class = shift;
        my $self = {};
	my $arg = shift;
	if (-e $arg) {
		my $json;
		open TO, "$arg";
		$json.=$_ while(<TO>);
		close TO;
		$self->{json} = $json;
	}
	else {
		$self->{json} = $arg;
	}
        $self->{'decode'} = JSON->new()->decode($self->{json});
        bless $self, $class;
}

sub dump {
	my $self = shift;
	print Dumper $self->{'decode'};
}

