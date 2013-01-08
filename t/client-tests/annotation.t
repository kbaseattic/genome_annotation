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
	'debug'	=>	\$debug,
	'localServer'	=>	\$localServer,
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

#  Test 3 - Can the object do all of the annotation related
#           methods that take genomeTO and returns genomeTO
my @annotation_methods = qw(
        assign_functions_to_CDSs
	annotate_genome
);
can_ok($obj, @annotation_methods);    

#  Test 4 - Can a new object be created with valid parameter? 
my $annotation_server = Bio::KBase::GenomeAnnotation::Client->new($url);
ok( defined $annotation_server, "Did an object get defined" );               

#  Test 5 - Is the object in the right class?
isa_ok( $annotation_server, 'Bio::KBase::GenomeAnnotation::Client', "Is it in the right class" );   

#  Test 6 - Download test data
unlink "MIT9313.genomeTO" if -e "MIT9313.genomeTO";

my $ua = LWP::UserAgent->new();
my $res = $ua->get("http://www.kbase.us/docs/build/MIT9313.genomeTO", 
		   ":content_file" => "MIT9313.genomeTO");

ok($res->is_success, "Downloaded test data");

# Create a genome typed object
my $genome_to = GenomeTO->new("MIT9313.genomeTO");
my %results;

note("Test the happy cases for annotation methods");

# Call genes
if ($debug)
{
	ok($genome_to->{decode} = $annotation_server->call_CDSs($genome_to->{decode}),
		'test call_CDSs');
} else {
	eval { $genome_to->{decode} = $annotation_server->call_CDSs($genome_to->{decode}); };
	ok(!$@, "Test call_CDSs");
}


foreach my $method (@annotation_methods) {
	if ($debug)
	{
		ok($results{$method} = $annotation_server->$method($genome_to->{decode}),
			"Test $method");

	} else {
		eval {$results{$method} = $annotation_server->$method($genome_to->{decode}); };
		ok(!$@, "Test $method" );
	}
}

done_testing();
Server::stop($pid) if $pid;
unlink "MIT9313.genomeTO" if -e "MIT9313.genomeTO";


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
     die "could not find ./lib/$service.psgi" unless -e "lib/$service.psgi";
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

package GenomeTO;
#    typedef structure {
#       genome_id id;
#       string scientific_name;
#       string domain;
#       int genetic_code;
#       string source;
#       string source_id;
#       list<contig> contigs;
#       list<feature> features;
#    } genomeTO;
#
# From the type spec
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

# dumps the internally stored perl data structure
sub dump {
	my $self = shift;
	print Dumper $self->{'decode'};
}

# decodes the internally stored json document
sub decode {
  my $self = shift;
  $self->{decode} = JSON->new()->decode($self->{json});
}

# encode the internally stored perl data structure
sub encode {
  my $self = shift;
  $self->{json} = JSON->new()->encode($self->{decode});

}

# stores a perl data structure
sub set_decode {
  my $self = shift;
  $self->{decode} = shift;
  $self->encode();
}

# stores a json document
sub set_json {
  my $self = shift;
  $self->{json} = shift;
  $self->decode();
}

# returns a perl data structure
sub get_decode {
  $_[0]->{decode};
}

# returns the encoded jason document
sub get_encode {
  $_[0]->{encode};
}
