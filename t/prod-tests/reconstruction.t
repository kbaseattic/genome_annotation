use strict;
use warnings;

use Test::More tests => 14;
use Data::Dumper;
use lib "lib";
use lib "t/prod-tests";
use AnnotationTestConfig qw(getHost getPort getURL);

use Bio::KBase::GenomeAnnotation::Client;

# Start a server on localhost
#my ($pid, $url) = Server::start('GenomeAnnotation');
my $obj;
# MAKE A CONNECTION (DETERMINE THE URL TO USE BASED ON THE CONFIG MODULE)
my $host=getHost(); my $port=getPort();  my $url=getURL();
print "-> attempting to connect to:'".$url."'\n";

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
my $command = "wget --quiet http://www.kbase.us/docs/build/MIT9313.genome.annotated.reconstructionTO";
eval { !system($command) or die $!; };
ok(!$@, "Downloaded test data");
diag("unable to run $command") if $@;

# Create a genome typed object
my $genome_to = ReconstructionTO->new("MIT9313.genome.annotated.reconstructionTO");

my $genomeTO;
my $reconstructionTO;
my %results;

note("Test the happy cases for reconstruction methods");

foreach my $method (@reconstruction_methods) {

	eval {$results{$method} = $annotation_server->$method($genome_to->{decode}); };
	is($@,'', "Test $method");

	my @tmp_ary = @{$results{$method}};
	isnt($#tmp_ary,-1,"Return array is not empty");
	print "Return size = $#tmp_ary \n";
}

note("Test the unhappy cases for gene calling methods");

foreach my $method (@reconstruction_methods) {
        eval {$results{$method} = $annotation_server->$method($genome_to); };
        isnt($@,'', "Test $method Bad Inputs (must fail to pass the test)");
        eval {$results{$method} = $annotation_server->$method(); };
        isnt($@,'', "Test $method No Inputs (must fail to pass the test)");
}

my $empty_a = [];
my $empty_h = {};

note("Test using bad values and bad structures");

#        Create a BAD genome typed object
my  %bad_genome_to = %{$genome_to->{'decode'}};

#       If there is a bad value for this key, test it first (should pass)
#       Then test the bad structure (should fail)
#       Finally test for missing value (should pass)

foreach my $method (@reconstruction_methods) {
	last;  # Skip for now because structures are not enforced or documented as such
    $bad_genome_to{'bindings'} = $empty_a;
    eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
    is($@,'', "Test $method with null array (must be okay to pass the test) ");

    $bad_genome_to{'bindings'} = $empty_h;
    eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
    is($@,'', "Test $method with empty hash (must be okay to pass the test) ");

    delete $bad_genome_to{'bindings'};
    eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
    is($@,'', "Test $method with missing key bindings  (must be okay to pass the test)");
}




done_testing();
#Server::stop($pid);
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
=pod

=head1 NAME

gene_caller.t

=head1 DESCRIPTION

Test the following methods 

        reconstructionTO_to_roles
        reconstructionTO_to_subsystems

All methods take a Genome Token Object (genomeTO) as input and return a genomeTO

=head1 TEST PLAN

=over 4

=item 1. Get a valid Genome Token Object

=over 4

=item * http://www.kbase.us/docs/build/MIT9313.genomeTO

=item * decode from JSON to referenced hash

=back

=item 2. Make happy tests. 

=over 4

=item * Test all methods with the valid genomeTO

=item * For call_RNAs and call_CDSs, verify that a non-zero number of rows were returned.

=back

=item 3. Make unhappy tests.

=over 4

=item * Test all methods with a genomeTO which has not been decoded.  Test for failure.

=item * Test all methods with null data.  Test for failure.

=item * Test all methods with empty array for bindings

=item * Test all methods with the empty hash for bindings

=item * Test all methods with null values for bindings

=back
 
=cut


