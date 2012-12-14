use strict;
use warnings;

use Test::More tests => 15;
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

#  Test 3 - Can the object do all of the misc_methods
# misc_methods that take a genomeTO as input and return something else
my @misc_methods = qw(
        genomeTO_to_reconstructionTO
        genomeTO_to_feature_data
        find_close_neighbors
);
can_ok($obj, @misc_methods);    

#  Test 4 - Can a new object be created with valid parameter? 
my $annotation_server = Bio::KBase::GenomeAnnotation::Client->new($url);
ok( defined $annotation_server, "Did an object get defined" );               

#  Test 5 - Is the object in the right class?
isa_ok( $annotation_server, 'Bio::KBase::GenomeAnnotation::Client', "Is it in the right class" );   

#  Test 6 - Download test data
unlink "MIT9313.genomeTO.annotated" if -e "MIT9313.genomeTO.annotated";
eval { !system("wget --quiet http://www.kbase.us/docs/build/MIT9313.genomeTO.annotated") or die $!; };
ok(!$@, "Downloaded test data");

# Create a genome typed object
my $genome_to = GenomeTO->new("MIT9313.genomeTO.annotated");
my %results;
$genome_to->{'decode'}->{'contigs'}->[0]->{'dna'} = substr($genome_to->{'decode'}->{'contigs'}->[0]->{'dna'},0,10000);

note("Test the happy cases for misc misc_methods");

foreach my $method (@misc_methods) {

	eval {$results{$method} = $annotation_server->$method($genome_to->{decode}); };
	ok(!$@, "Test $method");
}

note("Test the unhappy cases for gene calling methods");

foreach my $method (@misc_methods) {
        eval {$results{$method} = $annotation_server->$method($genome_to); };
        isnt($@,'', "Test $method Bad Inputs (must fail to pass the test)");
        eval {$results{$method} = $annotation_server->$method(); };
        isnt($@,'', "Test $method No Inputs (must fail to pass the test)");
}


my $empty = [];
my %bad_value = (
        'contigs'      => $empty,
        'features'     => $empty,
);
my %bad_struct = (
        'contigs'      => 'A',
        'features'     => 'A',
);

$genome_to->{'decode'}->{'contigs'}->[0]->{'dna'} = substr($genome_to->{'decode'}->{'contigs'}->[0]->{'dna'},0,5000);

note("Test using bad values and bad structures");
my %bad_genome_to;

foreach my $key (keys(%bad_struct))
{
	last;  # Skip the bad structure tests because they are not documented this way
#        Create a BAD genome typed object
        %bad_genome_to = %{$genome_to->{'decode'}};

#       If there is a bad value for this key, test it first (should pass)
#       Then test the bad structure (should fail)
#       Finally test for missing value (should pass)

        foreach my $method (@misc_methods) {
                if (exists $bad_value{$key})
                {
                        note("\nTesting $key=$bad_value{$key} and method $method\n");
                        $bad_genome_to{$key} = $bad_value{$key};
                        eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                        is($@,'', "Test $method with bad value or null array for $key  (must be okay to pass the test) ");
                }
                note("\nTesting $key=$bad_struct{$key}  and method $method\n");
                $bad_genome_to{$key} = $bad_struct{$key};
                eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                is($@,'', "Test $method with bad structure for $key (must okay to pass the test) ");
                delete $bad_genome_to{$key};
                eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                is($@,'', "Test $method with missing key $key  (must be okay to pass the test)");
        }
}



done_testing();
#Server::stop($pid);
unlink "MIT9313.genomeTO.annotated";





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

sub dump {
	my $self = shift;
	print Dumper $self->{'decode'};
}


=pod

=head1 NAME

misc.t

=head1 DESCRIPTION

Test the following methods 

        genomeTO_to_reconstructionTO
        genomeTO_to_feature_data
        find_close_neighbors

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


=back

=item 3. Make unhappy tests.

=over 4

=item * Test all methods with a genomeTO which has not been decoded.  Test for failure.

=item * Test all methods with null data.  Test for failure.

=item * Test all methods with invalid IDs or empty array for the following fields, one at a time (expect no failures):

        'contigs'      => $empty,
        'features'     => $empty,

=item * Test all methods with the wrong field type (e.g. array for string and visa versa) for the following fields, one at a time (expect failures):

        'contigs'      => 'A',
        'features'     => 'A',

=item * Test all methods with null values for each of the fields above, one at a time.  Expect no failures.

=back
 
=cut


