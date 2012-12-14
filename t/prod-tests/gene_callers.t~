use strict;
use warnings;

use Test::More tests => 11;
use Data::Dumper;

use Bio::KBase::GenomeAnnotation::Client;

# Start a server on localhost
#my ($pid, $url) = Server::start('GenomeAnnotation');
my $obj;
my $url = "http://localhost:7050";

#  Test 1 - Can a new object be created without parameters? 
$obj = Bio::KBase::GenomeAnnotation::Client->new();
ok( defined $obj, "Did an object get defined" );               

#  Test 2 - Is the object in the right class?
isa_ok( $obj, 'Bio::KBase::GenomeAnnotation::Client', "Is it in the right class" );   

# Gene calling methods that takes a genomeTO as input
my @genecall_methods = qw(
	call_selenoproteins
	call_pyrrolysoproteins
	call_RNAs
	call_CDSs
	call_CDSs_by_projection
);

#  Test 3 - Can the object do all of the methods
can_ok($obj, @genecall_methods);    

#  Test 4 - Can a new object be created with valid parameter? 
my $annotation_server = Bio::KBase::GenomeAnnotation::Client->new($url);
ok( defined $annotation_server, "Did an object get defined" );               

#  Test 5 - Is the object in the right class?
isa_ok( $annotation_server, 'Bio::KBase::GenomeAnnotation::Client', "Is it in the right class" );   

#  Test 6 - Download test data
unlink "MIT9313.genomeTO" if -e "MIT9313.genomeTO";
eval { !system("wget --quiet http://www.kbase.us/docs/build/MIT9313.genomeTO") or die $!; };
ok(!$@, "Downloaded test data");

# Create a genome typed object
my $genome_to = GenomeTO->new("MIT9313.genomeTO");
my %results;

note("Test the happy cases for gene calling methods");

#
#	all_CDSs and all_RNAs must return features for this genome
#
foreach my $method (@genecall_methods) {
	eval {$results{$method} = $annotation_server->$method($genome_to->{decode}); };
	is($@,'', "Test $method $@");
	if ($method eq 'call_RNAs' || $method eq 'call_CDSs')
	{
		if (ref($results{$method}->{'features'}) eq 'ARRAY'  )
		{
			my @ary = @{$results{$method}->{'features'}};
			isnt($#ary,0,"Non-zero number of features returned for this genome");
		}
	}
}

note("Test the unhappy cases for gene calling methods");

foreach my $method (@genecall_methods) {
	eval {$results{$method} = $annotation_server->$method($genome_to); };
	isnt($@,'', "Test $method Bad Inputs (must fail to pass the test)");
	eval {$results{$method} = $annotation_server->$method(); };
	isnt($@,'', "Test $method No Inputs (must fail to pass the test)");
}

my $empty = [];
my %bad_value = (
	'genetic_code' => 'A',
	'baddata'      => 'A',
	'domain'       => 'A',
	'contigs'      => $empty,
	'features'     => $empty,
	'id'           => 'A',
);
my %bad_struct = (
	'genetic_code' => $empty,
	'domain'       => $empty,
	'contigs'      => 'A',
	'features'     => 'A',
	'id'           => $empty,
	'scientific_name'  => $empty,
	'source'       => $empty,
	'source_id'    => $empty,
);

$genome_to->{'decode'}->{'contigs'}->[0]->{'dna'} = substr($genome_to->{'decode'}->{'contigs'}->[0]->{'dna'},0,5000);

note("Test using bad values and bad structures");
my %bad_genome_to;

foreach my $key (keys(%bad_struct))
{
#        Create a BAD genome typed object
        %bad_genome_to = %{$genome_to->{'decode'}};

#       If there is a bad value for this key, test it first (should pass)
#       Then test the bad structure (should fail)
#       Finally test for missing value (should pass)

        foreach my $method (@genecall_methods) {
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
                isnt($@,'', "Test $method with bad structure for $key (must fail to pass the test) ");
                delete $bad_genome_to{$key};
                eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                is($@,'', "Test $method with missing key $key  (must be okay to pass the test)");
        }
}


done_testing();
#Server::stop($pid);
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

gene_caller.t

=head1 DESCRIPTION

Test the following methods 

	call_selenoproteins
	call_pyrrolysoproteins
	call_RNAs
	call_CDSs
	call_CDSs_by_projection

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

=item * Test all methods with invalid IDs or empty array for the following fields, one at a time (expect no failures):

	'genetic_code' => 'A'
	'baddata'      => 'A'   (not a valid field in a genomeTO struct)
	'domain'       => 'A'
	'contigs'      => $empty
	'features'     => $empty
	'id'           => 'A'

=item * Test all methods with the wrong field type (e.g. array for string and visa versa) for the following fields, one at a time (expect failures):

	'genetic_code' => $empty
	'domain'       => $empty
	'contigs'      => 'A'
	'features'     => 'A'
	'id'           => $empty
	'scientific_name'  => $empty
	'source'       => $empty
	'source_id'    => $empty

=item * Test all methods with null values for each of the fields above, one at a time.  Expect no failures.

=back
 
=cut

