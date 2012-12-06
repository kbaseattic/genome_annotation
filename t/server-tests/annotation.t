use strict;
use warnings;

use Test::More tests => 27;
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

#  Test 3 - Can the object do all of the annotation related
#           methods that take genomeTO and returns genomeTO
my @annotation_methods = qw(
        assign_functions_to_CDSs
	annotate_genome
        annotate_proteins
);

can_ok($obj, @annotation_methods);    

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
# Create shorter sequence to make it go faster
$genome_to->{'decode'}->{'contigs'}->[0]->{'dna'} = substr($genome_to->{'decode'}->{'contigs'}->[0]->{'dna'},0,10000);
my %results;

note("Test the happy cases for annotation methods");

# Call genes
#eval { $genome_to->{decode} = $annotation_server->call_CDSs($genome_to->{decode}); };
#is($@,'', "Test call_CDSs");

foreach my $method (@annotation_methods) {
  eval {$results{$method} = $annotation_server->$method($genome_to->{decode}); };
  is($@,'', "Test $method" );
}

note("Test the unhappy cases for gene calling methods");

foreach my $method (@annotation_methods) {
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


note("Test using bad values and bad structures");
my %bad_genome_to;

foreach my $key (keys(%bad_struct))
{
#        Create a BAD genome typed object
        %bad_genome_to = %{$genome_to->{'decode'}};

#       If there is a bad value for this key, test it first (should pass)
#       Then test the bad structure (should fail)
#       Finally test for missing value (should pass)

        foreach my $method (@annotation_methods) {
                note("\nTesting $key=$bad_struct{$key}  and method $method\n");
		if (  $method eq 'annotate_genome' && $key eq 'contigs'  )
		{
                	$bad_genome_to{$key} = $bad_value{$key};
                	eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                	isnt($@,'', "Test $method with bad value or null array for $key  (must fail to pass the test) ");
                	$bad_genome_to{$key} = $bad_struct{$key};
                	eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                	isnt($@,'', "Test $method with bad structure for $key (must fail to pass the test) ");
                	delete $bad_genome_to{$key};
                	eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                	isnt($@,'', "Test $method with missing key $key  (must fail to pass the test)");
		}
		elsif ($method eq 'assign_functions_to_CDSs' && $key eq 'contigs'  )
		{
                	$bad_genome_to{$key} = $bad_value{$key};
                	eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                	is($@,'', "Test $method with bad value or null array for $key  (must be okay to pass the test) ");
                	$bad_genome_to{$key} = $bad_struct{$key};
                	eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                	is($@,'', "Test $method with bad structure for $key (must be okay to pass the test) ");
                	delete $bad_genome_to{$key};
                	eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                	is($@,'', "Test $method with missing key $key  (must be okay to pass the test)");
		}
		elsif ($key eq 'features' && $method ne 'annotate_proteins')
		{
                	$bad_genome_to{$key} = $bad_value{$key};
                	eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                	is($@,'', "Test $method with bad value or null array for $key  (must be okay to pass the test) ");
                	$bad_genome_to{$key} = $bad_struct{$key};
                	eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                	isnt($@,'', "Test $method with bad structure for $key (must fail to pass the test) ");
                	delete $bad_genome_to{$key};
			if ($method eq 'assign_functions_to_CDSs')
			{ 
                		eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                		isnt($@,'', "Test $method with missing key $key  (must fail pass the test)");
			}
			elsif ($method eq 'annotate_genome')
			{ 
                		eval {$results{$method} = $annotation_server->$method(\%bad_genome_to); };
                		is($@,'', "Test $method with missing key $key  (must be okay to pass the test)");
			}
		}
		else
		{
			print "No tests at this time\n";
		}
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

=pod

=head1 NAME

gene_caller.t

=head1 DESCRIPTION

Test the following methods 

        assign_functions_to_CDSs
        annotate_genome
        annotate_proteins

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

=item * Test all methods with invalid IDs or empty array for the following fields, one at a time:

        'contigs'      => $empty
        'features'     => $empty

	annotate_genomes will fail on contigs but not features
	assign_function_to_CDSs will be okay
	annotate_proteins - unknown

=item * Test all methods with the wrong field type (e.g. array for string and visa versa) for the following fields, one at a time (expect failures):

        'contigs'      => 'A'
        'features'     => 'A'

	annotate_genomes will fail
	assign_function_to_CDSs will fail on features but not on contigs
	annotate_proteins - unknown

=item * Test all methods with null values for each of the fields above, one at a time.  Expect no failures.

	annotate_genomes will fail on contigs but not features
	assign_function_to_CDSs will fail on features but not on contigs
	annotate_proteins - unknown

=back
 
=cut

