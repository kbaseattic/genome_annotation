package main;

my $string;
open TO, "MIT9313.genomeTO";
$string.=$_ while(<TO>);
close TO;

my $gto = GenomeTO->new($string);
$gto->dump();

package GenomeTO;
#    typedef structure {
#       genome_id id;
#       string scientific_name;
#       string domain;
#       int genetic_code;
#       string source;
#       string source_id;
#       
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
        $self->{string} = shift;
        $self->{'decode'} = JSON->new()->decode($self->{string});
        bless $self, $class;
}

sub dump {
	my $self = shift;
	print Dumper $self->{'decode'};
}



1;
