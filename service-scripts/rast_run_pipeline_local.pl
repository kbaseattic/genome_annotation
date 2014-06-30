use strict;
use Data::Dumper;
use Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl;
use JSON::XS;
use File::Slurp qw(read_file write_file);

my $json = JSON::XS->new->pretty(1);

@ARGV == 2 or die "Usage: $0 genome-workflow-file output-file\n";

my $gwfile = shift;
my $out_file = shift;

my $impl = Bio::KBase::GenomeAnnotation::GenomeAnnotationImpl->new();

open(OF, ">", $out_file) or die "Cannot open $out_file: $!";

my($gobj, $wobj);
{
    my $gtext = read_file($gwfile);
    $gtext or die "Error reading $gwfile: $!";
    my $obj = $json->decode($gtext);
    ($gobj, $wobj) = @$obj;
}

print STDERR Dumper($wobj);

my $out = $impl->run_pipeline($gobj, $wobj);

print OF $json->encode($out);
close(OF);


