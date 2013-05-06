
=head1 NAME                                                                                                                                

genomeTO_to_feature_data                                                                                                                                 
=head1 SYNOPSIS                                                                                                                         

genomeTO_to_feature_data [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> reconstruction-file                                                                                                                           
=head1 DESCRIPTION                                                                                                                   

detailed_description_of_purpose_of_script                                                                                    

Example:                                                                                                                                   

   example_of_use                                                                                                                         

example_description                                                                                                                       

=head1 COMMAND-LINE OPTIONS                                                                                              

Usage: genomeTO_to_feature_data [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> reconstruction-file]                                                                                                               
=head1 AUTHORS                                                                                                                         

L<The SEED Project|http://www.theseed.org>

=cut

use Carp;
use strict;
use Bio::KBase::GenomeAnnotation::Client;
use JSON::XS;

use Getopt::Long;
use Data::Dumper;

my $input_file;
my $output_file;
#my $url = "http://bio-data-1.mcs.anl.gov/services/genome_annotation";
my $url = "http://kbase.us/services/genome_annotation";


$| = 1;

my $help;
my $rc = GetOptions('url=s'     => \$url,
		    'input=s' 	=> \$input_file,
		    'output=s'  => \$output_file,
                    'help'      => \$help,
		    );

if (!$rc || $help || @ARGV != 0) {
       seek(DATA, 0, 0);
      while (<DATA>) {
             last if /^=head1 COMMAND-LINE /;
      }
      while (<DATA>) {
          last if (/^=/);
           print $_;
       }
      exit($help ? 0 : 1);
}


my $usage = "genomeTO_to_feature_data [--input genome-file] [--output genome-file] [--url service-url] [< genome-file] [> reconstruction-file]";

@ARGV == 0 or die "Usage: $usage\n";

my $anno_server = Bio::KBase::GenomeAnnotation::Client->new($url);

my $in_fh;
if ($input_file)
{
    open($in_fh, "<", $input_file) or die "Cannot open $input_file: $!";
}
else
{
    $in_fh = \*STDIN;
}

my $out_fh;
if ($output_file)
{
    open($out_fh, ">", $output_file) or die "Cannot open $output_file: $!";
}
else
{
    $out_fh = \*STDOUT;
}
my $json = JSON::XS->new;

my $input_genome;
{
    local $/;
    undef $/;
    my $input_genome_txt = <$in_fh>;
    $input_genome = $json->decode($input_genome_txt);
}

my $feature_data = $anno_server->genomeTO_to_feature_data($input_genome);
foreach $_ (@$feature_data)
{
    print $out_fh join("\t",@$_),"\n";
}

__DATA__
