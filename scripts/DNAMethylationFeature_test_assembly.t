#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

 DNAMethylationFeature_test_assembly.pl
  
=head1 SYNOPSIS

 DNAMethylationFeature_test_assembly.t [options]
 
 Assembly test for DNAMethylationFeature ResultSets

=head1 OPTIONS

perl DNAMethylationFeature_test_assembly.t --dbname name --dbpass pass --dbport port --dbhost host --dbuser user --dnadb_user dnadb_user --dnadb_host dnadb_host --dnadb_port dnadb_port --dnadb_name dnadb_name --species homo_sapiens --data_root /path/to/data/root/ --result_set_name "H1ESC_5mC_Publication_Lister2009_PMID19829295 IMR90_5mC_Publication_Lister2009_PMID19829295" --chromosomes "1 2 3 7 9" --start_position 1000 --query_slice_length 10000000


 Mandatory params:

  --dbname          efgdb database name
  --dbuser          efgdb user name
  --dbpass          efgdb password
  --dbport          efgdb port number  
  --dbhost          efgdb host
  --dnadb_user      dnadb user name
  --dnadb_host      dnadb host
  --dnadb_port      dnadb port number
  --dnadb_name      dnadb database name
  --data_root       data_root directory for ResultSet objects where DNAMethylationFeature bidbed files are hosted. 
  
 Optional params:
  --result_set_name     The user defined ResultSet names passed on as space separated list enclosed in double quotes (").
                        If not defined all the ResultSets corresponding to DNAMethylationFeatures in database will be tested.
  
  --chromosomes         The user defined list of chromosome seq_region_names passed on as space separated list enclosed in 
                        double quotes ("). If not defined all the chromosomes for the species will be tested
   
  --start_position      Position on the chromosomes where the test should begin. If not defined the last 5% length of each
                        defined chromosome will be tested
  
  
  --query_slice_length  Length of the slices in base paris to query. Default is 5 Mb.


=head1 DESCRIPTION

This program checks the assembly for DNAMethylationFeature ResultSets

=cut

use strict;
use warnings;
use Test::More;
use Bio::EnsEMBL::Funcgen::DNAMethylationFeature;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;

my (
    $dbuser,      $dbpass,          $dbname,
    $dbport,      $dbhost,          $dnadb_user,
    $dnadb_host,  $dnadb_port,      $dnadb_name,
    $species,     $result_set_name, $data_root,
    $chromosomes, $start_position,  $query_slice_length
);
my @tmp_args = @ARGV;

GetOptions(

    #DB connection
    "dbname=s" => \$dbname,
    "dbpass=s" => \$dbpass,
    "dbport=i" => \$dbport,
    "dbhost=s" => \$dbhost,
    "dbuser=s" => \$dbuser,

    #DNADB connection

    "dnadb_user=s" => \$dnadb_user,
    "dnadb_host=s" => \$dnadb_host,
    "dnadb_port=i" => \$dnadb_port,
    "dnadb_name=s" => \$dnadb_name,

    'species=s'         => \$species,
    'data_root=s'       => \$data_root,
    'result_set_name:s' => \$result_set_name,
    'chromosomes:s'     => \$chromosomes,

    'start_position:i'     => \$start_position,
    'query_slice_length:i' => \$query_slice_length

) or die("Specified invalid params: @tmp_args");

print "DNAMethylationFeature_test_assembly.t @tmp_args\n";

#obtain funcgen database adaptor

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
    -user       => $dbuser,
    -pass       => $dbpass,
    -driver     => 'mysql',
    -DNADB_USER => $dnadb_user,
    -DNADB_PORT => $dnadb_port,

    -species    => $species,
    -dbname     => $dbname,
    -host       => $dbhost,
    -DNADB_HOST => $dnadb_host,
    -DNADB_NAME => $dnadb_name,

);

########
##
#obtain Adaptors
my $rsa = $efgdba->get_adaptor("resultset");

#my $dnaa          = $efgdba->get_adaptor("DNAMethylationFeature");
my $slice_adaptor = $efgdba->get_adaptor("slice");

#set data root
$rsa->dbfile_data_root($data_root);

my ( @chromos, @result_sets );

if ($chromosomes) {
    my @defined_chrs = split /\s+/, $chromosomes;

    map {
        push( @chromos, $slice_adaptor->fetch_by_region( 'chromosome', $_ ) )
    } @defined_chrs;
}
else {
    @chromos = @{ $slice_adaptor->fetch_all('chromosome') };
}

if ( defined $result_set_name ) {
    my @result_set_names = split /\s+/, $result_set_name;
    map { push( @result_sets, @{ $rsa->fetch_all_by_name($_) } ) }
      @result_set_names;
}
else {
    @result_sets = @{
        $rsa->fetch_all_by_feature_class( 'dna_methylation',
            { status => 'DISPLAYABLE' } )
    };

}

my $slice_length = defined $query_slice_length ? $query_slice_length : 5000000;

#plan tests => 1;

my ($slice, $slice_seq, $dna_meth_features, $dnaa);

foreach my $chr (@chromos) {

    ok( $chr, 'defined chr ' . $chr->name );
    print "\nTesting :\t" . $chr->name . "\n";

    #set default 5% chromosome ends if start_position is not defined
    my $fetch_start =
      defined $start_position ? $start_position : int( 0.95 * $chr->end );

    my $fetch_end = $fetch_start + $slice_length;

    while ( $fetch_start < $fetch_end ) {

        #$slice = $chr->expand( -$fetch_start, -$chr->end + $fetch_end );
        $slice =
          $slice_adaptor->fetch_by_region( 'chromosome', $chr->seq_region_name,
            $fetch_start, $fetch_end );
        $slice_seq = $slice->seq;

        #my $slice     = $chr;

        foreach my $resultset (@result_sets) {
            $dnaa = $efgdba->get_adaptor("DNAMethylationFeature");
            $dna_meth_features =
              $dnaa->fetch_all_by_Slice_ResultSet( $slice, $resultset );

            foreach my $df ( @{$dna_meth_features} ) {

                #is($df->seq_region_name,$chr->seq_region_name);
                my $seq = substr( $slice_seq, $df->start + 1, 1 );
                if ( $seq eq 'C' ) {
                    is( $df->strand, 1,
                            'strand for '
                          . $df->set->name . "\t"
                          . $df->seq_region_name . ":"
                          . ( $df->seq_region_start + 2 )
                          . "\tstrand:"
                          . $df->strand
                          . "\tGenomic base:"
                          . $seq );
                }
                elsif ( $seq eq 'G' ) {
                    is( $df->strand, -1,
                            'strand for '
                          . $df->set->name . "\t"
                          . $df->seq_region_name . ":"
                          . ( $df->seq_region_start + 2 )
                          . "\tstrand:"
                          . $df->strand
                          . "\tGenomic base:"
                          . $seq );
                }
                else {
                    fail(   'sequence incorrect for '
                          . $df->set->name . "\t"
                          . $df->seq_region_name . ":"
                          . ( $df->seq_region_start + 2 )
                          . "\tstrand:"
                          . $df->strand
                          . "\tGenomic base:"
                          . $seq );
                }
            }
            $dna_meth_features = undef;
            $dnaa              = undef;

        }

        $fetch_start = $fetch_end;
        $fetch_end += $slice_length;
        $fetch_end = $chr->end if ( $fetch_end > $chr->end );
        $slice     = undef;
        $slice_seq = undef;

        #print STDERR "$fetch_start\t$fetch_end\n";

    }

}
