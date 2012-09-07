#!usr/bin/env/perl
#
# Ensembl module for Bio::EnsEMBL::Funcgen::DNAMethylationFeatureAdaptor
#

=head1 LICENSE

  Copyright (c) 1999-2012 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor - Adaptor to fetch

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -user       => 'user',
   -pass       => 'pass',
   -dbname     => 'dbname',
   -host       => 'host',
   -species    => 'species',
   -DNADB_USER => 'dnadb_user',
   -DNADB_HOST => 'dnadb_host',
   -DNADB_NAME => 'dnadb_name',
   -DNADB_PORT => dnadb_port
);



my $rsa = $efgdba->get_ResultSetAdaptor;
my @a   = @{ $rsa->fetch_all_by_name('ES_5mC_Stadler2011_PMID22170606') };

my $dnaa = $efgdba->get_DNAMethylationFeatureAdaptor;
#$dnaa->load_resultset( $a[0] );


my $slice_adaptor = $efgdba->get_adaptor("slice");
my $slice =
  $slice_adaptor->fetch_by_region( 'chromosome', 1, 3010493, 3011550 );


#my $dna_meth_features = $dnaa->get_DNAMethylationFeatures( -SLICE => $slice );


foreach my $df ( @{$dna_meth_features} ) {
    print "Methylated reads" . $df->methylated_reads . "\n";
    print "Total reads" .$df->total_reads . "\n";
    print "Percent Methylation" . $df->percent_methylation . "\n";
    print "Context" . $df->context . "\n";
    print "Display label" . $df->display_label . "\n";
    print "Cell Type" . $df->cell_type->name . "\n";
    print "Feature Type" . $df->feature_type->name . "\n";
    print "Analysis Method" . $df->analysis->logic_name . "\n";
    # . . .;
}


=head1 DESCRIPTION

The Bio::EnsEMBL::Funcgen::DNAMethylationFeatureAdaptor uses Lincoln Stein's Bio::DB::BigFile interface to
BigBed files and creates Bio::EnsEMBL::Funcgen::DNAMethylationFeature objects. For information about 
BigBed files please see http://genome.ucsc.edu/FAQ/FAQformat.html. The fourth field of the BigBed file is expected
to contain information about cytosine context and total reads. Fifth field in the file is score which is percentage methylation multiplied by 10 (ranges from 0-1000)
An example line of the bed file that represents a cytosine in CG context with 33 methylated reads out of a total of 37 reads would be as under:

chr1    3010492 3010493 CG/37   892     +


=head1 SEE ALSO

Bio::EnsEMBL::Funcgen::DNAMethylationFeature

Bio::DB::BigFile


=cut

package Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor;

use strict;
use warnings;

use Bio::DB::BigFile;
use Bio::DB::BigFile::Constants;
use Bio::EnsEMBL::Funcgen::DNAMethylationFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor);


my %strands = (
               '.' => 0,
               '+' => 1,
               '-' => -1,
              );

#currently assumes uncoverted data


sub get_DNAMethylationFeature_params_from_file_data {
  my ($self, $bed_data, $slice) = @_;

  if(! ($bed_data && ref($bed_data) ne 'ARRAY') ){
    throw('You must pass an ARRAYREF of bed data');
  }
  
  if(! ($slice && $slice->isa('Bio::EnsEMBL::SLice') )){
    throw('You must pass a Bio::EnsEMBL::Slice');
  }
  
  my $slice_start = $slice->start;
  my $slice_end   = $slice->end;
  my $total_reads;
  my ($start, $end, $context, $percent_methylation, $strand) = @$bed_data;
  $start       = $start - $slice_start;       #half open
  $end         = $end   -  $slice_start + 1 ;
  
  ( $context, $total_reads ) = split /\//, $context;
  $percent_methylation /= 10;  
  my $methylated_reads  = sprintf "%d", ( ( $percent_methylation / 100 ) * $total_reads );
  $strand = $strands{$strand} if defined $strand; #undef strand is not valid bed!

  return ( $start, $end, $strand, $methylated_reads,
           $total_reads, $percent_methylation, $context );
}


1;
