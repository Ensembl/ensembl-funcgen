#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor
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
Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor
Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFeatureAdaptor;

use strict;
use warnings;

use Bio::EnsEMBL::Funcgen::DNAMethylationFeature;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor;
use Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::DBSQL::BaseFeatureAdaptor); #Should also inherit from Funcgen BigBedAdaptor

#This 'mix-in' inheritance should be well defined, similar to a Moose 'has' (role) style relation ship.
#It should provide a well defined interface, which does not corrupt the interface of the inheriting class
#Any method naming clashes as these will never be executed in the mix-in, for example, should not contain 
#a new method.  But this behaviour could be employed to redefine a mix-in method in the inheritor

# To Do
#
# 1 Add file test and eval to Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor::bigbed_open
#
# 2 Migrate fetch_features (with coderef support) back to Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor
#
# 3 Create funcgen BigBedAdaptor which caches multiple Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptors

# Caveats
#
# This implementation adds extra dependancies on:
#    Bio::EnsEMBL::ExternalData
#     which in turn depends on
#    Bio::EnsEMBL::Web !!!
#    Bio::DB::BigFile



#All calls in this method should be generic wrt file format

sub fetch_all_by_Slice_ResultSet{
  my ($self, $slice, $rset, $params) = @_;

  #calling this on self right now
  #but fetch_features should probably move back to Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor
  
  return $self->fetch_features($slice->seq_region_name,
                               $slice->start,
                               $slice->end,
                               $rset->dbfile_data_dir,
                               $self->can('_create_DNAMethylationFeature_from_file_data'),#return a code ref
                               $self,
                               #list ref rather than hash ref
                               #as we deref as array on calling 
                               #$constructor_ref method
                               [
                                -slice => $slice,
                                -set   => $rset
                               ]
                              );
}




#All the following methods should ultimately be moved elsewhere

#Move the following to funcgen BigBedAdaptor
#What should the namespace of this be in the funcgen API?
#Bio::EnsEMBL::Funcgen::ExternalData::BigFile::BigBedAdaptor
#Bio::EnsEMBL::Funcgen::DBFile::BigFile::BigBedAdaptor
#which caches and delegates to separate Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptors

my %strands = (
               '.' => 0,
               '+' => 1,
               '-' => -1,
              );

sub get_DNAMethylationFeature_params_from_file_data {
  my ($self, $bed_data, $slice) = @_;

  if(! (defined $bed_data && ref($bed_data) ne 'ARRAYREF') ){
    throw("bed_data argument in not an ARRAYREF:\t$bed_data");
  }
  
  if(! (defined $slice && $slice->isa('Bio::EnsEMBL::Slice') )){
    throw("slice argument is not a valid Bio::EnsEMBL::Slice:\t$slice");
  }
  
  my $slice_start = $slice->start;
  my $slice_end   = $slice->end;
  my $total_reads;
  my ($chr, $start, $end, $context, $percent_methylation, $strand) = @$bed_data;
  #$chr is actually never used in the Ensembl code as we already have the query slice

  $start       = $start - $slice_start;       #half open
  $end         = $end   -  $slice_start + 1 ;
  
  ( $context, $total_reads ) = split /\//, $context;
  $percent_methylation /= 10;  
  my $methylated_reads  = sprintf "%d", ( ( $percent_methylation / 100 ) * $total_reads );
  $strand = $strands{$strand} if defined $strand; #undef strand is not valid bed!

  return ( $chr, $start, $end, $strand, $methylated_reads,
           $total_reads, $percent_methylation, $context );
}


# This will ultimately move to funcgen BigBedAdaptor
# This should use Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor


#This one should probably be in a base class

sub get_file_adaptor{
  my ($self, $path) = @_;

  #This doesn't need to be generically named as it will always be called from
  #funcgen BigBedAdaptor, not a feature adaptor


  #This is similar to some to core FileAdaptor::get_filehandle code

  #my $bb_path = $rset->dbfile_data_dir;
  #my $rset_id = $rset->dbID;


  if(! exists $self->{result_set_file_adaptors}{$path}){

    #Removed a lot of this as it should be in bigbed_open

    #my $bb_path = $rset->dbfile_data_dir;

    #unless ( ( -e $bb_path ) && 
    #         ( -r $bb_path ) 
    #       ) {
    #  throw("The ResultSet dbfile ' $bb_path ' is not a valid path or the file is not readable");
    #}

    #eval to catch croak from BigFile
    #eval { $self->{result_set_file_adaptors}{$rset_id} 
    #         = Bio::DB::BigFile->bigBedFileOpen($bb_path) };

    #if($@){
    #  throw("Failed to open bigBedFile:\t$bb_path\n$@");
    #}

    $self->{result_set_file_adaptors}{$path} = 
      Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor->new($path);

  }
  
  if(! defined $self->{result_set_file_adaptors}{$path}){
    #only warn once!
    warn "Encountered error when opening bigBedFile\t".$path;
  }
  
  return $self->{result_set_file_adaptors}{$path};
}


sub _create_DNAMethylationFeature_from_file_data{
  my ($self, $bed_data, @args) = @_;
  
  #Would do read count filtering here if required?
 
  return Bio::EnsEMBL::Funcgen::DNAMethylationFeature->new
    (
     -FILE_DATA => $bed_data, 
     -ADAPTOR   => $self, 
     @args);# args are slice and result_set
}



#This should be completely agnostic to calling adaptor
#And should really be migrated back to Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor
#as it is generic due to the passing of a coderef

sub fetch_features {
  my ($self, $chr_id, $start, $end, $bb_path, $constructor_ref, $constructor_obj, $constructor_args) = @_;

  my $bb_adaptor = $self->get_file_adaptor($bb_path);
  my $bb         = $bb_adaptor->bigbed_open;

  warn "Failed to open BigBed file" . $bb_path unless $bb;
  return [] unless $bb;
  
  #  Maybe need to add 'chr' 
  my $seq_id = $bb_adaptor->munge_chr_id($chr_id);
  
  return [] if ! defined $seq_id;

  # Remember this method takes half-open coords (subtract 1 from start)
  my $list_head = $bb->bigBedIntervalQuery($seq_id, $start-1, $end);
  my @features;
 

  for ( my $i=$list_head->head; $i; $i=$i->next ) {

    #my @bedline = ($chr_id,$i->start,$i->end,split(/\t/,$i->rest));
    #my $bed = EnsEMBL::Web::Text::Feature::BED->new(\@bedline);
    #$bed->coords([$chr_id,$i->start,$i->end]);
    ### Set score to undef if missing to distinguish it from a genuine present but zero score
    #$bed->score(undef) if @bedline < 5;

    
    ##This assumes the coderef is an object method and not a sub, 
    #we could easily reset $constructor_args[$bed_data_index]
    #and set bed_data_index to 0 or 1 for sub or method respectively
    #with $constructor_args[0] set to $constructor for method based refs. 
    
    #Passing $chr_id is redundant as we already know this from the query slice

    my $feature = $constructor_ref->($constructor_obj, 
                                     [ $chr_id, $i->start, $i->end, split(/\t/,$i->rest) ], #bed data
                                     @{$constructor_args} ); #other args e.g. {-slice=> $slice, , -set => $rset}

    #do we need to handle under score here?


    #This defined check is to allow the constructor to return undef
    #if there is any sort of filtering going on
    push @features, $feature if defined $feature;
  }
  
  return \@features;
}




1;
