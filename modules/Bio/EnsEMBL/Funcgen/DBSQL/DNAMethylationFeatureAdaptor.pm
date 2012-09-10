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
#Only needs this for non-DB based methods, e.g. projection, new_fast etc.
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


=head2 fetch_all_by_Slice_ResultSet

  Arg[1]     : Bio::EnsEMBL::Slice
  Arg[2]     : Bio::EnsEMBL::Funcgen::ResultSet
  Arg[3]     : HASHREF (optional) - contraint params e.g.
                 {
                  min_read_depth => 5,
                  context        => 'CG',
                 }
              
  Example    : my @dna_meth_feats = @{$dna_mf_adaptor->fetch_all_by_Slice_ResultSet($slice, $result_set)};
  Description: Fetches DNAMethylationFeatures for a given Slice and ResultSet
  Returntype : ARRAYREF of Bio::EnsEMBL::Funcgen::DNAMethylationFeature objects
  Exceptions : Throws if args not valid and if ResultSet feature class is not dna_methylation
  Caller     : General
  Status     : Stable

=cut

#All calls in this method should be generic wrt file format

sub fetch_all_by_Slice_ResultSet{
  my ($self, $slice, $rset, $params) = @_;

  #Validate contraints?
  #context CG, CHG or CZ is valid?
  #i.e. support ambiguity codes? and throw if invalid
  #_validate_$contraint

  # VALIDATE ARGS
  if(! (defined $slice && $slice->isa('Bio::EnsEMBL::Slice') )){
    throw("Slice argument is not a valid Bio::EnsEMBL::Slice:\t$slice");
  }
  
  if(! (defined $rset && $rset->isa('Bio::EnsEMBL::Funcgen::ResultSet') )){
    #Doesn't need to be stored
    throw("ResultSet argument is not a valid Bio::EnsEMBL::Funcgen::ResultSet:\t$rset");
  }

  if($rset->feature_class ne 'dna_methylation'){
    throw("ResultSet ".$rset->name." feature class is not dna_methylation:\t".$rset->feature_class);
  }

  my $constructor_wrapper = (defined $params) ? 
    '_create_DNAMethylationFeature_from_file_data_constraints' :
      '_create_DNAMethylationFeature_from_file_data';
  
  #calling fetch_features on self but should move back to
  #Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor
  
  return $self->fetch_features
    (
     $slice->seq_region_name,
     $slice->start,
     $slice->end,
     $rset->dbfile_data_dir,
     $self->can($constructor_wrapper),#return a code ref
     #Constructor wrapper args
     $self,
     {
      new_args => {
                   slice => $slice,
                   set   => $rset
                  },
      constraints => $params,
      #defined contraints like this prevent ne_'new_args' when constraining on each feature
     }
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
  
  #Most likely already validated in fetch method
  #Remove here for speed?
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

  return [ $chr, $start, $end, $strand, $methylated_reads,
           $total_reads, $percent_methylation, $context ];
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
    warn "Encountered error when opening bigBedFile\t".$path;
  }
  
  return $self->{result_set_file_adaptors}{$path};
}


sub _create_DNAMethylationFeature_from_file_data{
  my ($self, $bed_data, $args) = @_;
    
  $bed_data = $self->get_DNAMethylationFeature_params_from_file_data
    (
     $bed_data,
     $args->{new_args}{slice}
    );

  return Bio::EnsEMBL::Funcgen::DNAMethylationFeature->new_fast
    ({
      start               => $bed_data->[1],
      end                 => $bed_data->[2],
      strand              => $bed_data->[3],
      methylated_reads    => $bed_data->[4],
      percent_methylation => $bed_data->[5],
      total_reads         => $bed_data->[6],
      context             => $bed_data->[7],
      adaptor             => $self,
      %{$args->{new_args}}                     # args are slice and result_set
     });
}


sub _create_DNAMethylationFeature_from_file_data_constraints{
  my ($self, $bed_data, $args) = @_;
  
  #loop through the contraints
  #calling methods based on the contraint key, passing the constraint value.
  #This will mean we won't have to do this check here for every single feature if we don't want to filter
  #just once in the caller!

  #$bed_data is not yet processed, so get_DNAMethylationFeature_params_from_file_data
  #and pass attrs directly here

  $bed_data = $self->get_DNAMethylationFeature_params_from_file_data
    (
     $bed_data,
     $args->{new_args}{slice}
    );
  

  my $constrain_method;
  
  foreach my $constraint(keys %{$args->{constraints}} ){
    $constrain_method = '_constrain_'.$constraint;
    
    if( $self->$constrain_method($args->{constraints}{$constraint}, $bed_data) ){
      warn "constraining $constraint";
      return;
    }
  }
  
  return Bio::EnsEMBL::Funcgen::DNAMethylationFeature->new_fast
    ({
      start               => $bed_data->[1],
      end                 => $bed_data->[2],
      strand              => $bed_data->[3],
      methylated_reads    => $bed_data->[4],
      percent_methylation => $bed_data->[5],
      total_reads         => $bed_data->[6],
      context             => $bed_data->[7],
      adaptor             => $self,
      %{$args->{new_args}}                     # args are slice and result_set
     });
}

#Have to be class methods rather than subs
#when using strict refs. subs would probably be a little faster

sub _constrain_min_read_depth{
  my ($self, $min_read_depth, $bed_data) = @_;
  return ($bed_data->[6] < $min_read_depth) ? 1 : 0;
}

sub _constrain_context{
  my ($self, $context, $bed_data) = @_;
  return ($bed_data->[7] ne $context) ? 1 : 0;
}

sub _constrain_percent_methylation{
  my ($self, $perc, $bed_data) = @_;
  return ($bed_data->[5] < $perc) ? 1 : 0;
}





#This should be completely agnostic to calling adaptor
#And should really be migrated back to Bio::EnsEMBL::ExternalData::BigFile::BigBedAdaptor
#as it is generic due to the passing of a coderef

#Will still need a wrapper method in the Funcgen BigBed adaptor
#which should probably be renamed and standardised across all
#file adaptor wrapper

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
                                     $constructor_args ); 
    #other for use in constructor wrapper e.g. {new_args => {-slice=> $slice, -set => $rset}}
    
    #Do we have to maintain constructor arg as array?
    #This is to support direct constructor(new) usage
    #rather than need for wrapper method
    #Will always need a wrapper method to translate the bed data into the appropriate constructor params


    #do we need to handle under score here?


    #This defined check is to allow the constructor to return undef
    #if there is any sort of filtering going on
    push @features, $feature if defined $feature;
  }
  
  return \@features;
}




1;
