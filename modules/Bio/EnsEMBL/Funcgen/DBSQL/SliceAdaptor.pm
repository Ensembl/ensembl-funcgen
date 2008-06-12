
#
# Ensembl module for Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor
#
#
# Copyright (c) 2004 Ensembl
#
# You may distribute this module under the same terms as perl itself

# POD documentation - main docs before the code

=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor - A database aware adaptor responsible for
the creation of Slices in the context of eFG objects.

=head1 SYNOPSIS

  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
 

  $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);

  $slice_adaptor = $db->get_SliceAdaptor();

  # get a slice on the entire chromosome X
  $gene_regulation_slice = $slice_adaptor->fetch_by_gene($gene);


=head1 DESCRIPTION

This module is simple wrapper class for the core SliceAdaptor, extending new
methods to generate Slices for eFG features associated with a given gene or transcript.

=head1 CONTACT

This module is part of the Ensembl project http://www.ensembl.org

For more information email <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut


package Bio::EnsEMBL::Funcgen::DBSQL::SliceAdaptor;
use vars qw(@ISA);
use strict;


use Bio::EnsEMBL::DBSQL::SliceAdaptor;
use Bio::EnsEMBL::Slice;


use Bio::EnsEMBL::Utils::Exception qw(throw deprecate warning stack_trace_dump);


@ISA = ('Bio::EnsEMBL::DBSQL::SliceAdaptor');


=head2 fetch_by_gene

  Arg [1]    : Bio::EnsEMBL::Gene - $gene
  Example    : $name  = 'chromosome:NCBI34:X:1000000:2000000:1';
               $slice = $slice_adaptor->fetch_by_name($name);
               $slice2 = $slice_adaptor->fetch_by_name($slice3->name());
  Description: Fetches a slice using a slice name (i.e. the value returned by
               the Slice::name method).  This is useful if you wish to 
               store a unique identifier for a slice in a file or database or
               pass a slice over a network.
               Slice::name allows you to serialise/marshall a slice and this
               method allows you to deserialise/unmarshal it.

               Returns undef if no seq_region with the provided name exists in
               the database.

  Returntype : Bio::EnsEMBL::Slice or undef
  Exceptions : throw if incorrent arg provided
  Caller     : Pipeline
  Status     : Stable

=cut

sub fetch_by_gene {
  my $self = shift;
  my $gene = shift;

  if(! ( ref($gene) && $gene->isa('Bio::EnsEMBL::Gene'))){
    throw("You must pass a valid Bio::EnsEMBL::Gene object");
  }
  
  #Now we want to get all the eFG feature associated with this gene and build the slice
  #resitrict to RegFeats for now.
  my $cs_name      = $gene->slice->coord_system_name;
  my $gene_chr     = $gene->seq_region_name;
  my $start        = $gene->seq_region_start;
  my $end          = $gene->seq_region_end;
  

  ($start, $end) = $self->_set_bounds_by_regulatory_feature_xref($gene_chr, $start, $end, $gene);

  #Now need to do this for all transcripts and proteins too?
  
  foreach my $trans(@{$gene->get_all_Transcripts}){
	($start, $end) = $self->_set_bounds_by_regulatory_feature_xref($gene_chr, $start, $end, $trans);

	foreach my $prot(@{$trans->get_all_Translations}){
	  ($start, $end) = $self->_set_bounds_by_regulatory_feature_xref($gene_chr, $start, $end, $prot);
	}
  }

  return $self->fetch_by_region($cs_name, $gene_chr, $start, $end);
}

sub _set_bounds_by_regulatory_feature_xref{
  my ($self, $chr, $start, $end, $obj) = @_;

  my ($regf, $extdb_name);
  my $dbe_adaptor = $self->db->get_DBEntryAdaptor;
  my $regf_adaptor = $self->db->get_RegulatoryFeatureAdaptor;

  #Do we need a central store for ensembl db names?

  if($obj->is('Bio::EnsEMBL::Gene')){
	$extdb_name = 'ensembl_core_Gene';
  }
  elsif($obj->is('Bio::EnsEMBL::Transcript')){
	$extdb_name = 'ensembl_core_Transcript';
  }
  elsif($obj->is('Bio::EnsEMBL::Translation')){
	$extdb_name = 'ensembl_core_Translation';
  }
  else{
	throw('Currently only handles Ensembl Gene, Transcript and Translation xrefs');
  }


  foreach my $rf_id($dbe_adaptor->list_regulatory_feature_ids_by_extid
					($obj->stable_id, $extdb_name)){
	
	$regf = $regf_adaptor->fetch_by_dbID($rf_id);
	
	next if $regf->seq_region_name ne $chr;
	
	$start = $regf->seq_region_start if $regf->seq_region_start < $start;
	$end   = $regf->seq_region_end   if $regf->seq_region_end   > $end;
  }

  return ($start, $end);

}


1;
