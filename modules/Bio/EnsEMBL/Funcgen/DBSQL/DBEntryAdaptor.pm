
=head1 NAME

Bio::EnsEMBL::Funcgen::DBSQL::DBEntryAdaptor -
MySQL Database queries to load and store external object references.

=head1 SYNOPSIS

$db_entry_adaptor = $db_adaptor->get_DBEntryAdaptor();
$db_entry = $db_entry_adaptor->fetch_by_dbID($id);

my $gene = $db_adaptor->get_GeneAdaptor->fetch_by_stable_id('ENSG00000101367');
@db_entries = @{$db_entry_adaptor->fetch_all_by_Gene($gene)};
@gene_ids = $db_entry_adaptor->list_gene_ids_by_extids('BAB15482');


=head1 CONTACT

Post questions to the EnsEMBL developer list <ensembl-dev@ebi.ac.uk>

=head1 METHODS

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::DBEntryAdaptor;

use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Bio::EnsEMBL::DBSQL::BaseAdaptor;

use Bio::EnsEMBL::DBEntry;
#use Bio::EnsEMBL::IdentityXref;
#use Bio::EnsEMBL::OntologyXref;

use Bio::EnsEMBL::Utils::Exception qw(deprecate throw warning);

use vars qw(@ISA @EXPORT);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::DBEntryAdaptor Bio::EnsEMBL::DBSQL::BaseAdaptor);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});


=head2 fetch_all_by_FeatureType


  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureType $feature_type 
               (The feature type to retrieve DBEntries for)
  Arg [2]    : optional external database name
  Arg [3]    : optional external_db type 
  Example    : @db_entries = @{$db_entry_adaptor->fetch_by_FeatureType($feature_type)};
  Description: This returns a list of DBEntries associated with this gene.
               Note that this method was changed in release 15.  Previously
               it set the DBLinks attribute of the gene passed in to contain
               all of the gene, transcript, and translation xrefs associated
               with this gene.
  Returntype : listref of Bio::EnsEMBL::DBEntries; may be of type IdentityXref if
               there is mapping data, or GoXref if there is linkage data.
  Exceptions : throws if feature type object not passed
  Caller     : Bio::EnsEMBL::Funcgen::FeatureType
  Status     : At Risk

=cut

sub fetch_all_by_FeatureType {
  my ( $self, $feature_type, $ex_db_reg, $exdb_type ) = @_;

  if(!ref($feature_type) || !$feature_type->isa('Bio::EnsEMBL::Funcgen::FeatureType')) {
    throw("Bio::EnsEMBL::Funcgen::FeatureType argument expected.");
  }

  return $self->_fetch_by_object_type($feature_type->dbID(), 'FeatureType', $ex_db_reg, $exdb_type);
}


#Change these to take arrayref of external_names?
#This would allow gene based queries to retrieve all
#associated transcript xref'd DBEntries at the same time


=head2 list_feature_type_ids_by_extid

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_feature_type_ids_by_extid('BEAF-32');
  Description: Gets a list of regulatory_feature IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_feature_type_ids_by_extid {
  my ( $self, $external_name, $external_db_name ) = @_;

  return $self->_type_by_external_id( $external_name, 'FeatureType', 
									  undef, $external_db_name );
}



=head2 list_regulatory_feature_ids_by_extid

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_regulatory_feature_ids_by_extid('GO:0004835');
  Description: Gets a list of regulatory_feature IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_regulatory_feature_ids_by_extid {
  my ( $self, $external_name, $external_db_name ) = @_;

 
  return $self->_type_by_external_id( $external_name, 'RegulatoryFeature', 
									  undef, $external_db_name );
}

=head2 list_external_feature_ids_by_extid

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_external_feature_ids_by_extid('GO:0004835');
  Description: Gets a list of external_feature IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_external_feature_ids_by_extid {
  my ( $self, $external_name, $external_db_name ) = @_;

  return
    $self->_type_by_external_id( $external_name, 'ExternalFeature', undef,
                                 $external_db_name );
}

=head2 list_annotated_feature_ids_by_extid

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_annotated_feature_ids_by_extid('GO:0004835');
  Description: Gets a list of annotated_feature IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_annotated_feature_ids_by_extid {
  my ( $self, $external_name, $external_db_name ) = @_;

  return
    $self->_type_by_external_id( $external_name, 'AnnotatedFeature', undef,
                                 $external_db_name );
}

=head2 list_probe_feature_ids_by_extid

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_annotated_feature_ids_by_extid('ENST000000000001');
  Description: Gets a list of annotated_feature IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_probe_feature_ids_by_extid {
  my ( $self, $external_name, $external_db_name ) = @_;

  return
    $self->_type_by_external_id( $external_name, 'ProbeFeature', undef,
                                 $external_db_name );
}

=head2 list_probe_ids_by_extid

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_probe_id_by_extid('ENST000000000001');
  Description: Gets a list of probe IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_probe_ids_by_extid {
  my ( $self, $external_name, $external_db_name ) = @_;

  return
    $self->_type_by_external_id( $external_name, 'Probe', undef,
                                 $external_db_name );
}


=head2 list_probeset_ids_by_extid

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_probeset_ids_by_extid('ENST000000000001');
  Description: Gets a list of probeset IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : At risk

=cut

sub list_probeset_ids_by_extid {
  my ( $self, $external_name, $external_db_name ) = @_;

  return
    $self->_type_by_external_id( $external_name, 'ProbeSet', undef,
                                 $external_db_name );
}


=head2 list_regulatory_feature_ids_by_external_db_id

  Arg [1]    : string $external_id
  Example    : @gene_ids = $dbea->list_regulatory_feature_ids_by_external_db_id(1020);
  Description: Retrieve a list of regulatory_feature ids by an external identifier that is
               linked to  any of the genes transcripts, translations or the
               gene itself. NOTE: if more than one external identifier has the
               same primary accession then genes for each of these is returned.
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_regulatory_feature_ids_by_external_db_id{
   my ($self,$external_db_id) = @_;

   my %T = map { ($_, 1) }
            $self->_type_by_external_db_id( $external_db_id, 'RegulatoryFeature' );
   return keys %T;
}





=head2 _type_by_external_id

  Arg [1]    : string $name - dbprimary_acc
  Arg [2]    : string $ensType - ensembl_object_type
  Arg [3]    : (optional) string $extraType
  Arg [4]    : (optional) string $external_db_name
  	       other object type to be returned
  Example    : $self->_type_by_external_id($name, 'regulatory_feature');
  Description: Gets
  Returntype : list of dbIDs (regulatory_feature, external_feature )
  Exceptions : none
  Caller     : list_regulatory/external_feature_ids_by_extid
  Status     : Stable

=cut

sub _type_by_external_id {
  my ( $self, $name, $ensType, $extraType, $external_db_name ) = @_;

  my $from_sql  = '';
  my $where_sql = '';
  my $ID_sql    = "oxr.ensembl_id";

  if ( defined $extraType ) {

	
	#Use the Core behavior in this case? No, because we are not on the core DB!!
	#return $self->SUPER::_type_by_external_id($name, $ensType, $extraType, $external_db_name);
	throw('Extra types not accomodated in eFG xref schema');

    if ( lc($extraType) eq 'translation' ) {
      $ID_sql = "tl.translation_id";
    } else {
      $ID_sql = "t.${extraType}_id";
    }

    if ( lc($ensType) eq 'translation' ) {
      $from_sql  = 'transcript t, translation tl, ';
      $where_sql = qq(
          t.transcript_id = tl.transcript_id AND
          tl.translation_id = oxr.ensembl_id AND
          t.is_current = 1 AND
      );
    } else {
      $from_sql  = 'transcript t, ';
      $where_sql = 't.'
        . lc($ensType)
        . '_id = oxr.ensembl_id AND '
        . 't.is_current = 1 AND ';
    }
  }

  #if ( lc($ensType) eq 'gene' ) {
  #  $from_sql  = 'gene g, ';
  #  $where_sql = 'g.gene_id = oxr.ensembl_id AND g.is_current = 1 AND ';
  #} elsif ( lc($ensType) eq 'transcript' ) {
  #  $from_sql = 'transcript t, ';
  #  $where_sql =
  #    't.transcript_id = oxr.ensembl_id AND t.is_current = 1 AND ';
  #} elsif ( lc($ensType) eq 'translation' ) {
  #  $from_sql  = 'transcript t, translation tl, ';
  #  $where_sql = qq(
  #      t.transcript_id = tl.transcript_id AND
  #      tl.translation_id = oxr.ensembl_id AND
  #      t.is_current = 1 AND
  #  );
  #}
  #


  if(lc($ensType) eq 'regulatoryfeature'){
	$from_sql  = 'regulatory_feature rf, ';
	$where_sql = qq( rf.regulatory_feature_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'externalfeature'){
	$from_sql  = 'external_feature ef, ';
	$where_sql = qq( ef.external_feature_id = oxr.ensembl_id AND );
  } 
  elsif(lc($ensType) eq 'annotatedfeature'){
	$from_sql  = 'annotated_feature af, ';
	$where_sql = qq( af.annotated_feature_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'featuretype'){
	$from_sql  = 'feature_type ft, ';
	$where_sql = qq( ft.feature_type_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'probefeature'){
	$from_sql  = 'probe_feature pf, ';
	$where_sql = qq( pf.probe_feature_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'probe'){
	$from_sql  = 'probe p, ';
	$where_sql = qq( p.probe_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'probeset'){
	$from_sql  = 'probe_set ps, ';
	$where_sql = qq( ps.probe_set_id = oxr.ensembl_id AND );
  }
  else{
	throw("ensembl_object_type $ensType is not accommodated");
  }
  

  if ( defined($external_db_name) ) {
    # Involve the 'external_db' table to limit the hits to a particular
    # external database.

    $from_sql .= 'external_db xdb, ';
    $where_sql .=
        'xdb.db_name LIKE '
      . $self->dbc()->db_handle()->quote( $external_db_name . '%' )
      . ' AND xdb.external_db_id = x.external_db_id AND';
  }

  my @queries = (
    "SELECT $ID_sql
       FROM $from_sql xref x, object_xref oxr
      WHERE $where_sql x.dbprimary_acc = ? AND
             x.xref_id = oxr.xref_id AND
             oxr.ensembl_object_type= ?",
    "SELECT $ID_sql 
       FROM $from_sql xref x, object_xref oxr
      WHERE $where_sql x.display_label = ? AND
             x.xref_id = oxr.xref_id AND
             oxr.ensembl_object_type= ?"
  );

  if ( defined $external_db_name ) {
    # If we are given the name of an external database, we need to join
    # between the 'xref' and the 'object_xref' tables on 'xref_id'.

    push @queries, "SELECT $ID_sql
       FROM $from_sql xref x, object_xref oxr, external_synonym syn
      WHERE $where_sql syn.synonym = ? AND
             x.xref_id = oxr.xref_id AND
             oxr.ensembl_object_type= ? AND
             syn.xref_id = oxr.xref_id";
  } else {
    # If we weren't given an external database name, we can get away
    # with less joins here.

    push @queries, "SELECT $ID_sql
       FROM $from_sql object_xref oxr, external_synonym syn
      WHERE $where_sql syn.synonym = ? AND
             oxr.ensembl_object_type= ? AND
             syn.xref_id = oxr.xref_id";
  }

  # Increase speed of query by splitting the OR in query into three
  # separate queries.  This is because the 'or' statments render the
  # index useless because MySQL can't use any fields in it.

  my %hash   = ();
  my @result = ();

  foreach (@queries) {
    my $sth = $self->prepare($_);
    $sth->bind_param( 1, "$name",  SQL_VARCHAR );
    $sth->bind_param( 2, $ensType, SQL_VARCHAR );
    $sth->execute();

    while ( my $r = $sth->fetchrow_array() ) {
      if ( !exists $hash{$r} ) {
        $hash{$r} = 1;
        push( @result, $r );
      }
    }
  }

  return @result;
} ## end sub _type_by_external_id

=head2 _type_by_external_db_id

  Arg [1]    : string $type - external_db type
  Arg [2]    : string $ensType - ensembl_object_type
  Arg [3]    : (optional) string $extraType
  	       other object type to be returned
  Example    : $self->_type_by_external_id(1030, 'Translation');
  Description: Gets
  Returntype : list of dbIDs (gene_id, transcript_id, etc.)
  Exceptions : none
  Caller     : list_translation_ids_by_extids
               translationids_by_extids
  	       geneids_by_extids
  Status     : Stable

=cut

sub _type_by_external_db_id{
  my ($self, $external_db_id, $ensType, $extraType) = @_;

  my $from_sql = '';
  my $where_sql = '';
  my $ID_sql = "oxr.ensembl_id";

  if (defined $extraType) {

	throw('Extra types not accomodated in eFG xref schema');

    if (lc($extraType) eq 'translation') {
      $ID_sql = "tl.translation_id";
    } else {
      $ID_sql = "t.${extraType}_id";
    }

    if (lc($ensType) eq 'translation') {
      $from_sql = 'transcript t, translation tl, ';
      $where_sql = qq(
          t.transcript_id = tl.transcript_id AND
          tl.translation_id = oxr.ensembl_id AND
          t.is_current = 1 AND
      );
    } else {
      $from_sql = 'transcript t, ';
      $where_sql = 't.'.lc($ensType).'_id = oxr.ensembl_id AND '.
          't.is_current = 1 AND ';
    }
  }

 # if (lc($ensType) eq 'gene') {
 #   $from_sql = 'gene g, ';
 #   $where_sql = 'g.gene_id = oxr.ensembl_id AND g.is_current = 1 AND ';
 # } elsif (lc($ensType) eq 'transcript') {
 #   $from_sql = 'transcript t, ';
 #   $where_sql = 't.transcript_id = oxr.ensembl_id AND t.is_current = 1 AND ';
 # } elsif (lc($ensType) eq 'translation') {
 #   $from_sql = 'transcript t, translation tl, ';
  #   $where_sql = qq(
  #       t.transcript_id = tl.transcript_id AND
  #       tl.translation_id = oxr.ensembl_id AND
  #       t.is_current = 1 AND
  #   );
 # }els

  if(lc($ensType) eq 'regulatoryfeature'){
	$from_sql  = 'regulatory_feature rf, ';
	$where_sql = qq( rf.regulatory_feature_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'externalfeature'){
	$from_sql  = 'external_feature ef, ';
	$where_sql = qq( ef.external_feature_id = oxr.ensembl_id AND );
  } 
  elsif(lc($ensType) eq 'annotatedfeature'){
	$from_sql  = 'annotated_feature af, ';
	$where_sql = qq( af.annotated_feature_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'featuretype'){
	$from_sql  = 'featuretype ft, ';
	$where_sql = qq( ft.feature_type_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'probefeature'){
	$from_sql  = 'probe_feature pf, ';
	$where_sql = qq( pf.probe_feature_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'probe'){
	$from_sql  = 'probe p, ';
	$where_sql = qq( p.probe_id = oxr.ensembl_id AND );
  }
  elsif(lc($ensType) eq 'probeset'){
	$from_sql  = 'probe_set ps, ';
	$where_sql = qq( pf.probe_set_id = oxr.ensembl_id AND );
  }
  else{
	throw("ensembl_object_type $ensType is not accommodated");
  }

  my $query = 
    "SELECT $ID_sql
       FROM $from_sql xref x, object_xref oxr
      WHERE $where_sql x.external_db_id = ? AND
  	     x.xref_id = oxr.xref_id AND oxr.ensembl_object_type= ?";

# Increase speed of query by splitting the OR in query into three separate 
# queries. This is because the 'or' statments render the index useless 
# because MySQL can't use any fields in the index.

  my %hash = ();
  my @result = ();


  my $sth = $self->prepare( $query );
  $sth->bind_param(1, "$external_db_id", SQL_VARCHAR);
  $sth->bind_param(2, $ensType, SQL_VARCHAR);
  $sth->execute();
  while( my $r = $sth->fetchrow_array() ) {
    if( !exists $hash{$r} ) {
      $hash{$r} = 1;
      push( @result, $r );
    }
  }
  return @result;
}

#Placeholders to catch error from inherited methods
#These now work in reverse as the Gene/Transcript/Translation 
#is the xref not the ensembl_object as with the core code

sub fetch_all_by_Gene {
  my ( $self, $gene) = @_;
 
  if(! (ref($gene) && $gene->isa('Bio::EnsEMBL::Gene'))) {
    throw("Bio::EnsEMBL::Gene argument expected.");
  }

  throw('Not yet implemented for eFG,  maybe you want the core DBEntryAdaptor?');
  

  #This is going to be a bit of a work around as we should really have a separate fetch method
  #fetch_all_by_external_name_object_type?
  #No!! Because this simply pulls back the xrefs, not the object xrefs!!
  #This is the same for the fetch_by_dbID method???

  #_fetch_by_external_id
  #The problem here is that we want to return ox info aswell.
  #Just rewrite _fetch_by_object_type?
}


sub fetch_all_by_Transcript {
  my ( $self, $trans) = @_;

  throw('Not implemented in eFG, maybe you want the core DBEntryAdaptor?');

  return;
}
sub fetch_all_by_Translation {
  my ( $self, $trans) = @_;
  
  throw('Not implemented in eFG, maybe you want the core DBEntryAdaptor?');
  
  return;
}


1;

