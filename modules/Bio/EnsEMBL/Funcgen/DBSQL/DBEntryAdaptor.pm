# EnsEMBL External object reference reading writing adaptor for mySQL
#
# Copyright EMBL-EBI 2001
#
# Author: Arne Stabenau
# 
# Date : 06.03.2001
#

=head1 NAME

Bio::EnsEMBL::DBSQL::DBEntryAdaptor -
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
use Bio::EnsEMBL::IdentityXref;
use Bio::EnsEMBL::GoXref;

use Bio::EnsEMBL::Utils::Exception qw(deprecate throw warning);

use vars qw(@ISA @EXPORT);
use strict;

@ISA = qw( Bio::EnsEMBL::DBSQL::DBEntryAdaptor Bio::EnsEMBL::DBSQL::BaseAdaptor);
@EXPORT = (@{$DBI::EXPORT_TAGS{'sql_types'}});

=head2 fetch_by_db_accession

  Arg [1]    : string $dbname - The name of the database which the provided
               accession is for.
  Arg [2]    : string $accession - The accesion of the external reference to
               retrieve.
  Example    : my $xref = $dbea->fetch_by_db_accession('Interpro','IPR003439');
               print $xref->description(), "\n" if($xref);
  Description: Retrieves a DBEntry (xref) via the name of the database it is
               from and its primary accession in that database. Undef is
               returned if the xref cannot be found in the database.
  Returntype : Bio::EnsEMBL::DBSQL::DBEntry
  Exceptions : thrown if arguments are incorrect
  Caller     : general, domainview
  Status     : Stable

=cut



sub fetch_by_db_accession {
  my $self = shift;
  my $dbname = shift;
  my $accession = shift;

  my $sth = $self->prepare(
   "SELECT xref.xref_id, xref.dbprimary_acc, xref.display_label,
           xref.version, xref.description,
           exDB.dbprimary_acc_linkable, exDB.display_label_linkable, exDB.priority,
           exDB.db_name, exDB.db_display_name, exDB.db_release, es.synonym,
           xref.info_type, xref.info_text, exDB.type, exDB.secondary_db_name,
           exDB.secondary_db_table
    FROM   (xref, external_db exDB)
    LEFT JOIN external_synonym es on es.xref_id = xref.xref_id
    WHERE  xref.dbprimary_acc = ?
    AND    exDB.db_name = ?
    AND    xref.external_db_id = exDB.external_db_id");

  $sth->bind_param(1,$accession,SQL_VARCHAR);
  $sth->bind_param(2,$dbname,SQL_VARCHAR);
  $sth->execute();

  if(!$sth->rows() && lc($dbname) eq 'interpro') {
    #
    # This is a minor hack that means that results still come back even
    # when a mistake was made and no interpro accessions were loaded into
    # the xref table.  This has happened in the past and had the result of
    # breaking domainview
    #
    $sth->finish();
    $sth = $self->prepare
      ("SELECT null, i.interpro_ac, i.id, null, null, 'Interpro', null, null ".
       "FROM interpro i where i.interpro_ac = ?");
    $sth->bind_param(1,$accession,SQL_VARCHAR);
    $sth->execute();
  }

  my $exDB;

  while ( my $arrayref = $sth->fetchrow_arrayref()){
    my ( $dbID, $dbprimaryId, $displayid, $version, $desc, 
	 $primary_id_linkable, $display_id_linkable, $priority, $dbname, $db_display_name,
         $release, $synonym, $info_type, $info_text, $type, $secondary_db_name,
	 $secondary_db_table) = @$arrayref;

    if(!$exDB) {
      $exDB = Bio::EnsEMBL::DBEntry->new
        ( -adaptor => $self,
          -dbID => $dbID,
          -primary_id => $dbprimaryId,
          -display_id => $displayid,
          -version => $version,
          -release => $release,
          -dbname => $dbname,
	  -primary_id_linkable => $primary_id_linkable,
	  -display_id_linkable => $display_id_linkable,
	  -priority => $priority,
	  -db_display_name=>$db_display_name,
	  -info_type => $info_type,
	  -info_text => $info_text,
	  -type => $type,
	  -secondary_db_name => $secondary_db_name,
	  -secondary_db_table => $secondary_db_table);

      $exDB->description( $desc ) if ( $desc );
    }

    $exDB->add_synonym( $synonym )  if ($synonym);
  }

  $sth->finish();

  return $exDB;
}



sub fetch_all_by_Gene {
  my ( $self, $gene) = @_;
  throw('Not implemented in eFG, maybe you want the core DBEntryAdaptor?');
}

sub fetch_all_by_Transcript {
  my ( $self, $trans) = @_;
  throw('Not implemented in eFG, maybe you want the core DBEntryAdaptor?');
}


sub fetch_all_by_Translation {
  my ( $self, $trans) = @_;
  throw('Not implemented in eFG, maybe you want the core DBEntryAdaptor?');  
}



=head2 list_gene_ids_by_external_db_id

  Arg [1]    : string $external_id
  Example    : @gene_ids = $dbea->list_gene_ids_by_external_db_id(1020);
  Description: Retrieve a list of geneid by an external identifier that is
               linked to  any of the genes transcripts, translations or the
               gene itself. NOTE: if more than one external identifier has the
               same primary accession then genes for each of these is returned.
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_gene_ids_by_external_db_id{
   my ($self,$external_db_id) = @_;

throw('Not implemented in eFG, maybe you want the core DBEntryAdaptor?');



   my %T = map { ($_, 1) }
       $self->_type_by_external_db_id( $external_db_id, 'Translation', 'gene' ),
       $self->_type_by_external_db_id( $external_db_id, 'Transcript',  'gene' ),
       $self->_type_by_external_db_id( $external_db_id, 'Gene' );
   return keys %T;
}

=head2 list_gene_ids_by_extids

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @gene_ids = $dbea->list_gene_ids_by_extids('ARSE');
  Description: Retrieve a list of geneid by an external identifier that is 
               linked to  any of the genes transcripts, translations or the 
               gene itself 
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_gene_ids_by_extids {
  my ( $self, $external_name, $external_db_name ) = @_;

throw('Not implemented in eFG, maybe you want the core DBEntryAdaptor?');



  my %T = map { ( $_, 1 ) }
    $self->_type_by_external_id( $external_name, 'Translation', 'gene',
                                 $external_db_name ),
    $self->_type_by_external_id( $external_name, 'Transcript', 'gene',
                                 $external_db_name ),
    $self->_type_by_external_id( $external_name, 'Gene', undef,
                                 $external_db_name );

  return keys %T;
}


=head2 list_transcript_ids_by_extids

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_gene_ids_by_extids('BCRA2');
  Description: Retrieve a list transcript ids by an external identifier that 
               is linked to any of the genes transcripts, translations or the 
               gene itself 
  Returntype : list of ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_transcript_ids_by_extids {
  my ( $self, $external_name, $external_db_name ) = @_;

throw('Not implemented in eFG, maybe you want the core DBEntryAdaptor?');



  my %T = map { ( $_, 1 ) }
    $self->_type_by_external_id( $external_name, 'Translation',
                                 'transcript',   $external_db_name
    ),
    $self->_type_by_external_id( $external_name, 'Transcript', undef,
                                 $external_db_name );

  return keys %T;
}


=head2 list_translation_ids_by_extids

  Arg [1]    : string $external_name
  Arg [2]    : (optional) string $external_db_name
  Example    : @tr_ids = $dbea->list_gene_ids_by_extids('GO:0004835');
  Description: Gets a list of translation IDs by external display IDs
  Returntype : list of Ints
  Exceptions : none
  Caller     : unknown
  Status     : Stable

=cut

sub list_translation_ids_by_extids {
  my ( $self, $external_name, $external_db_name ) = @_;

  throw('Not implemented in eFG, maybe you want the core DBEntryAdaptor?');


  return
    $self->_type_by_external_id( $external_name, 'Translation', undef,
                                 $external_db_name );
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

  warn "in list with $external_name, $external_db_name ";

  return $self->_type_by_external_id( $external_name, 'regulatory_feature', 
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
    $self->_type_by_external_id( $external_name, 'external_feature', undef,
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
            $self->_type_by_external_db_id( $external_db_id, 'regulatory_feature' );
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

  if ( lc($ensType) eq 'gene' ) {
    $from_sql  = 'gene g, ';
    $where_sql = 'g.gene_id = oxr.ensembl_id AND g.is_current = 1 AND ';
  } elsif ( lc($ensType) eq 'transcript' ) {
    $from_sql = 'transcript t, ';
    $where_sql =
      't.transcript_id = oxr.ensembl_id AND t.is_current = 1 AND ';
  } elsif ( lc($ensType) eq 'translation' ) {
    $from_sql  = 'transcript t, translation tl, ';
    $where_sql = qq(
        t.transcript_id = tl.transcript_id AND
        tl.translation_id = oxr.ensembl_id AND
        t.is_current = 1 AND
    );
  }
  elsif(lc($ensType) eq 'regulatory_feature'){
	$from_sql  = 'regulatory_feature rf, ';
	$where_sql = qq( rf.regulatory_feature_id = oxr.ensembl_id AND );
	#rf.is_current = 1 AND );
  }
  elsif(lc($ensType) eq 'external_feature'){
	$from_sql  = 'external_feature ef, ';
	$where_sql = qq( ef.external_feature_id = oxr.ensembl_id AND );
	#rf.is_current = 1 AND );
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

  if (lc($ensType) eq 'gene') {
    $from_sql = 'gene g, ';
    $where_sql = 'g.gene_id = oxr.ensembl_id AND g.is_current = 1 AND ';
  } elsif (lc($ensType) eq 'transcript') {
    $from_sql = 'transcript t, ';
    $where_sql = 't.transcript_id = oxr.ensembl_id AND t.is_current = 1 AND ';
  } elsif (lc($ensType) eq 'translation') {
    $from_sql = 'transcript t, translation tl, ';
    $where_sql = qq(
        t.transcript_id = tl.transcript_id AND
        tl.translation_id = oxr.ensembl_id AND
        t.is_current = 1 AND
    );
  }elsif(lc($ensType) eq 'regulatory_feature'){
	$from_sql  = 'regulatory_feature rf, ';
	$where_sql = qq( rf.regulatory_feature_id = oxr.ensembl_id AND );
	#rf.is_current = 1 AND );
  }
  elsif(lc($ensType) eq 'external_feature'){
	$from_sql  = 'external_feature ef, ';
	$where_sql = qq( ef.external_feature_id = oxr.ensembl_id AND );
	#rf.is_current = 1 AND );
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


1;

