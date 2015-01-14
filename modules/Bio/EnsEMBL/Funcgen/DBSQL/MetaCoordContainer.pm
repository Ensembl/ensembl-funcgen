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

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::MetaCoordContainer;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Exception;

use base qw(Bio::EnsEMBL::DBSQL::BaseAdaptor);

#Can remove this whole class from the API and use the core class if we can resolve problems below
#remove new and inherit from Bio::EnsEMBL::DBSQL::MetaCoordContainer?

sub new {
  my $class = shift;

  my $self = $class->SUPER::new(@_);

  #
  # Retrieve the list of the coordinate systems that features are stored in
  # and cache them
  #

	my $sql = 'SELECT mc.table_name, mc.coord_system_id, mc.max_length  FROM meta_coord mc';
	my @args;
	if($self->is_multispecies()) {
		$sql .= ' join coord_system cs using (coord_system_id) where cs.species_id =?';
		push(@args, $self->species_id());
	}
  
  my $sth = $self->prepare($sql);
  $sth->execute(@args);

  while(my ($table_name, $cs_id, $max_length) = $sth->fetchrow_array()) {
    $self->{'_feature_cache'}->{lc($table_name)} ||= [];
    push @{$self->{'_feature_cache'}->{lc($table_name)}}, $cs_id;
    $self->{'_max_len_cache'}->{$cs_id}->{lc($table_name)} = $max_length;
  }
  $sth->finish();

  return $self;
}




=head2 fetch_all_CoordSystems_by_feature_type

  Arg [1]    : string $table - the name of the table to retrieve coord systems
               for.  E.g. 'gene', 'exon', 'dna_align_feature'
  Example    : @css = @{$mcc->fetch_all_CoordSystems_by_feature_type('gene')};
  Description: This retrieves the list of coordinate systems that features
               in a particular table are stored.  It is used internally by
               the API to perform queries to these tables and to ensure that
               features are only stored in appropriate coordinate systems.
  Returntype : listref of Bio::EnsEMBL::Funcgen::CoordSystem objects
  Exceptions : throw if name argument not provided
  Caller     : BaseFeatureAdaptor
  Status     : At risk

=cut

# can remove this if we can get get_CoordSystemAdaptor to return Funcgen rather than core

sub fetch_all_CoordSystems_by_feature_type {
  my $self = shift;
  my $table = lc(shift); #case insensitive matching

  throw('Name argument is required') unless $table;

  if(!$self->{'_feature_cache'}->{$table}) {
    return [];
  }

  my @cs_ids = @{$self->{'_feature_cache'}->{$table}};
  my @coord_systems;

  my $csa = $self->db->get_FGCoordSystemAdaptor();

  foreach my $cs_id (@cs_ids) {
    my $cs = $csa->fetch_by_dbID($cs_id);

    if(!$cs) {
      throw("meta_coord table refers to non-existant coord_system $cs_id");
    }

    push @coord_systems, $cs;
  }

  return \@coord_systems;
}



=head2 fetch_max_length_by_CoordSystem_feature_type

  Arg [1]    : Bio::EnsEMBL::Funcgen::CoordSystem $cs
  Arg [2]    : string $table
  Example    : $max_len = 
                $mcc->fetch_max_length_by_CoordSystem_feature_type($cs,'gene');
  Description: Returns the maximum length of features of a given type in
               a given coordinate system.
  Returntype : int or undef
  Exceptions : throw on incorrect argument
  Caller     : BaseFeatureAdaptor
  Status     : At risk

=cut

#can remove this if we can get Funcgen::Coordsystem to inherit from core CoordSystem

sub fetch_max_length_by_CoordSystem_feature_type {
  my $self = shift;
  my $cs = shift;
  my $table = shift;

  if(!ref($cs) || !$cs->isa('Bio::EnsEMBL::Funcgen::CoordSystem')) {
    throw('Bio::EnsEMBL::Funcgen::CoordSystem argument expected');
  }

  throw("Table name argument is required") unless $table;

  return $self->{'_max_len_cache'}->{$cs->dbID()}->{lc($table)};
}



=head2 add_feature_type

  Arg [1]    : Bio::EnsEMBL::Funcgen::CoordSystem $cs
               The coordinate system to associate with a feature table
  Arg [2]    : string $table - the name of the table in which features of
               a given coordinate system will be stored in
  Arg [3]    : int $length
               This length is used to update the max_length in the database
               and the internal cache. 
  Example    : $csa->add_feature_table($chr_coord_system, 'gene');
  Description: This function tells the coordinate system adaptor that
               features from a specified table will be stored in a certain
               coordinate system.  If this information is not already stored
               in the database it will be added.
  Returntype : none
  Exceptions : none
  Caller     : BaseFeatureAdaptor
  Status     : At risk

=cut


#Can also be removed if inheritance/get_CoordSystemAdaptor issues resolved

sub add_feature_type {
  my $self = shift;
  my $cs   = shift;
  my $table = lc(shift);
  my $length = shift;
  if(!ref($cs) || !$cs->isa('Bio::EnsEMBL::Funcgen::CoordSystem')) {
    throw('CoordSystem argument is required.');
  }

  if(!$table) {
    throw('Table argument is required.');
  }

  my $cs_ids = $self->{'_feature_cache'}->{$table} || [];

  my ($exists) = grep {$cs->dbID() == $_} @$cs_ids;
  if( $exists ) {
    if( !$self->{'_max_len_cache'}->{$cs->dbID()}->{$table} ||
        $self->{'_max_len_cache'}->{$cs->dbID()}->{$table} < $length ) {
      my $sth = $self->prepare('UPDATE meta_coord ' .
                               "SET max_length = $length " .
                               'WHERE coord_system_id = ? ' .
                               'AND table_name = ? '.
                               "AND (max_length<$length ".
                               "OR max_length is null)");
      $sth->execute( $cs->dbID(), $table );
      $self->{'_max_len_cache'}->{$cs->dbID()}->{$table} = $length;
    }
    return;
  }

  #store the new tablename -> coord system relationship in the db
  #ignore failures b/c during the pipeline multiple processes may try
  #to update this table and only the first will be successful
  my $sth = $self->prepare('INSERT IGNORE INTO meta_coord ' .
                              'SET coord_system_id = ?, ' .
                                  'table_name = ?, ' .
			   'max_length = ? ' 
			  );

  $sth->execute($cs->dbID, $table, $length );

  #update the internal cache
  $self->{'_feature_cache'}->{$table} ||= [];
  push @{$self->{'_feature_cache'}->{$table}}, $cs->dbID();
  $self->{'_max_len_cache'}->{$cs->dbID()}->{$table} = $length;

  return;
}


1;
