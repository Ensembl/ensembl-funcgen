=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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
=head1 DESCRIPTION
=cut

package Bio::EnsEMBL::Funcgen::Ftp::FetchQCRelatedData;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = bless {}, $class;

  my (
    $db_connection,
  )
    = rearrange([ 'db_connection' ], @_);

  $self->db_connection($db_connection);
  
  my $sth_fetch_input_subsets = $db_connection->prepare(
    qq(
      select 
	input_subset.name,
	input_subset.biological_replicate,
	input_subset.technical_replicate,
	input_subset.is_control
      from 
	result_set 
	join result_set_input using (result_set_id) 
	join input_subset on (table_id=input_subset_id and table_name="input_subset") 
	join experiment on (result_set.experiment_id=experiment.experiment_id)
      where 
	result_set.name=?
	and experiment.is_control=input_subset.is_control
    )
  );
  $self->_sth_fetch_input_subsets($sth_fetch_input_subsets);
  
  my $sth_fetch_xrefs = $db_connection->prepare(
    qq(
      select 
	external_db.db_display_name external_database,
	xref.dbprimary_acc accession
      from epigenome 
	join object_xref on (ensembl_id=epigenome_id and ensembl_object_type="epigenome")
	join xref using (xref_id)
	join external_db using (external_db_id)
      where 
	epigenome.display_label = ?
    )
  );
  $self->_sth_fetch_xrefs($sth_fetch_xrefs);

  return $self;
}

sub _sth_fetch_input_subsets { return shift->_generic_get_or_set('_sth_fetch_input_subsets', @_) }
sub _sth_fetch_xrefs         { return shift->_generic_get_or_set('_sth_fetch_xrefs',         @_) }
sub db_connection            { return shift->_generic_get_or_set('db_connection',            @_) }

sub fetch_input_subset_data_for_result_set {

  my $self            = shift;
  my $result_set_name = shift;
  
  my $sth = $self->_sth_fetch_input_subsets;
  
  $sth->bind_param(1, $result_set_name);
  $sth->execute;
  my $current_input_subset = $sth->fetchall_hashref("name");
  my @input_subsets = values %$current_input_subset;
  return \@input_subsets;
}

sub fetch_xrefs_for_epigenome {

  my $self                    = shift;
  my $epigenome_display_label = shift;
  
  my $sth = $self->_sth_fetch_xrefs;
  
  $sth->bind_param(1, $epigenome_display_label);
  $sth->execute;
  my $current_xrefs = $sth->fetchall_hashref("accession");
  my @xrefs = values %$current_xrefs;
  
  use Hash::Util qw( lock_hash );
  my $xref_hash = {};
  foreach my $current_xref (@xrefs) {  
    lock_hash(%$current_xref);
    $xref_hash->{$current_xref->{external_database}} = $current_xref->{accession};
  }
  
  my $epigenome_hash = {
    name                => $epigenome_display_label,
    external_references => $xref_hash
  };
  return $epigenome_hash;
}

sub _generic_get_or_set {
  my $self  = shift;
  my $name  = shift;
  my $value = shift;

  if(defined $value) {
    $self->{$name}  = $value;
  }
  return $self->{$name};
}

1;


