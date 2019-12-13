=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::AlignmentAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

use Bio::EnsEMBL::Utils::Exception qw( throw );

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::Alignment';
}

sub _tables {
  return ['alignment', 'a']
}

sub _load_dependencies {
    my $self = shift;
    my $alignment = shift;
    
    my @linked_read_file_ids;
    
    $self->sql_helper->execute_no_return(
      -SQL          => 'select read_file_id from alignment_read_file where alignment_id = ?',
      -PARAMS       => [ $alignment->dbID ],
      -USE_HASHREFS => 1,
      -CALLBACK     => sub {
          my $row = shift;
          my $read_file_id = $row->{read_file_id};
          push @linked_read_file_ids, $read_file_id;
          return;
        },
    );
    $alignment->read_file_ids(\@linked_read_file_ids);
    return $alignment;
}

sub _store_dependencies {
    my $self = shift;
    my $alignment = shift;

    my $alignment_id  = $alignment->dbID;
    my $read_file_ids = $alignment->read_file_ids;
    
    $self->sql_helper->batch(
      -SQL      => 'INSERT INTO alignment_read_file (alignment_id, read_file_id) values (?, ?)',
      -CALLBACk => sub {
        my $sth = shift;
        foreach my $read_file_id (@$read_file_ids) {
            $sth->execute( $alignment_id, $read_file_id );
        }
      }
    );
    return;
}

sub fetch_all_deduplicated_replicates_by_Alignment {
    my $self      = shift;
    my $alignment = shift;
    
    if (! defined $alignment) {
      throw("Alignment was undefined");
    }
    
    my $constraint = join ' and ', (
        'experiment_id  = ' . $alignment->experiment_id,
        'is_complete    = False',
        'has_duplicates = False'
    );
    return $self->fetch_all($constraint);
}

sub fetch_all_by_ReadFileExperimentalConfiguration {
  my $self = shift;
  my $read_file_experimental_configuration = shift;
 
  my @all_alignments;
  $self->sql_helper->execute_no_return(
    -SQL          => '
      select 
        distinct alignment_id
      from 
        alignment_read_file 
        join read_file_experimental_configuration using (read_file_id)
      where 
        read_file_experimental_configuration_id = ?
    ',
    -PARAMS       => [ $read_file_experimental_configuration->dbID ],
    -USE_HASHREFS => 1,
    -CALLBACK     => sub {
        my $row = shift;
        my $alignment_id = $row->{alignment_id};
        my $alignment = $self->fetch_by_dbID($alignment_id);
        
        push @all_alignments, $alignment;
        return;
      },
  );
  return \@all_alignments;
}

sub fetch_complete_deduplicated_by_Experiment {
    my $self       = shift;
    my $experiment = shift;
    
    if (! defined $experiment) {
      throw("Experiment was undefined");
    }
    
    my $constraint = join ' and ', (
        'is_complete    = 1',
        'has_duplicates = 0',
        'experiment_id  = ' . $experiment->dbID
    );
    return $self->fetch_single_object($constraint);
}

sub fetch_all_with_duplicates_by_Experiment {
    my $self       = shift;
    my $experiment = shift;
    
    if (! defined $experiment) {
      throw("Experiment was undefined");
    }
    
    my $constraint = join ' and ', (
        'has_duplicates = 1',
        'experiment_id  = ' . $experiment->dbID
    );
    return $self->fetch_all($constraint);
}

1;
