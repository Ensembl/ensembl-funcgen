#
# Ensembl module for Bio::EnsEMBL::DBSQL::Funcgen::AnnotatedFeatureAdaptor
#

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

=head1 SYNOPSIS

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DBSQL::PeakCallingAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericAdaptor';

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::PeakCalling';
}

sub _tables {
  return ['peak_calling', 'pc']
}

sub fetch_all_by_Epigenome {
    my $self      = shift;
    my $epigenome = shift;
    
    if (! defined $epigenome) {
      throw("Epigenome was undefined");
    }
    my $constraint = 'epigenome_id  = ' . $epigenome->dbID;
    return $self->fetch_all($constraint);
}

sub _fetch_all_by_constraints {
    my $self            = shift;
    my $constraint_hash = shift;
    
    my @constraints;
    
    if (exists $constraint_hash->{epigenomes}) {
      for my $epigenome (@{$constraint_hash->{epigenomes}}) {
        push @constraints, 'epigenome_id  = ' . $epigenome->dbID;
      }
    }
    if (exists $constraint_hash->{feature_types}) {
      for my $feature_type (@{$constraint_hash->{feature_types}}) {
        push @constraints, 'feature_type_id  = ' . $feature_type->dbID;
      }
    }
    if (exists $constraint_hash->{projects}) {
      for my $project (@{$constraint_hash->{projects}}) {
      
        my $experiment_group_name = $project->name;
      
        push @constraints, "experiment_id in (select experiment_id from experiment join experimental_group using (experimental_group_id) where experimental_group.name = '".$experiment_group_name."')";
      }
    }
    if (exists $constraint_hash->{evidence_types}) {
      for my $evidence_type (@{$constraint_hash->{evidence_types}}) {
        if ($evidence_type eq 'DNase1 & TFBS') {
          push @constraints, "feature_type_id in (select feature_type_id from feature_type where class in ('Open Chromatin', 'Transcription Factor'))";
        }
        if ($evidence_type eq 'Hists & Pols') {
          push @constraints, "feature_type_id in (select feature_type_id from feature_type where class in ('Histone', 'Polymerase'))";
        }
      }
    }
    my $constraint = join ' and ', @constraints;
    return $self->fetch_all($constraint);
}

# Legacy method that had to be carried over, because web relies on it to 
# display sources.
#
sub _fetch_feature_set_filter_counts {
  my $self = shift;

   my $sql =<<SQL
    SELECT 
      count(*), 
      eg.name, 
      eg.description, 
      eg.is_project, 
      ft.class, 
      epi.name, 
      epi.description 
    FROM 
      experimental_group eg, 
      experiment e, 
      peak_calling fs, 
      feature_type ft, 
      epigenome epi 
    WHERE 
      fs.experiment_id = e.experiment_id 
      AND e.experimental_group_id = eg.experimental_group_id 
      AND fs.feature_type_id = ft.feature_type_id 
      AND fs.epigenome_id = epi.epigenome_id 
    GROUP BY 
      eg.name, 
      eg.is_project, 
      ft.class, 
      epi.name
SQL
;

  #warn $sql;
  #Need to write HC around this as we sometimes get less than expect.


  my @rows = @{$self->db->dbc->db_handle->selectall_arrayref($sql)};
  my $ft_a = $self->db->get_FeatureTypeAdaptor;
  my $ftype_info = $ft_a->get_regulatory_evidence_info;

  my %filter_info =
    (
     #Project=> {},
     #'Cell/Tissue' => {},
     All =>
     {
      All =>{ count       => 0,
              description => 'All experiments',
            }
     }

    );

  foreach my $row(@rows){

    my ($count, $project, $proj_desc, $is_proj, $ft_class, $epigenome_name, $epigenome_desc) = @$row;

    #All counts
    $filter_info{All}{All}{count} += $count;

    #Project counts
    if($is_proj){

      if(! exists $filter_info{Project}{$project}){
        $filter_info{Project}{$project} =
          { count       => 0,
            description => $proj_desc,
          };
      }

      $filter_info{Project}{$project}{count} += $count;
    }

    #Cell/Tissue counts
    if(! exists $filter_info{'Cell/Tissue'}{$epigenome_name}){
      $filter_info{'Cell/Tissue'}{$epigenome_name} =
        { count       => 0,
          description => $epigenome_desc,
        };
    }
    $filter_info{'Cell/Tissue'}{$epigenome_name}{count} += $count;

    #Evidence class counts
    #Do we want to split this into ft.class
    #i.e. split 'DNase1 & TFBS'
    my $evidence_type  = $ft_a->get_regulatory_evidence_type($ft_class);
    my $ft_class_label =  $ftype_info->{$evidence_type}{label};

    if(! exists $filter_info{'Evidence type'}{$ft_class_label}){
      $filter_info{'Evidence type'}{$ft_class_label} =
        { count       => 0,
          description => $ftype_info->{$evidence_type}{long_name},
        };
    }
    $filter_info{'Evidence type'}{$ft_class_label}{count} += $count;
  }

  return \%filter_info;

  #Do we need to add an 'in_build' filter /data field?

}

1;
