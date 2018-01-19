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
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::DefineDataSet

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::Hive::DefineDataSet;

use warnings;
use strict;
 
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( scalars_to_objects 
                                               get_set_prefix_from_Set );

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub fetch_input {
  my $self = shift @_;
  $self->SUPER::fetch_input;  
  my $set = $self->fetch_Set_input('ResultSet');#defines result_set method

  my $anal_type = $self->param_required('feature_set_analysis_type'); 
  my $set_lname = $self->get_param_method($anal_type.'_analysis', 'silent');
  
  
  ### Set post IDR FeatureSet analysis 
  if($self->is_idr_ResultSet($set)) {
    if(my $ppeak_anal = $self->get_param_method('permissive_peaks') ) {
      $set_lname = $ppeak_anal.'_IDR';
    }
  }

  if(! defined $set_lname) {
  
    my %default_analyses;
    
    $default_analyses{feature_set} = $self->param_silent('default_feature_set_analyses');     
    
    if (! defined $default_analyses{feature_set}) {
      throw("Please define -feature_set_analysis or add to default_feature_set_analyses".
        " in the analysis config");
    }
    
    if (exists $default_analyses{feature_set}{$set->feature_type->name}) {
    
      $set_lname = $default_analyses{feature_set}{$set->feature_type->name}; 
      
    } elsif (exists $default_analyses{feature_set}{$set->feature_type->class}) {
    
      $set_lname = $default_analyses{feature_set}{$set->feature_type->class};
      
    } else {
    
      throw("No default feature_set analysis available for ".$set->feature_type->name.
        '('.$set->feature_type->class.
        ").\n Please add FeatureType name or class to  default_feature_set_analyses ".
        "in analysis config or specify -feature_set_analysis"); 
    }
  }

  # Catch undefs in config
  if(! defined $set_lname) {
    throw('Unable to identify defined feature_set analysis in config for '.$set->feature_type->name.
    ".\nPlease define in -default_feature_set_analyses config or specify -feature_set_analysis");     
  }
  
  $self->param('feature_set_analysis', 
               &scalars_to_objects($self->out_db,
                                   'Analysis', 
                                   'fetch_by_logic_name',
                                   [$set_lname])->[0]);
  return;
}

sub run {   # Check parameters and do appropriate database/file operations... 
  my $self   = shift;
  my $helper = $self->helper;  
  my $rset   = $self->ResultSet;
  my $fset_anal = $self->param('feature_set_analysis');
  my $set_prefix = get_set_prefix_from_Set($rset);
  my $set;    
  
    warn(
      "result set name:" . $rset->name . "\n"
      . "set_prefix:" . $set_prefix
    );
#    die($set_prefix);
  
  eval {
    # Never set -FULL_DELETE here!
    # It is unwise to do this in a pipeline and should be handled
    # on a case by case basis using a separate rollback script    
    # Should also never really specify recover here either?
    # This bascailly ignores the fact that a ResultSet may be linked to other DataSets
    $set = $helper->define_DataSet
    (-NAME                 => $set_prefix.'_'.$fset_anal->logic_name,
      -FEATURE_CLASS        => 'annotated', #Is there overlap with rset feature_class here?
      -SUPPORTING_SETS      => [$rset],
      -DBADAPTOR            => $self->out_db,
      -FEATURE_SET_ANALYSIS => $fset_anal,
      -RESULT_SET_ANALYSIS  => $self->param('result_set_analysis'),
      -RESULT_SET_MODE      => $self->param('result_set_mode'),
      -ROLLBACK             => $self->param('rollback'),
      -RECOVER              => $self->param('recover'),
      -SLICES               => $self->slices,    
      -EPIGENOME            => $rset->epigenome,
      -FEATURE_TYPE         => $rset->feature_type                    );  
  };
  if ($@) {
    $self->throw($@);
  }
  $self->param('output_id', 
               {%{$self->batch_params},
                dbID       => $set->dbID, 
                set_name   => $set->name,
                set_type   => 'DataSet'  });
  return ;
}

sub write_output {  # Create the relevant jobs
  my $self = shift;
  $self->helper->debug(1, 'DefineDataSet data flowing:', $self->param('output_id'));
  $self->dataflow_output_id($self->param('output_id'), 2);
  return;
}

1;
