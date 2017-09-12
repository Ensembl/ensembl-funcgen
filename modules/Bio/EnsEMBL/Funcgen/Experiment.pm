=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::Experiment

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Experiment;

my $exp = Bio::EnsEMBL::Funcgen::Experiment->new
 (-ADAPTOR             => $self,
  -NAME                => $name,
  -EXPERIMENTAL_GROUP  => $experimental_group,
  -PRIMARY_DESIGN_TYPE => 'binding_site_indentification',
  -DESCRIPTION         => $description,
  -ARCHIVE_ID          => $archive_id,
 );

my $db_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(...);
my $exp_adaptor = $db_adaptor->get_ExperimentAdaptor();
my $exp = $exp_adaptor->fetch_by_name($exp_name)

=head1 DESCRIPTION

The Experiment class represents an instance of an experiment i.e. a discrete set of data

=cut


################################################################################

package Bio::EnsEMBL::Funcgen::Experiment;

use warnings;
use strict;

use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref check_ref );

use base qw( Bio::EnsEMBL::Funcgen::Storable );

=head2 new

  Arg [-NAME]                : String - experiment name
  Arg [-EXPERIMENTAL_GROUP]  : Bio::EnsEMBL::Funcgen ExperimentalGroup associated with this experiment
  Arg [-CONTROL]             : Bio::EnsEMBL::Funcgen::Experiment object which is used as control for this experiment,
  Arg [-IS_CONTROL]          : Boolean - defines whether this experiment is control or not
  Arg [-DATE]                : String - Date of the experiment (YYYY-MM-DD)
  Arg [-PRIMARY_DESIGN_TYPE] : String - MGED term for the primary design of teh experiment e.g. binding_site_identification
  Arg [-DESCRIPTION]         : String

  Example    : my $array = Bio::EnsEMBL::Funcgen::Experiment->new
                (-NAME                => $name,
                 -EXPERIMENTAL_GROUP  => $group,
                 -CONTROL             => $control,
                 -IS_CONTROL          => 1,
                 -PRIMARY_DESIGN_TYPE => $p_design_type,
                 -DESCRIPTION         => $description,
                 -ARCHIVE_ID          => 'SRX000000',
                 -DISPLAY_URL         => $non_archive_url);

  Description: Creates a new Bio::EnsEMBL::Funcgen::Experiment object.
  Returntype : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions : Throws if name not defined
               Throws if ExperimentalGroup not valid
  Caller     : General
  Status     : Medium Risk

=cut

sub new {
  my $caller = shift;
  my $class  = ref($caller) || $caller;
  my $self   = $class->SUPER::new(@_);

  my ( $name, $group, $control, $is_control, $desc, $ftype, $epigenome,
    $archive_id, )
    = rearrange(
    [   'NAME',        'EXPERIMENTAL_GROUP',
        'CONTROL',     'IS_CONTROL',
        'DESCRIPTION', 'FEATURE_TYPE',
        'EPIGENOME',   'ARCHIVE_ID',
    ],
    @_
    );

  # Mandatory attr checks
  throw('You must provide a name parameter') if ! defined $name;
  throw('You must provide a is_control parameter') if ! defined $is_control;

  assert_ref( $group,     'Bio::EnsEMBL::Funcgen::ExperimentalGroup' );
  if (defined $epigenome) {
    assert_ref( $epigenome, 'Bio::EnsEMBL::Funcgen::Epigenome' );
  }
  assert_ref( $ftype,     'Bio::EnsEMBL::Funcgen::FeatureType' );
  assert_ref( $control,   'Bio::EnsEMBL::Funcgen::Experiment') if defined $control;


  #Direct assignment here so we avoid setter test in methods
  $self->{name}                = $name;
  $self->{group}               = $group;
  $self->{control}             = $control if defined $control;
  $self->{is_control}          = $is_control;
  $self->{description}         = $desc       if defined $desc;
  $self->{epigenome}           = $epigenome;
  $self->{feature_type}        = $ftype;
  $self->{archive_id}          = $archive_id;

  return $self;
}

sub count_biological_replicates {

    my $self = shift;

    my $read_file_experimental_configuration_adaptor = $self
        ->adaptor
        ->db
        ->get_ReadFileExperimentalConfigurationAdaptor;
    
    my $count_biological_replicates 
        = $read_file_experimental_configuration_adaptor
            ->count_biological_replicates_from_Experiment($self);
    
    return $count_biological_replicates;
}

sub count_technical_replicates {

    my $self = shift;

    my $read_file_experimental_configuration_adaptor = $self
        ->adaptor
        ->db
        ->get_ReadFileExperimentalConfigurationAdaptor;
    
    my $count_technical_replicates 
        = $read_file_experimental_configuration_adaptor
            ->count_technical_replicates_from_Experiment($self);
    
    return $count_technical_replicates;
}

=head2 epigenome

  Example    : my $epigenome_name = $exp->epigenome->name;
  Description: Getter for the Epigenome
  Returntype : Bio::EnsEMBL::Funcgen::Epigenome
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub epigenome { return shift->{epigenome}; }


=head2 feature_type

  Example    : my $ftype_name = $exp->feature_type->name;
  Description: Getter for the FeatureType.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : None
  Caller     : General
  Status     : Stable

=cut

sub feature_type { return shift->{feature_type}; }




=head2 name

  Example     : my $exp_name = $exp->name;
  Description : Getter for the experiment name
  Returntype  : String
  Exceptions  : None
  Caller      : General
  Status      : Stable

=cut

sub name{
  return shift->{name};
}


=head2 experimental_group

  Example     : my $exp_group_name = $exp->experimental_group()->name();
  Description : Getter for the experimental group
  Returntype  : Bio::EnsEMBL::Funcgen::ExperimentalGroup
  Exceptions  : None
  Caller      : General
  Status      : At risk - to be deprecated

=cut

#change in webcode before deprecating

sub experimental_group{
  return shift->{group};
}


=head2 get_ExperimentalGroup

  Example     : my $exp_group_name = $exp->experimental_group()->name();
  Description : Getter for the experimental group
  Returntype  : Bio::EnsEMBL::Funcgen::ExperimentalGroup
  Exceptions  : None
  Caller      : General
  Status      : Stable

=cut

sub get_ExperimentalGroup{ return shift->{group}; }


=head2 is_control

  Example     : my $is_control = $exp->is_control();
  Description : Getter for the is_control attribute
  Returntype  : Boolean
  Exceptions  : None
  Caller      : General
  Status      : Stable

=cut

sub is_control{
  return shift->{is_control};
}


=head2 get_control

  Example     : my $control_exp = $exp->get_control();
  Description : Getter for the experiment which is used as control
  Returntype  : Bio::EnsEMBL::Funcgen::Experiment
  Exceptions  : None
  Caller      : General
  Status      : Stable

=cut

sub get_control{
  return shift->{control};
}


=head2 get_InputSubsets

  Example     : my @issets = @{$exp->getInputSubsets()}
  Description : Retrieves all InputSubsets associated with this experiment.
  Returntype  : Arrayref of Bio::EnsEMBL::Funcgen::InputSubsets
  Exceptions  : None
  Caller      : General
  Status      : At risk

=cut

sub get_InputSubsets{
  my $self = shift;

  if(! exists $self->{'input_subsets'}){
    $self->{'input_subsets'} = {};

    foreach my $isset(@{$self->adaptor->db->get_InputSubsetAdaptor->
                          fetch_all_by_Experiments([$self])}){
      $self->{'input_subsets'}->{$isset->dbID} = $isset;
    }
  }

  return [values %{$self->{'input_subsets'}}];
}

=head2 archive_id

  Example     : my $archive_id = $exp->archive;
  Description : Getter for the archive_id
  Returntype  : String
  Exceptions  : None
  Caller      : General
  Status      : At risk

=cut

sub archive_id { return shift->{archive_id};  }

=head2 source_info

  Example    : my $source_info = $exp->source_info;
  Description: Getter for the experiment source info i.e. [[ $label, $url ]]
  Returntype : Arrayref
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

sub source_info{
  my $self = shift;

  if(! defined $self->{source_info}){
    #could have data_url as highest priority here
    #but we need to ensure removal when adding archive ids
    #so we link to the archive and not the old data url

    my $exp_group = $self->experimental_group;
    my @source_info;
    my ($proj_name, $proj_link);

    if($exp_group->is_project){
      $proj_name = $exp_group->name;
      $proj_link = $exp_group->url;
    }


    if(defined $self->archive_id ){
      #Need to handled comma separated values in here

      foreach my $archive_id(split/,/, $self->archive_id){

        push @source_info, [$archive_id, $self->display_url];
        #source_link can be undef here as archive_id overrides display url
        #undef links will automatically go to the SRA
        #If multiple IDs exists, they will all use the same display_url
      }
    }
    elsif(defined $proj_name){
      push @source_info, [$proj_name, ($self->display_url || $proj_link)];
    }

    $self->{source_info} = \@source_info;
  }

  return $self->{source_info};
}

1;

