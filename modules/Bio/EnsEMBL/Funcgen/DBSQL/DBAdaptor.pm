=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor

=head1 SYNOPSIS

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host => "ensembldb.ensembl.org",
   -dbname => "mus_musculus_funcgen_67_37",
   -species => "Mus_musculus",
   -user => "anonymous",
   -port => '3307',
  );

my $experiment_adaptor = $db->get_ExperimentAdaptor();

=back

=head1 DESCRIPTION

An adaptor to access the funcgen database and expose other available adaptors.

=cut

################################################################################


package Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

use strict;
use warnings;
use DBI; #for resolving core DB
use Bio::EnsEMBL::Utils::Exception         qw( throw deprecate warning );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( assert_ref_do );
use Bio::EnsEMBL::Registry;

use base qw(Bio::EnsEMBL::DBSQL::DBAdaptor);

=head2 new
  Description: Constructor for DBAdaptor.
  Returntype : Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor
  Caller     : general
  Status     : Stable
=cut
sub new {
  my ($class, @args) = @_;
  my $self = $class->SUPER::new(@args, '-group', 'funcgen');
  return $self;
}

=head2 dbfile_data_root

  Arg[1]     : Optional String: Root path of dbfile data directory
  Example    : $rset_adaptor->dbfile_data_root('/data/root/dir/);
  Description: This allows the root path to be defined. If an adaptor uses
               files, it will use this to find its data.
  Returntype : String
  Exceptions : None
  Caller     : Bio::EnsEMBL::Funcgen::DBAdaptor::ResultSet
  Status     : at risk - move this to SetAdaptor/FileAdaptor?

=cut

sub dbfile_data_root {
  my ($self, $root) = @_;

  if($root){
    $root =~ s/\/$//o;  # strip off trailing /, as this is present in dbfile_registry.path
    $self->{dbfile_data_root} = $root;
  }

  return $self->{dbfile_data_root} || '';  # Avoids concat warning
}

=head2 is_stored_and_valid

  Arg [1]    : String - class namespace
  Arg [2]    : Bio::EnsEMBL::Funcgen::Storable e.g. ResultSet etc.
  Arg [3]    : String (optional) - Name of variable to use in error output (for use with assert_ref)
  Example    : $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);
  DESCRIPTION: Validates object class and stored status
  Returntype : None
  Exceptions : Throws if Storable is not stored
  Caller     : General
  Status     : At risk

=cut

sub is_stored_and_valid {
  my ($self, $class, $obj, $name) = @_;
  
  use Carp;
  
  if (! defined $class) {
    confess("class is not defined!");
  }
  if (! defined $obj) {
    confess("obj is not defined!");
  }
  
  assert_ref($obj, $class, $name);

  if (! $obj->can('is_stored')) {
    return defined $obj->dbID
  }

  throw("$obj is not stored") if ! $obj->is_stored($self);
  return;
}


=head2 are_stored_and_valid

  Arg [1]    : String - Namespace of class
  Arg [2]    : ARRAYREF os Bio::EnsEMBL::Funcgen::Storable objects e.g. ResultSet
  Arg [3]    : String : return value method name
  Example    : $db->are_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', \@rsets);
  DESCRIPTION: Wrapper for is_stored_and_valid. Will optionally return array of values
               defined by calling method name arg on each object passed
  Returntype : ARRAYREF - contents defined by optional method name arg
  Exceptions : Throws if object list is not an ARRAY with at least one element
  Caller     : general
  Status     : At risk

=cut

sub are_stored_and_valid {
  my ($self, $class, $obj_list, $method_name) = @_;
  assert_ref($obj_list, 'ARRAY', 'object list');

  if(scalar(@$obj_list) == 0){
   throw('Objects Arrayref is empty');
  }

  my @return_vals;

  foreach my $obj (@$obj_list) {
    #$self->is_stored_and_valid($class, $obj);

    if(! $method_name){
      assert_ref($obj, $class, 'object');
    }
    else{
      push @return_vals, assert_ref_do($obj, $class, $method_name, 'object');
    }
  }
  return \@return_vals;
}

=head2 get_available_adaptors

  Example    : my %pairs = %{$dba->get_available_adaptors()};
  Description: gets a hash of the available adaptors
  ReturnType : reference to a hash
  Exceptions : none
  Caller     : Bio::EnsEMBL::Utils::ConfigRegistry
  Status     : Stable

=cut
sub get_available_adaptors{
    my $self = shift;

    my %pairs = (
        'Alignment'                         => 'Bio::EnsEMBL::Funcgen::DBSQL::AlignmentAdaptor',
        'Analysis'                          => 'Bio::EnsEMBL::DBSQL::AnalysisAdaptor',
        'Array'                             => 'Bio::EnsEMBL::Funcgen::DBSQL::ArrayAdaptor',
        'ArrayChip'                         => 'Bio::EnsEMBL::Funcgen::DBSQL::ArrayChipAdaptor',
        'BindingMatrix'                     => 'Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixAdaptor',
        'BindingMatrixFrequencies'          => 'Bio::EnsEMBL::Funcgen::DBSQL::BindingMatrixFrequenciesAdaptor',
        'Chance'                            => 'Bio::EnsEMBL::Funcgen::DBSQL::ChanceAdaptor',
        'CrisprSitesFile'                   => 'Bio::EnsEMBL::Funcgen::DBSQL::CrisprSitesFileAdaptor',
        'DataFile'                          => 'Bio::EnsEMBL::Funcgen::DBSQL::DataFileAdaptor',
        'DBEntry'                           => 'Bio::EnsEMBL::Funcgen::DBSQL::DBEntryAdaptor',
        'DNAMethylationFile'                => 'Bio::EnsEMBL::Funcgen::DBSQL::DNAMethylationFileAdaptor',
        'Epigenome'                         => 'Bio::EnsEMBL::Funcgen::DBSQL::EpigenomeAdaptor',
        'ExampleFeature'                    => 'Bio::EnsEMBL::Funcgen::DBSQL::ExampleFeatureAdaptor',
        'ExecutionPlan'                     => 'Bio::EnsEMBL::Funcgen::DBSQL::ExecutionPlanAdaptor',
        'Experiment'                        => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentAdaptor',
        'ExperimentalGroup'                 => 'Bio::EnsEMBL::Funcgen::DBSQL::ExperimentalGroupAdaptor',
        'ExternalFeature'                   => 'Bio::EnsEMBL::Funcgen::DBSQL::ExternalFeatureAdaptor',
        'FastQC'                            => 'Bio::EnsEMBL::Funcgen::DBSQL::FastQCAdaptor',
        'FeatureSet'                        => 'Bio::EnsEMBL::Funcgen::DBSQL::FeatureSetAdaptor',
        'FeatureType'                       => 'Bio::EnsEMBL::Funcgen::DBSQL::FeatureTypeAdaptor',
        'Frip'                              => 'Bio::EnsEMBL::Funcgen::DBSQL::FripAdaptor',
        'Idr'                               => 'Bio::EnsEMBL::Funcgen::DBSQL::IdrAdaptor',
        'MetaContainer'                     => 'Bio::EnsEMBL::DBSQL::MetaContainer',
        'MetaCoordContainer'                => 'Bio::EnsEMBL::DBSQL::MetaCoordContainer',
        'MirnaTargetFeature'                => 'Bio::EnsEMBL::Funcgen::DBSQL::MirnaTargetFeatureAdaptor',
        'MotifFeature'                      => 'Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureAdaptor',
        'MotifFeatureFile'                  => 'Bio::EnsEMBL::Funcgen::DBSQL::MotifFeatureFileAdaptor',
        'Peak'                              => 'Bio::EnsEMBL::Funcgen::DBSQL::PeakAdaptor',
        'PeakCalling'                       => 'Bio::EnsEMBL::Funcgen::DBSQL::PeakCallingAdaptor',
        'PeakCallingStatistic'              => 'Bio::EnsEMBL::Funcgen::DBSQL::PeakCallingStatisticAdaptor',
        'PhantomPeak'                       => 'Bio::EnsEMBL::Funcgen::DBSQL::PhantomPeakAdaptor',
        'Probe'                             => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeAdaptor',
        'ProbeFeature'                      => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureAdaptor',
        'ProbeFeatureTranscriptMapping'     => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeFeatureTranscriptMappingAdaptor',
        'ProbeMapping'                      => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeMappingAdaptor',
        'ProbeMappingStatistic'             => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeMappingStatisticAdaptor',
        'ProbeSequence'                     => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeSequenceAdaptor',
        'ProbeSet'                          => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeSetAdaptor',
        'ProbeSetTranscriptMapping'         => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeSetTranscriptMappingAdaptor',
        'ProbeTranscriptMapping'            => 'Bio::EnsEMBL::Funcgen::DBSQL::ProbeTranscriptMappingAdaptor',
        'ReadFile'                          => 'Bio::EnsEMBL::Funcgen::DBSQL::ReadFileAdaptor',
        'ReadFileExperimentalConfiguration' => 'Bio::EnsEMBL::Funcgen::DBSQL::ReadFileExperimentalConfigurationAdaptor',
        'RegulatoryActivity'                => 'Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryActivityAdaptor',
        'RegulatoryBuild'                   => 'Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildAdaptor',
        'RegulatoryBuildStatistic'          => 'Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryBuildStatisticAdaptor',
        'RegulatoryEvidenceLink'            => 'Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryEvidenceLinkAdaptor',
        'RegulatoryFeature'                 => 'Bio::EnsEMBL::Funcgen::DBSQL::RegulatoryFeatureAdaptor',
        'Segmentation'                      => 'Bio::EnsEMBL::Funcgen::DBSQL::SegmentationAdaptor',
        'SegmentationFile'                  => 'Bio::EnsEMBL::Funcgen::DBSQL::SegmentationFileAdaptor',
        'SegmentationStateAssignment'       => 'Bio::EnsEMBL::Funcgen::DBSQL::SegmentationStateAssignmentAdaptor',
        'SegmentationStateEmission'         => 'Bio::EnsEMBL::Funcgen::DBSQL::SegmentationStateEmissionAdaptor',
        'SegmentationStatistic'             => 'Bio::EnsEMBL::Funcgen::DBSQL::SegmentationStatisticAdaptor',
        'TranscriptionFactor'               => 'Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorAdaptor',
        'TranscriptionFactorComplex'        => 'Bio::EnsEMBL::Funcgen::DBSQL::TranscriptionFactorComplexAdaptor',
        'UnmappedObject'                    => 'Bio::EnsEMBL::DBSQL::UnmappedObjectAdaptor',
    );
    
    return (\%pairs);
}

1;
