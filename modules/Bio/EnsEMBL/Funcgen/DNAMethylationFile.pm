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

  Bio::EnsEMBL::Funcgen::DNAMethylationFile

=head1 SYNOPSIS

  use Bio::EnsEMBL::Registry;
  use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

  Bio::EnsEMBL::Registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org', # alternatively 'useastdb.ensembl.org'
    -user => 'anonymous'
  );

  my $species = 'homo_sapiens';

  my $dna_methylation_file_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'funcgen', 'DNAMethylationFile');
  my $dna_methylation_file = $dna_methylation_file_adaptor->fetch_all;

  print "\nThere are " . @$dna_methylation_file . " dna methylation files in the database for $species.\n";

  foreach my $current_dna_methylation_file (@$dna_methylation_file) {
    print to_string($current_dna_methylation_file);
    print "\n\n";
  }

  sub to_string {
    my $dna_methylation_file = shift;

    return join "\n", 
      "The name of the dna methylation track is: " . $dna_methylation_file->name,
      "The description is:                       " . $dna_methylation_file->get_Analysis->description,
      "The name of the dna methylation file is:  " . $dna_methylation_file->file;
  }

=head1 DESCRIPTION

  DNA methylation sites are not stored in the database, but in an external file. The
  DNAMethylationFile object helps users locate the file.

=cut

package Bio::EnsEMBL::Funcgen::DNAMethylationFile;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

=head2 new
=cut
sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self = bless {}, $class;

  my @field = qw(
    dbID
    name
    file
    file_type
  );
  
  my (
    $dbID,
    $name,
    $file,
    $file_type,
  )
    = rearrange([ @field ], @_);

  $self->dbID            ($dbID);
  $self->name            ($name);
  $self->file            ($file);
  $self->file_type       ($file_type);

  return $self;
}

sub dbID           { return shift->_generic_get_or_set('dbID',            @_) }
sub name           { return shift->_generic_get_or_set('name',            @_) }
sub file           { return shift->_generic_get_or_set('file',            @_) }
sub file_type      { return shift->_generic_get_or_set('file_type',       @_) }
sub _analysis      { return shift->_generic_get_or_set('_analysis',       @_) }
sub _feature_type  { return shift->_generic_get_or_set('_feature_type',   @_) }
sub adaptor        { return shift->_generic_get_or_set('adaptor',         @_) }

=head2 get_Analysis

  Description: Fetches the analysis used for DNAMethylation files.
  Returntype : Bio::EnsEMBL::Analysis
  Exceptions : none
  Status     : At Risk

=cut

sub get_Analysis {
  my $self = shift;
  return $self->_analysis
}

=head2 get_FeatureType

  Description: Fetches the feature type used for DNAMethylation files.
  Returntype : Bio::EnsEMBL::Funcgen::FeatureType
  Exceptions : none
  Status     : At Risk

=cut

sub get_FeatureType {
  my $self = shift;
  return $self->_feature_type
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


