=head1 LICENSE

Copyright [1999-2016] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

Bio::EnsEMBL::Funcgen::RegulatoryActivity

=head1 SYNOPSIS
=head1 DESCRIPTION
=head1 SEE ALSO
=cut

package Bio::EnsEMBL::Funcgen::RegulatoryActivity;

use strict;
use warnings;
use Bio::EnsEMBL::Utils::Argument  qw( rearrange );
use Bio::EnsEMBL::Utils::Exception qw( throw deprecate );

sub new {
  my $caller = shift;
  my $class = ref($caller) || $caller;
  my $self =  bless({}, $class);
  return $self;
}

sub db {
  my $self = shift;
  my $db   = shift;

  if ($db) {
    $self->{'_db'}  = $db;
  }
  return $self->{'_db'};
}

=head2 epigenome_id
=cut
sub epigenome_id {
  my $self = shift;
  my $epigenome_id = shift;

  if(defined $epigenome_id) {
    $self->{'_epigenome_id'}  = $epigenome_id;
  }
  return $self->{'_epigenome_id'};
}

=head2 epigenome
=cut
sub epigenome {
  my $self = shift;

  if(! defined $self->{'_epigenome'}) {
    $self->{'_epigenome'} = $self
      ->db
      ->get_EpigenomeAdaptor()
      ->fetch_by_dbID(
	$self->epigenome_id
      );
  }
  return $self->{'_epigenome'};
}

=head2 _is_multicell
=cut
sub _is_multicell {
  my $self = shift;
  my $is_multicell = shift;

  if(defined $is_multicell) {
    $self->{'_is_multicell'}  = $is_multicell;
  }
  return $self->{'_is_multicell'};
}

=head2 activity
=cut
sub activity {
  my $self = shift;
  my $activity = shift;

  if(defined $activity) {
    $self->{'_activity'}  = $activity;
  }
  return $self->{'_activity'};
}

=head2 regulatory_feature
=cut
sub regulatory_feature {
  my $self = shift;
  my $regulatory_feature = shift;

  if(defined $regulatory_feature) {
    $self->{'_regulatory_feature'}  = $regulatory_feature;
  }
  return $self->{'_regulatory_feature'};
}

=head2 regulatory_evidence
=cut
sub regulatory_evidence {
  my $self = shift;
  my $regulatory_evidence = shift;

  if ($regulatory_evidence) {
    $self->{'_regulatory_evidence'} = $regulatory_evidence;
  }
  if (! $self->{'_regulatory_evidence'}) {
    use Bio::EnsEMBL::Funcgen::RegulatoryEvidence;
    my $regulatory_evidence = Bio::EnsEMBL::Funcgen::RegulatoryEvidence->new;
    $regulatory_evidence->db($self->db);
    $self->{'_regulatory_evidence'} = $regulatory_evidence;
  }
  return $self->{'_regulatory_evidence'};
}

# sub SO_term {
#   # http://www.sequenceontology.org/browser/current_svn/term/SO:0005836
#   return 'SO:0005836';
# }
sub SO_term {
  my $self = shift;
  return $self->regulatory_feature->feature_type->so_accession;
}

sub summary_as_hash {
  my $self = shift;
  
  my $regulatory_feature = $self->regulatory_feature;
  my $epigenome          = $self->epigenome;

  return {
    regulatory_feature_stable_id => $regulatory_feature->stable_id,
    epigenome           => $epigenome->display_label,
#     source              => 'Regulatory_Build',
    source              => $regulatory_feature->analysis->logic_name,
    bound_start         => $regulatory_feature->bound_seq_region_start,
    bound_end           => $regulatory_feature->bound_seq_region_end,
    start               => $regulatory_feature->seq_region_start,
    end                 => $regulatory_feature->seq_region_end,
    strand              => $regulatory_feature->strand,
    seq_region_name     => $regulatory_feature->seq_region_name,
    activity            => $self->activity,
    description         => 'Activity of ' . $regulatory_feature->stable_id . ' in ' . $epigenome->display_label,
    feature_type        => "regulatory activity",
  };
}

1;
