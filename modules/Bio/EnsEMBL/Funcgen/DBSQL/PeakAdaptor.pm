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

package Bio::EnsEMBL::Funcgen::DBSQL::PeakAdaptor;

use strict;
use base 'Bio::EnsEMBL::Funcgen::DBSQL::GenericFeatureAdaptor';
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Utils::Scalar qw( assert_ref );

sub object_class {
    return 'Bio::EnsEMBL::Funcgen::Peak';
}

sub _tables {
  return ['peak', 'p']
}

sub insertion_method {
    return 'insert ignore'
}

sub _columns {
  return qw(
    p.peak_id               p.seq_region_id
    p.seq_region_start      p.seq_region_end
    p.seq_region_strand     p.peak_calling_id
    p.score                 p.summit
  );
}

=head2 fetch_all_by_PeakCalling

  Arg [1]    : Bio::EnsEMBL::Funcgen::PeakCalling
  Example    : None
  Description: Fetches a list of Peak objects by PeakCalling
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::Peak objects
  Exceptions : Throws if PeakCalling parameter is not specified
  Caller     : Internal
  Status     : At Risk

=cut

sub fetch_all_by_PeakCalling {

    my $self         = shift;
    my $peak_calling = shift;

    assert_ref($peak_calling, 'Bio::EnsEMBL::Funcgen::PeakCalling');

    my $constraint = 'peak_calling_id = ' . $peak_calling->dbID;

    return $self->fetch_all($constraint);
}

=head2 fetch_all_by_Slice_PeakCalling

  Arg [1]    : Bio::EnsEMBL::Slice
  Arg [2]    : Bio::EnsEMBL::Funcgen::PeakCalling
  Example    : None
  Description: Fetches a list of Peak objects by PeakCalling on a given Slice
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::Peak objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub fetch_all_by_Slice_PeakCalling {

  my $self         = shift;
  my $slice        = shift;
  my $peak_calling = shift;
  
  my $peak_calling_id = $peak_calling->dbID;
  my $features = $self->fetch_all_by_Slice_constraint(
    $slice, 
    "peak_calling_id = " . $peak_calling_id
  );
  return $features;
}

=head2 _fetch_overlapping_MotifFeatures

  Arg [1]    : Bio::EnsEMBL::Funcgen::Peak
  Example    : None
  Description: Fetches a list of MotifFeature objects that overlap
               with the Peak object
  Returntype : Arrayref of Bio::EnsEMBL::Funcgen::MotifFeature objects
  Exceptions : None
  Caller     : Internal
  Status     : At Risk

=cut

sub _fetch_overlapping_MotifFeatures {
    my ( $self, $peak ) = @_;

    if (! defined $peak){
      throw('Must provide a Peak parameter');
    }

    my $sth = $self->prepare( "
      SELECT DISTINCT motif_feature_id FROM 
      motif_feature_peak
      WHERE peak_id=?
      " );

    $sth->execute( $peak->dbID() );

    my @motif_feature_ids;

    while ( my @row = $sth->fetchrow_array ) {
        push @motif_feature_ids, $row[0];
    }

    my $motif_feature_adaptor
        = $self->db->get_adaptor('MotifFeature');

    my @motif_features;

    for my $id (@motif_feature_ids) {
        push @motif_features, $motif_feature_adaptor->fetch_by_dbID($id);
    }

    return \@motif_features;
}

sub _parse_bed_line {

  my $self     = shift;
  my $bed_line = shift;
  
  (
    my $seq_region_id,
    my $seq_region_start,
    my $seq_region_end,
    my $peak_id,
    my $score,
    my $strand,
  ) = split "\t", $bed_line;

  my $hash = {
    seq_region_name  => $seq_region_id,
    seq_region_start => $seq_region_start,
    seq_region_end   => $seq_region_end,
    peak_id          => $peak_id,
    score            => $score,
    strand           => $strand,
  };
  return $hash;
}

sub _bulk_export_to_bed_by_PeakCalling {

  my $self = shift;
  my $peak_calling = shift;
  #my $bed_fh       = shift;
  my @bed_fh       = @_;
  
  if (! defined $peak_calling) {
    throw("Peak calling parameter was undefined!");
  }
  
  my $species = $self->db->species;
  
  if ($species eq 'DEFAULT') {
    die;
  }
  
  my $slice_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'core', 
        'Slice'
    );

  my %seq_region_id_to_name_cache;

  $self->sql_helper->execute_no_return(
    -SQL          => 'select seq_region_id, seq_region_start, seq_region_end, peak_id, score, 0 from peak where peak_calling_id = ?',
    -PARAMS       => [ $peak_calling->dbID ],
    -USE_HASHREFS => 0,
    -CALLBACK     => sub {
        my $row = shift;
        
        my $seq_region_id = $row->[0];
        
        my $seq_region_name;
        if (exists $seq_region_id_to_name_cache{$seq_region_id}) {
          $seq_region_name = $seq_region_id_to_name_cache{$seq_region_id};
        } else {
          my $current_slice = $slice_adaptor->fetch_by_seq_region_id($seq_region_id);
          $seq_region_name = $current_slice->seq_region_name;
          $seq_region_id_to_name_cache{$seq_region_id} = $seq_region_name;
        }
        
        my $bed_line = join "\t", (
          $seq_region_name,
          $row->[1],
          $row->[2],
          $row->[3],
          $row->[4],
          $row->[5],
        );
        
        map { $_->print($bed_line . "\n"); } @bed_fh;
        return;
      },
  );
}

1;
