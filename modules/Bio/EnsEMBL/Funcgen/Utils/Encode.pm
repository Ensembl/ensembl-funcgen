# Personal Ensembl module for Bio::EnsEMBL::Funcgen::Utils::Encode
#
# Copyright (c) 2007 Ensembl
#

=head1 NAME

Bio::EnsEMBL::Funcgen::Utils::Encode - provides some handy methods to deal 
with Encode regions.

=head1 SYNOPSIS

    use Bio::EnsEMBL::Funcgen::Utils::Encode qw(get_encode_regions);

    my $encode_regions = &get_encode_regions($dnadb, $assembly_version);

=head1 DESCRIPTION

=head1 AUTHOR

Stefan Graf <graef@ebi.ac.uk>, EMBL-EBI, Ensembl Functional Genomics
This module is part of the Ensembl project: http://www.ensembl.org

=head1 CONTACT

Post questions to the Ensembl development list: ensembl-dev@ebi.ac.uk

=cut

use strict;
use warnings;

package Bio::EnsEMBL::Funcgen::Utils::Encode;

require Exporter;
use vars qw(@ISA @EXPORT_OK);
@ISA = qw(Exporter);
@EXPORT_OK = qw(get_encode_regions);

use Bio::EnsEMBL::Utils::Exception qw(throw);
use POSIX;

=head2 get_encode_regions

  Arg [1]     : Bio::EnsEMBL::DBSQL::DBAdaptor
  Arg [2]     : assembly version string, e.g. NCBI35 (optional)
  Description : fetch encode region coordinates from core db
  Returns     : reference to hash with Encode region names as keys
  Exceptions  : 
  Example     : 

=cut

sub get_encode_regions {

    my ($db, $version) = @_;
    
    $version = undef if (! @_);

    my @chromosomes=qw(1 2 3 4 5 6 7 8 9 10 11 12 13 14 
                       15 16 17 18 19 20 21 22 X Y MT);

    my $sa = $db->get_SliceAdaptor();
    my $mfa = $db->get_MiscFeatureAdaptor();

    my @mf;
    foreach my $c (@chromosomes)
    {
        my $s = $sa->fetch_by_region
            ('chromosome', $c, undef, undef, undef, $version);
        push @mf, @{$mfa->fetch_all_by_Slice_and_set_code($s, 'encode')};
    }
    
    my %encode_regions;
    foreach my $mf (@mf) 
    {
        $encode_regions{$mf->display_id} = sprintf
            ("%s:%s:%s:%d:%d:%d",
             $mf->slice->coord_system_name, 
             $mf->slice->coord_system()->version(),
             $mf->slice->seq_region_name,
             $mf->start, $mf->end, $mf->strand);
    }
    
    return \%encode_regions;

}

1;
