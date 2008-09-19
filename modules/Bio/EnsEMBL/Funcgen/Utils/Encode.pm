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
use Data::Dumper;

=head2 get_encode_regions

  Arg [1]     : Bio::EnsEMBL::DBSQL::DBAdaptor
  Description : fetch encode region coordinates from core db
  Returns     : reference to list of Bio::EnsEMBL::MiscFeatures Encode region slices
  Exceptions  : Throws if no Bio::EnsEMBL::DBSQL::DBAdaptor defined
  Example     : get_encode_regions($dnadb);

=cut

sub get_encode_regions {

    my $db = shift;
    
    throw("Need to pass a valid Bio::EnsEMBL::DBSQL::DBAdaptor") 
        if (! ($db  && $db->isa("Bio::EnsEMBL::DBSQL::DBAdaptor")));

    my $sa = $db->get_SliceAdaptor();
    my $tls = $sa->fetch_all('toplevel');
    #map { print Dumper $_->name } @$tls;

    my $mfa = $db->get_MiscFeatureAdaptor();

    my @encode_regions;
    map { 
        push @encode_regions, @{$mfa->fetch_all_by_Slice_and_set_code($_, 'encode')};
    } @$tls;

    return \@encode_regions;


    #my %encode_regions;
    #map { $encode_regions{$_->display_id} = sprintf
    #          ("%s:%s:%s:%d:%d:%d",
    #           $_->slice->coord_system_name, 
    #           $_->slice->coord_system()->version(),
    #           $_->slice->seq_region_name,
    #           $_->start, $_->end, $_->strand);
    #  } @encode_regions;
    #
    #return \%encode_regions;

}

1;
