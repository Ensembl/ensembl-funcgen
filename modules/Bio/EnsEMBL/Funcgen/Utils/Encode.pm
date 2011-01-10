# Ensembl module for Bio::EnsEMBL::Funcgen::Utils::Encode
#

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.


=head1 NAME

Bio::EnsEMBL::Funcgen::Utils::Encode - provides some handy methods to deal 
with Encode regions.

=head1 SYNOPSIS

    use Bio::EnsEMBL::Funcgen::Utils::Encode qw(get_encode_regions);

    my $encode_regions = &get_encode_regions($dnadb, $assembly_version);

=head1 DESCRIPTION

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
