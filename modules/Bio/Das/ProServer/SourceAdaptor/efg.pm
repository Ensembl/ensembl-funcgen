#########
# Author:  Stefan Graf
# Maintainer: sg550@cam.ac.uk
# Created: 2009-05-11
# Last Modified: 2009-05-11
# Builds DAS features from eFG databases
#
package Bio::Das::ProServer::SourceAdaptor::efg;

=head1 AUTHOR

Stefan Graf <sg550@cam.ac.uk>

Copyright (c) 2009 Department of Oncology, University of Cambridge

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

use strict;
use vars qw(@ISA);
use Data::Dumper;
use Bio::Das::ProServer::SourceAdaptor;
use Bio::EnsEMBL::Utils::Exception qw(stack_trace_dump);
@ISA = qw(Bio::Das::ProServer::SourceAdaptor);

sub init {
    my $self                = shift;
	
	#warn "<<< adaptor::init >>>";

    $self->{'capabilities'} = {
        'features' => '1.0',
        'stylesheet' => '1.0',
    };
    #$self->{'dsnversion'} = '1.1';
    #$self->{'coordinates'} = {
    #	'http://www.dasregistry.org/dasregistry/coordsys/CS_DS40' => '7:140079754,140272033'
    #};
}

sub title {
    my $self = shift;
	my $hydraname = $self->{'config'}->{'hydraname'};
    my ($title) = $self->dsn =~ m/${hydraname}_(.+)/;
    return "$title";
}

sub build_features {
    my ($self, $args) = @_;
	my @features;
	
	print Dumper $args if $self->{'debug'};
	my $segment = $args->{'segment'};
    my $start   = $args->{'start'};
    my $end     = $args->{'end'};
	my $slice   = $self->transport->slice_adaptor->fetch_by_region('chromosome', $segment, $start, $end);

    ### using max bins
    #if( $opts->{'maxbins'} && $gStart && $gEnd) {
    #    return $self->merge_features($opts);
    #}

	my $hydraname = $self->{'config'}->{'hydraname'};
    my $dsn       = $self->{'dsn'};
	(my $fset_name = $dsn) =~ s/^${hydraname}_//;

	print 'FeatureSet: ', Dumper $fset_name if $self->{'debug'};
	
	

	my $feature_set = $self->transport->feature_set_adaptor->fetch_by_name($fset_name);
	#print Dumper $feature_set if $self->{'debug'};

	my $id=1;
	foreach (@{$feature_set->get_Features_by_Slice($slice)}) {

		my $ft = $_->transform('chromosome');

        push @features, {
            'id'      => $id++,
            'label'   => $ft->display_label,
            #'method'  => $ft->feature_type->name,
            'type'    => 'default',
            'typecategory'=> 'default',
            'start'  => $ft->start,
            'end'    => $ft->end,
            'ori' => $ft->strand,
            'score' => $ft->score,
            #'note' => [],
        };

	}
	print Dumper @features;

    print "No. of features: ".scalar(@features)."\n";
    return @features;
}

sub das_stylesheet
{
    my ($self) = @_;
	
	my $feature_set = $self->transport->feature_set_adaptor->fetch_by_name($self->{'dsn'});
	my $feature_type = $feature_set->feature_type->name;

	if ($feature_type eq "Segment") {

        return <<EOT;
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
    <STYLESHEET version="1.0">
        <CATEGORY id="default">
            <TYPE id="default">
                <GLYPH>
                    <BOX>
                        <HEIGHT>30</HEIGHT>
                        <BGCOLOR>black</BGCOLOR>
                        <FGCOLOR>red</FGCOLOR>
                    </BOX>
                </GLYPH>
            </TYPE>
        </CATEGORY>
    </STYLESHEET>
</DASSTYLE>
EOT

	}

}


1;

=head1 NAME

Bio::Das::ProServer::SourceAdaptor::Transport::efg

=head1 AUTHOR

Stefan Graf <graef@ebi.ac.uk or sg550@cam.ac.uk>

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2009 EMBL-EBI and CRUK-CRI

=cut
