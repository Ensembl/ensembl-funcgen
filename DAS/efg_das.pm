#########
# Author:  Stefan Graf
# Maintainer: graef@ebi.ac.uk
# Created: 2008-03-13
# Last Modified: 2008-03-13
# Builds DAS features for sequencing reads
#
package Bio::Das::ProServer::SourceAdaptor::efg_das;

=head1 AUTHOR

Stefab Graf <graef@ebi.ac.uk>

Copyright (c) 2008 EBI

This library is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.  See DISCLAIMER.txt for
disclaimers of warranty.

=cut

use strict;
use vars qw(@ISA);
use Bio::Das::ProServer::SourceAdaptor;
@ISA = qw(Bio::Das::ProServer::SourceAdaptor);

sub init {
    my $self                = shift;
    $self->{'capabilities'} = {
        'features' => '1.0',
        'stylesheet' => '1.0',
    };
    $self->{'dsnversion'} = '1.1';
    #$self->{'coordinates'} = {
    #	'http://www.dasregistry.org/dasregistry/coordsys/CS_DS40' => '7:140079754,140272033'
    #};
}

sub title {
    my $self = shift;
    my ($title) = $self->dsn =~ m/heroic_(.+)/;
    return "$title";
}

sub build_features {
    my ($self, $opts) = @_;
    
    my $segment       = $opts->{'segment'};
    my $gStart        = $opts->{'start'};
    my $gEnd          = $opts->{'end'};

    ### using max bins
    #if( $opts->{'maxbins'} && $gStart && $gEnd) {
    #    return $self->merge_features($opts);
    #}

    my $dsn           = $self->{'dsn'};
    my $dbtable       = $dsn;
    
    #########
    # if this is a hydra-based source the table name contains the hydra name and needs to be switched out
    #
    #my $hydraname     = $self->config->{'hydraname'};
    #if ($hydraname) {
    #    my $basename = $self->config->{'basename'};
    #    $dbtable =~ s/$hydraname/$basename/;
    #}
    
    my $qsegment      = $self->transport->dbh->quote($segment);
    my $qbounds       = qq(AND start <= '$gEnd' AND end >= '$gStart') if($gStart && $gEnd);
    my $query         = qq(SELECT * FROM $dbtable WHERE  seq_region = $qsegment $qbounds); # ORDER BY start);
    my $ref           = $self->transport->query($query);
    my @features      = ();
    
    for my $row (@{$ref}) {

        my ($start, $end, $score) = ($row->{'start'}, $row->{'end'}, $row->{'score'});

        my $type;
        if ($self->{'dsn'} =~ m/_profile/) {
            $type = 'default';
        } elsif ($self->{'dsn'} =~ m/_reads/) {
            $type = $row->{'strand'} eq '+' ? 'forward' : 'reverse';
            $row->{'strand'} = ''
        } else {
            $type = $score > 0? 'uniq' : 'non-uniq';
        }

        push @features, {
            'id'      => $row->{'feature_id'},
            'label'   => $row->{'name'} !~ m/^.$/ ? $row->{'name'} : $row->{'feature_id'},
            'method'  => $self->config()->{'method'},
            'type'    => $type,
            'typecategory'=> $self->config()->{'category'},
            'start'  => $row->{'start'},
            'end'    => $row->{'end'},
            'ori' => $row->{'strand'},
            'score' => $row->{'score'},
            #'note' => [],
        };

    }       
    print "No. of features: ".scalar(@features)."\n";
    return @features;
}

sub das_stylesheet
{
    my $self = shift;

    if ($self->{'dsn'} =~ m/_profile/) {
        return <<EOT;
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
    <STYLESHEET version="1.0">
        <CATEGORY id="sequencing">
            <TYPE id="default">
                <GLYPH>
                    <TILING>
                        <HEIGHT>30</HEIGHT>
                        <COLOR1>brown4</COLOR1>
                        <MIN>0</MIN>
                        <MAX>15</MAX>
                    </TILING>
                </GLYPH>
            </TYPE>
        </CATEGORY>
    </STYLESHEET>
</DASSTYLE>
EOT

    } else {
        return <<EOT;
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
    <STYLESHEET version="1.0">
        <CATEGORY id="sequencing">
            <TYPE id="forward">
                <GLYPH>
                    <ARROW>
                        <HEIGHT>10</HEIGHT>
                        <FGCOLOR>black</FGCOLOR>
                        <BUMP>1</BUMP>
                        <BAR_STYLE>full</BAR_STYLE>
                        <SOUTHWEST>no</SOUTHWEST>
                    </ARROW>
                </GLYPH>
            </TYPE>
            <TYPE id="reverse">
                <GLYPH>
                    <ARROW>
                        <HEIGHT>10</HEIGHT>
                        <FGCOLOR>darkgreen</FGCOLOR>
                        <BUMP>1</BUMP>
                        <BAR_STYLE>full</BAR_STYLE>
                        <NORTHEAST>no</NORTHEAST>
                    </ARROW>
                </GLYPH>
            </TYPE>
        </CATEGORY>
    </STYLESHEET>
</DASSTYLE>
EOT
    }
}

1;

