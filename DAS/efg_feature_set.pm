=head1 NAME

efg_feature_set.pm

=head1 SYNOPSIS

Proserver module to run DAS server off the efg FeatureSet table

=head1 DESCRIPTION

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf <graef@ebi.ac.uk>, Ensembl Functional Genomics

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

package Bio::Das::ProServer::SourceAdaptor::efg_feature_set;

use warnings;
use strict;

use Data::Dumper;
use IO::File;

use base qw(Bio::Das::ProServer::SourceAdaptor);

sub init
{
    my ($self) = @_;

    $self->{'capabilities'} = {
        'features' => '1.0',
        'stylesheet' => '1.0'
        };



}

sub title {
    my $self = shift;
    my ($title) = $self->dsn;

    $title =~ s/\./_/;

    return "$title";
}

sub build_features
{
    my ( $self, $args ) = @_;
    print Dumper $args if ($self->{'debug'});
    
    my $segment = $args->{'segment'} || return ();
    my $start = $args->{'start'} || return ();
    my $end   = $args->{'end'} || return ();


    my $type = $self->config()->{'feature_set_type'};
    my $feature_set_id = $self->config()->{'feature_set_id'};
    print Dumper $feature_set_id if ($self->{'debug'});

    # get data
    
    my $qbounds = ($start && $end)?qq(AND seq_region_start <= $end AND seq_region_end >= $start):"";
    my $version = $self->config()->{'data_version'};
    print Dumper $version if ($self->{'debug'});
    my $sql = "SELECT *
              FROM ${type}_feature f, seq_region sr 
             WHERE sr.seq_region_id=f.seq_region_id
               AND sr.name=\"$segment\"
               AND sr.schema_build=\"$version\"
               AND f.feature_set_id=$feature_set_id
               $qbounds";
    
    print Dumper $sql if ($self->{'debug'});

    my $features = $self->transport->query($sql);
    
    print Dumper $features if ($self->{'debug'});
    #print "Number of features: ", scalar(@{$features}), "\n";

    my @features;
    foreach my $ft (@{$features}) {

        my $id = sprintf( "%s:%s,%s",
                          $segment,
                          $ft->{'seq_region_start'},
                          $ft->{'seq_region_end'}
                          );

        my $note = $ft->{'display_label'};
        
        push @features,
        {
            'id'          => $id,
            'label'       => $id,
            'start'       => $ft->{'seq_region_start'},
            'end'         => $ft->{'seq_region_end'},
            'ori'         => $ft->{'seq_region_strand'},
            'score'       => $ft->{'score'},
#            #'method'      => $self->config()->{'source'},
            'type'        => $type,
#            #'typecategory'=> $self->config()->{'category'},
            'note'        => $note,
#            'link'        => 'http://www.sanger.ac.uk/PostGenomics/epigenome/',
#            'linktxt'     => 'Human Epigenome Project (HEP)'
        };
    }

    print 'NoF: '.scalar(@features)."\n";
    return @features;
}

sub das_stylesheet
{
    my $self = shift;

    print Dumper $self->{'dsn'};

    if ($self->{'dsn'} =~ m/_BR/) {
        return <<EOT;
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
    <STYLESHEET version="0.01">
        <CATEGORY id="default">
            <TYPE id="default">
                <GLYPH>
                    <BOX>
                        <FGCOLOR>black</FGCOLOR>
                        <BGCOLOR>brown3</BGCOLOR>
                        <HEIGHT>5</HEIGHT>
                    </BOX>
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
    <STYLESHEET version="0.01">
        <CATEGORY id="default">
            <TYPE id="default">
                <GLYPH>
                    <BOX>
                        <FGCOLOR>red</FGCOLOR>
                        <BGCOLOR>black</BGCOLOR>
                    </BOX>
                </GLYPH>
            </TYPE>
            <TYPE id="annotated">
                <GLYPH>
                    <BOX>
                        <FGCOLOR>black</FGCOLOR>
                        <BGCOLOR>brown4</BGCOLOR>
                        <HEIGHT>5</HEIGHT>
                    </BOX>
                </GLYPH>
            </TYPE>
            <TYPE id="DNA methylation">
                <GLYPH>
                    <BOX>
                        <FGCOLOR>blue</FGCOLOR>
                        <BGCOLOR>black</BGCOLOR>
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
