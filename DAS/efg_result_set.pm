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

package Bio::Das::ProServer::SourceAdaptor::efg_result_set;

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

sub build_features
{
    my ( $self, $args ) = @_;
    #print Dumper $args;
    
    my $segment = $args->{'segment'} || return ();
    my $start = $args->{'start'} || return ();
    my $end   = $args->{'end'} || return ();


    my $dba = $self->transport->adaptor();
    my $slice = $dba->get_SliceAdaptor->fetch_by_region
        ('chromosome', $segment, $start, $end);
    print Dumper $slice if ($self->{'debug'});

    my $rset = $dba->get_ResultSetAdaptor->fetch_by_dbID
        ($self->config()->{'result_set_id'});

    my @features;

    map {
        my $ft_start = $start+$_->start;
        my $ft_end = $start+$_->end;
        my $id = sprintf( "%s:%s,%s",
                          $segment,
                          $ft_start,
                          $ft_end);
        print Dumper $id if ($self->{'debug'});
        push @features, {
            
            'id'          => $id,
            'label'       => $id,
            'start'       => $ft_start,
            'end'         => $ft_end,
            #'ori'         => 1,
            'score'       => $_->score,
            'method'      => $self->config()->{'source'},
            'type'        => $self->config()->{'type'},
            'typecategory'=> $self->config()->{'category'},
            #'note'        => $note,
            #'link'        => 'http://www.sanger.ac.uk/PostGenomics/epigenome/',
            #'linktxt'     => 'Human Epigenome Project (HEP)'
            
        }

    } @{$rset->get_ResultFeatures_by_Slice($slice)};
    print 'NoRF: '.scalar(@features)."\n";

    return @features;



#    my $dset = $dba->get_DataSetAdaptor->fetch_by_dbID
#        ($self->config()->{'data_set_id'});
#    #print Dumper $dset->get_supporting_sets;
#
#    my @features;
#
#    foreach my $rset (@{$dset->get_supporting_sets()}) 
#    {
#        
#        #print Dumper $rset->get_ResultFeatures_by_Slice($slice);
#        map {
#            #print Dumper $_;
#            my $ft_start = $start+$_->start;
#            my $ft_end = $start+$_->end;
#            my $id = sprintf( "%s:%s,%s",
#                              $segment,
#                              $ft_start,
#                              $ft_end);
#            #print Dumper $id;
#            push @features, {
#              
#                'id'          => $id,
#                'label'       => $id,
#                'start'       => $ft_start,
#                'end'         => $ft_end,
#                #'ori'         => 1,
#                'score'       => $_->score,
#                'method'      => $self->config()->{'source'},
#                'type'        => $self->config()->{'type'},
#                'typecategory'=> $self->config()->{'category'},
#                #'note'        => $note,
#                #'link'        => 'http://www.sanger.ac.uk/PostGenomics/epigenome/',
#                #'linktxt'     => 'Human Epigenome Project (HEP)'
#
#            } } @{$rset->get_ResultFeatures_by_Slice($slice)};
#
#        return @features;
#        
#    }


#    my $type = $self->config()->{'type'};
#    my $feature_set_id = $self->config()->{'feature_set_id'};
#    print Dumper $feature_set_id if ($self->{'debug'});
#
#    # get data
#    
#    my $qbounds = ($start && $end)?qq(AND seq_region_start <= $end AND seq_region_end >= $start):"";
#    my $version = $self->config()->{'version'};
#    my $sql = "SELECT *
#              FROM ${type}_feature f, seq_region sr 
#             WHERE sr.seq_region_id=f.seq_region_id
#               AND sr.name=\"$segment\"
#               AND sr.schema_build=\"$version\"
#               AND f.feature_set_id=$feature_set_id
#               $qbounds";
#    
#    print Dumper $sql if ($self->{'debug'});
#
#    my $features = $self->transport->query($sql);
#    
#    print Dumper $features if ($self->{'debug'});
#    #print "Number of features: ", scalar(@{$features}), "\n";
#
#    foreach my $ft (@{$features}) {
#
#        my $id = sprintf( "%s:%s,%s",
#                          $segment,
#                          $ft->{'seq_region_start'},
#                          $ft->{'seq_region_end'}
#                          );
#
#        my $note = $ft->{'display_label'};
#        
#        push @features,
#        {
#            'id'          => $id,
#            'label'       => $id,
#            'start'       => $ft->{'seq_region_start'},
#            'end'         => $ft->{'seq_region_end'},
#            'ori'         => $ft->{'seq_region_strand'},
#            'score'       => $ft->{'score'},
##            #'method'      => $self->config()->{'source'},
#            'type'        => $type,
##            #'typecategory'=> $self->config()->{'category'},
#            'note'        => $note,
##            'link'        => 'http://www.sanger.ac.uk/PostGenomics/epigenome/',
##            'linktxt'     => 'Human Epigenome Project (HEP)'
#        };
#    }
#
}

sub das_stylesheet
{
    my $self = shift;

    return <<EOT;
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
    <STYLESHEET version="0.01">
        <CATEGORY id="default">
            <TYPE id="default">
                <GLYPH>
                    <BOX>
                        <FGCOLOR>brown3</FGCOLOR>
                        <BGCOLOR>brown3</BGCOLOR>
                        <HEIGHT>30</HEIGHT>
                    </BOX>
                </GLYPH>
            </TYPE>
        </CATEGORY>
    </STYLESHEET>
</DASSTYLE>
EOT
}

1;
