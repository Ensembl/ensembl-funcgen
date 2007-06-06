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

sub build_features
{
    my ( $self, $args ) = @_;
    #print Dumper $args;
    
    my $segment = $args->{'segment'} || return ();
    my $start = $args->{'start'} || return ();
    my $end   = $args->{'end'} || return ();


    use Bio::EnsEMBL::DBSQL::DBAdaptor;
    my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
        (
         -user   => 'ensro',
         -dbname => 'homo_sapiens_core_44_36f',
         -host   => 'ens-livemirror',
         );
    my $sa = $cdb->get_SliceAdaptor();
    my $slice = $sa->fetch_by_region('chromosome', "$segment");

    my $seq_region_id = $slice->get_seq_region_id();
    print Dumper $seq_region_id if ($self->{'debug'});
    my ($name, $version) = (split(/:/, $slice->name()))[0,1];
    print Dumper ($name, $version) if ($self->{'debug'});
    (my $schema_build = $self->config()->{'dbname'}) =~ s/.*_(\d+_\d+[a-z])$/$1/;
    print Dumper ($schema_build) if ($self->{'debug'});

    my $sql = "SELECT coord_system_id
                 FROM coord_system
                WHERE name='$name'
                  AND version='$version'
                  AND schema_build='$schema_build'
                  AND attrib='default_version'";
    my $coord_system_id = $self->transport->query($sql)->[0]->{coord_system_id};
    print Dumper $coord_system_id if ($self->{'debug'});

    my $feature_set = $self->config()->{'feature_set'};
    print Dumper $feature_set if ($self->{'debug'});

    # get data
    
    my $qbounds = ($start && $end)?qq(AND seq_region_start <= $end AND seq_region_end >= $start):"";

    $sql = "SELECT *
              FROM feature_set fs, predicted_feature pf 
             WHERE fs.feature_set_id=pf.feature_set_id 
               AND fs.name=\"$feature_set\"
               AND coord_system_id=$coord_system_id
               AND seq_region_id=$seq_region_id $qbounds";
    
    #print Dumper $sql;

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

        my $note = $feature_set;
        if ($note =~ m/^(Overlap_|Regulatory)/) {
            $note = $ft->{'display_label'};
        }
        
        push @features,
        {
            'id'          => $id,
            'label'       => $id,
            'start'       => $ft->{'seq_region_start'},
            'end'         => $ft->{'seq_region_end'},
            'ori'         => $ft->{'seq_region_strand'},
            'score'       => $ft->{'score'},
#            #'method'      => $self->config()->{'source'},
#            #'type'        => $ft->{'type'},
#            #'typecategory'=> $self->config()->{'category'},
            'note'        => $note,
#            'link'        => 'http://www.sanger.ac.uk/PostGenomics/epigenome/',
#            'linktxt'     => 'Human Epigenome Project (HEP)'
        };
    }

    return @features;
}

sub das_stylesheet
{
    my $self = shift;

    return <<EOT;
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
    <STYLESHEET version="0.01">
        <CATEGORY id="epigenomic modification">
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

1;
