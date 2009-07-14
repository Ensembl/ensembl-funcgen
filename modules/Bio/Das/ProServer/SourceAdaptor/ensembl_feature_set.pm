=head1 NAME

ensembl_feature_set.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

Proserver module to access a DAS source from an efg FeatureSet

=head1 LICENSE

  Copyright (c) 1999-2009 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=cut

package Bio::Das::ProServer::SourceAdaptor::ensembl_feature_set;

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

	my $segment = $args->{'segment'};
    my $start = $args->{'start'};
    my $end   = $args->{'end'};

	if(! ($segment && $start && $end)){
	  warn "Insufficient loci params:\t${segment}:${start}:${end}";
	  return;
	}

    my $type = $self->config()->{'type'};
    #my $feature_set_id = $self->config()->{'set_id'};
	my $fset = $self->transport->adaptor->get_FeatureSetAdaptor->fetch_by_dbID($self->config()->{'set_id'});
    #print Dumper $feature_set_id if ($self->{'debug'});

    # get data
    


	#Data version/schema_build?
	#Was defined by param when writing config
	#We don't necessarily need this
	#But we do need to know which assembly to use?
	#Just use the FeatureSet directly.
	#Can we get the Adaptor to generate lihgt weight DAS features directly?
	#Therefore removing the duplication of feature generation?
	#This may be very handy with the read data. Let's test the implementation there.
	


    #my $qbounds = ($start && $end)?qq(AND seq_region_start <= $end AND seq_region_end >= $start):"";
    #my $version = $self->config()->{'data_version'};
    #print Dumper $version if ($self->{'debug'});
    #my $sql = "SELECT *
     #         FROM ${type}_feature f, seq_region sr 
     #        WHERE sr.seq_region_id=f.seq_region_id
     #          AND sr.name=\"$segment\"
     #          AND sr.schema_build=\"$version\"
     #          AND f.feature_set_id=$feature_set_id
     #          $qbounds";
    
    #print Dumper $sql if ($self->{'debug'});

    #my $features = $self->transport->query($sql);
    
    #print Dumper $features if ($self->{'debug'});
    #print "Number of features: ", scalar(@{$features}), "\n";



	my $slice    = $self->transport->chromosome_by_region($segment, $start, $end);

	#Could change method here to reflect other params
	#e.g. filter by type etc.
	my $features = $fset->get_Features_by_Slice($slice);
	print "Number of features: ".scalar(@{$features})."\n" if ($self->{'debug'});

	

    my ($id, $label, %score, @features);

    foreach my $ft (@{$features}) {

	  $id = sprintf( "%s:%s,%s",
					 $segment,
					 $ft->seq_region_start,
					 $ft->seq_region_end,
				   );
	  

	  #Need to build label and score dpendant on type
	  $label = ($type eq 'regulatory') ? $ft->stable_id : $id;
	  #elsif($type eq 'external'){
	  #}
	  %score = ($type eq 'annotated') ? (score => $ft->score) : ();

	  push @features,
        {
		 'id'          => $id,
		 'label'       => $label,
		 #'segment'     => $segment,#Does not return to DAS features test page???
		 'start'       => $ft->seq_region_start,
		 'end'         => $ft->seq_region_end,
		 'ori'         => $ft->seq_region_strand,
		 #            #'method'      => $self->config()->{'source'},

		 #Can change this to be more informative? Have display label/feature_type here
		 #And extended info(reg atts) in note?
		 'type'        => $type,

		 #            #'typecategory'=> $self->config()->{'category'},
		 'note'        => $ft->display_label,
		 #            'link'        => 'http://www.sanger.ac.uk/PostGenomics/epigenome/',
		 #            'linktxt'     => 'Human Epigenome Project (HEP)'
		 %score
        };
    }

    return @features;
}

sub das_stylesheet
{
    my $self = shift;


	#need to define this based on config
	#And some defaults
	#Will this be possible with hydra config?


    #print Dumper $self->{'set_name'};

    if ($self->{'set_name'} =~ m/_BR/) {
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
