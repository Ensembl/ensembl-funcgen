=head1 NAME

ensembl_funcgen_set.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

Proserver module to access an ensembl funcgen Result/FeatureSet as a DAS source.
Supports both standard and hydra implementation using the GenerateDASConfig function of the eFG environment.

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


#To do
# Move transport stuff to transport? This module should only be for config parsing and display functions
# Integrate colour ini file parsing from webcode to use same default colours as web display for different
# feature types.
#Add some debug info

package Bio::Das::ProServer::SourceAdaptor::ensembl_funcgen_set;

use warnings;
use strict;
use Data::Dumper;

use base qw(Bio::Das::ProServer::SourceAdaptor);

sub init
{
    my ($self) = @_;

    $self->{'capabilities'} = {
        'features' => '1.0',
        'stylesheet' => '1.0'
        };

	
	#Set the set name here if hydra
	my $set_name;

	if($self->hydra){
	  my $hydraname = $self->config->{'hydraname'};
	  ($set_name = $self->dsn) =~ s/${hydraname}_//;
	}

	#Then explicitly fetch the set using the name
	$self->{'set'} = $self->transport->fetch_set($set_name);


	$self->{'title'} = 	$self->config->{'title'} || $self->set->name;

	#Will this be valid for result_sets?

	$self->{'description'} = $self->config->{'description'} || 
	  ($self->config->{'set_type'} eq 'result') ? $self->set->name.' tiling' : $self->set->display_label;

	return;
}

sub title {
    my $self = shift;
  
    return $self->{'title'};
}

sub set{
  my $self = shift;
  
  #Never required as we don this in init.
  #if(! defined $self->{'set'}){
#	$self->{'set'} = $self->transport->fetch_set;
#  }

  return $self->{'set'};
}


sub build_features{
  my ( $self, $args ) = @_;
  print Dumper $args if ($self->{'debug'});
  
  my $segment = $args->{'segment'};
  my $start   = $args->{'start'};
  my $end     = $args->{'end'};
  
  if(! ($segment && $start && $end)){
	warn "Insufficient loci params:\t${segment}:${start}:${end}";
	return;
  }
  
  my $build_method = 'build_'.$self->config->{'set_type'}.'_set_features';
  return $self->$build_method($segment, $start, $end);
}


sub build_result_set_features{
  my ($self, $segment, $start, $end) = @_;

  my ($id, $label, %score, @features);
  my $slice        = $self->transport->chromosome_by_region($segment, $start, $end);
  print Dumper $slice if ($self->{'debug'});
  my $features     = $self->set->get_ResultFeatures_by_Slice($slice);
  print "Number of features: ".scalar(@{$features})."\n" if ($self->{'debug'});

  my $type     = $self->config()->{'type'} || 'default';
  my $source   =  $self->config()->{'source'} || $self->set->analysis->display_label;
  my $type_cat =  $self->config()->{'typecategory'} || 'result_set';

  foreach my $ft(@{$features}){
 
        my $ft_start = $start+$ft->start();
        my $ft_end = $start+$ft->end();
		
        my $id = sprintf( "%s:%s,%s",
                          $segment,
                          $ft_start,
                          $ft_end);
       
        push @features, {
            
						 'id'          => $id,
						 'label'       => $id,
						 'start'       => $ft_start,
						 'end'         => $ft_end,
						 #'ori'         => 1,
						 'score'       => $ft->score,
						 'method'      => $source,
						 'type'        => $type,
						 'typecategory'=> $type_cat,
						 #'note'        => $note,
						 #'link'        => '',
						 #'linktxt'     => '',
						 
						}

	  };

  return @features;
}



sub build_feature_set_features{
  my ($self, $segment, $start, $end) = @_;

  my ($id, $label, %score, @features);
  my $slice        = $self->transport->chromosome_by_region($segment, $start, $end);
  print Dumper $slice if ($self->{'debug'});
  my $features     = $self->set->get_Features_by_Slice($slice);
  print "Number of features: ".scalar(@{$features})."\n" if ($self->{'debug'});
  my $set_type = $self->config->{'set_type'};


  my $type = $self->config->{'type'}         || 'peak';
  #Need a way of setting this to DNA met, or even reg feats!

  my $cat  = $self->config->{'typecategory'} ||  'feature_set';
  my $label_method = ($set_type eq 'regulatory') ? 'stable_id' : 'display_label';
  my $source   =  $self->config()->{'source'} || $self->set->analysis->display_label;


  foreach my $ft (@{$features}) {

	$id = sprintf( "%s:%s,%s",
				   $segment,
				   $ft->seq_region_start,
				   $ft->seq_region_end,
				 );
	  
	#because we don't want to pass empty score?
	%score = ($set_type eq 'annotated') ? (score => $ft->score) : ();


	push @features,
	  {
	   'id'          => $id,
	   'label'       => $ft->$label_method,
	   #'segment'     => $segment,#Does not return to DAS features test page???
	   'start'       => $ft->seq_region_start,
	   'end'         => $ft->seq_region_end,
	   'ori'         => $ft->seq_region_strand,
	   'method'      => $source,
	   #Can change this to be more informative? Have display label/feature_type here
	   #And extended info(reg atts) in note?
	   'type'        => $type,
	   'typecategory'=> $cat,
	   #'note'        => $ft->display_label,
	   #            'link'        => 'http://www.sanger.ac.uk/PostGenomics/epigenome/',
	   #            'linktxt'     => 'Human Epigenome Project (HEP)'
	   %score
	  };
  }
  
  return @features;
}



sub das_stylesheet{
  my $self = shift;

  
  #need to define this based on config
  #And some defaults
  #Will this be possible with hydra config?
  
  
  #print Dumper $self->{'set_name'};
  
  if ($self->set->isa('Bio::EnsEMBL::Funcgen::ResultSet')){
	return <<EOT;
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
    <STYLESHEET version="0.01">
        <CATEGORY id="result_set">
            <TYPE id="default">
                <GLYPH>
                    <TILING>
<LABEL>no</LABEL>
                        <HEIGHT>30</HEIGHT>
                        <COLOR1>brown3</COLOR1>
                    </TILING>
                </GLYPH>
            </TYPE>
        </CATEGORY>
    </STYLESHEET>
</DASSTYLE>
EOT
  }
  elsif ($self->set->name =~ m/_BR/) {
	return <<EOT;
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
    <STYLESHEET version="0.01">
        <CATEGORY id="feature_set">
            <TYPE id="default">
                <GLYPH>
                    <BOX>
 <LABEL>no</LABEL>
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
        <CATEGORY id="feature_set">
            <TYPE id="default">
                <GLYPH>
                    <BOX>
                        <LABEL>no</LABEL>
                        <FGCOLOR>red</FGCOLOR>
                        <BGCOLOR>black</BGCOLOR>
                    </BOX>
                </GLYPH>
            </TYPE>
            <TYPE id="peak">
                <GLYPH>
                    <BOX>
                        <LABEL>no</LABEL>
                        <FGCOLOR>black</FGCOLOR>
                        <BGCOLOR>brown4</BGCOLOR>
                        <HEIGHT>5</HEIGHT>
                    </BOX>
                </GLYPH>
            </TYPE>
            <TYPE id="DNA methylation">
                <GLYPH>
                    <BOX>
                        <LABEL>no</LABEL>
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

