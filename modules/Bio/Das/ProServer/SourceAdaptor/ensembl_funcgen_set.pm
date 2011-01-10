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

ensembl_funcgen_set.pm

=head1 SYNOPSIS


=head1 DESCRIPTION

Proserver module to access an ensembl funcgen Result/FeatureSet as a DAS source.
Supports both standard and hydra implementation using the GenerateDASConfig function of the eFG environment.

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
use Carp;

use base qw(Bio::Das::ProServer::SourceAdaptor);

sub init
{
    my ($self) = @_;

    $self->{'capabilities'} = {
        'features' => '1.0',
        'stylesheet' => '1.0'
        };
	
	#Set the set name here if hydra
	my $set_name =  $self->config->{set_name};

	if($self->hydra){

	  my $hydraname = $self->config->{'hydraname'};
	  ($set_name = $self->dsn) =~ s/${hydraname}://;
	  my ($cs_level, $set_coords, $cs_version);
	  ($set_name, $cs_version) = split ':', $set_name;
	  #At present we don't strictly need this as dnadb is automatically
	  #set to db with cs_version as default
	  $self->coord_system_version($cs_version);
	  
	  #Set the coords for this particular source
	  my $cs;
	  my %css = %{$self->coordinates};

	  foreach my $coord_sys(keys %css){
		($cs = $coord_sys) =~ s/_//;

		if($cs =~ /${cs_version}\s*,/){
		  $self->{'coordinates'} = { $coord_sys => $css{$coord_sys} };
		  last;
		}
	  }
	}
	  


	#warn Data::Dumper::Dumper($self->coordinates);

	#No coord system stuff for non hydra as we expect dnadb to be configured correctly
	
	#Then explicitly fetch the set using the name
	$self->{'set'} = $self->transport->fetch_set($set_name);

	#Assembly version used in full source name
	#Title is used for display
	#Preferable that start of title should be enough to discriminate between sets.
	$self->{'title'} = $self->config->{'title'} || $set_name;
	#This is not picked up by SourceAdaptor::das_sourcedata hence we get warning about undef in sprintf!

	#More word humanly readable description
	#only shown on sources page
	$self->{'description'} = $self->config->{'description'} ||  $self->set->display_label;

	#warn "config coord are ".join(', ', keys(%{$self->config->{'coordinates'}}));
	warn $self->{'description'}." coords are ".join(', ', keys(%{$self->{'coordinates'}}))."\n" if $self->{debug};

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


sub coord_system_version{
  my ($self, $version) = @_;

  if (defined $version){
	$self->{'_coord_system_version'} = $version;
  }

  return $self->{'_coord_system_version'};
}



sub build_features{
  my ( $self, $args ) = @_;
  warn "build_features args are:\n".Dumper($args) if ($self->{'debug'});
  
  my $segment = $args->{'segment'};
  my $start   = $args->{'start'};
  my $end     = $args->{'end'};
  
  #Need to implement maxbins here for result_set feature

  if(! ($segment && $start && $end)){
	warn "Insufficient loci params:\t${segment}:${start}:${end}";
	return;
  }
  
  my $build_method = 'build_'.$self->config->{'set_type'}.'_set_features';

  my $dnadb = $self->transport->adaptor->dnadb;
  my $slice = $dnadb->get_SliceAdaptor->fetch_by_region('chromosome', $segment, $start, $end, undef, $self->coord_system_version);
  warn "Slice is:\t".$slice->name."\n" if ($self->{'debug'});

  return $self->$build_method($slice, $args);
}


sub build_result_set_features{
  my ($self, $slice, $args) = @_;

  my ($id, $label, %score, @features);
  my $features = $self->set->get_ResultFeatures_by_Slice($slice, undef, undef, $args->{'maxbins'});

  

  warn "Max bins:\t".$args->{'maxbins'}."\nNumber of features: ".scalar(@{$features}) if $self->{debug};


  if(@$features){
	my $bin_size = $features->[0]->window_size;
	warn "Using bin size:\t$bin_size"  if $self->{'debug'};


	my $type     = $self->config()->{'type'} || 'default';
	my $source   = $self->config()->{'source'} || $self->set->analysis->display_label;
	my $type_cat = $self->config()->{'typecategory'} || 'result_set';
	my $start    = $slice->start; 
	my $segment  = $slice->seq_region_name;
	my $method   = $self->config()->{'source'} || $self->set->analysis->display_label;
	
	if(! defined $method){
	  croak('Cannot determine mandatory \'method\' attribute. Please set \'source\' in DAS config or analysis display_label');
	}


	#We are currently storing everything as 1(+)
	#We should store actual orientation but
	#set to 0 here to get the display on the same strand as the peaks?
	
	my %ori = (
			   -1 => '-',
			   0  => '0',
			   1  => '+',
			  );

		
	foreach my $ft(@{$features}){
	  
	  #Set seq_region_start/end as we don't have direct access 
	  #using the current ResultFeature class
	  
	  my $ft_start = $start + $ft->start();
	  my $true_end = $start + $ft->end();
	  my $ft_end   = $ft_start - 1;
	  
	  #Can we use seq_region_start/end here or do they have to be local to the slice?
	  #If not they we can change as ResultFeature is now hash based
	  #Can only do this for 0 wsize! So not much point
		
	  foreach my $score(@{$ft->scores}){

		if($bin_size == 0){
		  $ft_end = $true_end;
		}
		else{
		  $ft_end += $bin_size;
		}
		
		
		my $id = sprintf( "%s:%s,%s",
						  $segment,
						  $ft_start,
						  $ft_end);
		
		#warn "Got score $score";
		
		
		push @features, {
						 
						 'id'          => $id,
						 'label'       => $id,
						 'start'       => $ft_start,
						 'end'         => $ft_end,
						 'ori'         => $ori{$ft->strand},#change this to seq_region_strand when we cahnge to hash feature?
						 'score'       => $score,
						 'method'      => $source,
						 'type'        => $type,
						 'typecategory'=> $type_cat,
						 #'note'        => $note,
						 #'link'        => '',
						 #'linktxt'     => '',
						 
						 #Mandatory attrs
						 'method' => $method,
						 'phase'  => '-', 
						};					   
		
		$ft_start += $bin_size;
	  }
	}
  }

  print "Returning ".scalar(@features).' '.$self->set->name." ResultFeatures\n" if $self->{debug};

  return @features;
}



sub build_feature_set_features{
  my ($self, $slice) = @_;

  my ($id, $label, %score, @features);
  my $features     = $self->set->get_Features_by_Slice($slice);
  print "Number of features: ".scalar(@{$features})."\n" if ($self->{'debug'});
  my $set_type = $self->config->{'set_type'};
  my $start    = $slice->start;
  my $segment  = $slice->seq_region_name;

  my $type = $self->config->{'type'}         || 'peak';
  #Need a way of setting this to DNA met, or even reg feats!

  my $cat  = $self->config->{'typecategory'} ||  'feature_set';
  my $label_method = ($set_type eq 'regulatory') ? 'stable_id' : 'display_label';
  my $source   =  $self->config()->{'source'} || $self->set->analysis->display_label;  


  if(! defined $source){
	croak('Cannot determine mandatory \'method\' attribute. Please set \'source\' in DAS config or analysis display_label');
  }


  my %ori = (
			 -1 => '-',
			 0  => '0',
			 1  => '+',
			);

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
	   'ori'         => $ori{$ft->seq_region_strand},
	   'method'      => $source,
	   'phase'       => '-',
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
  
  warn "Returning ".scalar(@features).' '.$self->set->name." features\n" if $self->{debug};

  return @features;
}



sub das_stylesheet{
  my $self = shift;

  
  #need to define this based on config
  #And some defaults
  #Will this be possible with hydra config?
  
  
  #print Dumper $self->{'set_name'};
  
  #LABEL here is not in spec for HISTOGRAM/colour gradient
  #Is used as Ensembl was adding labels by default which was not the correct behaviour
  #hence this prevents DAS validation
  #<LABEL>no</LABEL>
  if ($self->set->isa('Bio::EnsEMBL::Funcgen::ResultSet')){
	return <<EOT;
<!DOCTYPE DASSTYLE SYSTEM "http://www.biodas.org/dtd/dasstyle.dtd">
<DASSTYLE>
    <STYLESHEET version="0.01">
        <CATEGORY id="result_set">
            <TYPE id="default">
                <GLYPH>
                    <HISTOGRAM>
                        <HEIGHT>30</HEIGHT>
                        <COLOR1>brown3</COLOR1>
                    </HISTOGRAM>
                </GLYPH>
            </TYPE>
        </CATEGORY>
    </STYLESHEET>
</DASSTYLE>
EOT
  }
  elsif ($self->set->name =~ m/_BR/) {#? Needs updating? to match FeatureSet?
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

