#
# EnsEMBL module for Bio::EnsEMBL::Funcgen::Parsers::Sanger
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

Bio::EnsEMBL::Funcgen::Parsers::Sanger

=head1 SYNOPSIS

  my $parser_type = "Bio::EnsEMBL::Funcgen::Parsers::Sanger";
  push @INC, $parser_type;
  my $imp = $class->SUPER::new(@_);


=head1 DESCRIPTION

This is a definitions class which should not be instatiated directly, it 
normally inherited from the Importer.  Sanger contains meta data and methods 
specific to Sanger PCR arrays to aid parsing and importing of experimental data.

=cut

package Bio::EnsEMBL::Funcgen::Parsers::Sanger;

use Bio::EnsEMBL::Funcgen::Array;
use Bio::EnsEMBL::Funcgen::ProbeSet;
use Bio::EnsEMBL::Funcgen::Probe;
use Bio::EnsEMBL::Funcgen::ProbeFeature;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::ExperimentalChip;
use Bio::EnsEMBL::Funcgen::ArrayChip;
use Bio::EnsEMBL::Funcgen::Channel;
use Bio::EnsEMBL::Utils::Exception qw( throw warning deprecate );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use strict;

use vars qw(@ISA);
@ISA = qw(Bio::EnsEMBL::Funcgen::Utils::Helper);

=head2 new

  Example    : my $self = $class->SUPER::new(@_);
  Description: Constructor method for Sanger class
  Returntype : Bio::EnsEMBL::Funcgen::Parsers::Sanger
  Exceptions : throws if Experiment name not defined or if caller is not Importer
  Caller     : Bio::EnsEMBL::Funcgen::Importer
  Status     : at risk

=cut


sub new{
  my $caller = shift;

  my $class = ref($caller) || $caller;
  my $self = $class->SUPER::new();

  throw("This is a skeleton class for Bio::EnsEMBL::Importer, should not be used directly") if(! $self->isa("Bio::EnsEMBL::Funcgen::Importer"));
	
  $self->{'config'} =   {(
						#order of these data arrays is important!
						array_data   => [],	#["array_chip"],
						probe_data   => ["array_probe"],
						results_data => ["and_import_result"],
						#import_methods  => [],
						#data paths here?
						norm_method => undef,
						#is this disabling -input_dir override option?
					   )};



  return $self;
} 

=head2 set_config

  Example    : $imp->set_config();
  Description: Sets a attribute dependent variables
  Returntype : none
  Exceptions : None
  Caller     : Importer
  Status     : At risk 

=cut

sub set_config{
  my ($self) = @_;

  #placeholder method for setting any attr dependant vars e.g. file paths etc.


  return;
}


sub read_array_probe_data{
  my ($self, $array_file) = @_;

  warn("Remove hard coding for Sanger array import, and accomodate adf format");


  $array_file ||= $self->array_file();
  my ($line, $fh, @list, $array_file_format, $cmd);
  my ($op, $of, $imported, $fimported, $fanal);
  my $oa_adaptor = $self->db->get_ArrayAdaptor();
  my $op_adaptor = $self->db->get_ProbeAdaptor();
  my $of_adaptor = $self->db->get_ProbeFeatureAdaptor();
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $ac_adaptor = $self->db->get_ArrayChipAdaptor();
  my $slice_adaptor = $self->db->get_SliceAdaptor();
  my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("SangerPCR")->dbID();
  #have LiftOver? Could then use liftover in  pipeline to redo mappings

  #store now checks whether already stored and updates array chips accordingly
  my $array = Bio::EnsEMBL::Funcgen::Array->new
    (
     -NAME        => $self->array_name(),
     -FORMAT      => uc($self->format()),
     -VENDOR      => uc($self->vendor()),
     -TYPE        => 'PCR',
     -DESCRIPTION => "Sanger ENCODE PCR array 3.1.1",
    );

  ($array) = @{$oa_adaptor->store($array)};
         
  #This is treating each array chip as a separate array, unless arrayset is defined
  #AT present we have no way of differentiating between different array_chips on same array???!!!
  #Need to add functionality afterwards to collate array_chips into single array
  my $array_chip = Bio::EnsEMBL::Funcgen::ArrayChip->new(
														 -NAME      => $array->name(),
														 -DESIGN_ID => $array->name(),
														 -ARRAY_ID  =>$array->dbID(),
														);

  ($array_chip) = @{$ac_adaptor->store($array_chip)};
  $array->add_ArrayChip($array_chip);
  $self->add_Array($array);


  #we also need to test wether the array as been imported as well as the mappings
  #THis needs to use coord_sys-id not schema_build!!  Duplcaite entries for different schema_builds 
  #with same assembly

  my $dnadb_cs = $self->db->dnadb->get_CoordSystemAdaptor->fetch_by_name('chromosome');
  my $fg_cs = $self->db->get_FGCoordSystemAdaptor->validate_and_store_coord_system($dnadb_cs);


  #This fails if we're pointing to an old DB during the release cycle.  Will be fine if we manage to cs mapping dynamically


  if ($array_chip->has_status('IMPORTED')) {
    $imported = 1;
    $self->log("Skipping ArrayChip probe import (".$array_chip->name().") already fully imported");

	#need to build cache here, from file first else from DB????
	#This is required for feature only imports
	#as we won't have the probe dbID available

	if(! $self->get_probe_cache_by_Array($array)){
	  $self->get_probe_cache_by_Array($array, 1);
	}



  } elsif ($self->recovery()) {
    $self->log("Rolling back partially imported ArrayChip:\t".$array_chip->name());
    $self->db->rollback_ArrayChip([$array_chip]);	#This should really remove all CS imports too?
  }


  #should never really have CS imports if not IMPORTED
  #there is however the potential to trash a lot of data if we were to remove the CS importes by mistake
  #do we need to check whether any other sets are using the data?
  #we have to check for result using relevant cs_id and cc_id
  #no removal of probes is the key thing here as nothing is dependent on the feature_ids
  #get all result sets by array chip?  or get all ExperimentalChips by array chip
  #would have to be result set as we would find our own ecs.  May find our own rset
  
  
  throw('This needs updating');

  if ($array_chip->has_status('IMPORTED_CS_'.$fg_cs->dbID())) {
    $fimported = 1;
    $self->log("Skipping ArrayChip feature import (".$array_chip->name().") already fully imported for ".$self->data_version());
  } elsif ($self->recovery()) {
    $self->log("Rolling back partially imported ArrayChip features:\t".$array_chip->name());
    $self->db->rollback_ArrayChip_features($array_chip, $fg_cs);
  }


  #need to check whether already imported on specified schema_build
  #check for appropriate file given format in input dir or take path

  #if (! $fimported) {#now need to do this irrespective of import status due to x y requirements
  #need only do this once, i.e. if the cache isn't defined yet
  #this is assuming cache will be built properly
  #may cause problems if not cleaned up properly after use.

  #ignore xy requirements for now, these should be associated with results file



  #if (! defined $self->{'_probe_cache'}) {
  if (! $fimported) {

	

	if (! $array_file) {

	  if (! defined $self->get_dir('input')) {
		throw("No input_dir defined, if you are running in a non Experiment context please use -array_file");
      }
      
      #hacky ..do better?
      for my $suffix ("gff", "adf") {
		$cmd = $self->get_dir('input')."/".$self->array_name()."*".$suffix;
		@list = `ls $cmd 2>/dev/null`;
	
		if ((scalar(@list) == 1) && 
			($list[0] !~ /No such file or directory/o)) { ###this is only printed to STDERR?
	  
		  if (! defined $array_file) {
			$array_file = $list[0];
		  } else {
			throw("Found more than one array file : $array_file\t$list[0]\nSpecify one with -array_file");
		  }
		}
      }

      throw("Cannot find array file. Specify one with -array_file") if (! defined $array_file);
    }
    
    
    if ($array_file =~ /gff/io) {
      $array_file_format = "GFF";
    } elsif ($array_file =~ /adf/io) {
      $array_file_format = "ADF";
      throw("Does not yet accomodate Sanger adf format");
    } else {
      throw("Could not determine array file format: $array_file");
	}
	

	#if (! $fimported) {
	$fanal = $self->db->get_AnalysisAdaptor->fetch_by_logic_name(($array_file_format eq "ADF") ? "VendorMap" : "LiftOver");
	#}
   
	$self->log("Parsing ".$self->vendor()." array data (".localtime().")");
	$fh = open_file($array_file);
	my @lines = <$fh>;
	close($fh);

	

	my ($chr, $start, $end, $strand, $pid);#, $x, $y, $meta_x, $meta_y, @xy);
    
	#avoid mutliple calls for same array
	my $ac_dbid = $array->get_ArrayChip_by_design_id($array->name())->dbID();

	#sort file to enable probe cache method for new feature imports
	@lines = sort {(split/\t|\;/o, $a)[8] cmp (split/\t|\;/o, $b)[8]} @lines;

	#This is not sorting properly!!

	#my @tmp = map ((split/\t|\;/o, $_)[8], @lines);
	#@tmp = sort @tmp;


	#$self->log('Tmp sorted array is :\n'.join("\n", @tmp)."\n");




	foreach $line(@lines) {
	  $line =~ s/\r*\n//;
      
	  #($chr, $start, $end, $ratio, $pid) = split/\t/o, $line;
	  #($chr, undef, undef, $start, $end, undef, $strand, undef, $pid, $x, $y, $meta_x, $meta_y) = split/\t|\;/o, $line;
	  ($chr, undef, undef, $start, $end, undef, $strand, undef, $pid) = split/\t|\;/o, $line;
	  

	  if($self->ucsc_coords){
		$start += 1;
	  }


	  #$meta_x =~ s/META_X=//;
	  #$x =~ s/X=//;
	  #$x = $x + (($meta_x -1)*26);
	  #$meta_y =~ s/META_Y=//;
	  #$y =~ s/Y=//;
	  #$y = $y + (($meta_y -1)*25);
	  $pid =~ s/reporter_id=//o;
	  $chr  =~ s/chr//;
	  $strand = ($strand eq "+") ? 0 : 1;
	
	  #Hack!!!!!!  This is still maintaining the probe entry (and result?)
	  if (!  $self->cache_slice($chr)) {
		warn("-- Skipping non standard probe (${pid}) with location:\t${chr}:${start}-${end}\n");
		next;
	  }


	  #need to parse dependant on file format 
	  #also need to account for duplicate probes on grid

	  #need to test for imprted here for rebuilding the probe_info cache
	  #this will result in always using first x y for the inital import (i.e. skip any probe already in cache)
	  #or using last x y for previosuly imported as we can't check the cache as it will already be there
	  #could check for x y
	  #should always check x y as this will also implicitly check if it is in the cache
	
	  #if (! $self->get_probe_id_by_name($pid)) { #already present in cache
	  #if(! (@xy = @{$self->get_probe_x_y_by_name($pid)})){
	  
		#can we not use store_set_probes_features
		#would have to add x y to probe, which is not logical as probe can have many x y's
		#keep like this and just change cache_probe_info
	  
	  if (! $imported) {
		#when we utilise array coords, we need to look up probe cache and store again with new coords
		#we're currently storing duplicates i.e. different ids with for same probe
		#when we should be storing two records for the same probe/id
		#the criteria for this will be different for each vendor, may have to check container etc for NimbleGen
		
		$op = Bio::EnsEMBL::Funcgen::Probe->new(
												-NAME          => $pid,
												-LENGTH        => ($end - $start),
												-ARRAY         => $array,
												-ARRAY_CHIP_ID => $ac_dbid,
												-CLASS         => 'EXPERIMENTAL',
											   );
		
		($op) = @{$op_adaptor->store($op)};
		#$self->cache_probe_info($pid, $op->dbID, $x, $y);
	  } else {
		#update XY cache for previously imported array
		#$self->cache_probe_info($pid, $self->get_probe_id_by_name($pid), $x, $y);
	  }
	  
	  #if (! $fimported) {
	  $of = Bio::EnsEMBL::Funcgen::ProbeFeature->new(
													 -START         => $start,
													 -END           => $end,
													 -STRAND        => $strand,
													 -SLICE         => $self->cache_slice($chr),
													 -ANALYSIS      => $fanal,
													 -MISMATCHCOUNT => 0,
													 -PROBE_ID      => ($imported) ? 
													 $self->get_probe_id_by_name_Array($pid, $array) : $op->dbID(),
													);
	  
	  #get_probe_id will throw if not in cache, which means that we have an unimported probe 
	  #for an ArrayChip which is flagged as imported, must have been omitted from the import deisgn
	  #probably a manual fix required.  Can we log these and write an update/repair script.

	  $of_adaptor->store($of);
	  #}
	
	  #} else {
	  #warn("Sanger does not accomodate on plate duplicates yet, result are not linked to X Y coords, using first coords for probe if present in results for $pid\n");　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　　
	  #}
	}

	$array_chip->adaptor->set_status('IMPORTED_CS_'.$fg_cs->dbID(), $array_chip) if ! $fimported;
	$self->log("ArrayChip:\t".$array_chip->design_id()." has been IMPORTED_CS_".$fg_cs->dbID());
	
  }
  
  
  
  if (! $imported) {
    $array_chip->adaptor->set_status('IMPORTED', $array_chip);
    $self->log("ArrayChip:\t".$array_chip->design_id()." has been IMPORTED");
	$self->resolve_probe_data();
  }
  
  $self->log("Finished parsing ".$self->vendor()." array/probe data (".localtime().")");
  #warn("Finished parsing ".$self->vendor()." array/probe data (".localtime().")");  
  
  return;
}

=head2 read_and_import_result_data

  Example    : $imp->read_and_import_result_data();
  Description: Parses and imports result for the sanger PCR array platform
  Returntype : none
  Exceptions : none
  Caller     : Importer
  Status     : At risk

=cut

sub read_and_import_result_data{
  my $self = shift;

  #change this to read_gff_chip_results
  #as opposed to gff channel results
  #This should also use the default logic names for the Vendor, or take a user defined list 
  $self->log("Reading ".$self->vendor()." result data (".localtime().")");

  my ($file, $chip_uid, $line, $echip);
  my ($ratio, $pid, %chip_files, %roll_back);
  my $of_adaptor = $self->db->get_ProbeFeatureAdaptor();
  my $ec_adaptor = $self->db->get_ExperimentalChipAdaptor();
  my $chan_adaptor = $self->db->get_ChannelAdaptor();
  my $analysis = $self->db->get_AnalysisAdaptor->fetch_by_logic_name("SangerPCR");
  my $result_adaptor = $self->db->get_ResultSetAdaptor();
  #this is done to avoid having to self->array_name in loop, will make multiple array loop easier 
  my $array = ${$self->arrays()}[0];

  #This works a little differently as we're not parsing a meta file
  #so the echips haven't been added yet.
  #This is treating each array chip as a separate array, unless arrayset is defined
  #AT present we have no way of differentiating between different array_chips on same array???!!!
  #Need to add functionality afterwards to collate array_chips into single array

  #First add the echips to the Experiment
  
  if (! @{$self->result_files()}) {
    my $list = "ls ".$self->input_dir().'/[0-9]*-[0-9a-zA-Z]*\.all\.*';
    my @rfiles = `$list`;
    $self->result_files(\@rfiles);
  }

  
  foreach $file(@{$self->result_files()}) {
    chomp $file;
    ($chip_uid = $file) =~ s/.*\///;
    $chip_uid =~ s/\..*//;

	$self->log("Found SANGER results file for $chip_uid:\t$file");
    $chip_files{$chip_uid} = $file;
    

    $echip = $ec_adaptor->fetch_by_unique_id_vendor($chip_uid, 'SANGER');

    #this should throw if not recovery
    #Nee to check Nimbelgen methods

    if ($echip) {
  
      if (! $self->recovery()) {
		throw("ExperimentalChip(".$echip->unqiue_id().") already exists in the database\nMaybe you want to recover?");
      }else{
		#log pre-reg'd chips for rollback
		$roll_back{$echip->dbID()} = 1;
	  }
    } else {

      $echip =  Bio::EnsEMBL::Funcgen::ExperimentalChip->new
		(
		 -EXPERIMENT_ID  => $self->experiment->dbID(),
		 -ARRAY_CHIP_ID  => $self->arrays->[0]->get_ArrayChip_by_design_id($array->name())->dbID(),
		 -UNIQUE_ID      => $chip_uid,
		);
          
      ($echip) = @{$ec_adaptor->store($echip)};	
      $self->experiment->add_ExperimentalChip($echip); #if we need a contains method in  here , always add!!
    }

	#do we need DUMMY entries any more?

    #sub this passing the echip?
    foreach my $type ('DUMMY_TOTAL', 'DUMMY_EXPERIMENTAL') {

      my $channel = $chan_adaptor->fetch_by_type_experimental_chip_id($type, $echip->dbID());
      
      if ($channel) {
		if (! $self->recovery()) {
		  throw("Channel(".$echip->unique_id().":$type) already exists in the database\nMaybe you want to recover?");
		}
      } else {

		$channel =  Bio::EnsEMBL::Funcgen::Channel->new
		  (
		   -EXPERIMENTAL_CHIP_ID => $echip->dbID(),
		   -TYPE                 => $type,
		  );
	
		($channel) = @{$chan_adaptor->store($channel)};
      }
    }
  }


  
  #Now get rset using experiment echips
  my $rset = $self->get_import_ResultSet($analysis, 'experimental_chip');

  if ($rset) {					#we have some new data

    foreach my $echip (@{$self->experiment->get_ExperimentalChips()}) {
      
      if ($echip->has_status('IMPORTED_SangerPCR', $echip)) {
		$self->log("ExperimentalChip(".$echip->unique_id().") has already been imported");
      } else {

		my $cc_id = $rset->get_chip_channel_id($echip->dbID());
		
		if ($self->recovery() && $roll_back{$echip->dbID()}){
		  $self->log("Rolling back results for ExperimentalChip:\t".$echip->unique_id());
		  $self->rollback_results($cc_id);
		}

		$self->log("Reading SANGER result file for ".$echip->unique_id().":\t".$chip_files{$echip->unique_id()});
		$self->get_probe_cache_by_Array($array) || throw('Failed to reset probe cache handle');
		my $fh = open_file($chip_files{$echip->unique_id()});
		my @lines = <$fh>;
		close($fh);

		my $rfile_path = $self->get_dir("norm")."/result.SangerPCR.".$echip->unique_id().".txt";
		my $rfile = open_file($rfile_path, '>');
		my $r_string = "";
		
		
		@lines = sort {(split/\t|\:/o, $a)[5] cmp (split/\t|\:/o, $b)[5]} @lines;

		foreach my $line (@lines) {
		  $line =~ s/\r*\n//o;
		  
		  ($ratio, undef, $pid) = (split/\t|\:/o, $line)[3..5];
		  $pid =~ s/.*://o;
		  
		  $ratio = '\N' if $ratio eq 'NA'; #NULL is still useful info to store in result
		  #my ($x, $y) = @{$self->get_probe_x_y_by_name($pid)};
		  
		  #this is throwing away the encode region which could be used for the probeset/family?	
		  $r_string .= '\N'."\t".$self->get_probe_id_by_name_Array($pid, $array)."\t${ratio}\t${cc_id}\t".'\N'."\t".'\N'."\n";#${x}\t${y}\n";
		}
		
		print $rfile $r_string;
		close($rfile);

		$self->log("Importing:\t$rfile_path");
		$self->db->load_table_data("result",  $rfile_path);
		$self->log("Finished importing:\t$rfile_path");
		$echip->adaptor->set_status('IMPORTED_SangerPCR', $echip);
	  }
    }



  } else {
    $self->log("No new data, skipping result parse");
  }

  $self->log("Finished reading and importing ".$self->vendor()." result data (".localtime().")");
  return;
}



1;
