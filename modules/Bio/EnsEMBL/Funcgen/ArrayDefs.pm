package Bio::EnsEMBL::Funcgen::ArrayDefs;


use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use strict;


#define Array defs hashes here
#vendor specific paths etc
#method arrays for parse


sub set_defs{
	my ($self) = @_;

	#Is this only valid unless somethign inherits from from Experiment ?
	$self->throw("This is a skeleton class for Bio::EnsEMBL::Experiment, should not be used directly") if(! $self->isa("Bio::EnsEMBL::Funcgen::Experiment"));



	#Can we do some funky check on $self to see if it is Experiment
	#Abstract these to ${vendor}Defs.pm modules
	#we can't use methods to build vars :(
	
	#could tidy up below but would disable over-ride for inpur_dir option
	#$self->{'_input_dir'} = $self->get_dir("data")."/native/".$self->vendor()."/".$self->instance();

	my %array_defs = (
					  nimblegen => {
									#order of these data arrays is important!
									array_data   => ["array_chip"],
									probe_data   => ["probe"],
									results_data => ["results"],
									#import_methods  => [],
									#data paths here?
									
									#This input dir is being recast as _input_dit
									input_dir      => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->instance(),
									#but this is vendor specific and remains an array_def
									design_dir     => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->instance()."/DesignFiles",

									chip_file        => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->instance()."/SampleKey.txt",
									array_file       => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->instance()."/DesignNotes.txt",
									#feature_map_file => $self->get_dir("import")."/feature_map.tmp",
									results_file     => $self->get_dir("data")."/raw/".$self->vendor()."/".$self->instance()."/PairData/All_Pair.txt",

									dye_freqs => {(
												   Cy5 => 635,
												   Cy3 => 532,
												  )},

									#PROBE_DESIGN_ID	CONTAINER	DESIGN_NOTE	SELECTION_CRITERIA	SEQ_ID	PROBE_SEQUENCE	MISMATCH	MATCH_INDEX	FEATURE_ID	ROW_NUM	COL_NUM	PROBE_CLASS	PROBE_ID	POSITION	DESIGN_ID	X	Y
									#header_fields => [


									#this won't work as design_name is undef at this point
									#need to make this a ref so it will update the values?
									#gene_description_file => $self->get_array_def("design_dir")."/".$self->array_data("design_name").".ngd",
								   },
					 );

	#print "Setting ".$self->vendor()." defs to ".$array_defs{$self->vendor()}." ".Data::Dumper::Dumper($array_defs{$self->vendor()})."\n";

	$self->{'_array_defs'} = $array_defs{$self->vendor()};
	
	#Set mandatory defs here?



	return;

}



sub get_def{
	my ($self, $data_name) = @_;
	return $self->get_data('array_defs', $data_name);#will this cause undefs?
}





#How are we going to use different formats here?
#This is hardcoded for nimblegen at present
#Nimblgen data does not contain array data as such!

sub read_array_chip_data{
	my ($self) = shift;

	my ($design_desc);
	my $fh = open_file("<", $self->get_def("array_file"));
	($design_desc = do { local ($/); <$fh>;}) =~ s/\r*\n$//;
	close($fh);
	  
	#Currently 1 design = 1 chip = 1 array /DVD
	#Different designs are not currently collated into a chip_set/array in any ordered manner
	#Register each design as an array and an array_chip
	#May want to group array_chips into array/chip sets by association though the API

	$fh = open_file("<", $self->get_def("chip_file"));

	$self->log("Reading chip data");

	my ($line, $chip_uid, $dye, $design_name, $design_id, $sample_label, $species, $sample_type);
	my $tmp_uid = "FIRST";
	#Need some way of capturing other experimental variables?

	while ($line = <$fh>){
		next if $. == 1;#or should we have format check here for files

		$line =~ s/\r*\n//;#export from Helper::chump? changing $/ makes input go awry for some reason

		#Need to handle array class here i.e. two channel arrays will have two lines
		#ignore tissue treatment, too wordy with typos.
		#Need to map sample_label to dye and description?
		#Some duplicate sample labels do not have the same description > different channels
		
		(undef, $chip_uid, $dye, $design_name, $design_id, $sample_label, 
		 $species, undef, undef, $sample_type) = split/\t/o, $line;
		#Need to standardise species here, nimblegen has different names and typos
		
		if(! defined $self->{'_achips'}){
		
			#Need to capture array species here too! User define? Throw only if not already registered

			$self->{'_achips'}{$design_id} = {(
											   dbid        => undef,
											   design_name => $design_name,
											   description => $design_desc,
											  )};

			#Hack to enable generic array registration
			$self->{'_array'} = {(
								  dbid        => undef,
								  #design_id   => $design_id,#don't need this?
								  design_name => $design_name,
								  description => $design_desc,
								 )};

		}elsif(! exists $self->{'_achips'}{$design_id}){
			$self->throw("Found experiment with more than one design");
		}
	

		if ($tmp_uid eq "FIRST" || $chip_uid != $tmp_uid){
			$tmp_uid = $chip_uid;


			#These are experimental_chips
			$self->{'_echips'}{$chip_uid} = {(
											  design_id => $design_id,
											  dbid => undef,#This will be populated after import for use with results
											  channels  => {},
											)};
		}

		
		#Handles single/mutli

		#Only doing this method once at mo, but may need for other formats
		$self->set_channel($chip_uid, (
										dbid => undef,
										dye => $dye,
										sample_label => $sample_label,
										species => $species,#on channel/sample to enable multi-species chip/experiment
										type => uc($sample_type),
									   ));#do we need to inset {}?




	}

	close($fh);
	return;
}



sub array_data{
	my ($self, $data_name) = @_;
	return $self->get_data('array', $data_name);#will this cause undefs?
}

sub achip_data{
	my ($self, $design_id, $data_type, $value) = @_;

	$self->throw("Need to specify a design_id  and a data_type") if (! defined $data_type || ! defined $design_id);

	if(defined $value){
		${$self->get_data('achips', $design_id)}{$data_type} = $value;
	}
	else{
		return ${$self->get_data('achips', $design_id)}{$data_type};#will this cause undefs?
	}
}

sub echip_data{
	my ($self, $chip_uid, $data_type, $value) = @_;

	$self->throw("Need to specify a design_id  and a data_type") if (! defined $data_type || ! defined $chip_uid);
	
	if(defined $value){
		${$self->get_data('echips', $chip_uid)}{$data_type} = $value;
	}
	else{
		return ${$self->get_data('echips', $chip_uid)}{$data_type};#will this cause undefs?
	}
}

#have generic get_a/e/chips methods?

#sub get_chip{
#	my ($self, $type, $chip_uid) = @_;

#	#add s here
#	return $self->get_data("${type}s", $chip_uid);
#}



sub get_channels{
	my ($self, $chip_uid) = @_;
		
	return $self->echip_data($chip_uid, 'channels');
}

sub get_channel{
	my ($self, $chan_uid) = @_;
	
	#This may be different for other vendors!!		
	my ($chip_uid);
	($chip_uid = $chan_uid) =~ s/_[0-9]*//o;

	return ${$self->echip_data($chip_uid, 'channels')}{$chan_uid};
}


sub set_channel{
	my ($self, $chip_uid, %channel) = @_;

	$self->throw("Need to pass a chip_uid and a channel hash") if(! defined $chip_uid || ! %channel);

	#my $tmp = "${chip_uid}_".${$self->get_def("dye_freqs")}{$channel{'dye'}};


	$self->{'_echips'}{$chip_uid}{'channels'}{"${chip_uid}_".${$self->get_def("dye_freqs")}{$channel{'dye'}}} = {%channel};
	#$self->echip_data($chip_uid, 'channels'){$tmp} = {%channel};
	return;
}

sub channel_data{
	my ($self, $chan_uid, $data_type, $value) = @_;

	$self->throw("Need to provide a channel uid and a data_type") if (! defined $chan_uid || ! defined $data_type);

	if(defined $value){
		${$self->get_channel($chan_uid)}{$data_type} = $value;
	}
	else{
		return ${$self->get_channel($chan_uid)}{$data_type};
	}
}



#this handle probe_set, probe, probe_feature
#Assumes no chip_design per experimental set.
sub read_probe_data{
	my ($self) = shift;

	$self->log("Parsing probe data (".localtime().")...");

	### Read in
	#ndf file: probe_set, probe and probe_feature(.err contains multiple mappings)
	#pos file: probe locations(counts)
	#ndf file: build, chr, start stop,

	#Need to test whether array exists and is validated before importing probe data.
	#If already present, then can skip and just import results using probe_name to get probe_id and link to probe_feature.
	#just needs to create a mapping between probe_feature_ids and results
	#What if an array has been imported with only a few chips present?
	#Need check on size in import_array
	#e.g. someone has done an experiment with only a few chips from an array?
	#This would mean absent probes

	#also need to change how probe_names are generated for nimblegen?
	#native probe_ids may not be unique, but should be when combined with the seq_id which is currently being used as the xref_id
	#Handle with API!!

	#READ REGION POSITIONS
	my $fh = open_file("<", $self->get_def("design_dir")."/".$self->array_data("design_name").".ngd");

	my ($line, $start, $stop, @data, %hpos, %regions, %probe_pos);
	
	while ($line = <$fh>){
		$line =~ s/\r*\n//;
		@data =  split/\||\t/o, $line;


		#SEQ_ID	SEQ_UNIQUE|BUILD|CHROMOSOME|LOCATION|DESCRIPTION|DATE_ENTERED|SOURCE_DB
		if ($. == 1){
			%hpos = %{$self->set_header_hash(\@data)};
			next;
		}

		#What about strand!!!!!!!!!!!
		$data[$hpos{'CHROMOSOME'}] =~ s/chr//;
		($start, $stop) = split/-/o, $data[$hpos{'LOCATION'}];
		
		#Do we need seq_id check here for validity?
		#overkill?
		if(exists $regions{$data[$hpos{'SEQ_ID'}]}){
			croak("Duplicate regions\n");
		}else{
			#$data[$hpos{'CHROMOSOME'}] = species_chr_num($self->species(), 	$data[$hpos{'CHROMOSOME'}]);

			#Set region hash for SEQ_ID
			#Need to look up seq_region id here for given build
			#Build should be manually specified as we can't guarantee it will be in the correct format
			#or present at all

			$regions{$data[$hpos{'SEQ_ID'}]} = {
												start => $start,
												stop  => $stop,
												seq_region_id => $self->get_chr_seq_region_id($data[$hpos{'CHROMOSOME'}], $start, $stop),
												build => $data[$hpos{'BUILD'}],
											   };
		}
		
	}
	
	close($fh);
	

	#ONLY USE THIS FOR VALIDATION OF EXPERIMENTAL PROBES!!!! CAN REMOVE IF PROBE_CLASS POPULATED
	#SLURP PROBE POSITIONS
	$fh = open_file("<", $self->get_def("design_dir")."/".$self->array_data("design_name").".pos");

	#don't % = map ! Takes a lot longer than a while?!
	while($line = <$fh>){
		#$line =~ s/\r*\n//;#Not using last element
		@data =  split/\t/o, $line;

		#SEQ_ID	CHROMOSOME	PROBE_ID	POSITION	COUNT
		if ($. == 1){
			%hpos = %{$self->set_header_hash(\@data)};
			next;
		}
		#($seq_id, undef, $probe_id, $lstart) = split/\t/o, $line;

		#can we remove this?
		$self->throw("Found duplicate mapping for ".$data[$hpos{'PROBE_ID'}]) if(exists $probe_pos{$data[$hpos{'PROBE_ID'}]});
		
		$probe_pos{$data[$hpos{'PROBE_ID'}]} = {(
												 seq_id => $data[$hpos{'SEQ_ID'}],
												 lstart => $data[$hpos{'POSITION'}],
												)};
		
	}

	close($fh);
	

	#OPEN PROBE IN/OUT FILES
	$fh = open_file("<", $self->get_def("design_dir")."/".$self->array_data("design_name").".ndf");

	#Need to set these paths in each  achip hash, file names could be tablename.chip_id.txt
	my $p_out = open_file(">", $self->get_dir("import")."/probe.txt");
	my $ps_out = open_file(">", $self->get_dir("import")."/probe_set.txt");
	my $pf_out = open_file(">", $self->get_dir("import")."/probe_feature.txt");
	my $f_out = open_file(">", $self->get_dir("import")."/probe.fasta")	if($self->{'_dump_fasta'});
   	my ($loc, $probe_string, $ps_string, $pf_string, $f_string, $class, $length, %probe_set, %probe_map);
	my $size = 1;#needed for first set

	#If this chip has been pre-registered then get the the very first probe_feature_id and build the featuremap from that.
	#May change as changing link to probe rather than probe_feature table

	my $psid = $self->db_adaptor->get_last_table_id("probe_set");
	my $pid = $self->db_adaptor->get_last_table_id("probe");
	my $fid = $self->db_adaptor->get_last_table_id("probe_feature");
	
   	##HARDCODED!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
	my $strand = 0;
	my $anal_id = 1;#testing how are we going to define analysis methods?
	my $cig_line = "50M";
	
  	$probe_string = $ps_string = $pf_string = $f_string = "";
	#$self->Timer()->mark("Starting probe loop");
	
	while($line = <$fh>){
		#$line =~ s/\r*\n//;#don't need this as we're not using the last element?
		@data =  split/\t/o, $line;

		#PROBE_DESIGN_ID	CONTAINER	DESIGN_NOTE	SELECTION_CRITERIA	SEQ_ID	PROBE_SEQUENCE	MISMATCH	MATCH_INDEX	FEATURE_ID	ROW_NUM	COL_NUM	PROBE_CLASS	PROBE_ID	POSITION	DESIGN_ID	X	Y
		if ($. == 1){	
			%hpos = %{$self->set_header_hash(\@data)};
			next;
		}

		#These header fields should be stored in the nimblgen defs
		#overkill for now as this is nimblgen specific
	   
		$loc = "";
		$pid ++;
		
		###PROBE SETS
		#Shall we use x and y here to creete a temp cel file for checking in R?
		#	(undef, $family, undef, undef, $seq_id, $seq, $mismatch, $match_index, 
		#	 $probeset_id, $y, $x, undef, $probe_id, $lstart) = split/\t/o, $line;
		
		#THis probe_set business is a little backward as we need to count size before printing
		if($. == 2){#Capture first set
			$psid++;
			%probe_set = (
						  name    => $data[$hpos{'FEATURE_ID'}],
						  size    => $size,
						  family  => $data[$hpos{'CONTAINER'}],
						  xref_id => $data[$hpos{'SEQ_ID'}],
						  #design_id?
						 );
		}
	


	
		if($data[$hpos{'FEATURE_ID'}] ne $probe_set{'name'}){
			$psid++;

			#print "adding to ps string\n"."\t".$probe_set{'name'}."\t".$probe_set{'size'}."\t".$self->get_array_dbid()."\t".
			#  $probe_set{'family'}."\t".$probe_set{'xref_id'}."\n";
			

			$ps_string .= "\t".$probe_set{'name'}."\t".$probe_set{'size'}."\t".$self->get_array_dbid()."\t".
			  $probe_set{'family'}."\t".$probe_set{'xref_id'}."\n";

			%probe_set = (
						  name    => $data[$hpos{'FEATURE_ID'}],
						  size    => $size,
						  family  => $data[$hpos{'CONTAINER'}],
						  xref_id => $data[$hpos{'SEQ_ID'}],
						  #design_id?
						 );
			
			
			$size = 1;
		}else{
			$size ++;
		}
		
		
		###PROBES

		####THIS IS WRONG??  Actually looking at CONTAINER, not PROBE_CLASS?
		#Also does not catch non-experimental containers e.g. H_CODE?
		
		$class = "EXPERIMENTAL";
		
		if($data[$hpos{'CONTAINER'}] =~ /control/i){
			$class = "CONTROL";
		}
		elsif($data[$hpos{'CONTAINER'}] =~ /random/i){
			$class = "RANDOM";
		}
		elsif(! exists $probe_pos{$data[$hpos{'PROBE_ID'}]}){	#HACKY HACKY HACK HACK!! Needed for valid region retrival
			$class = "OTHER";
		}
	
		
		#should we cat $xref_id to $probe_id here to generate unique id?
		#would be messy to handle in the code, but would have to somewhere(in the retrieval code)
		
		$probe_string .= "\t${psid}\t".$data[$hpos{'PROBE_ID'}]."\t".$data[$hpos{'MATCH_INDEX'}]."\t${class}\n";
		push @{$self->{'_probe_map'}{$data[$hpos{'PROBE_ID'}]}}, $pid;
		
		
		###PROBE FEATURES
		#How can we be certain that we have the same mapping in the DB?
		#Put checks in here for build?
	
		
		#Need to handle controls/randoms here
		#won't have features but will have results!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


		#The format of the pos file looks like it should have all the data required, but 
		# chromsome is missing, first undef :(

		#Need to use $pos here instead of .pos file
		#However have problems defining probe class, as not populated in test set
		#derive from container! :(
		#Ignore controls/random as they won't have a region

		#Also need to handle multiple mappings? 
		#As results reference probes, no features, can have multiple features on different builds
	
		
		#if(exists $probe_pos{$probe_id}){
		if($class eq "EXPERIMENTAL"){
			$length = length($data[$hpos{'PROBE_SEQUENCE'}]);
			
			if(exists $regions{$data[$hpos{'SEQ_ID'}]}){
				$fid++;
				$pf_string .= "\t".$regions{$data[$hpos{'SEQ_ID'}]}{'seq_region_id'}."\t".$data[$hpos{'POSITION'}]."\t".
				  ($data[$hpos{'POSITION'}] + $length)."\t${strand}\t".$data[$hpos{'MISMATCH'}]."\t${pid}\t".
					$regions{$data[$hpos{'SEQ_ID'}]}{'build'}."\t${anal_id}\t${cig_line}\n";
				
				$loc .= $regions{$data[$hpos{'SEQ_ID'}]}{'seq_region_id'}.":".$data[$hpos{'POSITION'}].
				  "-".($data[$hpos{'POSITION'}] + $length).";" if ($self->{'_dump_fasta'});
			}
			else{ 
				die("No regions defined for ".$data[$hpos{'SEQ_ID'}]." ".$data[$hpos{'PROBE_ID'}].
					" with family ".$data[$hpos{'CONTAINER'}]);
			  }
		}
		
		if($self->{'_dump_fasta'}){			
			#filter controls/randoms?  Or would it be sensible to see where they map
			#wrap seq here?
			$f_string .= ">".$data[$hpos{'PROBE_ID'}]."\t".$data[$hpos{'SEQ_ID'}]."\t$loc\n".$data[$hpos{'PROBE_SEQUENCE'}]."\n";
		}
		
		
		#Should we build hash of probe_names:probe_feature_ids here for results import
		#Should we dump this as a lookup file for easier recoverability
		#This is the biggest step in the import
		#Building the hash would be fastest but least recoverable
		#Would have to write recover statments for each step i.e. build the most recent data structure required for the next import step
		#Individual queries for each result would take ages
		#This is all assuming there are no random records in the table i.e. ID series is linear with no gaps.
		
		
		#Could throw a random validation check in every X entries?
		#This would only work for non-parallel imports
		#periodic import when hit new probe_set, but no new data printed
	
	}
	
	#Need to print last probe_set here only if current and last probeset_id match
	if($probe_set{'name'} eq $data[$hpos{'FEATURE_ID'}]){
		$ps_string .= "\t".$probe_set{'name'}."\t".$probe_set{'size'}."\t".$self->get_array_dbid()."\t".
		  $probe_set{'family'}."\t".$probe_set{'xref_id'}."\n";
	}


	print $p_out $probe_string;
	print $ps_out $ps_string;
	print $pf_out $pf_string;
	print $f_out $f_string if($self->{'_dump_fasta'});
	#$self->Timer()->mark("End of probe loop");

	close($fh);
	close($ps_out);
	close($p_out);
	close($pf_out);
	close($f_out) if ($self->{'_dump_fasta'});

	$self->log("Finished parsing probe data\nTotal probe_sets:\t$psid\n".
			   "Total probes:\t$pid\nTotal probe_features:\t$fid");

	return;
}




sub get_probe_map{
	my ($self) = shift;

	my($fh, @tmp);

	if(! defined $self->{'_probe_map'}){
		
		#Do we really need this now, as we're mapping directly to probes
		if(-f $self->get_dir("import")."/probe_map.tmp"){
			$fh = open_file("<",  $self->get_dir("import")."/probe_map.tmp");#$self->get_file("feature_map"));
			
			$self->throw('Need to impement hash slurp');
			#%feature_map = map { chomp; ($_, @tmp) = split/\t/, $_; $feature_map{$_} = @tmp } %feature_map = <$fh>;
		}
		else{
			#Get feature_map from DB!!! dump also for recovery?
			


		}
	}

	return $self->{'_probe_map'};#should we be returning ref here?
}

sub read_results_data{
	my ($self) = shift;

	#TO DO
	#Recovery ? read feature_map.tmp into , query DB for last imported result(remove for safety?), restart from that point in file?
	#slurp here may require too much memory?
	#import files and do checks on ids to make sure they've imported properly
	#i.e. select the last entry based on the expected table id.
	
	$self->log("Parsing results(".localtime().")...");
	
	my $fh = open_file("<", $self->get_def("results_file"));

	my ($i, $line, $probe_elem, $first_result, @header, @data);
	my $r_string = "";

	#Hack var!
	my $metric_id = 1;
	
	while($line = <$fh>){
		$line =~ s/\r*\n//;
		
		if ($. == 1){
			#GENE_EXPR_OPTION	SEQ_ID	PROBE_ID	POSITION	43827_532	43827_635	43837_532	43837_635	46411_532	46411_635	46420_532	46420_635	47505_532	47505_635	47525_532	47525_635
			@header = split/\t/o, $line;

			for $i(0..$#header){
				$probe_elem = $i if($header[$i] eq "PROBE_ID");
				
				if($header[$i] =~ /[0-9]+_[0-9]+/){
					$first_result = $i;
					last;#assume all fields after first result, are result feilds
				}
			}
			next;
		}

		@data = split/\t/o, $line;
		
		#could validate here against reulsts vs. number of channels
		for $i($first_result..$#data){
			
			#multiple mappings
			foreach my $pid(@{${$self->get_probe_map()}{$data[$probe_elem]}}){
				#removed print to speed up io
				#print $r_out "\t${pid}\t".$data[$i]."\t$metric_id\t".$self->get_channel_dbid($header[$i])."\n";
				$r_string .= "\t${pid}\t".$data[$i]."\t$metric_id\t".$self->get_channel_dbid($header[$i])."\n";
			}
		}
	}
	my $r_out = open_file(">", $self->get_dir("import")."/result.txt");
	print $r_out $r_string;
	close($r_out);
	
	$self->log("Finished parsing and results");
	return;
}



1;
