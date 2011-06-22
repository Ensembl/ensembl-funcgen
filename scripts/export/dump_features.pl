#!/usr/bin/env perl

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

=cut

#This is generic format exporter which utilises the various FeatureAdaptor dump_FORMAT_by_Slice methods
#Or do we want an Exporter module which utilises various format parsers?
#Do we necessarily need to create the objects for this, or can we split _obj_from_sth
#into _data_from_sth and just have a separate wrapper for the dump
#This approach could also be used in dynamic collection generation
#quite often need object for associated object info.
#This may also double flow control loop unless we could stream or iterate this info?
#This would mean caching all the var in the method as attrs

#How an we over-ride calling _obj_from_sth? Basically we want all the underlying code, but the final method call to change 
#dependant on the output required. Could set a global var and use _obj_from_sth as wrapper. This would me one extra test per record?
#Or can we used a method var i.e. $self->$method_name where method_name can be _new_fast or dump_bed/gff or collect_feature

#Advantage here is that we don't double iterate.
#Disadvantage is obfucation _obj_from_sth really needs renaming
#Also method var performance vs hardcoded method name

#wrt to format dumping. Can we have base method for common feature attrs, then FeatureAdaptor specific methods for rest.
#Also probably need wrapper method in object.

#Need to pass hash to wrapped method. This will enable direct access or blessing of relevant attrs. No re-arrange!

#No way of using full namespace and method as var?
#my $data_method = 'Bio::EnsEMBL::Funcgen::ProbeFeature->new_fast';
#$data_method = '&dump_gff';

#Would have to wrap new_fast which would mean an extra method call for the core functionality :/
#Unless we put all the dump method in the Feature class and always use Bio::EnsEMBL::Funcgen::ProbeFeature->$data_method?
#other data methods would have to be subs 

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Utils::Iterator;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file strip_param_args strip_param_flags generate_slices_from_names);
use Cwd;
 
#To do
# 1 Integrate into Exporter module
# 2 Genericise this to dump_features, use various format parsers for output
# 3 enable different set dumps, i.e. incorporate get_data.pl
# 4 Add -contig_to_un mode which cats all contigs to one 'Un' file
# 5 Imlement Iterators for other features?
# 6 Add options to dump sets/arrays into separate files?
# Use logger?

# Bloat the feature class vs separate Exporter?
# which knows how to translate various features to each format. But is this necessary if we want to 
# avoid the extra iteration etc...

#For now have all subs here, in a way that can easily be moved to the relevant adaptors.

#force local (multi slice) doesn't work with -post_process yet
#Add merged dump mode? This already works for feature sets
#can we do a union query across different feature types filling adding the different fields as nulls
#to support a sorted query across all features? Would have to do a sub select and sort order union


$| =1;

my @tmp_args = @ARGV;
print"$0 @tmp_args\n";

my ($pass, @fset_names, @rset_names, $array, $vendor, $dump_name, $queue, $on_farm, $merged_dump);
my ($exp_name, $help, $pids, $man, @features, $chr, $not_status, $ensembl_coords, $farm, $cat_dumps);
my ($dbhost, $port, $user, $dbname, $cdbname, $species, @slices, @skip_slices, $out_root);
my ($dnadb_pass, $dnadb_user, $dnadb_name, $dnadb_port, $dnadb_host, $with_headers, $force_local);
my ($no_clean, $post_process, $window_size, $no_dups, $bin_dir);
my $inc_dups = 1;
my $format = 'GFF';
my $zip = 0;
my $out_dir = ".";
my $keep_colons   = 0;
my $write_file_header  = 1;
my $write_sr_header = 1;




GetOptions 
  (
   'feature_sets=s{,}' => \@fset_names,
   'result_sets=s{,}'  => \@rset_names,
   'window_size=s'     => \$window_size,
   #@rset_ids?
   #Change this to feature_sets to allow merged dumps from same feature class.
   'array=s'          => \$array,#Change this to \@arrays?
   'vendor=s'         => \$vendor,
   'format=s'         => \$format,
   'ensembl_coords=s' => \$ensembl_coords,#some UCSC formats are half open
   #shoudl really remove this and force format spec

   "pass=s"           => \$pass,
   'user=s'           => \$user,
   "port=s"           => \$port,
   'species=s'        => \$species,
   "dbname=s"         => \$dbname,
   "dbhost=s"         => \$dbhost,
   
   "dnadb_pass=s"     => \$dnadb_pass,
   'dnadb_user=s'     => \$dnadb_user,
   "dnadb_port=s"     => \$dnadb_port,
   "dnadb_name=s"     => \$dnadb_name,
   "dnadb_host=s"     => \$dnadb_host,
   
   "out_root=s"       => \$out_root,
   'dump_name=s'      => \$dump_name, #default to feature_set or array name, mandatory for multi set dump
   'slices=s{,}'      => \@slices,
   'skip_slices=s{,}' => \@skip_slices,
 
   'bin_dir=s'        => \$bin_dir,
   'no_dups'          => \$no_dups, #only works when -slices omitted
   'force_local'      => \$force_local,
   'farm'             => \$farm,  
   '_on_farm'         => \$on_farm,
   #'with_headers=i'   => \$with_headers, #This writes the file header when on farm?
   #default is to write on sr headers for slice farm jobs
   'merged_dump'      => \$merged_dump,
   'post_process'     => \$post_process,
   'no_clean'         => \$no_clean,#doen't remove files

   #merged dump i.e. sets or slices in same file
   #i.e. split into set or slice jobs but not both

   'queue=s'          => \$queue,
   'zip|z'            => \$zip,
   #'no_cat_un'        => \$no_cat_un,
   'keep_colons'      => \$keep_colons,
   

   "help|?"           => \$help,
   "man|m"            => \$man,
  ) 
  or pod2usage( -exitval => 1,
				-message => "Params are:\t@tmp_args"
			  );


$inc_dups = 0 if $no_dups;

pod2usage(0) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;




if(! $dbname ){
  throw("You must provide a funcgen -dbname paramter");
}

throw("Must define your funcgen -dbhost") if ! $dbhost;


my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host => $dbhost,
													  -dbname => $dbname,
													  -species => $species,
													  -user => $user,
													  -pass => $pass,
													  -port => $port,
													  -dnadb_pass => $dnadb_pass,
													  -dnadb_user => $dnadb_user,
													  -dnadb_port => $dnadb_port,
													  -dnadb_name => $dnadb_name,
													  -dnadb_host => $dnadb_host,
													 );


#Check DB connections
$db->dbc->db_handle;
$db->dnadb->dbc->db_handle;

#Check out_dir
#Use Helper!

my $slice_a  = $db->get_SliceAdaptor;
my $fset_a   = $db->get_FeatureSetAdaptor;
my $rset_a   = $db->get_ResultSetAdaptor;
my $array_a  = $db->get_ArrayAdaptor;

# VALIDATE input types
my $input_type_cnt = 0;
map {$input_type_cnt++ if $_ } (scalar(@fset_names), $array, scalar(@rset_names));

if( (! $input_type_cnt) || ($input_type_cnt > 1) ){
  die("You have specify dump type params from one of the following:\n".
	  "\t-array and -vendor\n\t-feature_sets\n\t-result_sets\n");
}


#Handle merged regulation evidence dumps
#These should have been specified: -feature_sets AnnotatedFeatures

if(@fset_names && 
   ($fset_names[0] eq 'AnnotatedFeatures')){
  
  #Could do this via meta strings?
  my %fset_names;
  my $dset_a = $db->get_DataSetAdaptor;

  foreach my $rf_fset(@{$fset_a->fetch_all_by_feature_class('regulatory')}){

	if($rf_fset->name =~ /v[0-9]+$/){
	  warn "Skipping archived set:\t".$rf_fset->name.
		"\nThis should be removed before release\n";
	  next;
	}

	map {$fset_names{$_->name} = undef} 
	  @{$dset_a->fetch_by_product_FeatureSet($rf_fset)->get_supporting_sets};
  }

  @fset_names = keys %fset_names;
}
  





### FORMAT CONFIG
#   Move to Exporter along with format methods?

$format      = uc($format); #validation done via introspection on Feature object
my $compress;

#handle compression formats
if($format eq 'BIGWIG'){
  $format = 'WIG';
  $compress = 1;
}
#add BIGBED etc

my %format_config = (
					 #Is no_dups format specific too?
					 

					 BED => {
							 has_header    => 1,
							 slice_headers => 0,
							 #allow_merged_dump => 1,
							 strand        => {
											   0  => '.',
											   1  => '+',
											   -1 => '-',
											  },
							},

					 WIG => {
							 has_header    => 1,
							 slice_headers => 1, #?
							 #allow_merged_dump => 0,
							 compression_script => 'wigToBigWig',

							 #add variable conf in ftype config							 
							},


					 GFF => {
							 has_header    => 0,
							 strand        => {
											   0  => '',
											   1  => '+',
											   -1 => '-',
											  },
							}
					);

my $conf_ref = $format_config{$format};




### FEATURE TYPE CONFIG ###
my ($feature_adaptor, $file_name, $compression_args);
my ($feature_class, $no_y, @fsets, @rsets, @fetch_params, @header_params, @output_names);

if($array){

  if(! ($array && $vendor)){
	die("Please choose either array and vendor or feature_sets parameters\n");
  }
  @output_names = ($array);
  my $array = $array_a->fetch_by_name_vendor($array, $vendor);

  if(! $array){
	die("Could not fetch array:\t$vendor\t$array\n");
  }

  $dump_name     ||= $array->vendor.'_'.$array->name;
  $file_name       = $dump_name;
  $feature_class   = 'ProbeFeature';
  @fetch_params    = ([$array]);
  @header_params   = ($vendor, $array);#arrays?
  $feature_adaptor = $db->get_ProbeFeatureAdaptor;
}
elsif(@fset_names){
 
  #default is merged dump here


  foreach my $fset_name(@fset_names){

	my $fset = $fset_a->fetch_by_name($fset_name);

	$feature_class ||= $fset->feature_class;
	
	if($feature_class ne $fset->feature_class){
	  die("Cannot dump mixed FeatureSet classes:\t".$feature_class."\t".$fset->feature_class);
	}

	push @fsets, $fset;
  }

  if(! $dump_name){

	if(scalar(@fsets) == 1){
	  $dump_name = $fset_names[0];
	}
	else{
	  die('You must specify a -dump_name to dump multiple FeatureSets');
	}
  }

  $file_name       = $dump_name;
  $feature_class    = ucfirst($feature_class).'Feature';
  @fetch_params     = (\@fsets);
  @header_params    = ();#?
  @output_names     = @fset_names;
  $feature_adaptor	= $fsets[0]->get_FeatureAdaptor;
}
else{ #@rset_names
  #can't do merged dump 
  #do this must farm by slice and result set!
  
  foreach my $rset_name(@rset_names){
	
	my @fetched_rsets = @{$rset_a->fetch_all_by_name($rset_name)};
	my $num_rsets = scalar(@fetched_rsets);

	if($num_rsets > 1){
	  die("Failed to fetch unique ResultSet using name:\t$rset_name\nNeed to implement -result_set_ids\n");
	}
	elsif(! $num_rsets){
	  die("Failed to fetch ResultSet:\t$rset_name\n");
	}
	
	push @rsets, $fetched_rsets[0];
  }



  if( ((! $farm) || $on_farm) &&
	  scalar(@rsets) != 1){
	
	#force local here?
	#This would require a fetch param loop in dump code below.
	die("Cannot run multiple ResultSet dumps in one dump on the farm. Please specify -farm\n");
  }


  $dump_name     ||= $rsets[0]->name;  
  $file_name       = $dump_name;
  $file_name      .= '.'.$window_size;# if ! $on_farm;

  $feature_class   = 'ResultFeature';
  @fetch_params    = ($rsets[0], undef, $window_size);
  @header_params   = ($rsets[0], $window_size);#?
  @output_names    = @rset_names;
  $feature_adaptor = $db->get_ResultFeatureAdaptor;   #Move get_FeatureAdaptor to Set.pm?
 
  #Validate window_size
  $feature_adaptor->set_collection_defs_by_ResultSet($rsets[0]);
  
  if(! (defined $window_size  &&
		$feature_adaptor->has_window_size($window_size) )){
	die("You must suppliy a valid -window_size from:\t".join(" ", @{$feature_adaptor->window_sizes}));
  }


  #This is an intersect between format and feature config
   #need to handle vars too:  -clip $wig_file $chrom_file $bw_file";
  #This is going to be set name and wsize specific
  
  #Assuming bigwig here
  #This also assumes we chdir so paths are local
  $compression_args = " -clip ${file_name}.wig ${file_name}.chrom.sizes ${file_name}.bigwig";
}


my $lsf_dir;

#This enables correct naming of dump dir and lsf files for slice jobs
if($on_farm){
  $out_dir = $out_root.'/';
}else{
  $out_dir = $out_root.'/'.$dump_name.'/';

  if(! -d $out_dir){
	system('mkdir -p '.$out_dir) == 0 or die("Could not make out dir:\t$out_dir");
  }

  if($farm){
	$lsf_dir  =  $out_dir.'LSF/';
	$lsf_dir .= $window_size.'/' if defined $window_size;#ResultFeature specific at present

	if(! -d $lsf_dir){
	  system('mkdir -p '.$lsf_dir) == 0 or die("Could not make LSF out dir:\t$lsf_dir");
	}
  }
}


#This ducktyping could be replaced with a Role if we were using a more recent version of PERL
print "Dumping $format ${feature_class}s for:\t".join("\t", @output_names)."\n";
print "Output directory:\t$out_dir\n";



if(! @slices){
  print "No slices defined defaulting to current toplevel\n";
}
else{
  print "Restricting to slices:\t\t\t".join(', ', @slices)."\n";
}

if(@skip_slices){
  print "Skipping slices:\t\t\t".join(', ', @skip_slices)."\n";
}


@slices = @{&generate_slices_from_names($slice_a, \@slices, \@skip_slices, 'toplevel', undef, $inc_dups)};

#Don't run multi slice in local mode unless we are cating?
#Need to hoik this out into EFGUtils::submit_Slice_jobs?



if(! $farm){
  
  if(scalar(@slices) !=1){
	  
	if(! ($force_local || $post_process)){
	  die("You are running with multiple Slices, maybe you want to specify -farm or -force_local\n");
	}
  }

  #Turn headers off here where appropriate
  $write_sr_header = 0 if ! $conf_ref->{slice_headers};

  if($on_farm){
	#Turn of global file header, this will be added with -post_process
	#do we want to over-ride this?
	$write_file_header = 0;
  }
}
else{ #Submit to farm
 
  if($post_process){
	die("-post_process should not be run in -farm mode. Please correct you parameters");
  }

  $queue ||= 'normal';
  print "Submitting parallel ".scalar(@slices)." parallel Slice jobs\n";
  #Need to handle feature_set and arrays here? wrt merged dumps
  #if(@fsets){
  #}
  #else{#array
  #}

  #$0 is name of script file
  my ($script_dir, $script_name) = ($0 =~ /(^.*\/)(.*$)/);
  my @args = @{&strip_param_args(\@tmp_args, ('slices', 'skip_slices', 'dump_name', 'out_dir'))};
  @args = @{&strip_param_flags(\@tmp_args, ('farm'))};
  my $cmd = " ${script_dir}/${script_name} @args -out_root $out_dir -_on_farm";
 
  foreach my $slice(@slices){
	#Count features first and only bsub jobs with featues
	#add sr_id to meta_coord?
	#This would also allow range queries to be catered for slice rather than whole genome
	#Also change this to job array?
	#would need to do this prior to this block so we have the correct slices for post porcessing
	
	my $job_name = $file_name.'.'.$slice->seq_region_name;	
	#my $bsub_cmd="bsub -q $queue -J \"${job_name}[1-${num_jobs}]\" -o ${output_dir}/${job_name}.".'%J.%I'.".out -e ${output_dir}/${job_name}.".'%J.%I'.".err";
	my $bsub_cmd="bsub -q $queue -J \"${job_name}\" -o ${lsf_dir}/${job_name}.".'%J'.".out -e ${lsf_dir}/${job_name}.".'%J'.".err";

	#dump name can't be job name here
	my $cmd_args = " -dump_name $dump_name -slices ".$slice->name;
	print $bsub_cmd.$cmd.$cmd_args."\n";
	system($bsub_cmd.$cmd.$cmd_args) == 0 or die("Failed to submit job:\t$job_name\n$bsub_cmd\n${cmd}${cmd_args}");
  }

  exit;
}

#Should not fail after here but before loop to avoid unnecessary fail on the farm


#Define code refs in hash
#Can remove this when we move this to the obj e.g. $feature->$format_method
#Should this be done on the adaptor instead to avoid creating the object if possible?
#Iterator prevents the need to override output of obj_from_sth method
#We are almost always going to want the object for dumping
#So the format methods can go in the objs not the adaptors
#Move these up to the config block and set all as generic ftype_conf

my %format_subs = (
				   'get_ProbeFeature_BED'      => \&get_ProbeFeature_BED,
				   'get_RegulatoryFeature_GFF' => \&get_RegulatoryFeature_GFF,
				   'get_AnnotatedFeature_GFF'  => \&get_AnnotatedFeature_GFF,
				   'get_AnnotatedFeature_BED'  => \&get_AnnotatedFeature_BED,
				   'get_ExternalFeature_GFF'   => \&get_ExternalFeature_GFF,
				   'get_ResultFeature_WIG'        => \&get_ResultFeature_WIG,
				   'get_ResultFeature_WIG_header' => \&get_ResultFeature_WIG_header,
				   'get_ProbeFeature_BED_header'      => \&get_ProbeFeature_BED_header,
				   #'get_RegulatoryFeature_GFF_header' => \&get_RegulatoryFeature_GFF_header,
				   #'get_AnnotatedFeature_GFF_header'  => \&get_AnnotatedFeature_GFF_header,
				   'get_AnnotatedFeature_BED_header'  => \&get_AnnotatedFeature_BED_header,
				   #'get_ExternalFeature_GFF_header'   => \&get_ExternalFeature_GFF_header,		   
				  );

#Change all of these to Slice Iterators to avoid high mem usage
#and allow parallel dumps

my %fetch_methods = (
					 'ProbeFeature'      => 'fetch_Iterator_by_Slice_Arrays',
					 'ResultFeature'     => 'fetch_Iterator_by_Slice_ResultSet',
					 'AnnotatedFeature'  => 'fetch_Iterator_by_Slice_FeatureSets',
					 'ExternalFeature'   => 'fetch_Iterator_by_Slice_FeatureSets',
					 'RegulatoryFeature' => 'fetch_Iterator_by_Slice_FeatureSets',
					);



### DUMP THE FEATURES
my $format_method = "get_${feature_class}_${format}";
my $fetch_method    = $fetch_methods{$feature_class}; #Can't use hashref as method name?
my ($seq_name, $feature, $outline, @output);
my (@dump_files, $ofile_name);


foreach my $slice(@slices){
  #Build file name
  $seq_name = $slice->seq_region_name();
  $ofile_name = $file_name.'.'.$seq_name.'.'.lc($format);

  if($slice->coord_system->name eq 'chromosome'){
	$seq_name = 'chr'.$seq_name;
  }


  if(! $keep_colons){
	#Remove colons from file path which can cause problems with scp
	$ofile_name =~ s/\:/_/go;
  }

  my $ofile_path = $out_dir.$ofile_name;
  my $ofile;
  
  #Cat files
  if($post_process){
	push @dump_files, $ofile_name;
	next;
  }

  #Dump to file
  print "\nDumping:\t\t\t\t".$slice->name."\n";
  my $cnt = 0;


  

  #Use Iterator method here to avoid chunking
  #warn "Fetching ".ref($adaptors{$feature_class})."->$fetch_method($slice, @fetch_params)\n";
  #This assumes Slice is always first
  #Set this dynamically using slice_idx?
  my $feats = $feature_adaptor->$fetch_method($slice, @fetch_params);


  #Validate format method
  #if(! $feature->can($format_method)){
  #die ref($feature)." does not support the $format_method method\n";
  #}
  
  #sub these out to avoid repeating header method calls
  my $seen_feats = 0;
  
  if(ref($feats) eq 'Bio::EnsEMBL::Utils::Iterator'){
	
	#Do this after feature fetch to avoid creating empty file
	if($feats->has_next){
	  $seen_feats = 1;
	  $ofile = open_file($ofile_path, '>', 0775);
  
	  if($conf_ref->{has_header}){
		print $ofile $format_subs{$format_method.'_header'}->(@header_params, $slice);
	  }
	}
	
	while( ($feature = $feats->next) &&
		   defined $feature){
	  &print_line;
	}
  }
  else{# Normal ARRAYREF of features
	
	#Do this after feature fetch to avoid creating empty file
	if(scalar(@$feats)){
	  $seen_feats = 1;
	  $ofile = open_file($ofile_path, '>', 0775);
	  
	  if($conf_ref->{has_header}){
		print $ofile $format_subs{$format_method.'_header'}->(@header_params, $slice);
	  }

	  foreach $feature(@$feats){
		&print_line;
	  }
	}
  }
	
  sub print_line{
	$cnt++;	
	#Move all this to the relevant Feature class method get_FORMAT
	#Utilise base Feature methods for common attrs i.e. $self->SUPER->get_FORMAT
	
	#push @output, $feature->$get_format_method;
	push @output, &{$format_subs{$format_method}}($feature, $seq_name, $dbname, $feature_class);
	
	if(scalar(@output) == 1000){#1000
	  print $ofile join('', @output);
	  @output = ();
	}
  }


  if($seen_feats){
	print $ofile join('', @output);
	@output = ();
	close($ofile);

	print "Features dumped:\t\t\t$cnt\nOutput can be found here:\t\t$ofile_name\n";
	

	if($zip){
	  system("gzip -f $ofile_name") == 0 or die("Could not gzip file:\t$ofile_name");
	}
  }
  else{
	print "No features found, file not written\n";
  }
}


# Cat file mode and tidy up
# change this to post_process?

if($post_process){
  #do all this in the out dir to avoid long cat cmd
  my $ofile_path =  $out_dir.$file_name.'.'.lc($format);
  my $cwd =  Cwd::cwd();
  chdir($out_dir);
  
  if(-e $ofile_path){
	die("Found existing output file, please remove to -post_process:\n\trm -f $ofile_path\n");
  }


  #Build wig chrom file
  my $chrom_fh;

  if($format eq 'WIG'){
	#this should actuall be the same for all window sizes
	my $chrom_file =  $out_dir.$file_name.'.chrom.sizes';
	open($chrom_fh, ">$chrom_file") or die("Failed to open file for writing:\t".$chrom_file);
  }


  #Only use non-empty files
  my @cat_files;

  for my $i(0..$#dump_files){
	  
	if(-e $dump_files[$i] &&
	   ! -z $dump_files[$i]){  #exists and is not empty
	  
	  if ($format eq 'WIG'){
		my $seq_name = $slices[$i]->seq_region_name;
		
		if($slices[$i]->coord_system->name eq 'chromosome'){
		  $seq_name = 'chr'.$seq_name;
		}
	

		
	
		print $chrom_fh $seq_name."\t".$slices[$i]->length."\n" if $format eq 'WIG';

		#The onverhanging bins exceed the seq length which results in:
		#chromosome chr18 has 78077248 bases, but item ends at 78077520
		#This can be ignored?
		#test chr 11, as this is near the end


		push @cat_files, $dump_files[$i];
	  }
	}
  }

  my $hfile = open_file($ofile_path.'.header', '>', 0775);

  #This need to include the global file header and omit slice header params
  $write_sr_header = 0;
  print $hfile $format_subs{$format_method.'_header'}->(@header_params);
  close($hfile);

  unshift @cat_files, $ofile_path.'.header';
  my $cmd = "cat @cat_files > $ofile_path";
  print "\n".$cmd."\n";
  system($cmd) == 0 or die("Cannot cat files:\n$cmd\n$!");


  #validate cat with du?

  #compress and zip should be mutually exclusive?
  
  if($compress){
	$cmd = $bin_dir.$conf_ref->{compression_script}.' '.$compression_args;
	print "\n".$cmd."\n";
	system($cmd) == 0 or die("Failed to compress compress output:\n$cmd\n$!");
	push @cat_files, $ofile_path;
	
  }


  if($zip){
	system("gzip -f $ofile_path") == 0 or die("Could not gzip file:\t$ofile_path");
	#gzip removes source file so no need to add to @dump_files
  }
  
  if(! $no_clean){
	$cmd = "rm -f @cat_files";
	#print $cmd."\n";
	system($cmd) == 0 or die("Failed to remove files:\n$cmd\n");
  }

  chdir($cwd);
}



#Move these to obj or possibly an Interface/Role?

#These methods rely on obj creation

sub get_RegulatoryFeature_GFF{
  #my ($feature, $seq_name, $dbname, $feature_class) = @_
  my ($feature) = @_;
  my $gff = &get_GFF(@_);


  #Including MotifFeatures?
  #Do we want this?
  my @attrs;
		
  foreach my $reg_attr(@{$feature->regulatory_attributes()}){
		
	if($reg_attr->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature')){

	  #Only need cell type here for MultiCell
	  my $attr_name = $reg_attr->feature_type->name;
	  $attr_name .= ':'.$reg_attr->cell_type->name if $feature->cell_type->name eq 'MultiCell';
	  my $mf_attrs = &get_associated_MotifFeature_GFF_attributes($reg_attr);
	  $attr_name .= '('.$mf_attrs.')' if $mf_attrs;

	  #We should add MotifFeatures here so they have there AF context.
	  push @attrs, $attr_name;
	}
	elsif(! $reg_attr->isa('Bio::EnsEMBL::Funcgen::MotifFeature')){
	  #warn we have an unsupported ftype
	  warn "Found unsupported RegulatoryFeature attribute Feature type:\t".ref($reg_attr);
	}
  }
  
  return join("\t", (@{$gff}, 
					 join('; ', ('Name='.$feature->feature_type->name, 'ID='.$feature->stable_id(), 
								 'bound_start='.$feature->bound_start, 'bound_end='.$feature->bound_end, 
								 'Note=Consists of following features: '.join(' ', @attrs)))))
	."\n";

}



sub get_AnnotatedFeature_GFF{
  my ($feature) = @_;

  my $gff = &get_GFF(@_);


  #Including MotifFeatures?
  #Do we want this?
  #Add peak summit here?

  #Need to add set name?

  #Need to add cell type for merged dumps, or is fset name enough?
  #Maybe concat with ftype name in Name? HUVEC:H3K4me3?

  my $attr_string = 'Name='.$feature->feature_type->name;
  $attr_string .= "; Alias=".$feature->feature_set->name;#Change this to data set acc?
  my $pwm_names = &get_associated_MotifFeature_GFF_attributes($feature);
  $attr_string .= "; Note=Contains the following PWMs:${pwm_names}" if $pwm_names;
  push @$gff, $attr_string;

  return join("\t", @$gff)."\n"; 
}

sub get_associated_MotifFeature_GFF_attributes{
  my ($mf_assocd_feature) = @_;

  my $pwm_names;
  my (@pwm_names);
  

  foreach my $assoc_mf(@{$mf_assocd_feature->get_associated_MotifFeatures()}){
	#mf display_label or binding_matrix id/name here?
	push @pwm_names, $assoc_mf->binding_matrix->name;
  }
  
  if(@pwm_names){
	$pwm_names = join(' ', @pwm_names);
  }

  #Just return arrayref?
  return $pwm_names;
}


sub get_ExternalFeature_GFF{
  my $gff = &get_GFF;

  #It may be more appropriate to have display_label as the Name
  #Now add attributes
  return join("\t", @$gff)."\tName=".$feature->feature_type->name.'; Alias='.$feature->display_label."\n";
}

sub get_AnnotatedFeature_BED{
  my ($feature) = @_;
   my $bed = &get_BED(@_);
  #two blocks with a 1bp intron at the peak summit?

  #chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
  $bed->[3] = $feature->feature_type->name;
  
  return join("\t", @$bed)."\n";
}


#Header methods should optionally take Slice arg
#so we can support slice dump and full genome dump

sub get_ProbeFeature_BED_header{
  my ($vendor, $array) = @_;
  #   track name=pairedReads description="Clone Paired Reads" useScore=1
 
  #use dump name instead of vendor array here?
  #this may have been polluted with seq_region_name
 
  return "track name=$dump_name description=\"Ensembl ".$array->vendor.':'.$array->name." mappings\" useScore=0\n";
}

sub get_ProbeFeature_BED{
   my ($feature) = @_;
   my $bed = &get_BED(@_);
   #Could do with enabling custom name score methods here

   #chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
   
   #Alter get_complete_probename to drop array
   #use complete name if we have merged arrays dumps in same file
   my $probeset = $feature->probeset->name;
   $probeset = defined($probeset) ? $probeset.':' : '';   
   $bed->[3] = $probeset.$feature->probe->get_probename($array);

   #$bed->[4] = we don't store the similarity score
   #encode gapped alignments in the blocks

   return join("\t", @$bed)."\n";
}


#Some commonality for SetFeature GFF i.e. feature_type name?


#Following should be moved to BaseFeature(Adaptor?)
#Shouldn't need to create the object for this info? 
#But may be redundant calcs if we need to for the more specific feature data anyway.


sub get_GFF{
  my ($feature, $seq_name, $dbname, $feature_class) = @_;
  
  #http://www.sequenceontology.org/gff3.shtml
  #seq_id, source, type, start, end, score, strand, phase, attributes
  
  ##gff-version 3
  #ctg123  .  exon  1300  1500  .  +  .  ID=exon00001

  #Change feature class here? to be fset display label? Or something more generic obviosu for AnnotatedFeatures?
  #DNA Binding site?
  #Histone Modification etc

  my @gff = ($seq_name, $dbname, $feature_class, $feature->start, $feature->end, 
			 '.', $conf_ref->{strand}->{$feature->seq_region_strand}, '.');

  #Set size to mandatory fields
  #$#gff = 7;


  return \@gff;
}




sub get_BED{
  #my ($feature, $seq_name, $dbname, $feature_class) = @_;
  my ($feature, $seq_name) = @_;
  

  #chrom, chromStart, chromEnd, name, score, strand, thickStart, thickEnd, itemRgb, blockCount, blockSizes, blockStarts
  #track name=pairedReads description="Clone Paired Reads" useScore=1
  #chr22 1000 5000 cloneA 960 + 1000 5000 0 2 567,488, 0,3512

  #Ensembl is 1 based, but BED is 0 based i.e.half open.
  #already prefixing with chr by default
  

  #Define known attrs and mandatory fields
  my @bed = ($seq_name, undef, $feature->seq_region_end, undef, $conf_ref->{strand}->{$feature->seq_region_strand});
  
  #Set size to mandatory fields
  # $#bed = 5;

  if($ensembl_coords){
	$bed[1] = $feature->seq_region_start;
  }
  else{
	$bed[1] = ($feature->seq_region_start -1);
  }
  

  return \@bed;
}

sub get_ResultFeature_WIG_header{
  my ($rset, $window_size, $slice) = @_;

  #http://genome.ucsc.edu/goldenPath/help/wiggle.html

  #browser position chr19:59304200-59310700
  #browser hide all

  #track type=wiggle_0 name=track_label description=center_label
  #     visibility=display_mode color=r,g,b altColor=r,g,b
  #     priority=priority autoScale=on|off alwaysZero=on|off
  #     gridDefault=on|off maxHeightPixels=max:default:min
  #     graphType=bar|points viewLimits=lower:upper
  #     yLineMark=real-value yLineOnOff=on|off
  #     windowingFunction=maximum|mean|minimum smoothingWindow=off|2-16
  #variableStep  chrom=chrN  [span=windowSize]
  #fixedStep  chrom=chrN  start=position  step=stepInterval  [span=windowSize]


#The track type with version is REQUIRED, and it currently must be wiggle_0:
#  type wiggle_0
#The remaining values are OPTIONAL:
#  name              trackLabel           # default is "User Track"
#  description       centerlabel          # default is "User Supplied Track"
#  visibility        full|dense|hide      # default is hide (will also take numeric values 2|1|0)
#  color             RRR,GGG,BBB          # default is 255,255,255
#  altColor          RRR,GGG,BBB          # default is 128,128,128
#  priority          N                    # default is 100
#  autoScale         on|off               # default is on
#  alwaysZero        on|off               # default is off
#  gridDefault       on|off               # default is off
#  maxHeightPixels   max:default:min      # default is 128:128:11
#  graphType         bar|points           # default is bar
#  viewLimits        lower:upper          # default is range found in data
#  yLineMark         real-value           # default is 0.0
#  yLineOnOff        on|off               # default is off
#  windowingFunction maximum|mean|minimum # default is maximum
#  smoothingWindow   off|[2-16]           # default is off
  



  #These are actually two headers:
  #Global file header - Should only be written in local mode i.e. single slice or cat_files
  #What about force local for multi slice dumps? $write_gheader = 1;
  #
  #Seq region header  - Always written


  my $header;

  #Move this to get_ResultFeature_WIG_file_header?
  #then remove need for $write_sr_header


  if($write_file_header){
	$header = 'track type=wiggle_0 name="'.$rset->name.'" description="'.$rset->display_label."\"\n"; 
	#Is windowing function max magnitude or just max? 'windowingFunction=maximum'
	$write_file_header = 0;
  }


  #ResultSet could do with a display name/label, but this is on the FeatureSet
  
  if($write_sr_header){
	$header .= "fixedStep chrom=".$seq_name." start=".$slice->start.
	  ' step='.$window_size.' span='.$window_size."\n";
  }

  return $header;

}

sub get_ResultFeature_WIG{
  my ($feature) = @_;
  #4 sig fig, as this just captures the noise
  #Maybe we don't want to capture this?
  #Maybe 3 sig fig will drop the noise in favour of speed?
  return join("\n", (map { sprintf("%.4f" , $_) } @{$feature->scores}))."\n";
}
