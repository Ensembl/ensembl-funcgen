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

# I guess the last issue is bloating the class rather than having a separate Exporter
# which knows how to translate various features to each format. But this is necessary if we want to 
# avoid the extra iteration etc...

#For now have all subs here, in a way that can easily be moved to the relevant adaptors.


$| =1;

my ($pass, @fset_names, @rset_names, $array, $vendor, $dump_name, $queue, $on_farm, $merged_dump);
my ($exp_name, $help, $pids, $man, @features, $chr, $not_status, $half_open, $farm, $cat_dumps);
my ($dbhost, $port, $user, $dbname, $cdbname, $species, @slices, @skip_slices, $out_root);
my ($dnadb_pass, $dnadb_user, $dnadb_name, $dnadb_port, $dnadb_host, $with_headers, $force_local);
my ($no_clean, $cat_files);
my $format = 'GFF';
my $zip = 0;
my $out_dir = ".";
my $keep_colons = 0;
my @tmp_args = @ARGV;

warn "$0 @tmp_args";


GetOptions 
  (
   'feature_sets=s{,}' => \@fset_names,
   'result_sets=s{,}'  => \@rset_names,
   #Change this to feature_sets to allow merged dumps from same feature class.
   'array=s'          => \$array,#Change this to \@arrays?
   'vendor=s'         => \$vendor,
   'format=s'         => \$format,
   'half_open=s'      => \$half_open,

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
 

   'force_local'      => \$force_local,
   'farm'             => \$farm,  
   '_on_farm'         => \$on_farm,
   'with_headers=i'   => \$with_headers,
   'merged_dump'      => \$merged_dump,
   'cat_files'        => \$cat_files,
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
my $array_a  = $db->get_ArrayAdaptor;
my ($feature_class, @fsets, @fetch_params, @header_params, @output_names);


if(@fset_names && 
   ($array || $vendor)){
  die("Please choose either array and vendor or feature_sets parameters\n");
}

if(! @fset_names){

  if(! ($array && $vendor)){
	die("Please choose either array and vendor or feature_sets parameters\n");
  }
  @output_names = ($array);
  my $array = $array_a->fetch_by_name_vendor($array, $vendor);

  if(! $array){
	die("Could not fetch array:\t$vendor\t$array\n");
  }

  $dump_name ||= $array->vendor.'_'.$array->name;
  $feature_class = 'ProbeFeature';
  @fetch_params = ([$array]);
  @header_params = ($vendor, $array);#arrays?
}
else{
 
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

  $feature_class = ucfirst($feature_class).'Feature';
  @fetch_params  = (\@fsets);
  @header_params = ();#?
  @output_names  = @fset_names;
}



#Need some format config here
#i.e. do we print headers for slice files or not?

#This enables correct naming of dump dir and lsf files for slice jobs

if($on_farm){
  $out_dir = $out_root;
}else{
  $out_dir = $out_root.'/'.$dump_name.'/';
  system('mkdir -p '.$out_dir) if(! -d $out_dir);
}


#This mean lsf out will be in the parent directory

$format = uc($format); #validation done via introspection on Feature object
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


@slices = @{&generate_slices_from_names($slice_a, \@slices, \@skip_slices, 1)};#toplevel flag




#Don't run multi slice in local mode unless we are cating?

#Need to hoik this out into Utils::submit_Slice_jobs?


#Add merge mode?
#Chrom files first in merge file?

#Need to turn headers on for local dump, but only if not farmed!!???
#This is format specific! So need config has for formats
#This should go in generic export config



my %format_config = (
					 bed => {
							 slice_headers => 0,
							},

					 wig => {
							 slice_headers => 1, #?
							},
					);


#Need to tidy this up to allow over-write of config
#i.e. undef params should use config
if(! $on_farm && ! defined $with_headers){
  $with_headers = 1;
}



if(! $farm){
  
  if(scalar(@slices) !=1){
	
	if(! $force_local){
	  die("You are running with multiple Slices, maybe you want to specify -farm or -force_local\n");
	}
	
	#Turn headers off here where appropriate
	$with_headers = 0 if ! $format_config{$format}{slice_headers};
  }
}
else{

  if($cat_files){
	die("-cat_files should not be run in -farm mode. Please correct you parameters");
  }

  $queue ||= 'normal';
  
  #need to strip params here

  #Do this as job array or with slice names?
  

  #"-job_name $job_name"
  print "Submitting parallel ".scalar(@slices)." parallel Slice jobs\n";

  #$0 is name of script file
  my ($script_dir, $script_name) = ($0 =~ /(^.*\/)(.*$)/);
  my @args = @{&strip_param_args(\@tmp_args, ('slices', 'skip_slices', 'dump_name', 'out_dir'))};
  @args = @{&strip_param_flags(\@tmp_args, ('farm'))};

   

  #Need to handle feature_set and arrays here? wrt merged dumps
  if(@fsets){
	
  }
  else{#array
	
  }


  my $cmd = " ${script_dir}/${script_name} @args -out_root $out_dir -_on_farm";
 
  foreach my $slice(@slices){

	
	my $job_name = $dump_name.':'.$slice->seq_region_name;
	#my $bsub_cmd="bsub -q $queue -J \"${job_name}[1-${num_jobs}]\" -o ${output_dir}/${job_name}.".'%J.%I'.".out -e ${output_dir}/${job_name}.".'%J.%I'.".err";
	my $bsub_cmd="bsub -q $queue -J \"${job_name}\" -o ${out_dir}/${job_name}.".'%J'.".out -e ${out_dir}/${job_name}.".'%J'.".err";

	my $cmd_args = " -dump_name $job_name -slices ".$slice->seq_region_name;

	print $bsub_cmd.$cmd.$cmd_args."\n";

	system($bsub_cmd.$cmd.$cmd_args) == 0 or die("Failed to submit job:\t$job_name\n$bsub_cmd\n${cmd}${cmd_args}");
  }

  exit;

}

#Should not fail after here but before loop to avoid unnecessary fail on the farm



my $format_method = "get_${feature_class}_${format}";
#my $header_method = $format_method.'_header';


#Define code refs in hash
#Can remove this when we move this to the obj e.g. $feature->$format_method
#Should this be done on the adaptor instead to avoid creating the object if possible?
#Iterator prevents the need to override output of obj_from_sth method
#We are almost always going to want the object for dumping
#So the format methods can go in the objs not the adaptors

#We still really need to wrap up Slice chunking in the Iterator
#As some queries are still too big?
#Can use code refs instead to pass the variables?
#i.e. pass code ref and params separately

my %format_subs = (
				   'get_ProbeFeature_BED'      => \&get_ProbeFeature_BED,
				   'get_RegulatoryFeature_GFF' => \&get_RegulatoryFeature_GFF,
				   'get_AnnotatedFeature_GFF'  => \&get_AnnotatedFeature_GFF,
				   'get_AnnotatedFeature_BED'  => \&get_AnnotatedFeature_BED,
				   'get_ExternalFeature_GFF'   => \&get_ExternalFeature_GFF,
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
					 #'ProbeFeature'      => 'fetch_all_by_Slice_Arrays',
					 'AnnotatedFeature'  => 'fetch_all_by_Slice_FeatureSets',
					 'ExternalFeature'   => 'fetch_all_by_Slice_FeatureSets',
					 'RegulatoryFeature' => 'fetch_all_by_Slice_FeatureSets',
					 
					);
my %adaptors;

if($feature_class eq 'ProbeFeature'){
  $adaptors{$feature_class}	= $db->get_ProbeFeatureAdaptor;
}else{
  $adaptors{$feature_class}	= $fsets[0]->get_FeatureAdaptor;
}

#my $feature_adaptor = $adaptors{$feature_class};
my $fetch_method    = $fetch_methods{$feature_class}; #Can't use hashref as method name?
my ($seq_name, $feature, $outline, @output);


my %gff_strand = ( 0  => '.',
				   1  => '+',
				   -1 => '-',
				 );

my %bed_strand = ( 0  => '',
				   1  => '+',
				   -1 => '-',
				 );


#How are we going to handle merged dumps with headers?
#-no_headers on farmed slice dumps?
#Then run in -merge mode to cat and create correct header for complete file?
#Make -no_headers over-ridable
my (@dump_files, $ofile_name);


foreach my $slice(@slices){

  #Build file name
  $seq_name = $slice->seq_region_name();

  if($on_farm){
	#We have already changed the dump name to include the slice
	$ofile_name = $dump_name.'.'.lc($format);
  }
  else{
	$ofile_name = $dump_name.'.'.$seq_name.'.'.lc($format);
  }


  if($slice->coord_system->name eq 'chromosome'){
	$seq_name = 'chr'.$seq_name;
  }


  if(! $keep_colons){
	#Remove colons from file patch which can cause problems
	#with scp
	$ofile_name =~ s/\:/_/go;
  }

  my $ofile_path = $out_dir.'/'.$ofile_name;

  
  #Cat files
  if($cat_files){


	if(! -e $ofile_path){
	  die("Could not find input for -cat_files merge:\t${ofile_path}");
	}
	elsif(-z $ofile_path){
	  print "Omitting empty file from file merge:\t${ofile_name}\n";
	  
	  if(! $no_clean){
		unlink $ofile_path;
	  }
	}
	else{
	  push @dump_files, $ofile_name;
	}
	
	next;
  }

  #Dump to file
  print "\nDumping:\t\t\t\t".$slice->name."\n";
  my $cnt = 0;



  my $ofile = open_file($ofile_path, '>', 0775);

  if($with_headers){
	print $ofile $format_subs{$format_method.'_header'}->(@header_params, $slice);
  }


  #Need to handle dumping header here
  

  #Use Iterator method here to avoid chunking
  #warn "Fetching ".ref($adaptors{$feature_class})."->$fetch_method($slice, @fetch_params)\n";
  my $feats = $adaptors{$feature_class}->$fetch_method($slice, @fetch_params);


  #Validate format method
  #if(! $feature->can($format_method)){
  #die ref($feature)." does not support the $format_method method\n";
  #}

  if(ref($feats) eq 'Bio::EnsEMBL::Utils::Iterator'){

	while($feature = $feats->next()){
	  &print_line;
	}
  }else{
 
	foreach $feature(@$feats){
	  &print_line;
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

  print $ofile join('', @output);
  @output = ();
  close($ofile);

  #Remove empty files and zip
  if(-z $ofile_path){
	print "No features found\n";
	#unlink $ofile_name; #Move this to cat_files?
  }
  else{
	print "Features dumped:\t\t\t$cnt\n";

	if($zip){
	  system("gzip -f $ofile_name") == 0 or die("Could not gzip file:\t$ofile_name");
	}
  }
}


# CAT File mode and tidy up

if($cat_files){
  #do all this in the out dir to avoid long cat cmd
  my $cwd =  Cwd::cwd();
  chdir($out_dir);
  my $ofile_path =  $out_dir.'/'.$dump_name.'.'.lc($format);
  
  if(-e $ofile_path){
	#fail?
  }

  my $hfile = open_file($ofile_path.'.header', '>', 0775);
  print $hfile $format_subs{$format_method.'_header'}->(@header_params);
  close($hfile);

  unshift @dump_files, $ofile_path.'.header';

  my $cmd = "cat @dump_files > $ofile_path";

  system($cmd) == 0 or die("Cannot cat files:\n$cmd\n$!");


  if($zip){
	system("gzip -f $ofile_path") == 0 or die("Could not gzip file:\t$ofile_path");
	#gzip removes source file so no need to add to @dump_files
  }
  
  if(! $no_clean){
	$cmd = "rm -f @dump_files";
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
	  
	  #We should add MotifFeatures here so they have there AF context.
	  push @attrs, $attr_name.&get_associated_MotifFeature_GFF_attributes($reg_attr);
	}
	elsif(! $reg_attr->isa('Bio::EnsEMBL::Funcgen::MotifFeature')){
	  #warn we have an unsupported ftype
	  warn "Found unsupported RegulatoryFeature attribute Feature type:\t".ref($reg_attr);
	}
  }
  
  return join("\t", (@{$gff}, 
					 join('; ', ('Name='.$feature->feature_type->name, 'ID='.$feature->stable_id(), 
								 'bound_start='.$feature->bound_start, 'bound_end='.$feature->bound_end, 
								 'Note=Consists of following features: '.join(',', @attrs)))))
	."\n";

}



sub get_AnnotatedFeature_GFF{
  my ($feature) = @_;

  my $gff = &get_GFF(@_);


  #Including MotifFeatures?
  #Do we want this?
  #Add peak summit here?

  my $attr_string = 'Name='.$feature->feature_type->name;
  my $pwm_names = &get_associated_MotifFeature_GFF_attributes($feature);
  $attr_string .= "; Note=Contains the following PWMs:${pwm_names}" if($pwm_names);
  push @$gff, $attr_string;

  return join("\t", @$gff)."\n"; 
}

sub get_associated_MotifFeature_GFF_attributes{
  my ($mf_assocd_feature) = @_;

  my $pwm_names = ''; #empty string as we cat in caller
  my (@pwm_names);
  

  foreach my $assoc_mf(@{$mf_assocd_feature->get_associated_MotifFeatures()}){
	#mf display_label or binding_matrix id/name here?
	push @pwm_names, $assoc_mf->binding_matrix->name;
  }
  
  if(@pwm_names){
	$pwm_names = '('.join(',', @pwm_names).')';
  }

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
			   '.', $gff_strand{$feature->seq_region_strand}, '.');

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
  my @bed = ($seq_name, undef, $feature->seq_region_end, undef, $bed_strand{$feature->seq_region_strand});
  
  #Set size to mandatory fields
  # $#bed = 5;

  if(! $half_open){
	$bed[1] = ($feature->seq_region_start -1);
  }
  else{
	$bed[1] = $feature->seq_region_start;
  }
  

  return \@bed;
}



