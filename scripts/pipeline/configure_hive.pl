#!/usr/bin/env perl

=head1 LICENSE


  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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



=head1 SYNOPSIS

  configure_hive.pl -host <string> -user <String> -pass <String> -dbname <String> -data_root <String> \
    -hive_script_dir <String> -configs <String> ... \
    [-port <Int> -dnadb_host <String> -dnadb_user <String> -dnadb_pass <String> -dnadb_name <String> \
     -dnadb_port <Int> -list -help -man]

=head1 PARAMETERS

  Mandatory:
    -data_root        <String>      Root data directory
    -hive_script_dir  <String>      Hive scripts directory
    -configs          <String> ...  Config module names e.g. Peaks (use -list for full lst)
    -host             <String>      Funcgen DB host
    -user             <String>      Funcgen DB user
    -pass             <String>      Funcgen DB pass
    -dbname           <String>      Fincgen DB name
    
  Optional:
    -port             <Int>         Funcgen database port
    -dnadb_host       <String>      Core DB host
    -dnadb_user       <String>      Core DB user
    -dnadb_pass       <String>      Core DB pass
    -dnadb_name       <String>      Core DB name
    -dnadb_port       <Int>         Core DB port
    -list             Print a list of the supported config module names
    -help             Prints a helpful message
    -man              Prints the man page

=head1 DESCRIPTION

This script is a wrapper for the main init_pipeline.pl hive script. Given a list
if config modules (-configs) it will initialise/create the pipeline DB if it 
does not exist and iteratively perform analysis_topup. It also inserts a hive_url 
meta_key into the output DB, and hive_conf meta keys into the pipeline DB.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;

use Bio::EnsEMBL::Funcgen::Utils::DBAdaptorHelper qw( process_DB_options
                                                      get_DB_options_config
                                                      create_Funcgen_DBAdaptor_from_options
                                                      create_DBAdaptor_from_params );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils        qw( run_system_cmd
                                                      url_from_DB_params
                                                      add_hive_url_to_meta );

#$| = 1;#for debug

&main();


sub main{
my (@configs, $species, $data_root, $hive_script_dir, $list);

#put this config in HiveHelper? Where else mightit be needed?
my %config_info = 
 (
  #ConfigModule => [priority, pre-req],
  #priorities should be unique!
  #although Peaks & Collections are effectively equal priority here
  #they need to be different numbers due to usage as indexes below
 
  'ReadAlignments'    => [5],#This will deal with InputSubsets
  #then we subsample(1/N reps), merge and rm dups again for final Peaks input

  #Where will the InputSets be defined?
  #after read alignment QC as we may want to drop some fastqs?
  #by default this will work on all subsets in an experiment?
  #what about control?
     
  #what was the argument that we should move InputSet registration to the pipeline?
  #Can we define a set of -controls_experiments which have useable controls for batch
  #i.e. same 'study' and match cell_type
  #this should throw if there is more than one experiment with control for a given cell type
  #as we don't know which one to pick
  #should throw if non-can be found
  #Experiments with no controls should be submitted with an -allow_no_controls flag
  #These options should be mutually exclusive to avoid danger of submitting
  #an experiment without control data by mistake
  #No longer subsampling!!!! Just merge! <---------------THIS
               
  #These can be linked by IdentifyInputSets
               
  #'DefineInputSets'   => [6],
  #Can't this be integrated into IDRPeaks?
  #This would well and truly make it impossible to 
  #run IDRPeaks an (non-IDR) Peaks in the same pipeline.
  
  #These can be linked by DefineOutputSets::IdentifyIDRInputSets
  #can't clash with IdentifyInputSets ->rename as IdentifyPeaksInputSets
  #Can we have two 'link in' analyses in the same config?
  #How do we protect again the same set details being flowed through both configs?
  #This may fail due to the presence of a duplicate input_id
  #when flowing to the next analysis
  #but not before we might attempt a rollback (if rollback/recover is set)
  #This may end up rolling back a set which is actively being processed.
  #Leave this for now due to safety concerns
               
               
  
  #Will we ever want to miss IDR out and link directly to Peaks?
  #We could do this by specifying a no_IDR batch param.
  #and dataflowing accordingly
  #DefineInputSets will create Merged and single rep IDR InputSets, dependat on -no_idr flag
  
  'IDRPeaks'          => [7],#Don't need to have a pre-req here
  #as we don't have a branch here, so this can integrate DefineReplicateOutputSets
  #This will have to create analyses for the FeatureSet
  #This will be flowed directly to DefineOutputSet(skipping IdentifyInputSets)
  #Where will the merged InputSet be created?
  #Must be in IDRPeaks
  #IDRPeaks will also create output sets, so we don't lose information
  #i.e. the IDR generated parameters need to be attached to the final FeatureSet
  
  
  #Let's define
  
  
  #The following currently require that the InputSets have already been defined
  'DefineOutputSets'      => [8], #A pre-req should not have a pre-req of it's own?
  'Peaks'                 => [9,     'DefineOutputSets'],
  'Collections'           => [10,    'DefineOutputSets'],
 );

#These are separate as we sort based on priority which gives an undef warning below
my %config_labels = 
 (
  'Peaks_and_Collections' => [undef, 'Peaks', 'Collections'],
  'UberPipe'              => [undef, 'Peaks', 'Collections', 'IDRPeaks',],
 );
                  
                  

#For now we are not storing the conf keys in the output DB
#but we may want to move them there, if we are to allow
#more than one hive to run on the same DB.


my @tmp_args = @ARGV;
my $db_opts  = get_DB_options_config();#This will get opts for funcgen, core and pipeline by default

GetOptions (
            #Mandatory
            %{$db_opts},
            'configs=s{,}'      => \@configs,
            'data_root=s'       => \$data_root,
            'hive_script_dir=s' => \$hive_script_dir,

            #Optional
            #'drop' => \$drop,#implement this when we use the updated version of init_pipeline?
            'species=s'   => \$species,
            'list'        => \$list,
            
            'help'             => sub { pod2usage(-exitval => 0); }, 
            #removed ~? frm here as we don't want to exit with 0 for ?
            
            'man|m'            => sub { pod2usage(-exitval => 0, -verbose => 2); },
           ) or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 
           
#Do we want to catch the rest of ARGV here by using -- in the input and assign these to init_pipeline as extra args?


### VALIDATE PARAMETERS ###

my @valid_confs = sort {$config_info{$a}->[0] <=> $config_info{$b}->[0]} keys %config_info;
#add some sprintf lpad action in here
push @valid_confs, map("$_ (= ".join(' ', @{$config_labels{$_}}[1..$#{$config_labels{$_}}]).')', 
                       keys %config_labels);
%config_info = (%config_info, %config_labels); #Add labels after sort

if($list){
  print "Valid config modules/labels are:\n\t".join("\n\t", @valid_confs)."\n";
  #This currently doesn't list those which are pre-reqs
  exit;  
}


if(! defined $hive_script_dir){
  throw('-hive_script_dir is a mandatory parameter');  
}      
elsif(! -d $hive_script_dir){
  throw("-hive_script_dir is not a valid directory:\t${hive_script_dir}");
}

if(! @configs){
  pod2usage(-exitval => 1, 
            -message => "Must provide at least one config module/label name. Valid names are:\n\t".
                        join("\n\t", @valid_confs));  
}

### HANDLE DB PARAMS ###

my $db_script_args = process_DB_options($db_opts, ['funcgen', 'core', 'pipeline'], undef, 'script');
my $pdb_params     = ${&process_DB_options($db_opts, ['pipeline'])}{pipeline};
my $db             = create_Funcgen_DBAdaptor_from_options($db_opts, 'pass');


### PRE-PROCESS CONFIGS ###

my (@confs, %abbrvs);

foreach my $conf(@configs){ 
  #Could remove this in favour of the throw in populate_config_array
  if(! exists $config_info{$conf}){
    pod2usage(-exitval => 1, 
              -message => "'$conf' is not a supported config module/label. Valid names are:\n\t".
                          join("\n\t", @valid_confs));  
  }  

  &populate_config_array($conf, \@confs, \%config_info);
}



#Now we have a pre-ordered (possibly gappy) array of confs to add
#Strip out any undef elements
my @tmp_confs;
map {push @tmp_confs, $_ if defined $_} @confs;
@confs = @tmp_confs;

### INITIALISE/CREATE PIPELINE DB ###

my $pipeline_params = "-data_root_dir $data_root -pipeline_name ".$pdb_params->{'-dbname'}.' '.
  $db_script_args->{funcgen}.' '.$db_script_args->{core}.' '.$db_script_args->{pipeline};
$pipeline_params .= " -species $species " if defined $species;  

my $pdb;  
eval { $pdb = create_DBAdaptor_from_params($pdb_params, 'core', 1); };

if($@){ #Assume the DB hasn't been created yet  
  #init the pipline with the first conf
  my $first_conf = shift @confs;
  my $init_cmd   = "perl $hive_script_dir/init_pipeline.pl ".
    "Bio::EnsEMBL::Funcgen::Hive::Config::${first_conf} $pipeline_params";
  
  print "\n\nINITIALISING DATABASE:\t".$pdb_params->{'-dbname'}."\n";
  run_system_cmd($init_cmd);  
  
  $pdb = create_DBAdaptor_from_params($pdb_params, 'core');
}

add_hive_url_to_meta(url_from_DB_params($pdb_params), $db);



### PERFORM ANALYSIS_TOPUP ###
my $mc         = $pdb->get_MetaContainer;
my $conf_key   = 'hive_conf';
my @meta_confs = @{$mc->list_value_by_key($conf_key)};

foreach my $conf(@confs){
 
  if( grep(/^$conf$/, @meta_confs) ){
    warn "Skipping hive -analysis_topup.  $conf config has already been added to the DB\n";  
  }
  else{
    my $topup_cmd = "perl $hive_script_dir/init_pipeline.pl Bio::EnsEMBL::Funcgen::Hive::Config::${conf} ".
      ' -analysis_topup '.$pipeline_params;  
    #warn $topup_cmd."\n";
    print "\n\nPERFORMING ANALYSIS TOPUP:\t".$conf."\n";
    run_system_cmd($topup_cmd);
    
     #Add key via API to store with appropriate species_id
    eval { $mc->store_key_value('hive_conf', $conf); };
    
    if($@){
      throw("Failed to store hive conf meta entry:\t$conf\n$@");  
    }
  }
}

}# end of main


sub populate_config_array{
  my ($conf_name, $confs, $config_info) = @_;  
  
  #Internal validation of config hash
  if(! exists $config_info->{$conf_name}){
    #Handle this in caller as this will also know the parent conf
    #No, need to handle here as this will recurse
    die("'$conf_name' is defined a pre-requisite, but is not defined in the config\nPlease amend config!"); 
  }
 
  #Must have at least a priority(can be undef) and a conf name
  if(scalar(@{$config_info->{$conf_name}}) < 1){
    die("'$conf_name' config is not valid, must have at least a priority and/or a pre-requisite module."); 
  }
 
  my $priority = $config_info->{$conf_name}->[0];
 
  if(defined $priority){ #We have a config rather than a label
  
    if(defined $confs->[$priority]  &&
      ($confs->[$priority] ne $conf_name)){
      die("Found a priority clash between configs:\t$conf_name and ".
            $confs->[$priority]."\nPlease correct config hash");    
    }
      
    $confs->[$priority] = $conf_name;
  }

  foreach my $conf_index(1..$#{$config_info->{$conf_name}}){
    &populate_config_array($config_info->{$conf_name}->[$conf_index], $confs, $config_info);
  }
  
  return $confs;
}# end of populate_config_array


### TA DAA! ###

1;
