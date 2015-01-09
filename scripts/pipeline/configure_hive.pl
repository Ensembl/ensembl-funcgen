#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

configure_hive.pl

=head1 SYNOPSIS

  configure_hive.pl [-configs <String> ...] \          #Normally all that is required when running within pipeline environent
    [-host <string> -user <String> -pass <String> -dbname <String> \
     -data_root <String> -hive_script_dir <String> ] \ #Mandatory params specified by pipeline enviornment
    [-port <Int> -dnadb_host <String> -dnadb_user <String> -dnadb_pass <String> 
     -dnadb_name <String> -dnadb_port <Int> ] \        #Optional params specified by the pipeline environment
    [ -list -help -man ] \                             # Truly optional params

=head1 PARAMETERS

  Mandatory:
    -configs          <String> ...  Config module names e.g. Peaks (use -list for full lst)

  Madatory if not run within pipeline environment:
    -data_root        <String>      Root data directory
    -hive_script_dir  <String>      Hive scripts directory
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
    -force            Forces analysis top up, even if config already exists in DB
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
                                                      add_DB_url_to_meta );
use Bio::EnsEMBL::Funcgen::Hive::Utils            qw( inject_DataflowRuleAdaptor_methods );

#TODO
# 1 Genericise pipeline param passing. Current allow_no_arch and archive_root are hardcoded.


my %config_info = 
 (
  #ConfigModule => [priority, pre-req],
  #ConfigLabel  => [undef, (List, Of, Consituent, Configs)]
  #Priorities should be unique and represent logical flow of pipeline!
  #Although Peaks & Collections are effectively equal priority here
  #they need to be different numbers due to usage as indexes below
  #Actual priority number are arbitrary, 
  #they just need to be in the right order
 
  'ReadAlignment'        => [5],
  'IDRPeaks'             => [7],
  'DefineMergedDataSets' => [8],
  'Peaks'                => [9,     'DefineMergedDataSets'],
  'Collections'          => [10,    'DefineMergedDataSets'],
 );

#These are separate as we sort based on priority which gives an undef warning below
my %config_labels = 
 (
  'Peaks_and_Collections' => [undef, 'Peaks', 'Collections'],
  'UberPipe'              => [undef, 'Peaks', 'Collections', 'IDRPeaks', 'ReadAlignment'],
 );

#Todo
#For now we are not storing the conf keys in the output DB
#but we may want to move them there, if we are to allow
#more than one hive to run on the same DB.

&main();


sub main{
  my @tmp_args = @ARGV;
  my (@configs, $species, $data_root, $hive_script_dir, $list);   
  my ($archive_root, $allow_no_arch, $force);            
  my $db_opts  = get_DB_options_config();#This will get opts for funcgen, core and pipeline by default
  
  GetOptions
   (#Mandatory
    %{$db_opts},
    'configs=s{,}'      => \@configs,
    'data_root=s'       => \$data_root,
    'hive_script_dir=s' => \$hive_script_dir,
  
    #Optional
    'archive_root=s'    => \$archive_root,
    'allow_no_archive'  => \$allow_no_arch,
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
  %config_info = (%config_info, %config_labels);
  
  if($list){
    print "Valid config modules/labels are:\n\t".join("\n\t", @valid_confs)."\n";
    #This currently doesn't list those which are just pre-reqs
    exit;  
  }
  
  
  if(! defined $hive_script_dir){
    throw('-hive_script_dir is a mandatory parameter');  
  }elsif(! -d $hive_script_dir){
    throw("-hive_script_dir is not a valid directory:\t${hive_script_dir}");
  }
  
  #Do some sensible pre-init dir checking, 
  #Thiss would be caught by the pipeline, but we don't want to 
  #have to tidy up if these are wrong
  
  if(! defined $data_root){
    throw('-data_root is a mandatory parameter');  
  }elsif(! -d $data_root){
    throw("-data_root is not a valid directory:\t${data_root}");
  }
  
  my $arch_params = '';
  
  if(defined $archive_root){
    
    if(! -d $archive_root){
      throw("-archive_root is not a valid directory:\t${archive_root}");  
    }
    $arch_params .= " -archive_root $archive_root ";  
  }
  
  $arch_params .= ' -allow_no_archive 1 ' if defined $allow_no_arch;
  
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
  
    &populate_config_array($conf, \@confs);
  }
  
  #Strip out any undef elements
  my @tmp_confs;
  map {push @tmp_confs, $_ if defined $_} @confs;
  @confs = @tmp_confs;
  
  
  ### INITIALISE/CREATE PIPELINE DB ###
  my $pipeline_params = "$arch_params -data_root_dir $data_root -pipeline_name ".$pdb_params->{'-dbname'}.' '.
    $db_script_args->{funcgen}.' '.$db_script_args->{core}.' '.$db_script_args->{pipeline};
  $pipeline_params .= " -species $species " if defined $species;  
  
  my ($cmd, $ntable_a, $pdb);  
  eval { $pdb = create_DBAdaptor_from_params($pdb_params, 'hive', 1); };
  
  if($@){ #Assume the DB hasn't been created yet  
    #init the pipline with the first conf
    my $first_conf = shift @confs;
    $cmd = "perl $hive_script_dir/init_pipeline.pl ".
      "Bio::EnsEMBL::Funcgen::Hive::Config::${first_conf} $pipeline_params";
    
    print "\n\nINITIALISING DATABASE:\t".$pdb_params->{'-dbname'}."\n";
    run_system_cmd($cmd);  
    
    $pdb      = create_DBAdaptor_from_params($pdb_params, 'hive');
    $ntable_a = $pdb->get_NakedTableAdaptor;
    $ntable_a->table_name('meta');
    _register_conf_in_meta($ntable_a, $first_conf);    
  }
  else{
    #$mc         = $pdb->get_MetaContainer; 
    $ntable_a = $pdb->get_NakedTableAdaptor;
    $ntable_a->table_name('meta');
  }
   
  my $hive_url = url_from_DB_params($pdb_params); 
  add_DB_url_to_meta('hive', $hive_url, $db);
  
  ### PERFORM ANALYSIS_TOPUP ###
  my $conf_key   = 'hive_conf';
  
  #my @meta_confs = @{$mc->list_value_by_key($conf_key)};#old method used this
  my @meta_confs = map {$_->{meta_value}} @{$ntable_a->fetch_all_by_meta_key($conf_key)};
  my $dfr_adaptor = $pdb->get_DataflowRuleAdaptor;
  inject_DataflowRuleAdaptor_methods($dfr_adaptor); #Injects get_semaphoring_analysis_ids
  
  foreach my $conf(@confs){
   
    if( grep(/^$conf$/, @meta_confs) ){
      warn "Skipping hive -analysis_topup.  $conf config has already been added to the DB\n";  
    }
    else{ #Add new config!
      
      #Handle potential resetting of pipeline wide 'can_run_AnalaysisLogicName' params 
      #These should be in the meta table and should be cached if they are 'true' 
      #as there is a danger that a subsequent top up of an upsteam conf may reset this to 0
      #This would result dataflow not occuring from the conf just added through the link 
      #analysis to the next conf(which has be added previously)
      #Put this method in HiveUtils? (with add_hive_url_to_meta?) where else would it be used?
                 
      my $meta_key        = 'can_%';
      my %meta_key_values = %{$ntable_a->fetch_all_like_meta_key_HASHED_FROM_meta_key_TO_meta_value($meta_key)};
      
      #now test failures
      #$ntable_a->fetch_like_meta_key_HASHED_FROM_meta_key_TO_meta_value($meta_key);
      #$ntable_a->fetch_all_like_meta_name_and_test_HASHED_FROM_meta_name_TO_meta_value($meta_key);
      #Both of these die nicely, but should probably throw?
      #die('should have failed by now');
            
      #Now do the top up
      $cmd = "perl $hive_script_dir/init_pipeline.pl Bio::EnsEMBL::Funcgen::Hive::Config::${conf} ".
        ' -analysis_topup '.$pipeline_params;  
      print "\n\nPERFORMING ANALYSIS TOPUP:\t".$conf."\n";
      
      #This will not catch non-fatal error output.
      run_system_cmd($cmd);
      
      #Reset can_run_AnalysisLogicName keys first, so we never assume that this has 
      #been done should things fail after adding the hive_conf key
      
      #We need to do these in order, with semaphored analyses reset last
      #As all link analyses are marked as DONE, when resetting a funnel job
      #before it's DONE fan jobs the beekeeper will look at the fan jobs 
      #and as they are all done will mark the funnel as READY insterad of SEMAPHORED
      #Even a 2nd attempt at resetting the funnel jobs will not work in this case
      #As they are seen as READY and effectively already reset :(
      #we may also need to run beekeeper -balance_semphores 
      #but this should be used with caution as it can go wrong.
      
      #Can we get the semaphore info from the DataFlowRuleAdaptor?
      
      
      #We need to pre-process these to order them such that the funnel jobs are reset last
      #We will need to bring back all the DataflowRules for this
      #as the funnel_dataflow_rule_id will likely be in a non-link analysis
      #let's do a direct SQL approach for this, rather than getting all the analyses
      
      #will -reset_all_bobs_for_analysis even reset to SEMAPHORED?
      #We could do this manually will an update and then -balance_semaphores
      
      
      my @can_run_keys;
      
      
      foreach my $can_run_key(keys %meta_key_values){
        (my $lname = $can_run_key) =~ s/^can_//o; 
        
        warn "push/unshift $lname, $can_run_key";
        
        if($dfr_adaptor->get_semaphoring_analysis_ids_by_logic_name($lname)){
          push @can_run_keys, $can_run_key;   
        }
        else{
          unshift @can_run_keys, $can_run_key;  
        }   
      }
      
      
      foreach my $can_run_key(@can_run_keys){
        
      
        if(scalar(@{$meta_key_values{$can_run_key}}) != 1){
          #this should not be possible as we are fetching a hash
          throw("Found multiple entries for meta_key $can_run_key:\t".join(' ', @{$meta_key_values{$can_run_key}}));  
        }
        
        #This is getting the old value! Not the new one!
        my $old_value = $meta_key_values{$can_run_key}->[0];
        my $new_value = $ntable_a->fetch_by_meta_key_TO_meta_value($can_run_key); 
        #warn "testing $can_run_key with values(new/old):\t".$new_value.' / '.$old_value;
        
        if($old_value){ #is defined and not 0
        
          if(! $new_value){
            die("Failed to process link analyses for $conf. $can_run_key meta_value has been reset from 1 to $new_value");  
          }
        }
        elsif($new_value){ #&& ! $old_value
       
          my $meta_id = $ntable_a->fetch_by_meta_key_TO_meta_id($can_run_key);            #PRIMARY KEY
          #$ntable_a->update_meta_value({meta_id=>$meta_id, meta_value =>$can_run_value}); #AUTOLOADED
       
          #Now reset the analysis if the value matches this config
          #i.e. we have just topped up with a downstream config, and want to reset and link
          #analyses which may have run.
          #$can_run_key will be double quoted as it is loaded from the config
          
          #Currently this will reset all DONE jobs. #This means that if we have 2 configs
          #which specify they're namsespace as the value, then truly DONE jobs will be reset
          #why was this change from 1 to config name?
          #The change from 0 to 1 enough here
          #Probably need to change this back?
          (my $can_run_analysis = $can_run_key) =~ s/^can_//o;
          
          #if($can_run_value eq "\"$conf\""){
            $cmd = "perl $hive_script_dir/beekeeper.pl -url $hive_url --reset_all_jobs_for_analysis $can_run_analysis";
            print "\nRESETTING LINK ANALYSIS JOBS FOR:\t$can_run_analysis\n";
            run_system_cmd($cmd);
            print "\nRESET LINK ANALYSIS JOBS FOR:\t$can_run_analysis\n";
          #}
        }
      }
      
      _register_conf_in_meta($ntable_a, $conf);    
    }
  }

}# end of main


#Add key via API to store with appropriate species_id i.e. 1
#shouldn't this be NULL? Probably but, core API only consideres the following meta_key type
#non-species specific: patch, schema_version, schema_type, ploidy
#Hive obviosuly does this with direct sql.

#This screws retrieval of the can_run_% meta values, as they are stored with species_id null
#but the BaseMetaContainer the expect a species ID when using any of the normal methods

#are these returned at all by the params?

#Work around with be to use direct mysql and do an update on them to reset the species ID to 1?
#Change this to use NakedTableAdaptor, as MetaContainer will disappear.

#Is this going to fail as don't we need the primary key defining?


sub _register_conf_in_meta{
  my $ntable_a   = shift;
  my $conf       = shift;
  
  eval { $ntable_a->store({meta_key => 'hive_conf', meta_value => $conf}); };
      
  if($@){
    throw("Failed to store hive conf meta entry:\t$conf\n$@");
  }
  
  return;  
}

sub populate_config_array{
  my ($conf_name, $confs) = @_;  
  
  #Internal validation of config hash
  if(! exists $config_info{$conf_name}){
    #Handle this in caller as this will also know the parent conf
    #No, need to handle here as this will recurse
    die("'$conf_name' is defined a pre-requisite, but is not defined in the config\nPlease amend config!"); 
  }
 
  #Must have at least a priority(can be undef) and a conf name
  if(scalar(@{$config_info{$conf_name}}) < 1){
    die("'$conf_name' config is not valid, must have at least a priority and/or a pre-requisite module."); 
  }
 
  my $priority = $config_info{$conf_name}->[0];
 
  if(defined $priority){ #We have a config rather than a label
  
    if(defined $confs->[$priority]  &&
      ($confs->[$priority] ne $conf_name)){
      die("Found a priority clash between configs:\t$conf_name and ".
            $confs->[$priority]."\nPlease correct config hash");    
    }
      
    $confs->[$priority] = $conf_name;
  }

  foreach my $conf_index(1..$#{$config_info{$conf_name}}){
    &populate_config_array($config_info{$conf_name}->[$conf_index], $confs);
  }
  
  return $confs;
}# end of populate_config_array


### TA DAA! ###

1;
