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

$| = 1;#for debug

my (@configs, $species, $data_root, $hive_script_dir, $list);

#put this config in HiveHelper?
my %config_info = (
                   #ConfigModule => [priority, abbreviation, pre-req],
                   #priorities should be unique!
                   #although Peaks & Collections are effectively equal priority here
                   #abbrevions are for use in the pipeline name/meta key
                   'DefineSets'  => [1], #A pre-req should not have a pre-req of it's own
                   'Peaks'       => [2, 'DefineSets'],
                   'Collections' => [3, 'DefineSets'],
                  );

#For now we are not storing the conf keys in the output DB
#but we may want to move them there, if we are to allow
#more than one hive to run on the same DB.


my @tmp_args = @ARGV;
my $db_opts = get_DB_options_config();#This will get opts for funcgen, core and pipeline by default

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


my @valid_confs = sort {$config_info{$a}->[0] cmp $config_info{$b}->[0]} keys %config_info;

if($list){
  print "Valid config modules are:\n\t".join("\n\t", @valid_confs)."\n";
  #This currently doesn't list those which are pre-reqs
  exit;  
}


if(! defined $hive_script_dir){
  throw('-hive_script_dir is a mandatory parameter');  
}      
elsif(! -d $hive_script_dir){
  throw("-hive_script_dir is not a valid directory:\t${hive_script_dir}");
}

my $db_script_args = process_DB_options($db_opts, ['funcgen', 'core', 'pipeline'], undef, 'script');
my $pdb_params     = ${&process_DB_options($db_opts, ['pipeline'])}{pipeline};
my $db             = create_Funcgen_DBAdaptor_from_options($db_opts, 'pass');

# Validate configs
if(! @configs){
  pod2usage(-exitval => 1, 
            -message => "Must provide at least one config module name. Valid modules are:\n\t".
                        join("\n\t", @valid_confs));  
}


my (@confs, %abbrvs);

foreach my $conf(@configs){
  
  if(! exists $config_info{$conf}){
    pod2usage(-exitval => 1, 
              -message => "'$conf' is not a supported config module. Valid modules are:\n\t".
                          join("\n\t", @valid_confs));  
  }  

  #todo validate Bio::EnsEMBL::Funcgen::Hive::Config::$conf exists

  #Populate conf array in order or priority
  $confs[$config_info{$conf}->[0]] = $conf;
  
  #Deal with pre-req 
   
  #foreach my $pre_req(1..$#{$config_info{$conf}}){
  #Now only supports 1
  
  if(scalar(@{$config_info{$conf}}) == 2){
    
    #assumes pre-reqs don't have pre-reqs
    #if(scalar(@{$config_info{$pre_req}}) != 1){
    #  die('Config pre-requisites cannot have pre-requisites themselves. Please check %config_info.');  
    #}
    
    $confs[$config_info{$config_info{$conf}->[1]}->[0]] = $config_info{$conf}->[1];  
  }
  
  #if(scalar(@{$config_info{$conf}}) == 2){
  #  $abbrvs{$config_info{$conf}->[0]} = $config_info{$conf}->[1]; 
  #}
}


#my $pipe_abbrv = join('_', (map { $abbrvs{$_} } sort keys %abbrvs));

#Now we have a pre-ordered (possibly gappy) array of confs to add
#Strip out any undef elements
my @tmp_confs;
map {push @tmp_confs, $_ if defined $_} @confs;
@confs = @tmp_confs;

my $pipeline_params = "-data_root_dir $data_root -pipeline_name ".$pdb_params->{'-dbname'}.' '.
  $db_script_args->{funcgen}.' '.$db_script_args->{core}.' '.$db_script_args->{pipeline};
$pipeline_params .= " -species $species " if defined $species;  


my $pdb;  
eval { $pdb = create_DBAdaptor_from_params($pdb_params, 'core', 1); };

if($@){ #Assume the DB hasn't been created yet
  
  #init the pipline with the first conf
  my $first_conf = shift @confs;
  my $init_cmd = "perl $hive_script_dir/init_pipeline.pl Bio::EnsEMBL::Funcgen::Hive::Config::${first_conf} ".
    $pipeline_params;
  
  print "\n\nINITIALISING DATABASE:\t".$pdb_params->{'-dbname'}."\n";
  run_system_cmd($init_cmd);  
  
  $pdb = create_DBAdaptor_from_params($pdb_params, 'core');
}

add_hive_url_to_meta(url_from_DB_params($pdb_params), $db);

#Now iterate through the confs using -analysis_top_up
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

1;
