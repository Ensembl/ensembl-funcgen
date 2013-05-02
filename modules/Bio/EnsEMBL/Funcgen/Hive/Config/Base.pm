
=pod

=head1 NAME

    Bio::EnsEMBL::Funcgen::HiveConfig::Alignment_conf;

=head1 SYNOPSIS

   # TODO: this could be easily merged with the Peaks pipeline...
   # TODO: allow subfolders which will represent replicates...
   # Allow semaphores so jobs can be run truly in parallel (see SemaStart and SemaLongMult_conf)

   # Example 1: specifying only the mandatory options (initial params are taken from defaults)
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Alignment_conf -password <mypass>

   # Example 2: specifying the mandatory options as well as setting initial params:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Alignment_conf -password <mypass> -p1name p1value -p2name p2value

   # Example 3: do not re-create the database, just load more tasks into an existing one:
init_pipeline.pl Bio::EnsEMBL::Funcgen::HiveConfig::Alignment_conf -job_topup -password <mypass> -p1name p1value -p2name p2value


=head1 DESCRIPTION

    This is the Config file for the Alignment Pipeline

    Please refer to Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf module to understand the interface implemented here.

    The Alignment pipeline consists of several "analysis":
        * SetupAlignmentPipeline verifies the existence of the files and creates alignment jobs ...
        * RunAlignment makes the alignment...
        * WrapUpAlignment merges the alignments, some QC and fills in the data tracking db

    Please see the implementation details in Runnable modules themselves.

=head1 CONTACT

    Please contact ensembl-dev@ebi.ac.uk mailing list with questions/suggestions.

=cut


package Bio::EnsEMBL::Funcgen::Hive::Config::Base;

use strict;
use warnings;
use Data::Dumper;
use base ('Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf');
# All Hive databases configuration files should inherit from HiveGeneric, directly or indirectly


=head2 default_options

    Description : Implements default_options() interface method of 
                  Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that is used
                  to initialize default options.
    Caller      : DependantOptions::process_options (init_pipeline.pl)

=cut

#These are used by init_pipeline.pl and according to the docs should the pipelinedb spec
#Anything set (e.g. test) in default_options will be accessible in other methods here 
#(e.g pipeline_wide_params, pipeline_create_commands) as $self->o(test), 
#unless specified on the cmdline in which case that will be the value of $self->o(test);
#This is basically mechanism to remove the needs for doing soemthing like this f
#or defaults in pipeline_wide_params: test => $self->o('test') || default_value
#or: test => (defined $self->o('test')) ? $self->o('test').'/ing' : '/default/ing'
#$self->o also captures undefined cmdline params, if they have not been set in defaults

#There seems to be no clear reason for not setting params in here, if we never want to set 
#a hardcoded default for them which we might want to over-write with a specific cmdlind opt
#in pipeline_wide_parameters

#It looks like it is possible to refer to a previously defined option within the same hash

sub default_options {
  my ($self) = @_;  
  
  return {
    #%{$self->SUPER::default_options},#Don't do this for Base
    #as we don't want any generic defs in there including the 
    #default ensembl_cvs_root_dir    
    ensembl_cvs_root_dir => $ENV{'SRC'},    
    bin_dir => '/software/ensembl/funcgen', 
 
    pipeline_name => $self->o('pipeline_name'),#Used for job names
    pipeline_db => 
	  {
	   -host   => $self->o('host'),
	   -port   => $self->o('port'),
	   -user   => $self->o('user'),
	   -pass   => $self->o('pass'),
	   -dbname => $self->o('pipeline_name'),
	   #todo deal with this in the env and don't corrupt pipeline_name add $ENV{USER}?
	  },
    
    #This could be inferred from the db, but it's probably safer(?) to pass as parameter...
    species      => $self->o('species'),
    #May pass this to input_id... to allow for files of different assemblies in the same pipeline run.
    assembly     => $self->o('assembly'),
    #add defaults for values for dir structure here
   
    data_dir   => $self->o('data_dir'), #'/lustre/scratch109/ensembl/funcgen',
    #The following should all be able to be redefined by the name sake cmdline opts
    #in pipeline_wide_params
   
   
    
	#These can't be modified by each conf, as they will be in the same DB things will get messy
	#todo These should be used as root dirs and specific confs should define new keys i.e. 'confname_data_dir'


	root_output_dir => $self->o('data_dir').'/output/',
	#need root_output_dir for non-DB specific analyses e.g. alignments
	db_output_dir   => $self->o('data_dir').'/output/'.$self->o('dbname'),
	hive_output_dir => $self->o('data_dir').'/output/'.$self->o('pipeline_name').'/hive_debug',
	use_tracking_db => 1,	
   };
}

=head2 resource_classes

    Description : Implements resource_classes() interface method of
      Bio::EnsEMBL::Hive::PipeConfig::HiveGeneric_conf that lists the LSF resource classes available

=cut


sub resource_classes {
  my ($self) = @_;
  return 
    {
     #Use this section when running on Sanger Farm
     default                 => { 'LSF' => '' },
     urgent                  => { 'LSF' => '-q yesterday' },
     normal_monitored        => { 'LSF' => " -R\"select[$ENV{LSF_RESOURCE_HOST}<1000] ".
                                            "rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1]\"" },
     normal_high_mem         => { 'LSF' => ' -M5000000 -R"select[mem>5000] rusage[mem=5000]"' },
     long_monitored          => { 'LSF' => "-q long -R\"select[$ENV{LSF_RESOURCE_HOST}<1000] ".
                                            "rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1]\"" },
     long_high_mem           => { 'LSF' => '-q long -M4000000 -R"select[mem>4000] rusage[mem=4000]"' },
     long_monitored_high_mem => { 'LSF' => "-q long -M4000000 -R\"select[$ENV{LSF_RESOURCE_HOST}<1000 && mem>4000]".
                                                " rusage[$ENV{LSF_RESOURCE_HOST}=10:duration=10:decay=1,mem=4000]\"" },

     #Use this section when running on EBI cluster???
     #    0 => { -desc => 'default',          'LSF' => '' },
     #    1 => { -desc => 'long_high_mem',      'LSF' => '-M5000 -R"select[mem>5000] rusage[mem=5000]"' },
     #    2 => { -desc => 'normal_high_memory',    'LSF' => '-M5000 -R"select[mem>5000] rusage[mem=5000]"' },
    };
}


=head2 pipeline_wide_parameters

    Description : Interface method that should return a hash of pipeline_wide_parameter_name->pipeline_wide_parameter_value pairs.
                  The value doesn't have to be a scalar, can be any Perl structure now (will be stringified and de-stringified automagically).
                  Please see existing PipeConfig modules for examples.
    Caller      : HiveGeneric_conf::run (init_pipeline.pl)

=cut

#$self->o here will take the cmdline option else default to what was set in default_options

#Will init_pipeline with analysis_topup, use previously stored opts so we don't have to define
#them again?


sub pipeline_wide_parameters {
    my ($self) = @_;
                            
    return {
      %{$self->SUPER::pipeline_wide_parameters},  # inheriting database and hive tables creation
    
      #Allow over-writing defaults set above
      bin_dir         => $self->o('bin_dir'),
      root_output_dir => $self->o('root_output_dir'),
	  db_output_dir   => $self->o('db_output_dir'),
  	  hive_output_dir => $self->o('hive_output_dir'),
  	  use_tracking_db => $self->o('use_tracking_db'),	  
    };
}



#default pipeline_create_commands defined in HiveGeneric_conf.pm
#specify more in sub-class conf if required but always do @{$self->SUPER::pipeline_create_commands},



1;

