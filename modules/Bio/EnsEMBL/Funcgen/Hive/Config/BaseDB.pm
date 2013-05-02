
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


package Bio::EnsEMBL::Funcgen::Hive::Config::BaseDB;

use strict;
use warnings;
use Data::Dumper;
use base ('Bio::EnsEMBL::Funcgen::Hive::Config::BaseDB');
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
    %{$self->SUPER::default_options},
    
    dnadb   => 
     {  
      -host   => $self->o('dnadb_host'),
      -pass   => $self->o('dnadb_pass'),
      -port   => $self->o('dnadb_port'),
      -user   => $self->o('dnadb_user'),
      -dbname => $self->o('dnadb_name'),
     },
        
    efgdb  => 
     {
      -host   => $self->o('dbhost'),
      -port   => $self->o('dbport'),
      -user   => $self->o('dbuser'),
      -pass   => $self->o('dbpass'),
      -dbname => $self->o('dbname'),
     },
    
    
    #'registry_host'       => 'ens-livemirror',
	#'registry_port'       => 3306,
	#'registry_user'       => 'ensro',	
	#'port'              => 3306,
	disconnect_when_inactive => 1, #Set on funcgen and core DBAdaptors
   };
}


__END__

#Generic resource_classes are in Base



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
      %{$self->SUPER::pipeline_wide_parameters}, 
    
      #todo move these to default_options?
      #and remove pipeline_wide_params
            
      dnadb   => 
        {
         -host   => $self->o('dnadb_host'),
         -pass   => $self->o('dnadb_pass'),
         -port   => $self->o('dnadb_port'),
         -user   => $self->o('dnadb_user'),
         -dbname => $self->o('dnadb_name'),
        },
        
      efgdb  => 
        {
         -host   => $self->o('dbhost'),
         -port   => $self->o('dbport'),
         -user   => $self->o('dbuser'),
         -pass   => $self->o('dbpass'),
         -dbname => $self->o('dbname'),
        },
    };
}


#pipeline_create_commands in Base and/or subclass config 
#pipeline_analyses in subclass config


1;

