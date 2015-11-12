#!/usr/local/bin/bash

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
#      http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
# CONTACT
#
#  Please email comments or questions to the public Ensembl
#  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.
#
#  Questions may also be sent to the Ensembl help desk at
#  <http://www.ensembl.org/Help/Contact>.

### GENERAL NOTES ###
# 1 Never exit! 
#   Execute has been changed to return the exit code rather than exiting.
#   Always catch and return 1, so we don't don't the use out of the subshell
#   Some of the older funcs from funcs.sh still exit, so these need updating to return or subshelling to avoid exit
#   As such you really don't want lot's of nested Execute via dependant functions calls
#   To avoid having to handle this, the depth of the function call stack should never exceed 2
#   i.e. Don't call other functions unless it is a very simple utility function, which you know does
#   does not call other functions. Should probably rename all these utility functions as internal prefixed with _
#   This would also facilate automatic listing of the env functions by ignore functions prefixed with _

# 2 Can it be moved to a perl script or the hive?
#   Always consider whether something is better written in perl with a simple wrapper here if required
#   Or even better, put it in the pipeline!


### TODOs ###
# 
# There are likely more specific TODOs inline. Search for TODO.
#
# 1 If we want the basic non-ERSA specific funcs & aliases available for other pipelines
#   We will have tosplit this back into hive_env.sh and ersa_env.sh
#   This is not-stricly necessary as non-ersa pipelines are unlikely to require the interactive
#   review and seeding process. i.e. they probably have a single seed point which can be initiated
#   without the need for the configure & seed functions. Not using the env means you won't get 
#   aliases, GetFailedJobs, ListJobs, DebugJobs etc.

# 2 _InitEnv CheckVariables  ARCHIVE_DIR and others
#
# 3 Merge ListJobs and GetFailedJobs
#
# 4 Review funcs.sh Re/Move all unused func elsewhere or delete
#
# 5 Add reseed support to DebugJob, such that we can add things like recover/rollback to a debug job. 

# 6 pass DB passwords via environment for security, currently done on cmdline

# 7 Add the config summary to Help -l

# 8 Change functions to return 0 for the Check/Validate style methods
#   Then we can handle returning from here, instead of exit the shell
#   This may break the array mapping env. In the interim, add the funcs here
#   until we have time to update/retro fit them into the array mapping env
#   We can get around this by subshelling the call e.g. if [ ! $(ValidateVariablesOrUsage) ]; then return; fi
#   This won't work as the value tested here is the string output, not the exit code 
#   which needs to be test with $?
#   And actually supresses the text output, hence we won't see usage output
#   CheckDirs currently exits the subshell horribly as will others.
#   But this is no longer subshelled as the subshell was the old pre-req of having the efg.env 
#   initialsed, so exits the terminal!!
#   The old exit version of the funcs are still used in arrays.env
#   Simpy Create new ones with _!

# 9 Either source efg.env (tidy up first) or move efg.env aliases/funcs to separate file, such that we can source them
#   Typing efg in ersa curretnly sources the efg env (as the alias hasn't been redefined by efg.env)
#   rather than simply cd'ing to the repo.
#   This is required to bring the git env config in! i.e. git prompt stuff.


ERSA_USAGE='>. ersa_env.sh /path/to/instance/config.sh [dbpass] [dnadb_pass]'
base_name=$(basename -- "$0")

if [[ "$base_name" != "bash" ]]; then
  echo -e "$base_name must be sourced not executed e.g.\n\t$ERSA_USAGE"
  # exit
fi

# Take existing $EFG_SRC or default to relative location of this file
export EFG_SRC=${EFG_SRC:=$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd | sed 's#/scripts/environments##')}
export SRC=$(echo $EFG_SRC | sed 's#/ensembl-funcgen##')


### Validate config files exist

INSTANCE_CONFIG=$1
if [[ ! $INSTANCE_CONFIG ]]; then
  echo -e "You must pass at least an instance config file path argument:\n\t$ERSA_USAGE\nSee:\t$EFG_SRC/scripts/environments/ersa_instance_config.sh.example"
  return 1
elif [[ ! -f $INSTANCE_CONFIG ]]; then
  echo -e "The instance config file does not exist:\n\t$INSTANCE_CONFIG\nSee:\t$EFG_SRC/scripts/environments/ersa_instance_config.sh.example"
  return 1
fi

PIPELINE_CONFIG=${PIPLINE_CONFIG:=$EFG_SRC/scripts/environments/ersa_site_config.sh}
if [[ ! -f $PIPELINE_CONFIG ]]; then
  echo -e "You must create a 'site' config file:\n\t$PIPELINE_CONFIG\nSee:\t$EFG_SRC/scripts/environments/ersa_site_config.sh.example"
  return 1
fi



echo ":::: Setting up the ERSA hive environment ::::" 
. $EFG_SRC/scripts/environments/funcs.sh
echo -e ":: Sourcing site config:\t$PIPELINE_CONFIG"
. $PIPELINE_CONFIG
echo -e ":: Sourcing instance config:\t$INSTANCE_CONFIG"
. $INSTANCE_CONFIG

# Have to do this after the sourcing the $PIPELINE_CONFIG
# As this resets all the config vars
export DB_PASS=$2
export DNADB_PASS=$3
echo $DB_PASS $DNADB_PASS
shift
shift
shift
# Just in case we want to pass any other arguments, but this is already getting close to needing GetOpts


### Some convenient aliases
alias datahome='cd $DATA_HOME'
alias workdir='cd $WORK_DIR' 
alias aligndir='cd $DATA_HOME/alignments/$SPECIES'
# Need to add default assembly here? This would be set or fetched from the DB.
# alias indexhome?
alias mysqlefg='mysql $DB_MYSQL_ARGS'
alias mysqlout='mysql $DB_MYSQL_ARGS'
alias mysqlcore='mysql $DNADB_MYSQL_ARGS'
alias mysqldna='mysql $DNADB_MYSQL_ARGS'
alias mysqlpipe='mysql $PDB_MYSQL_ARGS'
alias mysqlhive='mysql $PDB_MYSQL_ARGS'




################################################################################
# Func      : _InitHiveEnv(
# Desc      : Sets up all config defined variables. Should only be called from
#             instance config file (species_VERSION.arrays) e.g. mouse_51.arrays 
# Args [n]  : 
# Return    : none
# Exception : 
################################################################################



#TODO
# 1 CheckVariables $ARCHIVE_DIR? This is actually optional? 
# but would need to  conditionally pass it in ConfigurePipeline

#2 CheckVariables DNADB_SCRIPT_ARGS (or make it optional in ConfigurePipeline)


# This is now actually a merge of generic pipeline and hive specific environments

_InitErsaEnv(){
	#Set all generic pipeline stuff e.g. DB vars, PATH/PERL5LIB etc
	#_InitPipelineEnv || \
  #    ( echo "Failed to initialise pipeline environment" && return 1 )
  echo $DB_PORT $DB_NAME $DB_HOST $SPECIES $DB_PASS $PIPELINE_PACKAGE $DATA_HOME
  _CheckGlobalVariables DB_PORT DB_NAME DB_HOST SPECIES DB_PASS PIPELINE_PACKAGE DATA_HOME || return 1
  _CheckDirs $DATA_HOME || return 1

  export EHIVE_ROOT_DIR=$SRC/ensembl-hive
  export LC_SPECIES=$(echo $SPECIES | tr [A-Z] [a-z])

  #SCHEMA_BUILD now optional
  #Can grab from DBNAME
  #Maybe useful if using non-standard efg DBNAME
  
  if [ ! $SCHEMA_BUILD ]; then
    local schema_build
    schema_build=$(GetSchemaBuild $DB_NAME)

    if [ $? -ne 0 ]; then
      echo -e $schema_build
      echo  'Failed to set $SCHEMA_BUILD. Please define this in your config file or change your $DB_NAME'
      return 1
    fi

    SCHEMA_BUILD=$(echo $schema_build | sed 's/ /_/')

  fi
    
  # Not required at present
  # VERSION=$(echo $SCHEMA_BUILD | sed 's/_.*//')
  # BUILD=$(echo $SCHEMA_BUILD | sed 's/.*_//')

  export QUEUE_MANAGER=${QUEUE_MANAGER:=LSF}
  # Only required for bsubJob from GenerateIndexes. Can be set to LOCAL
  ValidateGlobalVariable QUEUE_MANAGER VALID_QUEUE_MANAGERS 1 || return 1

  # This was originally part of an attempt to allow both ensembl-pipeline
  # and ensembl-hive support from the same base pipeline.env
  # However, it also prevents eternal growing of PATH/PERL5LIB
  # if we re-source the hive.env

  if [[ $PIPELINE_ENV != $PIPELINE_PACKAGE  ]]; then
      
      if [[ ! $PIPELINE_ENV ]];then
          export ORIG_PATH=$PATH
          export ORIG_PERL5LIB=$PERL5LIB
      fi

      PATH=${ORIG_PATH:=$PATH}
      PERL5LIB=${ORIG_PERL5LIB:=$PERL5LIB}
      export EFG_SCRIPTS=${SRC}/ensembl-funcgen/scripts
      export PIPELINE_SCRIPTS=${SRC}/${PIPELINE_PACKAGE}/scripts
      _CheckDirs $PIPELINE_SCRIPTS $EFG_SCRIPTS || return 1
      export PATH=${PIPELINE_SCRIPTS}:${EFG_SCRIPTS}:$PATH
      
      # Test for BioPerl here and suggest soft_link
      # TODO Validate these aren't already in $PERL5LIB and skip
      local d
      for d in $SRC/ensembl-funcgen/modules $SRC/current_bioperl $SRC/ensembl/modules $SRC/$PIPELINE_PACKAGE/modules; do
          _CheckDirs $d || return 1
          export PERL5LIB=$PERL5LIB:$d
      done

      export PIPELINE_ENV=$PIPELINE_PACKAGE
  fi

  # Pipeline DB config
  export PDB_USER=${PDB_USER:=$DB_USER}
  export PDB_PASS=${PDB_PASS:=$DB_PASS}
  export PDB_HOST=${PDB_HOST:=$DB_HOST}
  export PDB_NAME=${PDB_NAME:=${ENV_NAME}_${DB_NAME}}
  export PDB_PORT=${PDB_PORT:=$DB_PORT}

  # DNADB/Core DB config
  export DNADB_USER=${DNADB_USER:=$DB_USER}
  export DNADB_HOST=${DNADB_HOST:=$DB_HOST}
  export DNADB_PORT=${DNADB_PORT:=$DB_PORT}

  # MySQL client args string
  export DB_MYSQL_ARGS="-h${DB_HOST} -P${DB_PORT} -u${DB_USER} -p${DB_PASS} $DB_NAME"
  export PDB_MYSQL_ARGS="-h${PDB_HOST} -P${PDB_PORT} -u${PDB_USER} -p${PDB_PASS} $PDB_NAME"

  if [[ $DNADB_NAME ]]; then
    local mysqlargs_pass
    local args_pass
    #export REGISTRY_HOST=${REGISTRY_HOST:=$DNADB_HOST}
    #export REGISTRY_VERSION=${REGISTRY_VERSION:=$(echo $DNADB_NAME | sed -r 's/.*_([0-9]+)_[0-9]+$/\1/')}
    
    if [[ $DNADB_PASS ]]; then
      args_pass="-dnadb_pass $DNADB_PASS"
      mysqlargs_pass="-p $DNADB_PASS"
    fi
    
    export DNADB_SCRIPT_ARGS="-dnadb_host $DNADB_HOST -dnadb_user $DNADB_USER $args_pass -dnadb_name $DNADB_NAME -dnadb_port $DNADB_PORT"
    export DNADB_MYSQL_ARGS="-h${DNADB_HOST} -P${DNADB_PORT} -u${DNADB_USER} $mysqlargs_pass $DNADB_NAME"
  fi

  # Sript args
  export PDB_SCRIPT_ARGS="-dbhost $PDB_HOST -dbuser $PDB_USER -dbpass $PDB_PASS -dbname $PDB_NAME -dbport $PDB_PORT"
  export DB_SCRIPT_ARGS="-host $DB_HOST -user $DB_USER -pass $DB_PASS -dbname $DB_NAME -port $DB_PORT"
  echo ""

  # TODO  Improve this a little with some directoried
  # Convert it to eval at runtime, such that it can be included in Help -l

  echo "DB:               ${DB_USER}@${DB_HOST}:${DB_NAME}:${DB_PORT}
DNADB:            ${DNADB_USER}@${DNADB_HOST}:${DNADB_NAME}:${DNADB_PORT}
PIPELINEDB:       ${PDB_USER}@${PDB_HOST}:${PDB_NAME}:${PDB_PORT}"

  #Will need this for array support?
	#if [ $MULTI_SPECIES ]; then
#
#		if [ $MULTI_SPECIES -ne 1 ]; then
#		echo 'MULTI_SPECIES is a boolean variable, please specify 1 or omit for default of 0'
#		else
#		#change to script paramter
#			export MULTI_SPECIES=' -multi_species '
#		fi
#	fi
	
  # Define directories 
	export WORK_DIR=${DATA_HOME}/output/${DB_NAME}
	export BACKUP_DIR=${WORK_DIR}/backup
	export PIPELINE_OUT=${WORK_DIR}/hive_debug

	MakeDirs $WORK_DIR

  #This is required as the config uses this to identify default config locations
  export ENSEMBL_CVS_ROOT_DIR=$SRC
  export HIVE_URL="mysql://${PDB_USER}:${PDB_PASS}@${PDB_HOST}:${PDB_PORT}/${PDB_NAME}"
  Execute perl $EFG_SCRIPTS/pipeline/add_hive_url_to_meta.pl $DB_SCRIPT_ARGS -url $HIVE_URL || return 1

  alias beekeeper="$PIPELINE_SCRIPTS/beekeeper.pl -url $HIVE_URL -hive_log_dir ${WORK_DIR}/debug"

  export PS1="\
\[\033[${PS1_COLOUR}m\]\
${ENV_NAME}:${DB_NAME}>\
\[\033[0m\]"
  Help
}


# TODO Change to perl script wrapper. No see ENSREGULATION-253

GenerateIndexes(){
	local usage
  local assembly
  local genders
  local force
  local valid_genders
  local genders
  valid_genders="male female"
  usage="Usage: GenerateHeaders  -a(ssembly) e.g. GRCh38 [ -g(ender) e.g. male|female (omit to generate both) -f(orce over-write, will also re-link fasta input) ]"

  ### PROCESS OPTIONS ###
	OPTIND=1
	while getopts ":a:g:fh" opt; do
		case $opt in 
          a  ) assembly=$OPTARG ;;
          g  ) genders=$OPTARG ;;
          f  ) force=1 ;;
          h  ) PrintColour -c blue -s "$usage"; return 0;;
          \? ) PrintColour -c blue -s "$usage"; return 1;;
    esac 
	done
  

  if [[ $genders ]]; then
      _ValidateVariable -v $gender -V "$valid_genders" -n gender || echo $usage && return 1
  else
      genders=$valid_genders
  fi

  if [[ ! $assembly ]]; then echo -e $usage && return 1; fi

  local index_root
  local fasta_root
  local sam_root
  local gender
  local f_fasta
  local i_fast
  local job_name
  local bsub_params
  local cmd  
  local s_faidx
  local f_faidx
  local s_header
  index_root=$DATA_HOME/bwa_indexes/$LC_SPECIES/
  fasta_root=$DATA_HOME/fasta/$LC_SPECIES/
  sam_root=$DATA_HOME/sam_header/$LC_SPECIES/

  ### GENERATE INDEXES FOR EACH GENDER ###
  for gender in $genders; do
      echo "Generating indexes for $gender"
      f_fasta="$fasta_root/${LC_SPECIES}_${gender}_${assembly}_unmasked.fasta"
         
    #Assume we already have indexes if rsa file is present as
    #this is last to be created
      if [[ -f "$index_root/${LC_SPECIES}_${gender}_${assembly}_unmasked.fasta.rsa"  && 
            ! $force ]]; then
          echo "Found previously existing BWA index files:"
          ls $index_root/${LC_SPECIES}_${gender}_${assembly}_unmasked.fasta.*
          echo "Specify -f(orce) flag to overwrite and re-link the fasta file"
      else
          i_fasta="$index_root/${LC_SPECIES}_${gender}_${assembly}_unmasked.fasta"
          
          if [[ ! -d $index_root ]]; then
              echo -e "Creating absent BWA index directory:\t$index_root";
              Execute "mkdir -p $index_root" || return 1
          fi   
          
          if [[ $force ||  ! -f $i_fasta ]]; then
              
              if [[ ! -f $f_fasta ]]; then
                  echo -e "Could not find fasta file:\t$f_fasta"
                  return 1
              fi
              
              #check index dir exists
              Execute "ln -fs $f_fasta $i_fasta" || return 1
          fi
          
          
          job_name="BWA_IDX-${LC_SPECIES}_${gender}_${assembly}_unmasked.fasta" 
          #4GB to be safe
          bsub_params='-M4000 -R"select[mem>4000] rusage[mem=4000]"'
          bsub_params=$bsub_params" -o ${WORK_DIR}/GenerateIndexes.${gender}.out -e ${WORK_DIR}/GenerateIndexes.${gender}.err"
          #-I'
          #interactive mode appears to hang when launched from here, 
          #but it is probably just that the STDOUT is not being 
          #being printed to the screen for some reason.
          #works fine on cmdline, but reverted to normal bsub as these jobs take ages.

          cmd="cd $index_root && bwa index -a bwtsw $i_fasta"
          #omission of -a btwsw causes seg fault!
          echo -e "Submitting bsub job:\t$bsub_params $cmd"
          bsubJob $job_name "$bsub_params" "$cmd" || return 1
        
          #This does not submit an interactive job!
          #bsub $bsub_params $cmd
          
          #samtools faidx some.fasta
          s_faidx="$sam_root/${LC_SPECIES}_${gender}_${assembly}_unmasked.fasta.fai"
          f_faidx="$fasta_root/${LC_SPECIES}_${gender}_${assembly}_unmasked.fasta.fai"
          
          if [[ -f $s_faidx && ! $force ]]; then
              echo -e "Found previously existing SAM index file:\n$s_faidx"
          else
              echo -e "Executing:\tsamtools faidx $f_fasta"
              Execute "samtools faidx $f_fasta" || return 1
              Execute "mv $f_faidx $s_faidx" || return 1

              #We actaully don't use the headers any more as 
              #the headers are already in the files and they can actaully be inserted 
              #with samtools view (although this is less efficient)
              #we can't generate the sam header here from samtools as it does not allow
              #viewing of the header only from a fai file (but it could!)
              #samtools view -t $fai_file -H
              #However, it's fairly easy to generate these from the fai
              s_header="$sam_root/${LC_SPECIES}_${gender}_${assembly}_unmasked.header.sam"
              echo -e "Creating sam header file:\n$s_header"
              awk '{print "@SQ\t" "SN:" $1 "\tLN:" $2;}' $s_faidx > $s_header
          fi          
      fi
  done
}


# These two are sequencing specific

# WE need to add similar functions to ConfigureMotifPipeline
# ConfigureArrayPipeline

ConfigurePipeline(){
  # Add the tracking procedures which maybe absent
  # This will fail if the tracking tables are not present
  echo ":: Loading tracking_procedures.sql"
  Execute mysql $DB_MYSQL_ARGS < $EFG_SRC/sql/tracking_procedures.sql || return 1
  # This is slightly redundant if we configure more than once, could put this in _InitErsaEnv
  # but that would prevent using the env with a non-tracking DB

  perl $EFG_SCRIPTS/pipeline/configure_hive.pl -pdb_user $PDB_USER -pdb_pass $PDB_PASS -PDB_NAME $PDB_NAME \
    -pdb_host $PDB_HOST -pdb_port $PDB_PORT -user $DB_USER -pass $DB_PASS -host $DB_HOST -port $DB_PORT -dbname $DB_NAME \
    $DNADB_SCRIPT_ARGS -hive_script_dir $PIPELINE_SCRIPTS -species $LC_SPECIES -data_root $DATA_HOME \
    -archive_root $ARCHIVE_DIR $@
  PrintColour -c blue -s "Type 'Help' for a list of aliases and commands"
}


SeedPipeline(){
  perl $EFG_SCRIPTS/pipeline/seed_hive.pl -url $HIVE_URL -hive_script_dir $PIPELINE_SCRIPTS $@
  PrintColour -c blue -s "Type 'Help' for a list of aliases and commands"
}


ReseedJobs(){
  perl $EFG_SCRIPTS/pipeline/reseed_jobs.pl -url $HIVE_URL $@
  PrintColour -c blue -s "Type 'Help' for a list of aliases and commands";
}



#add -can_respecialize 1 here? This might be on by default?

Beekeeper(){
  perl $PIPELINE_SCRIPTS/beekeeper.pl -url $HIVE_URL  -submit_workers_max 100 $@
#-hive_log_dir ${WORK_DIR}/debug
#-sleep 0.5
}

Beekeeper_with_logs(){
  perl $PIPELINE_SCRIPTS/beekeeper.pl -url $HIVE_URL -hive_log_dir ${WORK_DIR}/debug  -submit_workers_max 100 $@
}


AnalysisProgress(){
  mysqlpipe -e 'select * from progress';
}

# Currently this is list control experiments as unprocessed, as this is determined by the presence
# of a correponding result set. Need to update this to check for ALIGNED_CONTROL status instead
# Need to be able to filter this with additional sql

SegmentationProgress(){
  mysqlefg -e 'call SummariseSegBuild' > /dev/null;
  local cmd
  cmd="select * from progress_view $@"
  # Nesting $@ directly below causes interpolation issues
  mysqlefg -e"$cmd"
}

# TODO Add ReguBuildProgress (and maybe wrappers to views here too?)


SetSegmentationFeatureTypes(){
  local in_string

  for ft in "$@"; do
    in_string="$in_string '$ft',"
  done
  in_string=$(echo $in_string | sed 's/,$//')

  # Also need to validate feature types actually exist

  # Now insert ignore and then update just to force it in for now
  # todo incorporate check and force over-write behaviour  
  # Or keep is simply and just print out the old ones, just in case we want to reset them?
  mysqlefg -e "INSERT IGNORE INTO regbuild_string SELECT NULL, 'segmentation.feature_type_ids', 1,\
    group_concat(feature_type_id) from feature_type where name in($in_string)"

  mysqlefg -e "update regbuild_string set string=(select group_concat(feature_type_id) from \
    feature_type where name in($in_string)) where name='segmentation.feature_type_ids'"  
}


GetSegmentationFeatureTypes(){
  local ftype_ids
  ftype_ids=$(mysqlefg --skip-column-names -e "SELECT string from regbuild_string where name='segmentation.feature_type_ids'")
  mysqlefg --skip-column-names -e"select group_concat(name) from feature_type where feature_type_id in($ftype_ids)"
}


SetRegBuildCellTypes(){
  local in_string
  
  for ct in "$@"; do
    in_string="$in_string '$ct',"
  done
  in_string=$(echo $in_string | sed 's/,$//')

  # Also need to validate cell types actually exist

  # Now insert ignore and then update just to force it in for now
  # todo incorporate check and force over-write behaviour  
  # Or keep is simply and just print out the old ones, just in case we want to reset them?
  mysqlefg -e "INSERT IGNORE INTO regbuild_string SELECT NULL, 'regbuild.cell_type_ids', 1,\
    group_concat(cell_type_id) from cell_type where name in($in_string)"

  mysqlefg -e "update regbuild_string set string=(select group_concat(cell_type_id) from \
    cell_type where name in($in_string)) where name='regbuild.cell_type_ids'"  
}

GetRegBuildCellTypes(){
  local ctype_ids
  ctype_ids=$(mysqlefg --skip-column-names -e "SELECT string from regbuild_string where name='regbuild.cell_type_ids'")
  mysqlefg --skip-column-names -e "select group_concat(name) from cell_type where cell_type_id in($ctype_ids)"
}



# GetFailedJobs

# TODO Merge with ListJobs
#This is still getting retrying jobs, as FAILED can be assigned
#temporarily until the beekeeper detects whether it can be run again 
#These aren't actively being retried. But they are still unsafe
#can we omit them anyway without knowing th emax retry count permissable for that analysis?
#or are these fail_no_retry job?

#default behaviour now mirrors beekeeper, but with richer output
#although it will still return non-active retrying jobs see above

#todo add opt to list the worker.err too

#Having direct sql here instead of using the beekeeper allows us to run this on the farm
#when debugging large memory jobs

#TODO Add 1 line mode, which compressesjob_id line and throw MSG line onto one line
#GetFailedJobs -l WriteSignalCollection | grep -E '^(job|MSG)' | awk '{if(/^job_id/){ line = $0; }else{ print line, $0, "\n";}}'
#This is useful for harvesting jobs which have the same error, for resetting/retring or whatever.
# and even pass an optional grep string to filter on? e.g.
#GetFailedJobs -l WriteSignalCollection | grep -E '^(job|MSG)' | awk '{if(/^job_id/){ line = $0; }else{ print line, $0, "\n";}}' | grep ccat | sed -r 's/job_id=([0-9]+).*/\1/'
#and even specify whenther your wnt the prev_job_id?


GetFailedJobs(){
  #Declare vars local first!  
	local usage
  local limit_clause
  local lname_clause
  local input_id
  local status_clause
  local lm_join

  usage="Usage: GetFailedJobs  [ -l(ogic_name) <String>  -i(nclude input_id) -r(etrying jobs too) -L(imit) <Int> -h(elp))]"
  status_clause='="FAILED"' 
  lm_join='LEFT JOIN'
  
	OPTIND=1
	while getopts ":l:rihL:" opt; do
		case $opt in 
          l  ) lname_clause=" and a.logic_name='$OPTARG' " ;;
          r  ) status_clause='!="DONE"'; lm_join='JOIN' ;;
          i  ) input_id='concat("\nInput ID:\t", j.input_id),' ;;
          L  ) limit_clause=" limit $OPTARG" ;;
          h  ) PrintColour -c blue -s "$usage"; return 0;;
          \? ) PrintColour -c blue -s "$usage"; return 1;;
    esac 
	done

  local sql
  sql='SELECT concat("\n\njob_id=", j.job_id,"(", j.prev_job_id,") "), concat(a.logic_name, "(", a.analysis_id, ") "), concat("retry=", lm.retry, "(", j.retry_count, ")"), j.status, '$input_id' "\nErrorFile:\t", jf.stderr_file, concat("\n", lm.msg) 
FROM analysis_base a 
JOIN job j using(analysis_id)'$lm_join'(select lm2.job_id, max(lm2.retry) as retry, lm3.time from log_message lm2 
      JOIN (select job_id, max(time) as time from log_message where is_error=1 and job_id is not NULL group by job_id) lm3 \
      USING(job_id, time) 
GROUP BY lm2.job_id) l using (job_id) 
LEFT JOIN log_message lm using (job_id, time, retry) 
LEFT JOIN job_file jf on lm.job_id=jf.job_id and lm.worker_id=jf.worker_id and lm.retry=jf.retry 
WHERE j.status '$status_clause' '$lname_clause' order by lm.time '$limit_clause;

  #echoing through subshelled command removes all mysql table padding characters
  #but it does flatten every record onto 1 line, so we have to add formating
  #This is flattening the field header too, which we can't easily add, so  --skip-column-names 
  echo -e $(mysqlpipe --skip-column-names -e "$sql") 
}


# TODO Merge this with GetFailedJobs

ListJobs(){
  local job_id
  local sql
  job_id=$1

  if [[ $job_id ]]; then
    job_id=" and j.job_id=$job_id"
  fi

   sql='select concat("\n", j.job_id, "(", j.prev_job_id, ")", a.logic_name) as Job, j.input_id, j.worker_id, j.status, j.retry_count, j.completed from job j JOIN analysis_base a using (analysis_id) WHERE j.input_id not like "\_extended\_data\_id %"'$job_id'
UNION
select concat("\n", j.job_id, "(", j.prev_job_id, ")", a.logic_name) as Job, ad.data as input_id, j.worker_id, j.status, j.retry_count, j.completed from job j JOIN analysis_base a using (analysis_id) join analysis_data ad on ad.analysis_data_id=replace(j.input_id, "_extended_data_id ", "") WHERE j.input_id like "\_extended\_data\_id %"'$job_id
  echo "SQL is $sql"
  echo -e $(mysqlpipe -e "$sql")  
}




#Is the solution here to add functionality to reset the job locally then submit the interactive bsub job?
#this could only launch the env, as this will not run any further commands
#unless, can we take these as args to the script?
#This does work, but we'd have to send it through efg first, and pass through the path of this script to execute
#This would be a neat way to handle the bsub interactive Debug jobs

#


# It would be nice if we could -reset_job to a DEBUG_READY state, which can only be taken by 
# a worker in debug mode. This could not force reset if the job is already claimed.

#processing of getopts order reflects the order of the case statement,
#not as they are specified on the cmd line
#which means this shifting in these case statements doesn't work to
#leave the reamning arguments for passing onto runWorker.pl

#todo pick up resource / meadow from farm, and submit accordinly
#this is the use case for being able to bsub -Is directly

#todo delete all jobs with job_id as prev_job_id, else we may get erroneous jobs hanging around from the
#previous run, which will either keep failing or mess up the data if they are reset by accident.
#only do this if no_write is not set.
#as no_write is where the dataflow occurs, but we will still have normal output.

#todo add reseed support!

DebugJob(){
  local job_id
  local debug_level
  local force
  local no_write
	local usage
  usage="Usage: DebugJob -j(ob id) <Int> [ -d(ebug level) <Int> -f(orce DONE/SEMAPHORED jobs) -h(elp) -n(o_write, prevents deleting depedant jobs and shouldn't dataflow) ]"
  debug_level=1

	OPTIND=1  
	while getopts ":j:d:fn" opt; do
		case $opt in 
        j  ) job_id=$OPTARG  ;;
        d  ) debug_level=$OPTARG ;;
        f  ) force='-force 1' ;;
        n  ) no_write='-no_write' ;;
        h  ) PrintColour -c blue -s "$usage"; return 0;;
        \? ) PrintColour -c blue -s "$usage"; return 1;;
    esac 
	done

  #Removes the opts and opt args so we can process the path arguments
  shift $(($OPTIND - 1))
    
  local status
  local cmd

  _CheckVariablesOrUsage "-j(ob_id) not specified\n$usage" job_id || return 1
  status=$(mysqlpipe --skip-column-names -e "select status from job where job_id='$job_id'")

  echo "Job $job_id status is $status";

  if [[ $status = FAILED ]]; then 
      force='-force 1'
  elif [[ $status = DONE || $status = SEMAPHORED ]]; then
      
      #-z
      if [[ ! $force ]]; then
          echo "Job $job_id is the $status state. Please specify -f(orce) to over-ride"
          return 1
          #The Queen would actually output this sort of error anyway
      else
          echo "### WARNING: Forcing running of worker for $status job $job_id ###"
      fi
  else 
      #it is READY or some sort of run state, so let the worker handle that.
      #Just reset force to null in case it was specified as an opt
      force=

      #This could actually be a DebugJob which was killed and hence cannot be restarted
      #without a reset
      #the beekeeper -sync will not rectify this
      #but beekeeper -loop/run
      
      #persistant run state and so will require a beekeeper.pl -reset_job
  fi


  if [[ ! $no_write ]]; then
      #Do this so we don't have old potentially incorrect downstream jobs

      echo "Deleteing dependant jobs with prev_job_id=$job_id"
      mysqlpipe -e "DELETE from job where prev_job_id='$job_id'" || \
          ( echo "Failed to delete downstream jobs for job $job_id" && return 1 )
  fi

  #don't do this, as it may result in resetting a running job
  #which will then cause a race condition
  #perl $PIPELINE_SCRIPTS/beekeeper.pl -url $HIVE_URL -reset_job_id $job_id -job_id $job_id -debug $debug_level -meadow_type LOCAL $@
  
  cmd="time perl $PIPELINE_SCRIPTS/runWorker.pl -url $HIVE_URL -job_id $job_id \
-debug $debug_level $force $no_write $*"
  echo $cmd
  Execute $cmd
}


#Only for analysis, as for individual jobs, DebugJob does this
#will we ever want to delete a batch of job IDs but not the whole analysis?
#do not take analysis_id here as it is too risky

#Doing a simple job delete like this can be dangerous as the status entries will become
#incorrect e.g.
#When deleting all jobs downstream of Preprocess_bwa_samse_control we could be left
#with experiments which have the ALIGNED_CONTROL status
#This will need reverting to ALIGNING_CONTROL, otherwise subsequent depedant jobs
#could be submited resulting in a clash as they should really wait fo rthe control jobs
#to finish



DeleteDownstreamJobs(){
 #Declare vars local first!  
	local usage=
  local lname=
  local sql=
  local reset=

  usage="Usage: DeleteDownstreamJobs -l(ogic_name) <String> [ -r(eset analysis) -h(elp) ]"
  
  #This makes sure we reset the getopts ind if we have used it previously
	OPTIND=1

	while getopts ":l:rh" opt; do
		case $opt in 
          l  ) lname=$OPTARG ;;
          r  ) reset=1 ;;
          h  ) PrintColour -c blue -s "$usage"; return 0;;
          \? ) PrintColour -c blue -s "$usage"; return 1;;
    esac 
	done

  AskQuestion "Before commencing delete, have you considered the relevant status entries which may need patching/removing" "(y|c)"

  if [[ $REPLY == c ]]; then
      echo -e "Cancelling DeleteDownstreamJobs for $lname";
      return
  fi

  sql='DELETE j from analysis_base ab JOIN job j1 ON ab.analysis_id=j1.analysis_id '
  sql+="JOIN job j ON j1.job_id=j.prev_job_id WHERE ab.logic_name='$lname'"
  
  echo -e "Deleting jobs downstream of $lname using:\n$sql";
  mysqlpipe -e "$sql" || return 1


  if [[ $reset ]]; then
      echo "Resetting $lname jobs"
      Beekeeper -reset_all_jobs_for_analysis $lname
  else
      #Just do a sync so we can see that we have been successful
      Beekeeper -sync
  fi
}

GenerateGraph(){
  local format
  local no_display
  local usage
  format=pnge
  usage="GenerateGraph
  Description: Calls the generate_graph.pl script. Default ouput is:
                 ${WORK_DIR}/hive_diagram.${format}
  Usage:       GenerateGraph [ -f(ormat default=$format) -d(isplay, automatically opens graph in firefox) -h(elp prints this message) ]"

  OPTIND=1
  while getopts ":hf:d" opt; do
      case $opt in 
          h  ) PrintColour -c blue -s "$usage"; return 0;;
          f  ) format=$OPTARG ;;
          d  ) display=1 ;;  
          \? ) PrintColour -c blue -s "$usage"; return 1;;
        esac 
  done

  Execute perl $PIPELINE_SCRIPTS/generate_graph.pl -url $HIVE_URL -output ${WORK_DIR}/hive_diagram.${format}

  if [[ $display ]]; then
      Execute firefox ${WORK_DIR}/hive_diagram.${format} &
      #this will likely fail if you are using screen
  fi  

  echo -e "Your pipeline diagram has been written to:\n\t${WORK_DIR}/hive_diagram.${format}"  
}



GenerateLSFReport(){
  local usage
  usage="GenerateGraph
  Description: Calls the lsf_report.pl script. Default ouput is:
                 ${WORK_DIR}/lsf_report.txt
               Also see the lsf_report and lsf_usage tables for more info.
  Usage:       GenerateLSFReport [ -h(elp prints this message) ]"

  OPTIND=1
  while getopts ":h" opt; do
      case $opt in 
          h  ) PrintColour -c blue -s "$usage"; return 0;;
          \? ) PrintColour -c blue -s "$usage"; return 1;;
        esac 
  done

  Execute perl $PIPELINE_SCRIPTS/lsf_report.pl -url $HIVE_URL -output ${WORK_DIR}/lsf_report.txt
}

# TODO Will need updating
# Include eval'd version of config summary

Help(){
  #Declare vars local first!
  local long=
	local usage=
    
  #Take optional opts here to change regex i.e. include all functions

	usage='Help
  Description: Prints help documentation for this pipeline
  Usage:       Help [ -l(ists all functions and aliases) -h(elp prints this message) ]'

	OPTIND=1
	while getopts ":hl" opt; do
		case $opt in 
			h  ) PrintColour -c blue -s "$usage"; return 0;;
      l  ) long=1 ;;
      \? ) PrintColour -c blue -s "$usage"; return 1;;
    esac 
  done

  #if [[ $long ]]; then
  #    EFGHelp
  #fi


  PrintColour -c blue -s "

#################################################### 
###  Ensembl Regulation Analysis Pipeline Help  ####
####################################################

Useful aliases:

  mysqlpipe - Access the hive DB via mysql client
  mysqlefg  - Access the output funcgen DB via mysql client
  mysqlcore - Access the core DB via mysql client
  datahome  - cd to \$DATA_HOME i.e. $DATA_HOME
  workdir   - cd to the working directory i.e. $WORK_DIR
  aligndir  - cd to the alignments directory for this species i.e. $DATA_HOME/alignments/$SPECIES

Useful functions (use -h for help where specified):

  ConfigurePipeline    [optional init_pipeline.pl args] - Initialises and adds 'Configs' to the hive pipeline DB
  GenerateGraph        [-h] - Generates the pipeline graph
  GenerateLSFReport    [-h] - Generates lsf report
  SeedPipeline         [-h] - Seeds the pipeline with jobs or see what would be seeded
  ReseedPipeline       [optional reseed_pipeline.pl args]

  Beekeeper            [optional beekeeper.pl args] - Simply calls beekeeper.pl -url \$HIVE_URL 
  Beekeeper_with_logs  [optional beekeeper.pl args] - As above but with:\t-hive_log_dir \$WORK_DIR/debug
  DebugJob             [-h]
  DeleteDownstreamJobs [-h]

  AnalysisProgress             - Selects * from hive progress table
  SegBuildProgress             - Calls the SummariseSegBuild procedure. Useful selecting what to seed.
  ListJobs          [ job_id ] - Lists jobs (performing necessary join to analysis_data)
  GetFailedJobs     [ -h ] - Lists failing jobs
   
  Help          [-l(ong)]  - Prints this message, the -l(ong) option includes full list of aliases and functions
"

}

# DropPipeline      [-h] - Drops the hive pipeline tables or DB
# Left this in pipeline.env for now. Just do this manually, also safer.
# Funcs to incorporate/remove





# MOTIF PIPELINE METHODS #
# TODO Ensure that these work as is
# Update this pipeline to work as a fully formed hive
# Incorporate into same hive DB as the sequencing pipeline


# Dump peps using ensembl-analysis/scripts/protein/dump_translations.pl -dbhost $DNADB_HOST -dbuser $DNADB_USER -dbname $DNADB_NAME -file $EFG_DATA/fasta/$SPECIES/${DNADB_NAME}.pep.fasta
# format using formatdb -i $EFG_DATA/fasta/$SPECIES/${DNADB_NAME}.pep.fasta

#$EFG_SRC/scripts/pwm_mappings/load_jaspar_matrices.pl --user $DB_USER --host $DB_HOST --dbname $DB_NAME --pass $DB_PASS \
# -dnadb_user $DNADB_USER -dnadb_host $DNADB_HOST -dnadb_name $DNADB_NAME \
# -jdb_user $DB_USER -jdb_host $JDB_HOST --out_dir ~/tmp_jaspar -pep_fasta $EFG_DATA/fasta/$SPECIES/${DNADB_NAME}.pep.fasta  \
# -pfm_file ~/data/scratch109/binding_matrices/Jaspar_5.0/JASPAR_CORE/pfm/nonredundant/pfm_all.txt -pfm_file ~/data/scratch109/binding_matrices/Jaspar_5.0/JASPAR_CORE/pfm/redundant/pfm_all.txt -pfm_file ~/data/scratch109/binding_matrices/JASPAR2010/all_data/matrix_only/matrix.txt -pfm_file ~/data/scratch109/binding_matrices/Jaspar_pre_2010/all_data/matrix_only/matrix_only.txt


RunMotifMatch(){
  echo ":: RunMotifMatch $*"

  usage='usage: RunMotifMatch [-f(eature type) e.g. Max]+ [-h help]'

  folder=${DATA_DIR}/binding_matrices/Jaspar
  feature_list=
  OPTIND=1

  while getopts ":f:h" opt; do
      feature=
      case $opt in 
                f  ) feature=$OPTARG;;
          \? ) echo $usage; return 1;;
          h  ) echo $usage; return 0;;
       esac 
       
       if [[ $feature ]]; then
     feature_list="$feature_list $feature"
       fi
     
  done    
  
  map_dir=$WORK_DIR/pwm_mappings
  Execute "mkdir -p $map_dir"

  if [[ $feature_list ]]; then
      feature_list="-feature_type_list $feature_list"
  fi
       
  #At the moment this needs to be run in a specific folder. TODO update the scripts
  cd $map_dir

  _CheckGlobalVariables SCHEMA_BUILD
  

  #todo make this use $DB_READ_SCRIPT_ARGS

  #TODO extra job management...
  bsub -J pwm_mappings_${DB_NAME} -e $map_dir/pwm_mappings_${DB_NAME}.err -o $map_dir/pwm_mappings_${DB_NAME}.out -q long -M35000 -R"select[mem>35000] rusage[mem=35000]" perl $EFG_SRC/scripts/pwm_mappings/run_pwm_mapping_pipeline.pl $DB_SCRIPT_ARGS $DNADB_SCRIPT_ARGS -species $LC_SPECIES -assembly $ASSEMBLY -workdir $DATA_HOME -outputdir $map_dir $feature_list
  #echo $cmd
  #Execute  $cmd 

  echo "Now wait for the job to finish and check logs and output"

}


ImportMotifMatch(){
  echo ":: ImportMotifMatch $*"

  usage="usage: ImportMotifMatch [-s(slice) e.g. 22]+ [-h help]"
  #TODO do some tests when pipeline already exists... may need a DropPipeline before...

  slice_list=
  OPTIND=1

  while getopts ":s:h" opt; do
      slice=
      case $opt in 
                s  ) slice=$OPTARG;;
          \? ) echo $usage; return 1;;
          h  ) echo $usage; return 0;;
       esac 
       
       if [[ $slice ]]; then
         slice_list="$slice_list $slice"
       fi
     
  done    

  # TODO Make this use the same PDB as the reset
        #Make this more generic...
  PDB_NAME=mf_${DB_NAME}

 #todo make this use $MYSQL_READ_ARGS


  exists=$(echo "show databases like '$PDB_NAME'" | mysql $MYSQL_ARGS)
  if [[ $exists ]]; then
      echo "Pipeline DB $PDB_NAME already exists: Need to start from scratch. Use DropPipelineDB motif_import"
      echo "Aborting"
      return 0
  fi
  
  if [[ $slice_list ]]; then
      slice_list="-slices $slice_list"
  fi

  map_dir="$WORK_DIR/pwm_mappings/filtered"
    
  cmd="perl $EFG_SRC/scripts/pwm_mappings/run_binding_site_import_pipeline.pl $DB_SCRIPT_ARGS $DNADB_SCRIPT_ARGS -workdir $map_dir -output_dir $WORK_DIR $slice_list"
  echo $cmd
  Execute $cmd 

  cmd="beekeeper.pl -url mysql://${DB_USER}:${DB_PASS}@${DB_HOST}:${DB_PORT}/${PDB_NAME} -sync"
  Execute $cmd

  #Cannot use Execute in case we want to interrupt it
  echo "beekeeper.pl -url mysql://${DB_USER}:${DB_PASS}@${DB_HOST}:${DB_PORT}/${PDB_NAME} -loop"
  beekeeper.pl -url mysql://${DB_USER}:${DB_PASS}@${DB_HOST}:${DB_PORT}/${USER}_motif_import_${DB_NAME} -loop

  echo "Now check results in $WORK_DIR/motif_features/results"

}


### NEW FUNCS.SH METHODS ###
# These redefine those sourced from func.sh via the pipeline.env
# at least until we can update the array mapping env
# although these may break dependant functions in func.sh which are not over-ridden here
# do we even need funcs.sh for hive.env?
# Will these be enough?


# TODO
# Remove OrUsage functions, this should be easily/more appropriately handled in the caller


################################################################################
# Func      : _ValidateVariable()
# Desc      : Check that the passed variable is contained in pass valid variable
#             array
# Arg [1]   : Variable to validate 
# Arg [n+1] : List of valid variable values 
# Return    : 
# Exception : Exits if variable is not present in valid variable list
################################################################################

#warn if valid_var only contains 1 value?

_ValidateVariable(){
  #declare vars first before setting as they are only declared local after they are defined  
  local var_name= 
	local var=
	local valid_vars=
	local valid=
  local usage=
  usage="Usage:\tValidateVariable -v(ar to test) <String> -V(alid vars) \"<String ...>\" [ -n(ame of var/argument/parameter ]"
	
  OPTIND=1
	while getopts ":v:V:n:hl" opt; do
		case $opt in 
      v  ) var=$OPTARG ;;
      V  ) valid_vars=$OPTARG ;;
      n  ) var_name=$OPT_ARG ;;
			h  ) PrintColour -c blue -s "$usage"; return 0;;
      \? ) PrintColour -c blue -s "$usage"; return 1;;
    esac 
  done
  #Removes the opts and opt args so we can process the path arguments
  shift $(($OPTIND - 1))
    

  if [[ $* ]]; then
    echo "Found more arguments than expected number of parameters, perhaps you forgot to quote the -V(alid vars)?"
    echo -e "var = $var\tvalid_vars = $valid_vars"
    return 1   
  fi

  if [[ ! $var || ! $valid_vars ]]; then
      echo "Mandatory parameters not set:"
      echo -e "var = $var\tvalid_vars = $valid_vars\n$usage"
      return 1
  fi

  for valid_var in $valid_vars; do

    if [[ $valid_var = $var ]]; then 
      valid=1
    fi
  done

  if [[ ! $valid ]]; then
    echo -e "$var_name $var is not a valid, must be one of:\t$valid_vars"
    return 1
  fi
}


CheckFile(){
  if [[ ! -f $1 ]]; then
    echo "Error : file $1 does not exist"
    return 1
  fi
}

#todo make this bit shift to get the true exit code?

Execute(){
  #declare var first before setting as they are only declared local after they are defined  
  local rtn=
  $*
  rtn=$?

	if [[ $rtn != 0 ]]
	then 
		echo -e "Failed to:\t$*\nReturned exit code:\t$rtn"
		return $rtn 
	fi
}


bsubJob(){
  local job_name=
	local bsub_params=
	local job_cmd=
	local eval_flag=
  local usage=

	job_name=$1
	bsub_params=$2
	job_cmd=$3
	eval_flag=$4
  usage='Usage:\tbsubJob $job_name "$bsub_params" "$job_cmd" $eval_flag'

  if [[ ! $job_name  || ! $bsub_params  || ! $job_cmd ]]; then
      echo  "Mandatory arguments not set:"
      echo -e "\tjob_name = $job_name\n\tbsub_param = $bsub_params\n\tjob_cmd = $job_cmd"
      echo -e $usage
      return 1
  fi

	#Need to change this to take job_name bsub_params and job args

	if [[ "$QUEUE_MANAGER" = 'LSF' ]]; then

		checkJob $job_name
    #Sets global $JOB_ID based on $job_name
    # This can still exit if it fails to get the info 

	#We should test for more than one job here?
	#jid will currently catch all ids

	
		if [[ $JOB_ID ]]; then
			echo "Job($JOB_ID) $job_name already exists"
		else
		  #echo through bash to avoid weird resources parameter truncation

			if [[ $eval_flag ]]; then
				#We can get another level of interpolation here by evaling
				#But this may reveal more env vars
				#Trick is to not double quote things in the first place
				job_cmd=$(eval "echo \"$job_cmd\"")
			fi

			echo bsub -J $job_name $bsub_params $job_cmd

			#quote job cmd to avoid ; breaking bsub			
			#quote job name to handle job arrays
			JOB_ID=$(echo bsub -J "\"$job_name\"" $bsub_params "'$job_cmd'" | bash)
			
			#This is not catching some failures!!
      echo $JOB_ID

			if [[ $? -ne 0 ]]; then
				echo "Failed to submit job $job_name"
				return 1
			fi
							
			JOB_ID=$(echo $JOB_ID | sed 's/Job <//')
			JOB_ID=$(echo $JOB_ID | sed 's/>.*//')
		fi

	  #To let LSF process the job
		sleep 5
		echo "To monitor job: bjobs -lJ $job_name"
		bjobs -w $JOB_ID
		
	elif [[ "$QUEUE_MANAGER" = 'Local' ]]; then
		echo -e "Running Local job:\t$job_cmd"
		Execute $job_cmd

	else
		echo -e "bsubJob does not support QUEUE_MANAGER:\t$QUEUE_MANAGER"
		return 1
	fi
}


_InitErsaEnv
