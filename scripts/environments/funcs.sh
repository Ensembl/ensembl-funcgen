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



#echo ":::: Sourcing funcs.sh"
export PROMPT_COMMAND='echo -ne "\033]0;${USER}@${HOSTNAME}: ${PWD}\007"'

#To do
#Add more error echos when exiting
#Document error codes
#rename funcs lcfirst or _ first to enable us to filter public and 'private' methods in help functions

#We should add these to an array and then use the array values in ReportExit to list a nice error rather than a number
# File : funcs.sh
# Error Codes: #restrict these to 100s so we can use 200s in any environment which sources these?
# 205 - Incorrect user
# 100 - Variable not defined
# 201 - Binary not found
# 101 - Variable not valid
# 255 - Execution of script failed


# LIST OF RESERVED VARIABLES
# These may be set in one function to be used in the caller
# i.e. do not use/set these unless using the function which sets them
# Also check func.sh and efg.env for similar reserved variables e.g REPLY
# JOB_ID


_TRUE=1;
_FALSE=!$_



_setOptArgArray(){
    local array_var=
    local args=
    local nextarg=
   	local OptArgArray=

	#Set var in sub
	#Subshelling would cause problems with setting OPTIND 
	array_var=$1
	shift
	#Should this be $@ to avoid spliting 'words with' with spaces?
	args=($*)   
	nextarg=0


	#Need to check args here!
	#can I ever be less than 0?

	CheckVariables array_var
	#This is largely pointless as if it is not passed we just suck in the args
	#$i will be -1 if called outside of getopts?
	#Actaully Could be anything dependant on the last value of OPTIND
	#Is also susceptible of polluting of OPTIND if not reset to 0 before calling getopts

	i=$(($OPTIND - 2))

	while [[ $nextarg = 0 ]]; do
	  #echo "i is $i"
		#echo "arg is ${args[$i]}"
		#if [[ ("${args[$i]}" = -*) || (-z "${args[$i]}") ]]; then	
		if [[ ( ${args[$i]} = -* ) || ( ! ${args[$i]} ) ]]; then
			nextarg=1
			OPTIND=$(($i + 1))
			#echo "next arg is at $OPTIND"
				#do we need to set this to -1?
		else
			#echo "using arg ${args[$i]}"
			OptArgArray="$OptArgArray ${args[$i]}"		
		fi
		
		i=$(($i + 1))
	done

	eval "$array_var=($OptArgArray)"	
}



################################################################################
# Func      : jobWait()
# Desc      : Waits for lsf job
# Arg  [1]  : Job ID
# Arg  [2]  : Amount of time to wait before polling in seconds. Default is 280(5mins)
#             Remove this and allow a list of Job IDs?
# Return    : Exit code from lsf job
# Exception : Exits if JobID not defined
################################################################################

jobWait(){
  local job_id=
  local secs=
  local complete=
  local exited=
  local status=

  job_id=$1
	secs=$2
	secs=${secs:=300}
	CheckVariables job_id

	complete=
	exited=

	#Test for valid job first?
	echo -n "Checking for job every $secs secs"

	while [[ ! $complete ]]; do
		job_info=($(bjobs $job_id))
		status=${job_info[10]}
		#We probably need to catch lsf daemon not responding here
		#as well as bjobs failure

		if [[ $status = 'EXIT' ]]; then		
			#Get exit code
			complete=$(bjobs -l $job_id | grep 'Exited with exit code')
			complete=$(echo $complete | sed 's/.*code //')
			complete=$(echo $complete | sed 's/\. The CPU.*//')
			
		elif [[ $status = 'DONE' ]]; then
			complete=0
		else
			echo -n '.'
			sleep $secs
		fi
	done

	echo -e "Finished waiting for job:\t$job_id"
	return $complete
}



checkJob(){
  local job_name=
  local exit=
	JOB_ID=

	job_name=$1
	exit=$2

	CheckVariables job_name


	echo -e "Checking for job:\t$job_name"

	
	if [[ "$QUEUE_MANAGER" = 'LSF' ]]; then

	#This does not catch Job <gallus_gallus_transcripts.55_2m.fasta> is not found
		JOB_ID=$(bjobs -J $job_name | grep -e "^[0-9]" | sed -r 's/([0-9]+)[[:space:]].*/\1/')

		if [ $? -ne 0 ]; then
			echo "Failed to access job information"
			exit 1
		elif [ $JOB_ID ]; then
			echo "$job_name is still running"
		fi
		

	elif [[ "$QUEUE_MANAGER" = 'Local' ]]; then
		echo -e "Skipping job check for Local job:\t$job_name"
	else
		echo -e "funcs.sh does not support QUEUE_MANAGER:\t$QUEUE_MANAGER"
		exit
	fi


}

submitJob(){

	#Change this to take opts so we can wait for jobs?

	job_name=$1
	shift
	bsub_params=$1
	shift
	job_cmd=$1
	shift
	eval_flag=$1

	#default to LSF
	QUEUE_MANAGER=${QUEUE_MANAGER:=LSF}

	CheckVariables job_name bsub_params job_cmd


	#Need to change this to take job_name bsub_params and job args

	if [[ "$QUEUE_MANAGER" = 'LSF' ]]; then

		checkJob $job_name

	#We should test for more than one job here?
	#jid will currently catch all ids

	
		if [ $JOB_ID ]; then
			echo "Job($JOB_ID) $job_name already exists"
		else
		#echo through bash to avoid weird resources parameter truncation

			if [ $eval_flag ]; then
				#We can get another level of interpolation here by evaling
				#But this may reveal more env vars
				#Trick is to not double quote things in the first place
				job_cmd=$(eval "echo \"$job_cmd\"")
			fi

			echo bsub -J $job_name $bsub_params $job_cmd

			#exit;
			#quote job cmd to avoid ; breaking bsub			
			#quote job name to handle job arrays
			JOB_ID=$(echo bsub -J "\"$job_name\"" $bsub_params "'$job_cmd'" | bash)
			

			#This is not catching some failures!!

			if [ $? -ne 0 ]; then
				echo $JOB_ID
				echo "Failed to submit job $job_name"
				exit 1
			fi
			
			echo $JOB_ID
			
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
		echo -e "funcs.sh does not support QUEUE_MANAGER:\t$QUEUE_MANAGER"
		exit
	fi



}





################################################################################
# Func      : CheckGlobalVariables()
# Desc      : Checks that all passed environment variable names have defined values
# Args [n]  : List of variable names to check 
# Returns   : True is all are defined
#             False if string evaluated variable is of length zero
################################################################################

#This is really low utility now, apart from in the _InitEnv
#functions, and also in functions which require optional global vars

CheckGlobalVariables(){
  local line=

  if [[ ! $* ]]; then
		echo 'Must pass variable name(s) to check'
		return 1
	fi

  local var=

  for var in $*; do
    val=$(eval "echo \$$var")
    
    if [[ -z $val ]];  then 
      line="$line \$$var"
    fi
  done
  
  if [[ ! -z $line ]];  then
    echo "Environment variable(s) :$line are not set"
    return 1
  fi

}







################################################################################
# Func      : CheckVariables()
# Desc      : check that all passed environment variables have been assigned 
#             values
# Args [n]  : list of variables to check 
# Return    : #list of variables that cannot be found
# Exception : Exits if string evaluated variable is of length zero
################################################################################


# REMOVE THIS IN FAVOUR OF CheckGlobalVariables

CheckVariables(){
	line=

	if [ ! "$*" ]; then
		echo 'Must pass variable name(s) to check'
		exit 100
	fi


    for var in $*
    do
        val=$(eval "echo \$$var")

        if [ -z "$val" ]
		then 
			line="$line \$$var"
		fi
    done

    if [ ! -z "$line" ]
    then
        echo "Environment variable(s) :$line are not set"
		#This is not exiting if it is called from another environment
		exit 100
    fi
}

#No point in these OrUsage funcs
#Best to do this in the caller using a subshell error=$(func)
#As we are going to have to capture the error anyway!
#This will also stop the error getting flattened

CheckVariablesOrUsage(){
	usage=$1
	shift
	variable_names="$*"

	#Could CheckVariabl

	tmp=$(CheckGlobalVariables $variable_names)

	if [ $? != 0 ]; then
		echo -e "$tmp\n$usage"
		#This get's flattened into one line if we capture the output for returning rather than exit
		#So we don't get full error
		exit 1;
	fi
}

ValidateVariableOrUsage(){
	usage=$1
	shift

	CheckGlobalVariables usage

	tmp=$(ValidateGlobalVariable $*)

	if [ $? != 0 ]; then
		echo "$tmp"
		echo "$usage"
		exit 1;
	fi
}


################################################################################
# Func      : ValidateGlobalVariable()
# Desc      : Check that the passed variable is contained in pass valid variable
#             array
# Arg [1]   : Variable to validate 
# Arg [n+1] : List of valid variable values 
# Return    : 
# Exception : Exits if variable is not present in valid variable list
################################################################################

#Can specify no_exit here
#But probably better just to subshell, then we catch all the errors echoed

# eg. error=$(ValidateVariablesOrUsage "$usage" type valid_types)
# Then deal with error in context

ValidateGlobalVariable(){
	var_name=$1
	valid_vars_name=$2
	no_exit=$3
	#

	valid=

	CheckVariables $var_name $valid_vars_name

	var_value=$(eval "echo \$$var_name")
	valid_var_values=$(eval "echo \$$valid_vars_name")


	#Could we grep here instead?

    for valid_var_value in $valid_var_values
    do
        if [ $valid_var_value = $var_value ]
		then 
			valid=1
		fi
    done

    if [ ! $valid ]
    then
        echo "$var_value is not a valid $var e.g. $valid_var_values"

		if [ $no_exit ]; then
			return 101
		else
			exit 101
		fi
    fi
	
	return

}


################################################################################
# Func      : CheckBinaries() 
# Desc      : checks that binary files exist
# Args [n]  : list of binaries to check 
# Return    : list of binaries that cannot be found
# Exception : none
################################################################################

CheckBinaries(){
    for var in $*
    do
        which $var | grep ^no > /dev/null
        
        if [ $? -eq 0 ]; 
		then 
			#where is $line coming from?
			line="$line '$var'"
		fi
    done
    
    if [ ! -z "$line" ]
    then 
        echo "Binary(s) : $line cannot be found"
		exit 201
    fi
}

################################################################################
# Func      : CheckDirs()
# Desc      : checks to see whether passed parameter is a directory
# Arg [1]   : list of directory paths
# Return    : returns true if it is a directory
# Exception : exists if not a directory
################################################################################

CheckDirs(){

	for dir_name in $*
    do
		if [ ! -d $dir_name ]
		then
			echo "error : directory $dir_name does not exist"
			exit 203
		fi
	done
}



################################################################################
# Func      : MakeDirs() 
# Desc      : if a given directory does not exist then it will make it 
# Args [n]  : list of directory(s) to make 
# Return    : none
# Exception : exits if 'mkdir' fails 
################################################################################

MakeDirs(){
     for dir_name in $*
     do
         if [ ! -d $dir_name ]
         then
			mkdir -p $dir_name;
			rtn=$?; 

			if [ $rtn -ne 0 ]; 
			then 
				exit $rtn; 
			fi;
         fi
     done
}


################################################################################
# Func      : BackUpFile() 
# Desc      : Backs up file to $BACKUP_DIR if set
# Args [n]  : File name(not path, so have to back up from file directory)
# Return    : none
# Exception : 
################################################################################

BackUpFile(){

	#Can we add some file path parsing here to enable back up of file paths

	file=$(GetFilename $1)
	file_dir=$(GetDir $1)

	backup_dir=${BACKUP_DIR:=$file_dir}
	MakeDirs $backup_dir

	if [[ -f $1 ]]; 
	then 
		Execute cp $1 ${backup_dir}/${file}.$$
	fi
}

################################################################################
# Func      : ArchiveData() 
# Desc      : rsyncs(-essh -Wav) directory of files to archive dir. Mirrors 
#             directory structure between $DATA_DIR AND $ARCHIVE_DIR
# Args [n]  : List of files and/or directories
# Return    : none
# Exception : Returns if $ARCHIVE_DIR or $DATA_DIR not set
################################################################################

#Not handling hosts in these paths yet, so -essh is redundant at present
#Could do with lsa func, which lists filedir in the archive?
#Change this to ArchiveFileDir
#Could add $HOME support, but shouldn't have archiveable data in their



ArchiveData(){
	OPTIND=1

	compress=
	delete_source=
	aname=
	usage='usage: ArchiveData [ -c(compress, no name) -a(rchive name, compress) -d(elete source) -h(elp) ] FILES|DIRS'
	#add -l(ist) option here to see if we have already archived this in some form?
	#or we could just test for target file first and AskQuestion to confirm?
	#would need force flag to overwrite without AskQuestion

	while getopts ":dca:h" opt; do
		case $opt in 
	        d  ) delete_source=1 ;; 
            c  ) compress=1 ;;
            a  ) aname=$OPTARG ;;
            h  ) echo $usage; return 0;;
            \? ) echo $usage; exit 1;;
		esac 
	done

	if [ $aname ]; then
		compress=1
	fi

	i=1
	
	while [ $i -lt $OPTIND ]; do
		shift
		let i+=1
	done
	
	filedirs=$*
	TARGET_ROOT=
	TARGET_ROOT_NAME=
	SOURCE_ROOT=

	for filedir in $filedirs; do
    	#Get the full dereferenced path
		#Also strips trailing /
		#Need to capture readlink error here?
		
		_SetTargetAndSourceRoot $filedir
		#sets SOURCE_ROOT and return $derefd_filedir
   		retval=$?
		
		if [ $retval -ne 0 ]; then
			echo -e "Failed to archive:\t$filedir"
			return $retval
		fi
	
		filedir=$derefd_filedir
		
				
		#Need to match $DATA_DIR here or skip with warning
		if [[ ! -d $TARGET_ROOT ]]; then
			echo -e "Target archive dir not set or valid:$TARGET_ROOT_NAME\t=\t$TARGET_ROOT"
			return 1
		else
	        #Use ' instead of / due to / in path
			archive_filedir=$(echo $filedir | sed -r "s'^${SOURCE_ROOT}''")
			archive_filedir="$TARGET_ROOT/$archive_filedir"	
			
			#mkdir in archive incase we are using tar or rsync using a file
			#rsync will create dirs
			archive_dir=$(GetDir $archive_filedir)
			
			if [ ! -d $archive_dir ]; then
				echo -e "Making archive directory:\t$archive_dir"
				mkdir -p $archive_dir
			fi
			
			
			if [ -d $filedir ] && [ ! $compress ]; then
		    	#Remove the last dir name				
				archive_filedir=$(echo $archive_filedir | sed -r 's/\/$//')
				archive_filedir=$(echo $archive_filedir | sed -r 's/[^/]+$//')
			fi
			
			
			archive_cmd=
			
			if [ $compress ]; then
				
				if [ $aname ]; then
					aname=".${aname}"
				fi	
				#This overwrites existing files
				archive_cmd="tar -cvzf $archive_filedir${aname}.tar.gz $filedir"
			else
				archive_cmd="rsync -essh -Wavm $filedir $archive_filedir/"
	       		#-essh only necessary for remote archiving
	       		#-a archive mode; equals -rlptgoD (no -H,-A,-X)
	    		#-r recurse into directories
	       		#-l copy symlinks as symlinks
		     	#-p preserve permissions
			    #-t preserve modification times
    			#-g preserve group
    			#-o preserve owner (super-user only)
    			#-D same as:
	    		#    --devices     preserve device files (super-user only)
                #    --specials    preserve special files
                #-m prune emtpy dirs
                #-v verbose
			fi	
			
			echo $archive_cmd
			Execute $archive_cmd
			#assign Execute output here to catch exit and output, then return nicely?
			echo -e "Finished:\t$archive_cmd"
			
			if [[ $delete_source ]]; then
				echo -e "Removing source:\t $filedir"
				rm -rf $filedir
			fi
			
		fi
	done
}

#Need to test for aliases before defining these
#Was failing to compile as rm was already aliased

#rm(){
#	args=$*
#	seen_file=
#	i=0

#	files=$(echo $args | sed -r s'/^-[^[:space:]]+//') #deal with leading - i.e. not preceded
#	files=$(echo $files | sed -r s'/ -[^[:space:]]+//g') #deal with other opts
	#This will not restrict opts to start of args

	#echo files $files
#	opts=$(echo $args | sed "s'$files''");
# 	opts=$(echo $opts | sed "s'-'-o '");
 	   
	#echo del $opts $files
#	del -r $opts $files

	#Currently this is del'ing dirs from del roots
	#even if -r isn't specified
#}

#Enables fast removal of files by moving to .del folder in root of current path
#Post-pones actual rm to cron job, based on files last mod'd more than N days

#Was failing to compile as del was already aliased

#Probably want to hoik out delete > age code to separate func
#so we can use it separately (for logs)


#
#del(){
#	OPTIND=1
#	days=
#	del_verbose=
#	rm_opts=
#	rm_caller=
#	
#	cmd_line="del $*"
#	usage='usage: del  [ -o(pt for rm)+ -d(ays, purge .del of files older than this value, at root defined by) ] FILES|DIRS'
#
#	while getopts ":d:o:vrh" opt; do
#		case $opt in 
#	        d  ) days=$OPTARG ;;
#            o  ) rm_opts="$rm_opts -${OPTARG}" ;;
#			r  ) rm_caller=1 ;;
#			v  ) del_verbose=1 ;;
#			h  ) echo $usage; return 0;;
#			\? ) echo -r "$cmd_line\n$usage"; exit 1;;
#		esac 
#	done
#  
#	#Assume rm is okay if we have defined opt for rm
#	if [ $rm_opts ]; then
#		rm_caller=1
#	fi
#
#	i=1
#	while [ $i -lt $OPTIND ]; do
#		shift
#		let i+=1
#	done
#
#	filedirs=$*
#
#	if [ ! $filedirs ]; then
#		echo -e "You must define at least one file or directory\n$usage";
#	fi
#
#
#	#test $days is +ve int to avoid -gt test failure later
#	if [ $days ] &&
#		! [[ $days =~ ^[0-9]+$ ]]; then
#		echo -e "Parameter -d must be an integer\t$usage"
#		return 1
#	fi
#
#	#Build error log rather than bailing out asap
#	error_log=
#
#	for filedir in $filedirs; do
#		_SetTargetAndSourceRoot -n $filedir
#		#sets SOURCE_ROOT and derefd_filedir
#		retval=$?
#
#	
#		if [ ! $derefd_filedir ]; then
#			error_log="${error_log}File/directory does not exist:\t$filedir\n"
#			continue
#		fi
#
#	
#		if [ $retval -ne 0 ]; then
#				
#			if [ $rm_caller ]; then
#				#Only rm if we are calling from rm func
#			    #Not from del directly
#			    #-i is over-ridden by -f
#				rm_cmd="$(which rm) -i $rm_opts $filedir"
#
#				if [ $del_verbose ]; then
#					echo $rm_cmd;
#				fi
#
#				$rm_cmd
#
#				#capture $? here
#				#Need to capture error message too for summary report
#
#			else
#				error_log="${error_log}Failed to del as no .del dir available for:\t$filedir\nUse rm instead?\n"
#			fi
#
#			#Or enable a .del in /nfs home too?
#			#This would not work with _SetTargetAndSourceRoot
#			continue
#		fi
#
#		filedir=$derefd_filedir
#		del_dir="${SOURCE_ROOT}/.del"
#
#		if [ ! -d $del_dir ]; then
#			mkdir -p $del_dir
#	
#			if [ $? -ne 0 ]; then
#				echo -e "Failed create .del dir:\t$del_dir"
#				return $retval
#			fi
#		fi
#
#
#		if [ $days ]; then # PURGE!
#		
#			#Test valid .del $filedir to purge
#			#No need to match trailing as readlink strips this
#			if [[ $filedir != $del_dir ]]; then
#				echo -e "You must supply a valid .del dir to purge e.g.\n\t$del_dir\n$usage"
#				return 1
#			fi 
#
#			#Use find instead of ls to allow for many files
#			#-mindepth ignores $del_dir dir itself
#			#-depth processes dir contents before dir itself
#			#to prevent deleting dir before finding the next file which has already been deleted
#			#However, this means we are doing many deletes rather than one on the parent dir
#
#			find_cmd="find ${del_dir}/ -mindepth 1 -depth"
#	
#     		#could -delete here if we can implement an age test in find
#			#could then remove for loop below
#			#-ctime $days ? This seems to be an = rather than a -ge
#			#we want changed age, so we don't delete moved but unmodified files straight away.
#
#			
#			#This will currently delete files in subdirs
#			#and maintain the parent dir if another more recent file is del'd
#			#We want as the parent dir will keep getting refreshed and hence
#			#old data may accrue in subdirs
#
#			if [ $del_verbose ]; then
#				echo $find_cmd;
#			fi
#
#			for delfile in $($find_cmd); do
#
#				age=$(GetFileAge -c $delfile)
#					
#				if [ $? -ne 0 ]; then
#					#$age is error in this context
#					error_log="${error_log}${age}Failed to purge deleted file:\t$delfile\n"
#					continue
#				fi
#								
#
#				if [ $del_verbose ]; then
#					echo -e "$age days old:\t$delfile";
#				fi
#
#
#				if [ $age -ge $days ]; then
#					rm_cmd="$(which rm) -rf $delfile"
#
#					if [ $del_verbose ]; then
#						echo $rm_cmd;
#					fi
#
#					$rm_cmd
#
#					if [ $? -ne 0 ]; then
#						error_log="${error_log}Failed to purge deleted file:\t$delfile\n"
#						continue
#					fi
#				fi
#			done
#			
#		else               # MV ENTIRE PATH TO .DEL
#			del_path=$(echo $filedir | sed "s^${SOURCE_ROOT}^${del_dir}^")
#			#readlink already stripped trailing / for dir mv
#			del_path=$(GetDir $del_path)
#
#			if [ ! -d $del_path ]; then
#				mkdir -p $del_path
#				#catch error?
#			fi
#
#			#Do we want to have interactive by default here and override with -f?
#			mv_cmd="mv $filedir $del_path/$file"
#
#			if [ $del_verbose ]; then
#				echo $mv_cmd;
#			fi
#
#			$mv_cmd         	#mv'ing file updates last modified & changed date
#		
#			#Mirror the whole path under .del!
#			# - handle redundant naming
#			# - easier recovery/navigation
#			#This will require a recursive find when purging!
#
#			if [ $? -ne 0 ]; then
#				error_log="${error_log}Failed del file:\t$mv_cmd\n"
#				continue
#			fi
#		fi
#	done
#
#
#	if [ "$error_log" ]; then
#		echo -e "\nSummary of errors:\n$error_log"
#		return 1
#	fi
#
#
#	}


#This works slightly differently for ArchiveData, Distribute and del
#Do not subshell this as this will hide SOURCE/TARGET_ROOT/NAME
#e.g. error_of_fildir=$(_SetTargetAndSourceRoot $file)

_SetTargetAndSourceRoot(){
	OPTIND=1
	type=
	dir_txt=
	data_dir_txt="\n\t\$DATA_DIR\t=\t$DATA_DIR\nor\n\t\$GROUP_DATA_DIR\t=\t$GROUP_DATA_DIR"
	archive_dir_txt="\n\t\$ARCHIVE_DIR\t=\t$ARCHIVE_DIR\nor\n\t\$GROUP_ARCHIVE_DIR\t=\t$GROUP_ARCHIVE_DIR"
	derefd_filedir=
	no_warnings=

	#usage='usage: del  [ -d(ays, purge .del of files older than this value, at root defined by) ] FILES|DIRS'
	#add -l(ist) option here to see if we have already archived this in some form?
	#or we could just test for target file first and AskQuestion to confirm?
	#would need force flag to overwrite without AskQuestion

	while getopts ":adnh" opt; do
		case $opt in 
			a  ) type='ARCHIVE'; dir_text=$archive_dir_txt ;;
	        d  ) type='DATA'; dir_txt=$data_dir_txt ;;
            n  ) no_warnings=1 ;;
			h  ) echo $usage; return 0;;
			\? ) echo $usage; exit 1;;
		esac 
	done
  
	data_dir_txt="\n\t\$DATA_DIR\t=\t$DATA_DIR\nor\n\t\$GROUP_DATA_DIR\t=\t$GROUP_DATA_DIR"
	archive_dir_txt="\n\t\$ARCHIVE_DIR\t=\t$ARCHIVE_DIR\nor\n\t\$GROUP_ARCHIVE_DIR\t=\t$GROUP_ARCHIVE_DIR"


	if [ ! $type ]; then
		type='BOTH'
		dir_txt="${data_dir_txt}${archive_dir_txt}"
	fi
	
	i=1

	while [ $i -lt $OPTIND ]; do
		shift
		let i+=1
	done

	tmpfiledir1=$1

	#Check for more args here?
	TARGET_ROOT=
	TARGET_ROOT_NAME=
	SOURCE_ROOT=
	
	#text exists here
	derefd_filedir=$(readlink -e $tmpfiledir1)

	if [[ ! -e $derefd_filedir ]]; then
		
		if ! [ $no_warnings ]; then
			echo -e "File/dir argument does not exist:\t$tmpfiledir1"
		fi

		return 1
	fi

	#Could loop through array of var names and eval them to test

	if [[ -d $DATA_DIR ]] && [[ $derefd_filedir = $DATA_DIR* ]] && 
		( [ $type == 'BOTH' ] || [ $type == 'DATA' ] ); then
		TARGET_ROOT=$ARCHIVE_DIR
		TARGET_ROOT_NAME=ARCHIVE_DIR
		SOURCE_ROOT=$DATA_DIR
	elif [[ -d $GROUP_DATA_DIR ]] && [[ $derefd_filedir = $GROUP_DATA_DIR* ]] &&
		( [ $type == 'BOTH' ] || [ $type == 'DATA' ] ); then
		TARGET_ROOT=$GROUP_ARCHIVE_DIR
		TARGET_ROOT_NAME=GROUP_ARCHIVE_DIR
		SOURCE_ROOT=$GROUP_DATA_DIR
	elif [[ -d $GROUP_ARCHIVE_DIR ]] && [[ $derefd_filedir = $GROUP_ARCHIVE_DIR* ]] &&
		( [ $type == 'BOTH' ] || [ $type == 'ARCHIVE' ] ); then
		TARGET_ROOT=$GROUP_DATA_DIR
		TARGET_ROOT_NAME=GROUP_DATA_DIR
		SOURCE_ROOT=$GROUP_ARCHIVE_DIR
	elif [[ -d $ARCHIVE_DIR ]] && [[ $derefd_filedir = $ARCHIVE_DIR* ]] &&
		( [ $type == 'BOTH' ] || [ $type == 'ARCHIVE' ] ); then
		TARGET_ROOT=$DATA_DIR
		TARGET_ROOT_NAME=DATA_DIR
		SOURCE_ROOT=$ARCHIVE_DIR
	else
		if ! [ $no_warnings ]; then
			echo -e "Could not identify target/source root dir for:\t$tmpfiledir1"
			echo -e "Source needs to be in either:$dir_txt"
		fi

		return 1
	fi

	#Failure needs to be caught in caller using $?
}


GetFileAge(){
	OPTIND=1
	date_format=
	usage= 
	tmp_opts="$*"
	usage='usage:\t  GetFileAge -c(hanged age) || -m(odified age) file_or_dir_path'

	while getopts ":cmh" opt; do
		case $opt in 
			c  ) date_format='%Z' ;;
	        m  ) date_format='%Y' ;;
			h  ) echo -e $usage; return 0;;
			\? ) echo -e $usage; exit 1;;
		esac 
	done	
	
	i=1
	while [ $i -lt $OPTIND ]; do
		shift
		let i+=1
	done

	filename=$1


	if [ ! $date_format ] ||
		[ ! $filename ]; then
		echo -e "Invalid options\t${tmp_opts}\n${usage}"
		return 1
	fi


	#enable day and hours?
	#enable last modified?
	#apparently can't get created date from cmd line
	#would have to use perl or similar?

	day_factor=$((60 * 60 * 24))

	#seconds since Epoch 
	NOW=`date +%s`
	OLD=`stat -c $date_format $filename` 
	#Do we want modified or changed?
	#modified is actual content changed
	#change can be file moved or perms changed
	#So for del purge we want changed
	#other wise unmodified files would be removed straight aways
	


	if [ $? -ne 0 ]; then
		echo -e "Failed to GetFileAge for:\t$filename"
		return 1
	fi

	#This rounds down, which is what we want
	echo $(( ($NOW - $OLD) / $day_factor))
}

#remove this?

DistributeData(){
	#Keep this to data dir only
	#Or can we specify an optional destination dir?
	#Currently just mirror ArchiveData func




	error=$(CheckVariables ARCHIVE_DIR DATA_DIR)

	if [ $? -ne 0 ]; then
		echo $error
		return 1;
		#Should exit? If we ever use this in a script, or catch return in caller?
	fi


	for filedir in $*; do
	
		#Make sure we have a full path
		if [[ $filedir != /* ]]; then
			filedir=$(ls -d $PWD/$filedir)
		fi

		#This does not work if full path already defined!

		#Problem with this expanding links and not matching ARCHIVE_DIR
		#Could resolve by not using sym link for DIR vars?


		#Need to match $DATA_DIR here or skip with warning
		if [[ ! $filedir = $ARCHIVE_DIR* ]]; then
			echo -e "Skipping non archive file/directory:\t$filedir"
			echo -e "Source needs to be in \$ARCHIVE_DIR:\t$ARCHIVE_DIR"
			return 1
		else
			
			#Use ' instead of / due to / in path
			data_filedir=$(echo $filedir | sed -r "s'^${ARCHIVE_DIR}''")
			data_filedir="$DATA_DIR/$data_filedir"

		    #Need to mkdir in the arhcive?
			data_dir=$(GetDir $data_filedir)
			
			if [ ! -d $data_dir ]; then
				echo -e "Making data directory:\t$data_dir"
				mkdir -p $data_dir
			fi
			
			#-essh only necessary for remote archiving
			#This may cause data clashes, so we need to make sure paths are different?
			echo "rsync -essh -Wav $filedir $data_filedir"
			rsync -essh -Wav $filedir $data_filedir/
		fi
	done


}





################################################################################
# Func      : ChangeDir() 
# Desc      : if a given directory exists then it moves into it else exit 
# Args [1]  : path of directory to move into 
# Return    : none
# Exception : exits if 'cd' fails 
################################################################################

ChangeDir(){
	Execute cd $1
}

################################################################################
# Func      : ChangeMakeDir() 
# Desc      : if a given directory exists then it moves into it, else it creates
#             the directory and then moves into it 
# Args [1]  : path of directory to move into 
# Return    : none
# Exception : exits if 'cd' fails or MakeDir fails
################################################################################

ChangeMakeDir(){
    MakeDirs $1
	ChangeDir $1
}


################################################################################
# Func      : CheckFile() 
# Desc      : checks to see whether passed parameter is a file or a link to an 
#             existing file
# Arg [1]   : file name 
# Return    : none 
# Exception : exists if not a file 
################################################################################

CheckFile(){
    if [ ! -f $1 ]
    then
        echo "error : file $1 does not exist"
		exit 204
    fi
}

#Could change this to CheckFilesOrUsage?
#And just pass variable names like CheckVariables

CheckFilesOrUsage(){
	usage_string=$1
	shift
	file_variables=$*
	
	usage='usage: CheckFilesOrUsage "usage string" [/file/path]+ '
 	CheckVariablesOrUsage "$usage" usage_string $file_variables


	for file_var in $file_variables; do
		file=$(eval "echo \$$file_var")

		if [ ! -f $file ]; then
			echo "error : $file does not exist.
$usage_string"
			exit 204
		fi
	done
}




################################################################################
# Func      : Execute() 
# Desc      : executes unix command and tests return status 
# Arg [n]   : command
# Return    : none
# Exception : exits if command fails
################################################################################

# Can we trap errors here somehow, using signal handling?

# We're having problems with perl -e, mysql -e and mysql < file!
# Also haveing probles with cmd="echo 'sql statement' | mysql"
# Execute $cmd
# This just echos the whole of cmd and doesn't pipe!
# This works if we don't use $cmd and Execute directly. 
# Use $@ instead?

# No point in using this unless you want to exit
# Just test $? otherwise

Execute(){	
	#echo "Executing $*"	
  $*
	rtn=$? 

	if [ $rtn != 0 ]
	then 
		echo -e "Failed to:\t$*"
		exit $rtn 
		#This will exit the shell if we are running on cmdline
	fi

  return $rtn
}




################################################################################
# Func      : SedFile()
# Desc      : executes a  given (simple) sed command on a given file 
# Args [1]  : sed command 
# Args [2]  : file name 
# Return    : none
# Exception : exits if MoveFile fails 
################################################################################

SedFile(){
    _SED_CMD=$1
    _FILE=$2

    _TMP_FILE=/tmp/$(GetFilename $_FILE).$$
  
    sed "$_SED_CMD" $_FILE > $_TMP_FILE; _RTN=$?;

    if [ $_RTN -eq 0 ]
    then 
        MoveFile -f $_TMP_FILE $_FILE;
    else
        exit $_RTN
    fi
}

################################################################################
# Func      : GetDir() 
# Desc      : extracts the directory path from a path that includes a file name 
# Args [1]  : file name including path
# Return    : returns directory without file name  
# Exception : none 
################################################################################

GetDir(){
    echo ${1%/*}
}

################################################################################
# Func      : GetFilename() 
# Desc      : extracts the file name from a path that includes the directory 
# Args [1]  : file name including path
# Return    : returns file name without directory info 
# Exception : none 
################################################################################

GetFilename(){ 
    echo ${1##*/}
}


SetFileSize(){
	file=$1
	shift
	
	CheckFile $file
		
   	#Need to follow links here to get real size!
	file_size=($(ls -lLk $file))
	export file_size=${file_size[4]} 
}


################################################################################
# Func      : AskQuestion()
# Desc      : adds the appropriate switch (for Mac or OSF1) to prevent a new line
#             being added after a question.
# Args [1]  : question string
# Return    : #returns the given question with the appropriate switch ???
#             Return REPLY value. This is prevent having to access REPLY, which 
#             may persist from previous AskQuestion, hence preventing while [ $READ != $answer ]
#             Actually, this makes no difference as the caller still has to clean the reply var, but
#             it is more obvious what is going on.
# Exception : none 
################################################################################

#This causes a hang
#REPLY=$(AskQuestion "Would you like to use the following $seq_type file? [y|n] $file")


AskQuestion(){
    local question=
    local answer_regex=
    question=$1
    answer_regex=$2
      
    if [[ $answer_regex ]]; then
      question+=" $answer_regex";
      local first_try=
      echo -en "$question?" 
      read REPLY
    
      while [[ ! $REPLY =~ ^${answer_regex}$  ]]; do
          echo -e "$REPLY is not a recognized answer, please restrict to the following:  $answer_regex"
          read REPLY
      done

    else
        echo -en "$question?" 
        read REPLY
    fi
}


################################################################################
# Func      : ContinueOverride()
# Desc      : Executes AskQuestion unless a force flag has been specified in 
#             which case continues.  Probably need to set REPLY to Y here??????
# Args [1]  : question string
# Args [2]  : override flag
# Return    : 
# Exception : none 
################################################################################

ContinueOverride(){
    CheckVariables 1

    if [ ! $2 ]
    then
		AskQuestion "$1 [y|n]"

		if [[ $REPLY != [yY]* ]]; then
			echo "Exiting"
			exit
		fi
    else
		echo "Auto Continue"
    fi
}




################################################################################
# Func      : CheckCompressedFile()
# Desc      : checks whether a file is of a compressed format or not
# Args [1]  : file name
# Return    : none
# Exception : exits if not a compressed file and prints error msg
################################################################################

CheckCompressedFile(){
    _FILE=$1

    if [ $(isCompressedFile $_FILE) -eq $_FALSE ] 
    then 
        echo "Error: '$_FILE' is not a compressed file"
        exit 202 
    fi
}

################################################################################
# Func      : isCompressedFile()
# Desc      : tests whether a file is of a compressed format or not
# Args [1]  : file name 
# Return    : returns true if file is compressed else false  
# Exception : none
################################################################################

isCompressedFile(){
    _FILE2=$1

    file $_FILE2 | grep "compressed data" > /dev/null

    if [ $? -eq 1 ]
    then 
        echo $_FALSE
    else
        echo $_TRUE
    fi
}

################################################################################
# Func      : getDay()
# Desc      : gets the current day
# Args [0]  : none
# Return    : returns the current day in the format "1"
# Exception : none
################################################################################

getDay(){
    echo $(date '+%d') 
}

################################################################################
# Func      : getMonth()
# Desc      : gets the current month
# Args [0]  : none
# Return    : returns the current month in the format "1" 
# Exception : none
################################################################################

getMonth(){
    echo $(date '+%m') 
}

################################################################################
# Func      : getPreviousMonth()
# Desc      : gets the month previous to the current month
# Args [0]  : none
# Return    : returns the previous month in the format "1"
# Exception : none
################################################################################

getPreviousMonth(){
    if [ $(getMonth) -eq "01" ]
    then
        _PREV_MNTH=12
    else
        _PREV_MNTH=$(expr $(getMonth) - 1)
    fi
    echo $_PREV_MNTH
}

################################################################################
# Func      : getDaysInMonth()
# Desc      : gets the number of days in a given month
# Args [1]  : month
# Args [2]  : year
# Return    : returns the number of days in a month
# Exception : none
################################################################################

getDaysInMonth(){
    _MNTH=$1
    _YEAR=$2
    
    _DAYS=$(for i in `cal $_MNTH $_YEAR`; do echo; done; echo $i); _DAYS=$(expr $_DAYS);
    echo $_DAYS
}

################################################################################
# Func      : getYear()
# Desc      : gets the current day
# Args [0]  : none
# Return    : returns the current year in the format "2003"
# Exception : none
################################################################################

getYear(){
    echo $(date '+%Y')
}

################################################################################
# Func      : getPreviousYear()
# Desc      : gets the year previous to the current year
# Args [0]  : none
# Return    : returns the previous year in the format "2003"
# Exception : none
################################################################################

getPreviousYear(){
    echo $(expr $(getYear) - 1) 
}

################################################################################
# Func      : getPreviousRelativeYear()
# Desc      : gets the year relative to the previous month
#             eg  if the previous month was in the previous year it returns the
#             the previous year else it returns the current year.
# Args [0]  : none
# Return    : returns the and relevant year in the string format
#             "01 2003"
# Exception : none
################################################################################

getPreviousRelativeYear(){
    if [ $(getPreviousMonth) -eq "12" ]
    then
        _PREV_YEAR=$(getPreviousYear)
    else
        _PREV_YEAR=$(getYear)
    fi
    echo $_PREV_YEAR
}

################################################################################
# Func      : PadNumber()
# Desc      : reformats a number from "1" to "01"
# Args [1]  : number
# Return    : returns the given number in the format "01"
# Exception : none
################################################################################

padNumber(){
    _NUM=$1

    _PAD_NUM=$(echo $_NUM | awk '{printf("%02d",$1)}')
    echo $_PAD_NUM
}



isMac(){

	if [[ $(uname -s) = 'Darwin' ]]; then
		echo $_TRUE
	fi
}


PrintColour(){
	#Declare vars local first!
  local colour=
  local string=
  local usage=
  local nocolour=


  nocolour="\e[0m"
	usage="PrintColour\n
  Description:\tPrints text in colours\n
  Usage:\t\tPrintColour -c(olour) red|green|blue -s(tring) 'YOUR STRING TO PRINT']"

	OPTIND=1
	while getopts ":hc:s:" opt; do
		case $opt in 
			h  ) echo -e "$usage"; return 0;;
      c  ) colour=$OPTARG ;;
      s  ) string=$OPTARG ;;
      \? ) echo -e "$usage"; return 1;;
    esac 
  done

  [[ ! $string ]] && ( echo -e "Mandatory param -s(tring) not set\n$usage" && return 1)
  local colour_code=


  if [[ "$colour" = red ]];then
      colour_code="\e[0;31m"
  elif [[ "$colour" = blue ]];then
      colour_code="\e[0;34m"
  elif [[ "$colour" != green ]]; then
      colour_code="\e[0;32m"
  else
      echo -e "PrintColour: The colour $colour is not supported"
      return 1
  fi

  printf "$colour_code%s $nocolour\n" "$string"

}
