#!/usr/local/bin/bash


# File : funcs.sh

_TRUE=1;
_FALSE=!$_


################################################################################
# Func      : CheckUser()
# Desc      : checks the a user user name against the current user
# Args [1]  : user name
# Return    : none 
# Exception : exits if user names dont match 
################################################################################

Checkuser(){
    _REQUIRED_USER=$1
    _CURRENT_USER=$(whoami)

    if [ $_REQUIRED_USER != $_CURRENT_USER ]
    then
        echo "error : script must be run by '$_REQUIRED_USER' not '$_CURRENT_USER'";exit 205;
    fi
}

################################################################################
# Func      : CheckVariables()
# Desc      : check that all passed environment variables have been assigned 
#             values
# Args [n]  : list of variables to check 
# Return    : list of variables that cannot be found
# Exception : none
################################################################################

CheckVariables(){
    for _VAR in $*
    do
        _VAL=$(eval "echo \$$_VAR")

        if [ -z "$_VAL" ]; then _LINE="$_LINE \$$_VAR"; fi;
    done

    if [ ! -z "$_LINE" ]
    then
        echo "Environment variable(s) :$_LINE are not set"; exit 100;
    fi
}

################################################################################
# Func      : CheckBinaries() 
# Desc      : checks that binary files exist
# Args [n]  : list of binaries to check 
# Return    : list of binaries that cannot be found
# Exception : none
################################################################################

CheckBinaries(){
    for _VAR in $*
    do
        which $_VAR | grep ^no > /dev/null
        
        if [ $? -eq 0 ]; then _LINE="$_LINE '$_VAR'"; fi;
    done
    
    if [ ! -z "$_LINE" ]
    then 
        echo "Binary(s) : $_LINE cannot be found"; exit 201;
    fi
}

################################################################################
# Func      : CheckDIR()
# Desc      : checks to see whether passed parameter is a directory
# Arg [1]   : directory name
# Return    : returns true if it is a directory
# Exception : exists if not a directory
################################################################################

CheckDir(){
    _DIR=$1

    if [ ! -d $_DIR ]
    then
        echo "error : directory $_DIR does not exist"; exit 203;
    fi
}

################################################################################
# Func      : MakeDir() 
# Desc      : if a given directory does not exist then it will make it 
# Args [n]  : list of directory(s) to make 
# Return    : none
# Exception : exits if 'mkdir' fails 
################################################################################

MakeDir(){
     for _DIR_NAME in $*
     do
         if [ ! -d $_DIR_NAME ]
         then
             mkdir -p $_DIR_NAME; _RTN=$?; if [ $_RTN -ne 0 ]; then exit $_RTN; fi;
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
    _DIR_NAME=$1

    cd $_DIR_NAME; _RTN=$?; if [ $_RTN -ne 0 ]; then exit $_RTN; fi;
}

################################################################################
# Func      : ChangeDirMake() 
# Desc      : if a given directory exists then it moves into it, else it creates
#             the directory and then moves into it 
# Args [1]  : path of directory to move into 
# Return    : none
# Exception : exits if 'cd' fails or MakeDir fails
################################################################################

ChangeDirMake(){
    _DIR_NAME=$1

    MakeDir $_DIR_NAME

    cd $_DIR_NAME; _RTN=$?; if [ $_RTN -ne 0 ]; then exit $_RTN; fi;
}


################################################################################
# Func      : CheckFile() 
# Desc      : checks to see whether passed parameter is a file 
# Arg [1]   : file name 
# Return    : none 
# Exception : exists if not a file 
################################################################################

CheckFile(){
    _FILE=$1

    if [ ! -f $_FILE ]
    then
        echo "error : file $_FILE does not exist"; exit 204;
    fi
}

################################################################################
# Func      : CopyFile() 
# Desc      : executes unix 'cp' command and tests return status 
# Arg [n]   : [switches] source target (eg CopyFile -fr a b)
# Return    : none
# Exception : exits if 'cp' fails
################################################################################

CopyFile(){

    cp $*; _RTN=$?; if [ $_RTN -ne 0 ]; then exit $_RTN; fi;
}

################################################################################
# Func      : MoveFile()
# Desc      : executes unix 'mv' command and tests return status 
# Arg [n]   : [switches] source target (eg MoveFile -fi a b)
# Return    : none
# Exception : exits if 'mv' fails
################################################################################

MoveFile(){

    mv $*; _RTN=$?; if [ $_RTN -ne 0 ]; then exit $_RTN; fi;
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

################################################################################
# Func      : AskQuestion()
# Desc      : adds the appropriate switch (for Mac or OSF1) to prevent a new line
#             being added after a question 
# Args [1]  : question string
# Return    : returns the given question with the appropriate switch 
# Exception : none 
################################################################################

AskQuestion(){
   _QUESTION=$1

   #if [ $(isMac) -eq $_TRUE ]
   #then
   #    echo -n "$_QUESTION"
   #else
  #     echo -n "$_QUESTION" 
   #fi

    echo -n "$_QUESTION?  " 
    read REPLY
}

################################################################################
# Func      : isMac()
# Desc      : tests whether machine is a Mac or not 
# Args [0]  : 
# Return    : returns true if OS is Darwin else false 
# Exception : none 
################################################################################

isMac(){
    if [ $(uname) = "Darwin" ]
    then
        echo $_TRUE
    else
        echo $_FALSE
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

################################################################################
# Func      : toLowerCase()
# Desc      : converts a string to lower case
# Args [1]  : string
# Return    : returns lower case version of string
# Exception : none
################################################################################

toLowerCase(){
    echo $1  | tr '[A-Z]' '[a-z]'
}

################################################################################
# Func      : toUpperCase()
# Desc      : converts a string to upper case
# Args [1]  : string
# Return    : returns upper case version of string
# Exception : none
################################################################################

toUpperCase(){
    echo $1  | tr '[a-z]' '[A-Z]'
}

################################################################################
# Func      : Distribute()
# Desc      : securely copies a local directory to a remote directory
# Args [1]  : source directory
# Args [2]  : remote machine name
# Args [3]  : remote user name
# Args [4]  : remote directory
# Return    : scp exit status
# Exception : none
# Note      : Although -r is used in the scp command, which is for recursing
#             through directories, this function can also upload single and
#             multiple files as -r has no adverse effects when doing this.
################################################################################

Distribute(){
   CheckBinaries scp
   
   _LOCAL_DIR=$1
   _REMOTE_HOST=$2
   _REMOTE_USER=$3
   _REMOTE_DIR=$4

   _REMOTE_LOGIN="$_REMOTE_USER@$_REMOTE_HOST"
   scp -r  $_LOCAL_DIR $_REMOTE_LOGIN:$_REMOTE_DIR
}


