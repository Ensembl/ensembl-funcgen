
#$1 is dbanme $2 is password $3 is drop arg = 'drop'
#Set e.g. MYSQL_ARGS="-hecs2 -P3364 -uenadmin"
#else localhost, default port and -uensadmin are used

#MYSQL_ARGS=${MYSQL_ARGS:="-uensadmin"}" -p$2"

if [[ ! $2 ]]
then
    echo "Need to define dbname and password arguments"
    echo "e.g. create_db.sh homo_sapiens_funcgen_44_36f 'password'"
else

    dbname=$1
    pass=$2

    if [[ $3 ]]
    then 
    
	if [[ $3 != drop ]]
	then 
	    echo "Unrecognised argument, did you mean to specify 'drop' as the 3rd argument?"
	    exit
	else
	    drop=$3
	fi
    fi



    MYSQL_ARGS="-hens-genomics1 -uensadmin -P3306 -p${2}"

    echo "Connecting using $MYSQL_ARGS"
    
    dblike=$(echo "show databases like '$dbname'" | mysql $MYSQL_ARGS)
    
    dblike=$(echo $dbname | sed 's/Database (.*)//')

   
    if [[ $dbname = $dblike ]]
    then

	if [[ $drop ]]
	then
	    echo "Dropping DB $dbname"
	    echo "DROP DATABASE IF EXISTS \`$1\`;" |  mysql $MYSQL_ARGS
	else
	    echo "DB $dbname already exists, please drop the database manually or using the 3rd argument 'drop'"
	    echo "e.g. e.g. create_db.sh homo_sapiens_funcgen_44_36f 'password' drop"
	    exit;
	fi
    fi

    echo "Creating DB $dbname"
    echo "CREATE DATABASE \`$dbname\`;" | mysql $MYSQL_ARGS

    mysql $MYSQL_ARGS $dbname < efg.sql
fi
