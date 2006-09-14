
#$1 is dbanme $2 is password
#Set e.g. MYSQL_ARGS="-hecs2 -P3364 -uenadmin"
#else localhost, default port and -uensadmin are used

#MYSQL_ARGS=${MYSQL_ARGS:="-uensadmin"}" -p$2"

MYSQL_ARGS="-hecs2 -uensadmin -P3362 -p${2}"

echo "Connecting using $MYSQL_ARGS"

echo "DROP DATABASE IF EXISTS \`$1\`;	CREATE DATABASE \`$1\`;" | mysql $MYSQL_ARGS

mysql $MYSQL_ARGS $1 < efg_test.sql
