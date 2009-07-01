#!/software/bin/perl

####!/sw/arch/bin/perl

=head1 NAME

load_bed.pl -- 

=head1 SYNOPSIS

=head1 DESCRIPTION

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);

my ($pass,$port,$host,$user,$dbname,$prefix,$species,$help,$man, 
    $data_version,$debug);

$host='ens-genomics1';
$port=3306;
$user='ensadmin';
$pass='ensembl';
#$dbname='encode_das';
$dbname='eFG_DAS';

GetOptions (
            'host|h=s'         => \$host,
            'port:i'           => \$port,
            'user|u=s'         => \$user,
            'pass|p:s'         => \$pass,
            'dbname|d=s'       => \$dbname,
            'prefix=s'         => \$prefix,
            'data_version|v=s' => \$data_version,
            'species=s'        => \$species,
            'help|?'           => \$help,
            'man|m'            => \$man,
            'debug'            => \$debug,
            );


throw("No hydra source table prefix specified!") if (!$prefix);


# Connect to the database.
use DBI();
my $dbh = DBI->connect("DBI:mysql:database=$dbname;host=$host;port=$port",
                       "$user", "$pass",
                       {'RaiseError' => 1});


$ENV{LSB_JOBINDEX} = 1 unless (exists $ENV{LSB_JOBINDEX});
my $file = @ARGV[$ENV{LSB_JOBINDEX}-1];

throw("File $file does not exit") unless (-f $file);

print "loading $file...\n";

open(CMD, "file $file |")
    or die "Can't execute command: $!";
my $gzip = grep {/gzip compressed data/} (<CMD>);
close CMD;

if ($gzip) {
    print "decompressing file $file\n";
    system("gzip -d $file") == 0
        or die "Can't decompress file $file: $!";
    $file =~ s/\.gz$//;
    #print $file;
    throw("File $file does not exit") unless (-f $file);
}

(my $table_name=$file) =~ s,^(.*/)?(.+)\.bed$,$2,;
$table_name =~ s,\.,_,g;
$table_name = $prefix.'_'.$table_name;
print 'table name: ', $table_name, "\n";

my $sth = $dbh->prepare("DROP TABLE IF EXISTS `$table_name`;");
$sth->execute();

$sth = $dbh->prepare("CREATE TABLE `$table_name` (
    `feature_id`    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
    `seq_region`    VARCHAR(20) NOT NULL,
    `start`         INT(10) UNSIGNED NOT NULL DEFAULT '0',
    `end`           INT(10) UNSIGNED NOT NULL DEFAULT '0',
    `name`          VARCHAR(40) NOT NULL,
    `score`         FLOAT NOT NULL DEFAULT '0',
    `strand`        ENUM('0','+','-') DEFAULT '0',
    `note`          VARCHAR(40) DEFAULT NULL,
    PRIMARY KEY     (`feature_id`),
    KEY `seq_region_idx` (`seq_region`, `start`)
) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;");
$sth->execute();

if ($file =~ m/_reads/) {

    $sth = $dbh->prepare("LOAD DATA LOCAL INFILE '$file' INTO TABLE $table_name 
               (seq_region,start,end,name,\@mm,strand,score) 
               set seq_region=replace(seq_region, 'chr', ''), note=concat('mm=',\@mm);");

} elsif ($file =~ m/_profile/) {

    $sth = $dbh->prepare("LOAD DATA LOCAL INFILE '$file' INTO TABLE $table_name 
               (seq_region,start,end,name,score,strand) 
               set seq_region=replace(seq_region, 'chr', '');");

}

$sth->execute();
$sth->finish();
$dbh->disconnect();

if ($gzip) {
    print "compressing file $file...\n";
    system("gzip $file") == 0
        or die "Can't compress file $file: $!";
}

__END__




warn ($sql);

#open(SQL, "echo \'$sql\' | $MYSQL heroic_das |")
#    or throw "Can't execute SQL $sql: $!";
#print while (<SQL>);
#close SQL;


__END__

if [ $# -lt 1 ]; then
    echo "USAGE: $0 <file(s)>"
    exit 1
fi

# create array with file names to be used with LSB_JOBINDEX
n=0
for f in $@; do
    #echo $f;
    ((n++));
    file[$n]=$f;
done

f=${file[$LSB_JOBINDEX]}
echo $f

if [ ! -f $f ]; then
    echo "File $f does not exist!"
    exit 1;
fi

TABLE=heroic_$(echo $f | sed -e 's,.*/,,;s,\.bed$,,;s,\.,_,g')
#TABLE=heroic_$(echo $1 | sed -e 's,.*/,,;s,\.bed$,,;s,\.,_,g')

echo $TABLE

exit;

MYSQL='mysql -h mysql-heroic -u admin -P 4127 -pNqrC80eM heroic_das'

echo | $MYSQL <<EOF

DROP TABLE IF EXISTS \`$TABLE\`;
CREATE TABLE \`$TABLE\` (
                   \`feature_id\`    INT(10) UNSIGNED NOT NULL AUTO_INCREMENT,
                   \`seq_region\`    VARCHAR(30) NOT NULL,
                   \`start\`         INT(10) UNSIGNED NOT NULL DEFAULT '0',
                   \`end\`           INT(10) UNSIGNED NOT NULL DEFAULT '0',
                   \`name\`          VARCHAR(40) NOT NULL,
                   \`score\`         FLOAT NOT NULL DEFAULT '0',
                   \`strand\`        ENUM('0','+','-') DEFAULT '0',
                   \`note\`          VARCHAR(40) DEFAULT NULL,
                   PRIMARY KEY   (\`feature_id\`),
                   KEY \`seq_region_idx\` (\`seq_region\`, \`start\`)
               ) ENGINE=MyISAM DEFAULT CHARSET=latin1 COLLATE=latin1_bin;

EOF

case "$TABLE" in
    *_profile* )
        echo "LOAD DATA LOCAL INFILE '$f' INTO TABLE $TABLE (seq_region,start,end,name,score,strand) 
              set seq_region=replace(seq_region, 'chr', '');" | $MYSQL
        ;;
    *_reads )
        echo "LOAD DATA LOCAL INFILE '$f' INTO TABLE $TABLE (seq_region,start,end,name,@mm,strand,score) 
              set seq_region=replace(seq_region, 'chr', ''), note=concat('mm=',@mm);" | $MYSQL
    ;;
esac

# reads
#    printf \"INSERT INTO $TABLE (seq_region,start,end,name,score,strand) VALUES ('%s',%d,%d,'%s',%.2f,'%s');\n\",
#           \$1, \$2, \$3, \$4, \$7, \$6;
# profile
#    printf \"INSERT INTO $TABLE (seq_region,start,end,name,score,strand) VALUES ('%s',%d,%d,'%s',%.2f,'%s');\n\",
#           \$1, \$2, \$3, \$4, \$5, \$6;
