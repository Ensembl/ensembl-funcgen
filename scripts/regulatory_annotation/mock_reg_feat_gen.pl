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

=head1 DESCRIPTION

A mock feature creation system that allows us to specify where mock features cannot occur (eg in repeat regions) and where they must occur eg in introns

Used to create a set of randomly placed mock regulatory features which are used for chi square calculation in the regulatory feature classification script reg_feat_gen_feat_overlaps.pl

=head1 USAGE

run after reg_feats_4_classification.pl

=head1 EXAMPLES

  mock_reg_feat_gen.pl -e dk_genomic_features_36k -i regulatory_features_filtered -o mockreg_features_filtered


create features which lie inside introns but not in repeats

 mock_reg_feat_gen.pl -i regulatory_features_filtered -l havana_intron -r rmsk


=head1 SEE ALSO

=cut


use strict;
use DBI;
use Env;
use Getopt::Std;
use IO::Handle;
use IO::File;
use lib '/nfs/users/nfs_d/dkeefe/src/personal/ensembl-personal/dkeefe/perl/modules/';
use DBSQL::Utils;

use constant  NO_ROWS => '0E0';


my($user, $password, $driver, $host, $port);
my $verbose = 0;
my $intable='';
my $outtable = 'mockreg_features_filtered';
my $enc_db = 'dk_genomic_features_36k'; 
my $create_non_redundant = '';
my @oblig_tables = ();
my @avoid_tables = ();
my $oblig_region_table = 'mock_loc_table';
my $banned_region_table = 'banned_mock_table';

my %opt;
if ($ARGV[0]){
&Getopt::Std::getopts('u:p:H:h:e:no:i:P:l:r:', \%opt) || &help_text ;
}else{
&help_text; 
}

# get configuration from environment variables
&config; # this may fail but config can be on command line

&process_arguments;
my $originals_table = $intable;
if(@oblig_tables < 1){ # there is nowhere the mocks MUST lie
    $oblig_region_table = 'none';
}


# hook up with the server
my $dbh = &make_contact($enc_db);
my $dbu = DBSQL::Utils->new($dbh);

unless( $dbu->table_exists('seq_region')){die ("ERROR: You need to import the seq_region table into database $enc_db before running this script")}

#goto FRAG_GEN;

# add coords from other tables to banned regions
&prepare_banned_enc_table($dbh,$banned_region_table,@avoid_tables);

unless($oblig_region_table eq 'none'){
    &prepare_oblig_enc_table($dbh,$oblig_region_table,@oblig_tables);
}

FRAG_GEN:

my @sql;
my @res;
# first we need the lengths of the seq regions so we know the max 
# tfrag_enc_end value we can generate
my %enc_len;
my $q = "select distinct seq_region_name from $originals_table";
my $col_ref = $dbh->selectcol_arrayref($q) or die $dbh->errstr;
my $clause = " where name in('";
$clause .= join("','",@$col_ref);
$clause .= "')";

$q="select name, length from seq_region" .$clause;
my $aaref = $dbh->selectall_arrayref($q) or die $dbh->errstr;
foreach my $aref (@$aaref){
    $enc_len{$aref->[0]} = $aref->[1];
}


# seed random number generator
srand(time);


# next we need the seq_region and length of each original feature
my %frag_len;
my $q = "select seq_region_name, feature_end - feature_start +1 as length,feature_type,feature_strand,regulatory_feature_id from $originals_table order by (feature_end - feature_start +1) desc ";
my $aaref = $dbh->selectall_arrayref($q) or die "failed on \n$q\n".$dbh->errstr;

foreach my $aref (@$aaref){
    my $enc = $aref->[0];    
    my $len = $aref->[1];
    my $type = $aref->[2];
    my $ori = $aref->[3];
    my $id =  $aref->[4]; 
    &commentary( "length $len \n") if $verbose;

    my($start,$end) = &get_rand_region($dbh,$enc,$len,
                                       $banned_region_table,
                                       $oblig_region_table,
                                       %enc_len
                                       );

    #print "$start $end\n";


    my $q = "insert into temp_$$ values(";
    $q .= $id.",";
    $q .= "\'".$enc."\'," ;
    $q .= $start.",";
    $q .= $end.",";
    $q .= "\'".$ori."\'," ;
    $q .= "\'".$type."\'" ;
    $q .= ")";

    &commentary( $q."\n") if $verbose;

    push @res,$q;

    &add_frag_to_banned_list($dbh,$enc,$start,$end,$banned_region_table);
 
}


push @sql, "drop table if exists temp_$$";
push @sql, "create table temp_$$ (regulatory_feature_id int(10) unsigned,seq_region_name char(10),feature_start int(10) unsigned,feature_end int(10) unsigned,feature_strand tinyint(2),feature_type varchar(40))"; 

push @sql,@res;

push @sql, "drop table if exists $outtable";
push @sql, "create table $outtable select t.*,o.binary_string from temp_$$ t, $originals_table o where t.regulatory_feature_id = o.regulatory_feature_id"; # select f.seq_region_name, f.feature_type,f.feature_start,f.feature_end ,f.feature_strand from temp f";

push @sql, "alter table $outtable add index(regulatory_feature_id)";
push @sql, "alter table $outtable add index(seq_region_name)";
push @sql, "alter table $outtable add index(feature_start)";
push @sql, "alter table $outtable add index(feature_type)";
push @sql, "alter table $outtable add index(feature_strand)";
push @sql, "drop table if exists $oblig_region_table";
push @sql, "drop table if exists $banned_region_table";
#push @sql, "drop table if exists temp";
&execute($dbh,@sql) or die;




$dbh->disconnect;
exit;


 
###################################################################
sub add_frag_to_banned_list{
    my($dbh,$enc,$start,$end,$banned_region_table)= @_;

    my $q = "insert into $banned_region_table values('$enc',$start,$end) ";

    &execute($dbh,$q) or die("failed on \n$q\n".$dbh->errstr);


}


sub prepare_banned_enc_table{
    my $dbh = shift;
    my $banned_region_table = shift;
    my @tables = @_;

    my @sql;

    push @sql,"drop table if exists $banned_region_table";
    push @sql,"create table $banned_region_table ( seq_region_name varchar(10), feature_start int(11), feature_end int(11))";
    my $iter = 1;
    foreach my $table (@tables){

        push @sql,"insert into $banned_region_table select seq_region_name, feature_start, feature_end  from $table ";

    }

    push @sql,"alter table $banned_region_table add index(seq_region_name)";
    push @sql,"alter table $banned_region_table add index(feature_start)";
    push @sql,"alter table $banned_region_table add index(feature_end)";

    &execute($dbh,@sql) or die;

}


sub prepare_oblig_enc_table{
    my $dbh = shift;
    my $oblig_region_table = shift;
    my @tables = @_;

    my @sql;

    push @sql,"drop table if exists $oblig_region_table";

    my $iter = 1;
    foreach my $table (@tables){

	if($iter == 1){
            push @sql,"create table $oblig_region_table select seq_region_name, feature_start, feature_end  from $table ";
	}else{
            push @sql,"insert into $oblig_region_table select seq_region_name, feature_start, feature_end  from $table ";
	}
	$iter ++;
    }

    push @sql,"alter table $oblig_region_table add index(seq_region_name)";
    push @sql,"alter table $oblig_region_table add index(feature_start)";
    push @sql,"alter table $oblig_region_table add index(feature_end)";

    &execute($dbh,@sql) or die;

}



# returns start and end 
sub get_rand_region{
    my ($dbh,$enc,$len,$banned_regions_table,$good_regions_table,%enc_len) = @_;

    my $bad = 1;
    my $end;
    my $start;
    my $try = 1;
    while($bad){
	&commentary("try $try\n") if $verbose;
	$try++;
	if($try > 1000000){
	    &err("can't find a random region for seq_region $enc, length $len\n");
	    exit(1);
        }

    # get a random number which is within the encode region
	$end = int(rand($enc_len{$enc}));

    # the number must be higher than the length of the frag
        if($end > $len){
            $bad = 0;
	}
        $start = $end - $len +1;

        unless($good_regions_table eq 'none'){
	    if(&in_good_region($dbh,$good_regions_table,$enc,$start,$end)){
		$bad = 0;
	    }else{
		$bad = 1;
		&commentary( "not good region\n") if $verbose;
	    }
	}

        if($bad == 0){
	    if(&in_banned_region($enc,$start,$end,$banned_regions_table)){
		&commentary( "len $len   banned\n") if $verbose;
		$bad = 1;
	    } 
	}
    }

    return ($start,$end);

}

# inside the good regions
sub in_good_region{
    my($dbh,$good_regions_table,$enc,$start,$end)=@_;


    my $q = "select * from $good_regions_table where seq_region_name = '$enc' and $start >= feature_start and $end <= feature_end limit 1 ";

    my $aaref = $dbh->selectall_arrayref($q);

    unless(defined $aaref){die $dbh->errstr}

    #print scalar(@$aaref)."\n";
	
    if(scalar(@$aaref) == 0){return(0)}
   
    return(1);
}


# overlapping the bad regions
sub in_banned_region{
    my ($enc,$start,$end,$banned_regions_table) = @_;

    my $q = "select * from $banned_regions_table where seq_region_name = '$enc' and $start <= feature_end and $end >= feature_start limit 1 ";

    my $aaref = $dbh->selectall_arrayref($q);

    unless(defined $aaref){die $dbh->errstr}

    #print scalar(@$aaref)."\n";
	
    if(scalar(@$aaref) == 0){return(0)}
   
    return(1);

}





# execute
#
#  Arg [1]   : scalar database handle from DBI module
#  Arg [2]   : array of scalars containing lines of text in the format of SQL statements
#  Function  : submits each element of Arg2 to the $dbh->prepare and 
#              $dbh->execute methods and checks for successful execution. 
#              Returns 0 on failure, 1 on success. Emits error messages 
#              via &err.
#  Returntype: int
#  Exceptions: none
#  Example   : execute($dbh, @array);

sub execute{
    my $dbh = shift;
    my (@array)=@_;

    my $sth;

    &commentary("processing SQL.") if $verbose;

    foreach my $query(@array){
	
    	&commentary(".") if $verbose;
    	&commentary("$query\n") if $verbose;

	    unless($sth=$dbh->prepare($query)){
               &err("Error: preparation of statement failed on: $query\n");
               &err("database error message: ".$dbh->errstr."\n");
               return(0);
	    }

        unless($sth->execute){ # returns true on success
            &err("Error: statement execution failed on: $query\n");
            &err("statement handle error message:".$sth->errstr."\n");
            return(0);
	    }
    }

    &commentary("\n") if $verbose ;

    return(1);
}


# when using backticks to exec scripts the caller captures STDOUT
# its best therefore to have error on STDOUT and commentary on STDERR
sub commentary{
    print  STDERR "$_[0]";
}

   
sub err{
    print STDOUT "$_[0]\n";
}
  
sub errexit{

    err "$_[0]";
    exit(1);

}

sub config{
 
($user =     $ENV{'ENSFGUSER'}) or return(0); 
($password = $ENV{'ENSFGPWD'}) or return(0); 
($host   =   $ENV{'ENSFGHOST'}) or return(0); 
($port =     $ENV{'ENSFGPORT'}) or return(0); #3360
($driver  =  $ENV{'ENSFGDRIVER'}) or return(0); #mysql
 
}



sub make_contact{
    my $db = shift;

    unless($driver && $db && $host && $port && $user){
	&err("DB connection parameters not set");
        exit(1);
    }

    # hook up with the server
    my $dsn = "DBI:$driver:database=$db;host=$host;port=$port;mysql_local_infile=1";
    my $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 0});
    if ($dbh){
        print STDERR ("connected to $db\n");
    }else{
        &err("failed to connect to database $db");
        exit(1);
 
    }   

            
    return $dbh;
}



sub process_arguments{

    if ( exists $opt{'h'} ){ 
        &help_text;
    }


    if (exists $opt{H}){
        $host = $opt{H}; 
    }

    if (exists $opt{p}){
        $password = $opt{p}; 
    }

    if (exists $opt{P}){
        $port = $opt{P}; 
    }

    if (exists $opt{u}){
        $user = $opt{u}; 
    }


    if (exists $opt{o}){
        $outtable = $opt{o};
    }  


    if (exists $opt{i}){
        $intable = $opt{i};
    }else{
	&help_text("you must give a tablename for the original features");
    }  

    if (exists $opt{n}){
        $create_non_redundant = 1;
    }

    if (exists $opt{l}){
	@oblig_tables = split(/\s+/,$opt{l});
    }

    if (exists $opt{r}){
	@avoid_tables = split(/\s+/,$opt{r});
    }

    if  (exists $opt{e}){
        $enc_db = $opt{e};
    }else{
	#&help_text("you must give a database name");
    }


} 


sub help_text{
    my $msg=shift;

    if ($msg){
      print STDERR "\n".$msg."\n";
    }

    print STDERR <<"END_OF_TEXT";

   fg_mock_gen.pl [-h] for help
                   -e <database_name> name of database for encode work
                  [-i] <input_table> - name of a table with original features
                  [-l] <string> list of table names in quotes: regions where
                       mock features must lie
                  [-r]  <string> list of table names in quotes: regions where
                       mock features must NOT lie eg rmsk
                  [-H] <host machine> eg ecs2
                  [-u] <database user> 
                  [-o] <output_table> - name of table containing mock features
                  [-p] <mysql password> 
                  [-P] <mysql port> 



END_OF_TEXT


    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}
