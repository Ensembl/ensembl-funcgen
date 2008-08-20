#!/usr/local/ensembl/bin/perl -w

=head1 DESCRIPTION

Copies the reg_feature data needed for the overlap analysis into a
specified database. Does a bit of denormalising and quite a lot of
filtering/QC. Final product is a table in the same format as the
genomic features called regulatory_features_filtered.

=head1 AUTHOR(S)

dkeefe@ebi.ac.uk

=head1 USAGE

set the current regulatory database in file ~/dbs/current_funcgen

reg_feats_4_classification.pl -e dk_reg_feat_classify_49

=head1 EXAMPLES

 reg_feats_4_classification.pl -e dk_funcgen_classify_51_1  -v2 
          
=head1 SEE ALSO

=head1 TO DO

  tidy up
  POD

=head1 CVS

 $Log: not supported by cvs2svn $
 Revision 1.1  2008/04/11 10:55:11  dkeefe
 Copies the reg_feature data needed for the overlap analysis into a
 specified database. Does a bit of denormalising and quite a lot of
 filtering/QC. Final product is a table in the same format as the
 genomic features called regulatory_features_filtered.








=cut


use strict;
use DBI;
use Env;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Getopt::Std;
use IO::Handle;
use IO::File;
use lib '/nfs/acari/dkeefe/src/personal/ensembl-personal/dkeefe/perl/modules/';
use DBSQL::Utils;
use DBSQL::DBS; # to use databases listed in ~/dbs.

use constant  NO_ROWS => '0E0';


my($user, $password, $driver, $host, $port);
my @temp_tables;
my $dump_dir = '/lustre/scratch1/ensembl/dkeefe/reg_feats_junk/';
my $id_list;
my $sp;
my $verbose = 0;
my $jump_label = '';
my $filtered_features_table = 'regulatory_features_filtered';

my %opt;

if ($ARGV[0]){
&Getopt::Std::getopts('v:u:p:s:H:h:e:j:P:', \%opt) || die ;
}else{
&help_text; 
}

# get configuration from environment variables
&config; # this may fail but config can be on command line


my $enc_db = 'dk_genomic_features_36i';# default, can be overridden by args
&process_arguments;

# hook up with the server
my $dbh = &make_contact($enc_db);
my $dbu = DBSQL::Utils->new($dbh);


# we need access to the latest ensembl funcgen
my $source_dbs = DBSQL::DBS->new('current_funcgen');
my $source_dbh = $source_dbs->connect();
my $source_dbu = DBSQL::Utils->new($source_dbh);

&source_checks($dbh,$dbu,$source_dbs->name);

if($jump_label){
    goto $jump_label;
}

# we need a subset of the source tables so we dump them and read them in to
# the target db
&backtick("rm -rf $dump_dir");
&backtick("mkdir $dump_dir");
my $dump_templ = "mysqldump --opt -h ".$source_dbs->host.
                                 " -u ".$source_dbs->user.
                                 " -P ".$source_dbs->port.
                                 " -p".$source_dbs->pass.
                                 " ".$source_dbs->name.
                                 ' %s '.
                                 " > $dump_dir".'%s'.".dump";

my $load_templ = "mysql -h ".$host.
                 " -u ".$user.
                 " -P ".$port.
                 " -p".$password.
                 " ".$enc_db.
                 " < $dump_dir".'%s'.".dump";

my @core_tables =( 'regulatory_attribute',
                   'feature_type',
                   'regulatory_feature',
                   'feature_set',
                   'data_set',
                   'supporting_set',
                   'external_feature',
                   'meta',
                   'annotated_feature'
		   );



foreach my $table (@core_tables){
                   
    unless( $source_dbu->table_exists($table)){die ("ERROR: The core database does not contain table $table")}

    my $command = sprintf($dump_templ,$table,$table);
    &backtick($command);
    $command = sprintf($load_templ,$table);
    #print $command."\n";
    &backtick($command);
}



&copy_with_rename($source_dbh,'seq_region',$dbh,'func_seq_region');


POST_IMPORT:
    my @sql;
my $temp1 = &new_temp();
push @sql, "drop table if exists seq_name_lookup";
push @sql, "create table seq_name_lookup select distinct seq_region_id,name as seq_region_name from func_seq_region";
push @sql,"alter table seq_name_lookup add index(seq_region_id)";
push @sql,"alter table seq_name_lookup add index(seq_region_name)";
#add seq_region_name to regulatory_features
push @sql, "create table $temp1 select r.*,s.seq_region_name from regulatory_feature r,seq_name_lookup s where r.seq_region_id = s.seq_region_id";
push @sql, "drop table regulatory_feature";
push @sql, "alter table $temp1 rename as regulatory_feature";

push @sql, "alter table regulatory_feature add index(feature_set_id)";


&execute($dbh,@sql) or die;

@sql = ();

# remove all but the current set of features
my $feature_set_id = $dbu->get_count("select feature_set_id from feature_set where name = 'RegulatoryFeatures'") or die "failed to get feature_set_id for regulatory features"; 
push @sql,"delete  from regulatory_feature where feature_set_id != $feature_set_id";
push @sql, "alter table regulatory_feature add index(regulatory_feature_id)";
push @sql, "alter table regulatory_feature add index(seq_region_start)";
push @sql, "alter table regulatory_feature add index(seq_region_end)";
push @sql, "alter table regulatory_feature add index(seq_region_name)";
push @sql, "alter table regulatory_feature add index(feature_type_id)";
&execute($dbh,@sql) or die;

FILTER:

# we remove features which have very long 'whiskers'
# these tend to contain multiple focus features which results in
# a group of regulatory features all of which have the same attributes
my $temp2 = &new_temp();
@sql = ();

push @sql,"drop table if exists $filtered_features_table";
push @sql,"create table $filtered_features_table select *, 'regulatory_feature' as feature_type from regulatory_feature";
push @sql, "alter table $filtered_features_table add index(seq_region_name)";
push @sql, "alter table $filtered_features_table add index(feature_type_id)";
push @sql, "alter table $filtered_features_table add index(feature_set_id)";
push @sql, "alter table $filtered_features_table add index(regulatory_feature_id)";
push @sql, "alter table $filtered_features_table add index(seq_region_start)";
push @sql, "alter table $filtered_features_table add index(seq_region_end)";

my $temp3 =  &new_temp();

push @sql,"drop table if exists $temp3";

push @sql,"create table $temp3 select ra.regulatory_feature_id , af.* from regulatory_attribute ra, annotated_feature af where ra.attribute_feature_id = af.annotated_feature_id";

push @sql,"create table $temp2 select j.regulatory_feature_id,count(*) as n_attribs,count(distinct j.feature_set_id) as n_attrib_types,max(j.seq_region_end)-min(j.seq_region_start) as len,min(j.seq_region_start)as attribs_start,max(j.seq_region_end) as attribs_end,rf.display_label from $temp3 j, $filtered_features_table rf where j.regulatory_feature_id=rf.regulatory_feature_id group by regulatory_feature_id ";
push @sql,"alter table $temp2 add index(regulatory_feature_id)";
push @sql,"alter table $temp2 add index(attribs_start)";
push @sql,"alter table $temp2 add index(attribs_end)";


push @sql,"delete r from $filtered_features_table r, $temp2 t2 where r.regulatory_feature_id = t2.regulatory_feature_id and t2.len > 5000";

# we also remove any reg_features which only contain a single type of attribute
push @sql,"delete r from $filtered_features_table r, $temp2 t2 where r.regulatory_feature_id = t2.regulatory_feature_id and t2.n_attrib_types < 2";


push @sql, "alter table $filtered_features_table add index(seq_region_name)";
push @sql, "alter table $filtered_features_table add index(feature_type_id)";
push @sql, "alter table $filtered_features_table add index(feature_set_id)";
push @sql, "alter table $filtered_features_table add index(regulatory_feature_id)";
push @sql, "alter table $filtered_features_table add index(seq_region_start)";
push @sql, "alter table $filtered_features_table add index(seq_region_end)";
push @sql, "alter table $filtered_features_table add index(display_label)";
# we want to remove features which occur in the same place and have the same
# attributes as one another, leaving only one representative
push @sql, "delete b from $filtered_features_table a, $filtered_features_table b, $temp2 t2a, $temp2 t2b where a.seq_region_name=b.seq_region_name and b.display_label = a.display_label and a.regulatory_feature_id < b.regulatory_feature_id and a.regulatory_feature_id = t2a.regulatory_feature_id and b.regulatory_feature_id = t2b.regulatory_feature_id and t2a.attribs_start=t2b.attribs_start and t2a.attribs_end = t2b.attribs_end";




&execute($dbh,@sql) or die;

# reg  features which overlap centromeric repeats tend to contain too many marks# which suggests there is something wrong with the mappings.
# so we get rid of them

@sql = ();
push @sql,"drop table if exists $temp3";

push @sql,"create table $temp3 select f.regulatory_feature_id from $filtered_features_table f, Satellite_centr_repeat s where s.seq_region_name = f.seq_region_name and s.feature_end >= f.seq_region_start and s.feature_start <= f.seq_region_end";
push @sql,"delete f from $filtered_features_table f,$temp3 t where f.regulatory_feature_id = t.regulatory_feature_id";


push @sql,"drop table if exists $temp3";

push @sql,"create table $temp3 select f.regulatory_feature_id from $filtered_features_table f, Satellite_repeat s where s.seq_region_name = f.seq_region_name and s.feature_end >= f.seq_region_start and s.feature_start <= f.seq_region_end";
push @sql,"delete f from $filtered_features_table f,$temp3 t where f.regulatory_feature_id = t.regulatory_feature_id";




# reg feats on the mitochondrial DNA are unlikely to use the same
# 'histone code' as the rest of the genome. they might even be artefacts.
# so we remove them

push @sql, "delete from $filtered_features_table where seq_region_name = 'MT'";

# we make the col names compatible with the genomic features in the analysis
push @sql,"alter table $filtered_features_table change column seq_region_start feature_start int(10) unsigned";
push @sql,"alter table $filtered_features_table change column seq_region_end feature_end int(10) unsigned";
push @sql,"alter table $filtered_features_table change column seq_region_strand feature_strand tinyint(2)";
#push @sql,"alter table $filtered_features_table change column seq_region_strand feature_strand tinyint(2)";
push @sql,"alter table $filtered_features_table modify column feature_type varchar(45)";

&execute($dbh,@sql) or die;


&clean_temp();
unlink($dump_dir);
$dbh->disconnect;
$source_dbh->disconnect;
exit;


 
###################################################################
sub source_checks{
    my($dbh,$dbu,$func_db) = @_;

    unless($dbu->table_exists('sources')){
	die "you must run the script gen_feats_4_classification.pl before this one";
    }

    my $q = "select db_name from sources where db_type = 'core'";
    my $core_db = $dbh->selectrow_array($q);

    my ($core_release) = $core_db =~ /.*_([0-9][0-9]_.*)/;
    my ($func_release) = $func_db =~ /.*_([0-9][0-9]_.*)/;

    unless($core_release eq $func_release){
	die "database versions don't match $core_db $func_db";
    }

    &commentary("version $core_release\nupdating sources table\n") if $verbose > 1;
    my @sql;
    push @sql, "delete from sources where db_type = 'funcgen'";
    push @sql, "insert into sources values('funcgen','$func_db')";
    &execute($dbh,@sql) or die;
    
}



sub copy_with_rename{
    my($s_dbh,$source_table,$t_dbh,$targ_table)=@_;

    my $q = "desc $source_table";
    my $desc_aaref=$s_dbh->selectall_arrayref($q);
    unless(defined $desc_aaref){die "failed on:\n $q\n".$s_dbh->errstr}

    &execute($t_dbh,"drop table if exists $targ_table");
    $q = " create table $targ_table (";
    my $select;
    foreach my $aref (@$desc_aaref){
	$q  .= $aref->[0].' '.$aref->[1].',';
        $select .=  $aref->[0].',';
    }
    chop $q; #remove last comma
    chop $select;
    $q .= ")";
    print $q."\n";
    &execute($t_dbh,$q) or die;

    $q = "select $select from $source_table";
    my $aaref = $s_dbh->selectall_arrayref($q);
    unless(defined $aaref){die "failed on:\n $q\n".$s_dbh->errstr}
    my @sql;
    foreach my $aref (@$aaref){
	$q= "insert into $targ_table values('".join("','",@$aref);
	$q .= "')";
	push @sql,$q;

    }

    &execute($t_dbh,@sql) or die;

}



# uses global array @temp_tables
sub new_temp{

    my $nom = 'temp_'.$$.'_'.scalar(@temp_tables);
    push @temp_tables,$nom;
    return $nom;

}

# uses global array @temp_tables
sub clean_temp{

    my @sql;
    foreach my $table (@temp_tables){
	push @sql,"drop table if exists $table";
    }

    &execute($dbh,@sql) or die;
}

sub backtick{
    my $command = shift;

    &commentary( "executing $command \n") if $verbose >1;

    my $res = `$command`;
    if($?){
        warn "failed to execute $command \n";
        warn "output :-\n $res \n";
	die "exit code $?";
    }

    return $res;
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
	
    	&commentary("Executing  -  $query\n") if $verbose;

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

    &commentary("\n");

    return(1);
}


# when using backticks to exec scripts the caller captures STDOUT
# its best therefore to have error on STDOUT and commentary on STDERR
sub commentary{
    print  "$_[0]";
}



sub config{
 
($user =     $ENV{'ENSFGUSER'}) or return(0); # ecs1dadmin
($password = $ENV{'ENSFGPWD'}) or return(0); #
($host   =   $ENV{'ENSFGHOST'}) or return(0); #localhost
($port =     $ENV{'ENSFGPORT'}) or return(0); #3360
($driver  =  $ENV{'ENSFGDRIVER'}) or return(0); #mysql
 
}
   
sub err{
    print STDERR "$_[0]\n";
}
  


sub make_contact{
    my $db = shift;

    unless($driver && $db && $host && $port && $user){
	&err("DB connection parameters not set");
        exit(1);
    }

    # hook up with the server
    my $dsn = "DBI:$driver:database=$db;host=$host;port=$port;mysql_local_infile=1";
    my $dbh = DBI->connect("$dsn","$user",$password, {RaiseError => 0,PrintError=>0});
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

    $sp = 'homo_sapiens';
    if (exists $opt{s}){
        $sp = lc($opt{s});
    }

    if (exists $opt{v}){
        $verbose = $opt{v};
    }



    if (exists $opt{j}){
        $jump_label = $opt{j} ;
    }

    if  (exists $opt{e}){
        $enc_db = $opt{e};
    } 

} 


sub help_text{
    my $msg=shift;

    if ($msg){
      print STDERR "\n".$msg."\n";
    }

    print STDERR <<"END_OF_TEXT";

    .pl [-h] for help
                  [-e] use a particular local encode catalog
                  [-H] <host machine> eg ecs2
                  [-u] <database user>
                  [-j] <label> jump to label then start execution
                       POST_IMPORT, FILTER
                  [-p] <mysql password> 
                  [-P] <mysql port> 
                  [-s] <species> eg -smus_musculus, default = homo_sapiens
                  [-v] <integer> verbosity level 0,1 or 2 
                  [-] 
                  [-] 
                  [-] <> 
                  [-] <> 


END_OF_TEXT


    if($msg){
        exit(1);
    }else{
        exit(0);
    }
}
