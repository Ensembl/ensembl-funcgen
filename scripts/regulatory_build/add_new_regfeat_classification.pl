#!/software/bin/perl -w

=head1 CVS

  $Log: not supported by cvs2svn $

=head1 AUTHOR

  Damian Keefe dkeefe@ebi.ac.uk

=cut


###!/opt/local/bin/perl -w

BEGIN{
    if($ENV{SHELL} eq 'bash'){
	if(! defined $ENV{'EFG_DATA'}){
		if(-f "$ENV{'HOME'}/src/ensembl-functgenomics/scripts/.efg"){
			system (". ~/src/ensembl-functgenomics/scripts/.efg");
		}else{
			die ("This script requires ensembl-functgenomics/scripts/.efg\n".
				 "Please source it before running this script\n");
		}
	}
    }else{

        use lib '/nfs/acari/dkeefe/src/personal/ensembl-personal/dkeefe/perl/modules/';
        use lib '/nfs/acari/dkeefe/src/personal/ensembl-functgenomics/modules';
        use lib '/nfs/acari/dkeefe/src/personal/ensembl/modules';

    }

}

use Bio::EnsEMBL::Analysis;
#use Bio::EnsEMBL::Funcgen::Helper;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;#

use Bio::EnsEMBL::DBSQL::DBAdaptor;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(species_chr_num open_file);
use Bio::EnsEMBL::Utils::ConfigRegistry;
use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Funcgen::RegulatoryFeature;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use File::Path qw(mkpath);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(mean);
use Devel::Size qw(size total_size);
use Bio::EnsEMBL::Funcgen::ResultFeature;
use Tie::File;
use Data::Dumper;
use DBSQL::Utils;
#use Bio::MAGE::XMLUtils;
#use Bio::MAGE::Experiment::Experiment;

$| = 1; #no output buffer
use strict;

my @temp_tables;
my $quick = 0;
my $verbose = 1;

my $reg = "Bio::EnsEMBL::Registry";

# applied patch to classification db to remove mitochondrial reg. features
#mysql> delete rfc from regulatory_feature rf, regulatory_features_classified rfc  where rf.regulatory_feature_id = rfc.regulatory_feature_id and rf.feature_set_id = 57 and rf.seq_region_id =137;
#Query OK, 102 rows affected (3.42 sec)




my $cdb;
unless($quick){
  $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   #-host => 'ens-staging',#livemirror',#"ensdb-1-12",
   #-host => "ensdb-archive",
   -host => "ensembldb.ensembl.org",
   -port => '5306',
   -dbname => 'homo_sapiens_core_49_36k',
   #-dbname => 'homo_sapiens_funcgen_45_36g',
   -species => "Homo_sapiens",
   #-user => "ensro",
   -user => 'anonymous',
  #-pass => "",
   #-group => 'core',
   #-port => '3304',

  );
}

# new classification database
my $ncdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -host => 'ens-genomics2',#livemirror',#"ensdb-1-12",
   -dbname => 'dk_genomic_features_36k',
   -species => "Homo_sapiens",
   -user => "ensro",
   #-user => 'anonymous',
  #-pass => "",
   -group => '',
   #-port => '3304',
  );



# database to add new classification to
my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host => 'ens-genomics2',#livemirror',#"ensdb-1-12",
   #-host => "ensdb-archive",
   #-dbname => 'homo_sapiens_funcgen_46_36h',
   -dbname => 'dk_homo_sapiens_funcgen_49_36k',
   -species => "Homo_sapiens",
   -dnadb => $cdb,
   -user => "ensadmin",
   #-user => 'anonymous',
   -pass => "ensembl",
   #-group => 'core',
  # -port => '3307',
  );
my $retiring_feature_set_name = 'RegulatoryFeatures_v49';
my $dbu = DBSQL::Utils->new($db->dbc->db_handle);



#my $afa = $db->get_AnnotatedFeatureAdaptor();
#my $pfa = $db->get_PredictedFeatureAdaptor();
#my $slice_a = $db->get_SliceAdaptor();
#my $slice = $slice_a->fetch_by_region('chromosome', 14);#, 23083751, 23084751);#23118316);


# get the new classifications
my $q = "select regulatory_feature_id,feature_type_id from regulatory_features_classified ";
my $aaref = $ncdb->dbc->db_handle->selectall_arrayref($q);
unless(defined $aaref && @$aaref > 0){
    die "got no new classifications from ".$ncdb->dbname;
}

# QUICK
my $new_class_table;
if($quick){
    # put the new classifications in a temp table in the fg db
    $new_class_table = &new_temp();
    $q = "create table $new_class_table select regulatory_feature_id,feature_type_id from regulatory_feature where 1=0";
    &execute($db->dbc->db_handle,$q);
    
    $q = "insert into  $new_class_table values ";
    foreach my $aref (@$aaref){
	$q .= "(".join(",",@$aref)."),";
    }
    chop $q; # remove trailing ,
    &execute($db->dbc->db_handle,$q);
    $q = "alter table $new_class_table add index(regulatory_feature_id)";
    &execute($db->dbc->db_handle,$q);
}
# QUICK

# get access to existing reg feats
my $fsa = $db->get_FeatureSetAdaptor;
my $fset = $fsa->fetch_by_name('RegulatoryFeatures');

# have we tried this update before and failed? if so we need to clean up first
my $new_fset =  $fsa->fetch_by_name('RegulatoryFeatures_new');
if(defined $new_fset){ 
    &commentary( "rolling back feature set RegulatoryFeatures_new".
                 " with feature_set_id ".$new_fset->dbID."\n");


    if($quick){
	my @sql;

	push @sql,"delete  ra from regulatory_attribute ra, regulatory_feature rf where ra.regulatory_feature_id = rf.regulatory_feature_id and rf.feature_set_id = ".$new_fset->dbID;
	push @sql,"delete from regulatory_feature where feature_set_id = ".$new_fset->dbID;
        #we don't need the next line we can just re-use the existing entry
        #push @sql,"delete from feature_set where feature_set_id = ".$new_fset->dbID;
	&execute($db->dbc->db_handle,@sql) or die;
    }else{
	my $Helper=Bio::EnsEMBL::Funcgen::Utils::Helper->new(debug_level => 3);
	$Helper->rollback_FeatureSet($new_fset);
    }


}else{

    $new_fset =  Bio::EnsEMBL::Funcgen::FeatureSet->new(     
                                       #-dbid  => $dbid,
                                       -analysis    => $fset->analysis,
                                       -feature_type => $fset->feature_type,
                                       -cell_type => $fset->cell_type,
                                       -name =>'RegulatoryFeatures_new' ,
                                       -type => $fset->type
			                                    );


    my $ret = $fsa->store($new_fset);
    #print join("\t", @$ret)."\n";# returns an array ?? exit;
    $new_fset =  $fsa->fetch_by_name('RegulatoryFeatures_new');

    #$new_fset = $fsa->fetch_by_dbID($fsa->store($new_fset));
}

print  "existing feature_set_id ".$fset->dbID."\n";
print  "new feature_set_id ".$new_fset->dbID."\n";


#exit;

my @sql=();
#QUICK
if($quick){
    # we create a temp table containing all the current regulatory features
    # and the new classification(ie feature_type_id)
    my $temp_reg_feats = &new_temp();
    push @sql,"create table $temp_reg_feats select rf.*,nc.feature_type_id as new_type from regulatory_feature rf,$new_class_table nc where nc.regulatory_feature_id= rf.regulatory_feature_id and rf.feature_set_id = ".$fset->dbID;
    
    &execute($db->dbc->db_handle,@sql) or die;

    # check that the new classification table and the old have the same number
    # of rows
    my $old_count = $dbu->get_count("select count(distinct regulatory_feature_id)  from $temp_reg_feats");
    my $new_count = $dbu->get_count("select count(distinct regulatory_feature_id)  from $new_class_table");

    unless($old_count == $new_count){
        die "The number of regulatory features differs between the old and".
            "newly classified sets old = $old_count new = $new_count";
    }
    
    # this is a good time to analyse $temp_reg_feats table to see how the
    # classification has changed


    @sql = ();
    push @sql, "drop table $new_class_table";
    # give our newly classified reg_feats their new feature_set_id
    push @sql,"update $temp_reg_feats set feature_set_id = ".$new_fset->dbID;
    # load them without regulatory_feature_ids so they get new ones
    push @sql,"insert into regulatory_feature select NULL , seq_region_id, seq_region_start, seq_region_end, seq_region_strand, display_label, new_type, feature_set_id, stable_id from $temp_reg_feats";
    &execute($db->dbc->db_handle,@sql) or die;

    # now we need to know the mapping between the old and new 
    # regulatory_feature_ids so we can create new regulatory attribute entries
    my $old_new_mapping = &new_temp();
    @sql = ();
    push @sql,"alter table $temp_reg_feats add index(seq_region_id,seq_region_start,seq_region_end)";
    push @sql, "create table $old_new_mapping select o.regulatory_feature_id as old_id,n.regulatory_feature_id as new_id from $temp_reg_feats o, regulatory_feature n where n.feature_set_id = ".$new_fset->dbID." and o.seq_region_id = n.seq_region_id and o.seq_region_start=n.seq_region_start and o.seq_region_end = n.seq_region_end";
    push @sql, "alter table $old_new_mapping add index(old_id)";

    my $new_attributes = &new_temp();
    push @sql, "create table $new_attributes select m.new_id as regulatory_feature_id,ra.attribute_feature_id,ra.attribute_feature_table from $old_new_mapping m, regulatory_attribute ra where ra.regulatory_feature_id = m.old_id";
    push @sql, "insert into regulatory_attribute select regulatory_feature_id,attribute_feature_id, attribute_feature_table from $new_attributes";
    &execute($db->dbc->db_handle,@sql) or die;

    
    
    my $qc_ok = 0;

    $new_count = $dbu->get_count("select count(*) from $new_attributes");
    $old_count = $dbu->get_count("select count(*) from regulatory_attribute ra, regulatory_feature rf where ra.regulatory_feature_id = rf.regulatory_feature_id and rf.feature_set_id = ".$fset->dbID);
    if($new_count == $old_count){
	&commentary("adding $new_count new regulatory attributes\n")if $verbose;
        $qc_ok = 1;
    }else{
        die "The number of regulatory attributes differs between the old and".
            "newly classified sets old = $old_count new = $new_count";
    }


    # if QC looks OK archive the old set by changing its feature_set.name
    if( $qc_ok){
	@sql = ();
	push @sql,"update feature_set set name = '$retiring_feature_set_name' where name = 'RegulatoryFeatures'";
	push @sql,"update data_set set name = '$retiring_feature_set_name' where name = 'RegulatoryFeatures'";

	# set the new regulatory features to be the current set
	push @sql,"update feature_set set name = 'RegulatoryFeatures' where name ='RegulatoryFeatures_new'";
	# create a new entry in the data_set table
	push @sql,"insert into data_set values (NULL,".$new_fset->dbID.", 'RegulatoryFeatures') ";

        &execute($db->dbc->db_handle,@sql) or die;
        #&report_changes($db->dbc->db_handle,$dbu,$temp_reg_feats);
        # select a.feature_type_id,aft.name,b.feature_type_id,bft.name,count(*) from regulatory_feature a, regulatory_feature b, feature_type aft, feature_type bft where a.feature_set_id = 57 and b.feature_set_id = 67 and  a.seq_region_id = b.seq_region_id and a.seq_region_start=b.seq_region_start and a.seq_region_end = b.seq_region_end and a.feature_type_id = aft.feature_type_id and b.feature_type_id = bft.feature_type_id group by a.feature_type_id,b.feature_type_id;

        &clean_temp($db->dbc->db_handle);
    }
#QUICK
}else{
#use API

    # we want a lookup hash for the new feature_types
    # keyed on the regulatory feature_id
    my $fta = $db->get_FeatureTypeAdaptor();
    my %new_type_id;
    # and a hash containing feature_type objects keyed on feature_type_id
    my %ft_objs;
    foreach my $aref (@$aaref){
	$new_type_id{$aref->[0]} = $aref->[1];
	$ft_objs{$aref->[1]} = 1; # collect the type ids
    }
    foreach my $id (keys(%ft_objs)){$ft_objs{$id}=$fta->fetch_by_dbID($id)}


    # we create a new set of regulatory features with the feature_set.name 
    #'RegulatoryFeatures_new' from the existing features, keeping the existing 
    # attributes but changing the regulatory_feature_type to be as indicated in
    # the classification db 
    my $feat_adaptor = $fset->get_FeatureAdaptor();
    #print "got feature_adaptor\n";
    foreach my $dbID(@{$feat_adaptor->list_dbIDs()}){
    # the above gets all the regulatory features - there can be several sets

        my $regf = $feat_adaptor->fetch_by_dbID($dbID);
        #print "regf $regf dbID $dbID feature_set_id ".$regf->feature_set->dbID."\n";
        # only work with the feature_set 'RegulatoryFeatures'
        unless($regf->feature_set->dbID == $fset->dbID){next}

        $regf->feature_set($new_fset);

	#print join(' ',%ft_objs)."\n";
        #print $regf->dbID()."\n";
	#print join(' ',%new_type_id)."\n";
	my $ft = $ft_objs{$new_type_id{	$regf->dbID() }};
	unless(defined $ft){
            die "undefined new feature_type for feature_id".$regf->dbID()."\n";
        }
	#print "new_feature_type ".$ft->dbID."\n";
	$regf->feature_type($ft);
	$regf->dbID(''); 
	$regf->adaptor(''); 
	$feat_adaptor->store($regf);
	
    }

    
    
    my $qc_ok = 0;

    my $new_count =  $dbu->get_count("select count(*) from regulatory_attribute ra, regulatory_feature rf where ra.regulatory_feature_id = rf.regulatory_feature_id and rf.feature_set_id = ".$new_fset->dbID);
    my $old_count = $dbu->get_count("select count(*) from regulatory_attribute ra, regulatory_feature rf where ra.regulatory_feature_id = rf.regulatory_feature_id and rf.feature_set_id = ".$fset->dbID);
    if($new_count == $old_count){
	&commentary("adding $new_count new regulatory attributes\n")if $verbose;
        $qc_ok = 1;
    }else{
        die "The number of regulatory attributes differs between the old and".
            "newly classified sets old = $old_count new = $new_count";
    }


    # if QC looks OK archive the old set by changing its feature_set.name
    if( $qc_ok){
	@sql = ();
	push @sql,"update feature_set set name = '$retiring_feature_set_name' where name = 'RegulatoryFeatures'";
	push @sql,"update data_set set name = '$retiring_feature_set_name' where name = 'RegulatoryFeatures'";

	# set the new regulatory features to be the current set
	push @sql,"update feature_set set name = 'RegulatoryFeatures' where name ='RegulatoryFeatures_new'";
	# create a new entry in the data_set table
	push @sql,"insert into data_set values (NULL,".$new_fset->dbID.", 'RegulatoryFeatures') ";

        &execute($db->dbc->db_handle,@sql) or die;

    }


}

# we update the analysis.created field to reflect the current date.
@sql = ();
push @sql,"update analysis set created = now() where logic_name = 'RegulatoryRegion'";
&execute($db->dbc->db_handle,@sql) or die;

exit;


#######################################################################
# uses global array @temp_tables
sub new_temp{

    my $nom = 'temp_'.$$.'_'.scalar(@temp_tables);
    push @temp_tables,$nom;
    return $nom;

}

# uses global array @temp_tables
sub clean_temp{
    my $dbh = shift;

    my @sql;
    foreach my $table (@temp_tables){
	push @sql,"drop table if exists $table";
    }

    &execute($dbh,@sql) or die;
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
    print  STDERR "$_[0]";
}
