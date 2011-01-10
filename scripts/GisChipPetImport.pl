#!/usr/bin/env perl

=head1 LICENSE

  Copyright (c) 1999-2011 The European Bioinformatics Institute and
  Genome Research Limited.  All rights reserved.

  This software is distributed under a modified Apache license.
  For license details, please see

    http://www.ensembl.org/info/about/code_licence.html

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 CVS

  $Log: not supported by cvs2svn $
  Revision 1.2  2011-01-10 13:40:37  nj1
  updated boiler plate

  Revision 1.1  2008-07-01 13:39:43  dkeefe
  quick hack demonstrating the wrong way to add annotated features to
  the funcgen database.

=cut


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
#use Bio::EnsEMBL::Funcgen::ExperimentalGroup; not yet implemented
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::CellType;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::Utils::Helper;
use File::Path qw(mkpath);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw(mean);
use Devel::Size qw(size total_size);
use Bio::EnsEMBL::Funcgen::ResultFeature;
use Bio::EnsEMBL::Funcgen::ExperimentalSubset;
use Tie::File;
use Data::Dumper;
use DBSQL::Utils;
#use Bio::MAGE::XMLUtils;
#use Bio::MAGE::Experiment::Experiment;

$| = 1; #no output buffer
use strict;

my @temp_tables;
my $quick = 1;
my $verbose = 1;
my @sql;
my $analysis_logic_name = 'Chip_PET_cluster_depth'; 

my %new_analysis = ('-logic_name' => $analysis_logic_name,
		    '-description' => 'approximate sequence depth of PET tags thresholded',
                    '-parameters' => 'gap=100, min_length=34, min_depth=6',
                    '-program' => 'GisChipPetRegionCall.pl'
		    );

my $cell_type = 'hES3';
my $cell_type_description = 'Human embryonic stem cell line hES3 (46XX, Chinese)';

my $group_name = 'Genome Institute of Singapore';
my $group_location = 'Genome Technology and Biology Group, Genome Institute of Singapore, 138672, Singapore';
my $group_contact ='Chia-Lin Wei weicl@gis.a-star.edu.sg';


my $feature_type = 'H3K27me3';
#my $feature_type = 'H3K4me3';

my $experiment_description='http://dx.doi.org  doi:10.1016/j.stem.2007.08.004';
my $experiment_name = $cell_type.'_'.$feature_type;
my $primary_design_type = 'binding_site_identification';

my $experimental_set_format = 'ChIP_PET';
my $experimental_set_name = $experiment_name;
my $experimental_set_vendor = 'not obvious';

my $experimental_sub_set_name = $feature_type.'_hES3_GIS_ChIP_PET.bed';

my %new_feature_type = ();# dummy - the two types above exist

my $feature_set_name = $cell_type.'_'.$feature_type;
my $feature_set_description = $feature_type." histone modifications in  $cell_type cells detected by ChIP_PET";
my $source_table = lc("$feature_type".'_'."$cell_type".'_gis_chip_pet_merged_peaks');


my $reg = "Bio::EnsEMBL::Registry";


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

# chip_pet features database
my $ncdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
  (
   -host => 'ens-genomics2',#livemirror',#"ensdb-1-12",
   -dbname => 'dk_gis_chip_pet',
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
   -dbname => 'dk_homo_sapiens_funcgen_50_36l',
   -species => "Homo_sapiens",
   -dnadb => $cdb,
   -user => "ensadmin",
   #-user => 'anonymous',
   -pass => "ensembl",
   #-group => 'core',
  # -port => '3307',
  );

my $dbu = DBSQL::Utils->new($db->dbc->db_handle);



my $analysis_id = get_or_update_analysis($db,$analysis_logic_name,\%new_analysis);
print "analysis id = $analysis_id\n";
my $cell_type_id = get_or_update_cell_type($db,$cell_type,$cell_type_description);
print "cell type id = $cell_type_id\n";
my $feature_type_id = &get_or_update_feature_type($db,$feature_type,\%new_feature_type);

print "feature type id = $feature_type_id\n";

my $group_id = &get_or_update_experimental_group($db,$group_name,$group_location,$group_contact);

print "experimental_group_id  = $group_id\n";


my $clobber = 0;
# this removes associated experimental_set, experimental_sub_set 
# and supporting_set if $clobber is true
my $experiment_id = &get_or_update_experiment($db,
                                             $experiment_name,
                                             $group_id,
                                             $primary_design_type,
                                             $experiment_description
                                             );

print "experiment_id = $experiment_id\n";


my %new_experimental_set = (-experiment_id => $experiment_id,
                            -feature_type_id => $feature_type_id,
                            -cell_type_id => $cell_type_id,
                            -format =>$experimental_set_format,
                            -name => $experimental_set_name
			    );

my $experimental_set_id = &store_experimental_set($db,\%new_experimental_set);
print "experimental_set_id = $experimental_set_id\n";

my $experimental_sub_set_id =&store_experimental_sub_set($db,$experimental_set_id,$experimental_sub_set_name); 
print "experimental_sub_set_id = $experimental_sub_set_id\n";



if($quick){
    my $fset_id=$db->dbc->db_handle->selectrow_array("select feature_set_id from feature_set where name = '$feature_set_name'");
    if(defined $fset_id && $fset_id > 0){
        if($clobber){
	    @sql = ();
	    push @sql, "delete from data_set where feature_set_id = $fset_id";
	    push @sql, "delete from annotated_feature where feature_set_id = $fset_id";
            push @sql, "delete from feature_set where feature_set_id = $fset_id";
            &execute($db->dbc->db_handle,@sql) or die;
	}else{
            die "A feature set called $feature_set_name already exists ".
	        "with feature_set_id $fset_id\n".
                "set clobber = 1 to remove it\n" ;
        }
    }

    @sql = ();
    push @sql, "insert into feature_set values(NULL,$feature_type_id,$analysis_id,$cell_type_id,'$feature_set_name','annotated','$feature_set_description')";        
    &execute($db->dbc->db_handle,@sql) or die;
    $fset_id=$db->dbc->db_handle->selectrow_array("select feature_set_id from feature_set where name = '$feature_set_name'");
    unless(defined $fset_id && $fset_id > 0){
	die "failed to create entry for $feature_set_name in feature_set_table";
	

    }

    # get the new annotated features
    my $q = "select chr_name,feature_start,feature_end from $source_table";
    my $aaref = $ncdb->dbc->db_handle->selectall_arrayref($q);
    unless(defined $aaref && @$aaref > 0){
        die "got no new features from ".$ncdb->dbname." with query: $q ";
    }

    my $temp1 =  &new_temp();
    #my $new_features_table = &new_temp();
    $q = "create table $temp1 ( seq_region_id int(10),seq_region_start int(11),seq_region_end int(11)) ";
    &execute($db->dbc->db_handle,$q);

    my %seq_id = &get_fg_seq_region_lookup($db);    

    $q = "insert into  $temp1  values ";
    foreach my $aref (@$aaref){
	$q .= "(".$seq_id{ $aref->[0] }.",".$aref->[1].",".$aref->[2]."),";
    }
    chop $q;
    &execute($db->dbc->db_handle,$q) or die;

    @sql = ();
    #push @sql, "create table  $new_features_table select * from annotated_feature where 1=0";
    push @sql, "insert into annotated_feature select NULL,t1.seq_region_id,t1.seq_region_start,t1.seq_region_end,0,NULL,NULL,$fset_id from $temp1 t1";

    push @sql,"insert into data_set values(NULL,$fset_id,'$feature_set_name')";
    &execute($db->dbc->db_handle,@sql) or die;

    # now need to create a supporting_set to link the data_set with the 
    # experimental tables
    $q = "select data_set_id from data_set where feature_set_id = $fset_id";

    my $data_set_id = $db->dbc->db_handle->selectrow_array($q);
    unless(defined $data_set_id && $data_set_id > 0){
	die "failed to get a data set id with feature_set_id = $fset_id";
    }
    print "data_set_id = ". $data_set_id."\n";

    $q = "insert into supporting_set values ($data_set_id,$experimental_set_id,'experimental')";
    &execute($db->dbc->db_handle,$q) or die;
    

}

&clean_temp($db->dbc->db_handle);
exit;


# we update the analysis.created field to reflect the current date.
@sql = ();
push @sql,"update analysis set created = now() where logic_name = 'RegulatoryRegion'";
&execute($db->dbc->db_handle,@sql) or die;

exit;


#######################################################################

# API not yet implemented for experimental_group
sub  get_or_update_experimental_group{
    my($db,$group_name,$group_location,$group_contact)=@_;
    my $dbh = $db->dbc->db_handle;

    my $q = "select experimental_group_id from experimental_group where name = '$group_name'";
    my $aref = $dbh->selectcol_arrayref($q);
    unless(defined $aref){die "query failed : $q\n ".$dbh->errstr}
    if($aref->[0] >0){return $aref->[0]}

    $q = "insert into experimental_group values(NULL,'$group_name','$group_location','$group_contact')";
    &execute($dbh,$q) or die;
    $q = "select experimental_group_id from experimental_group where name = '$group_name'";
    my $id = $dbh->selectrow_array($q);
    unless(defined $id){die "query failed : $q\n ".$dbh->errstr}
    if($id >0){
        return $id;
    }else{
	die "ERROR: got experimental_group_id with value < 1";
    }
}


sub get_fg_seq_region_lookup{
    my($db) = @_;

    my $dbh = $db->dbc->db_handle;
    my $Helper=Bio::EnsEMBL::Funcgen::Utils::Helper->new();
    my $aref=$Helper->get_schema_and_build($db->dbc->dbname);
    my $schema_build = join("_",@$aref);

    my $q = "select sr.name,sr.seq_region_id from seq_region sr,coord_system cs where sr.coord_system_id = cs.coord_system_id and cs.rank = 1 and sr.schema_build = '$schema_build'  and cs.schema_build = '$schema_build' ";

    my $aaref = $dbh->selectall_arrayref($q);
    unless(defined $aaref && @$aaref > 0){
	die "Failed to get chromosome names using query:\n $q";
    }

    my %ret;
    foreach my $aref (@$aaref){
	$ret{$aref->[0]} = $aref->[1];
    }

    return %ret;

}


sub store_experimental_sub_set{
    my($db,$experimental_set_id,$experimental_sub_set_name)=@_;
    # subsets are stored through the ExperimentalSet ie they don't
    # have their own adaptor

    my $a = $db->get_ExperimentalSetAdaptor();
    die "didn't get ExperimentalSetAdaptor adaptor" unless $a;

    my $eset = $a->fetch_by_dbID($experimental_set_id);
    unless(defined $eset){
	die "failed to get experimental_set with dbID $experimental_set_id";
    }

    my $subset = Bio::EnsEMBL::Funcgen::ExperimentalSubset->new(
                 -experimental_set => $eset, 
	         -name => $experimental_sub_set_name
                                                                );
    $a->store_ExperimentalSubsets([$subset]);
    $eset = $a->fetch_by_dbID($experimental_set_id);
    $subset = $eset->get_subset_by_name($experimental_sub_set_name);
    return $subset->dbID();
}


sub store_experimental_set{
    my($db,$new_href)=@_;

    my $a = $db->get_ExperimentalSetAdaptor();
    die "didn't get ExperimentalSetAdaptor adaptor" unless $a;

    $new_href->{-experiment} = 
	$db->get_ExperimentAdaptor->fetch_by_dbID($new_href->{-experiment_id});

    $new_href->{-feature_type} = 
	$db->get_FeatureTypeAdaptor->
	fetch_by_dbID($new_href->{-feature_type_id});

    $new_href->{-cell_type} = 
	$db->get_CellTypeAdaptor->
	fetch_by_dbID($new_href->{-cell_type_id});


    my $new_obj = Bio::EnsEMBL::Funcgen::ExperimentalSet->new(%$new_href);
    $a->store($new_obj);
    $new_obj = $a->fetch_by_name($new_href->{-name});
    return $new_obj->dbID();

}


sub get_or_update_feature_type{
    my ($db,$name,$new_href) = @_;

    my $a = $db->get_FeatureTypeAdaptor();
    die "didn't get feature_type adaptor" unless $a;

    my $feature_type = $a->fetch_by_name($name);

    if(defined $feature_type){return $feature_type->dbID()}

    $feature_type = Bio::EnsEMBL::Funcgen::FeatureType->new(%$new_href);
    $a->store($feature_type);
    $feature_type = $a->fetch_by_name($name);
    return $feature_type->dbID();
}


# this is specific to this script ie not generic
sub get_or_update_experiment{
    my($db, $experiment_name, $group_id, $primary_design_type, 
       $experiment_description)=@_;

    my $a = $db->get_ExperimentAdaptor();
    die "didn't get experiment adaptor" unless $a;

    my $experiment = $a->fetch_by_name($experiment_name);
    if(defined $experiment){
    # we are adding new data so there should not be an experiment with this name    
	if($clobber){ #roll back experiment and any links to it
	    my @sql;
            # is there an experimental_set associated with this expt
	    my $esa = $db->get_ExperimentalSetAdaptor();
            die "didn't get experimental set adaptor" unless $esa;
            # this call returns an array THIS IS NOT AS DOCUMENTED
	    my $eset = $esa->fetch_all_by_Experiment($experiment);
	    print $eset."\n";
            warn "didn't get experimental set " unless $eset->[0];
	    if($eset->[0]){
		$eset=$eset->[0];
	        print "experimental_set_id = ".$eset->dbID."\n";
                # get the experimental subsets in the set
                my @subsets = @{$eset->get_subsets()};
                foreach my $subs (@subsets){
                    push @sql, "delete from experimental_subset where".
		       " experimental_subset_id =".$subs->dbID;
	        }
                # destroy any links in supporting_set
	        push @sql,"delete from supporting_set where type = 'experimental' and supporting_set_id = ".$eset->dbID;
	        push @sql,"delete from experimental_set where experimental_set_id = ".$eset->dbID
	    }

	    
	    push @sql,"delete from experiment where experiment_id = ".$experiment->dbID;
	    &execute( $db->dbc->db_handle, @sql) or die;

	}else{
            die "An entry with name $experiment_name already exists in the ".
                "experiment table. Set clobber = 1 to remove and overwrite it";
	}
    }

=head1 waiting for  GroupAdaptor.pm to be implemented
    $experiment = Bio::EnsEMBL::Funcgen::Experiment->new(
                  -name => $experiment_name,
                  -experimental_group_id => $group_id,
                  -primary_design_type => $primary_design_type,
                  -description => $experiment_description
							 );
    $a->store($experiment);
    $experiment= $a->fetch_by_name($experiment_name);
    return $experiment->dbID();
=cut

    my $q = "insert into experiment values(NULL,'$experiment_name',$group_id,NULL,'$primary_design_type','$experiment_description',NULL)";
    &execute($db->dbc->db_handle,$q) or die;
    $q = "select experiment_id from experiment where name= '$experiment_name'";
    my $id = $db->dbc->db_handle->selectrow_array($q);
    unless(defined $id){die "query failed : $q\n ".$db->dbc->db_handle->errstr}
    if($id >0){
        return $id;
    }else{
	die "ERROR: got experiment_id with value < 1";
    }


}

sub get_or_update_analysis{
    my ($db,$logic_name,$ana_href) = @_;

    my $aa = $db->get_AnalysisAdaptor();
    die "didn't get analysis adaptor" unless $aa;

    my $analysis = $aa->fetch_by_logic_name($logic_name);

    if(defined $analysis){return $analysis->dbID()}

    $analysis = Bio::EnsEMBL::Analysis->new(%$ana_href);
    my $analysis_id = $aa->store($analysis);

    if($analysis_id){
        return $analysis_id;

    }else{
	die "got invalid analysis_id = $analysis_id";
    }
}

#returns dbID of analysis
sub get_or_update_cell_type{
    my ($db,$cell_type_name,$description) = @_;

    my $cta = $db->get_CellTypeAdaptor();
    die "didn't get cell type adaptor" unless $cta;

    my $cell_type = $cta->fetch_by_name($cell_type_name);

    if(defined $cell_type){return $cell_type->dbID()}

    $cell_type = Bio::EnsEMBL::Funcgen::CellType->new(-name => $cell_type_name,
                                             -description => $description);
    $cta->store($cell_type); # doesn't return anything useful

    $cell_type = $cta->fetch_by_name($cell_type_name);

    return $cell_type->dbID();

}









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

    if(@sql >   0){&execute($dbh,@sql) or die}
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

    unless(@array >0){carp("no statements to execute")}

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
