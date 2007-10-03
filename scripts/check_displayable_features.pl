#!/usr/local/ensembl/bin/perl -w


use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;


my $reg = "Bio::EnsEMBL::Registry";

my $port = 3306;
my $user = 'ensro';
my $host = 'ens-genomics1';
my ($help, $pass, $species, $schema_build, $slice, $num_feats);
my $usage;

GetOptions( 'help'             => \$help,
            'host|h=s'         => \$host,
            'port=i'           => \$port,
            'user|u=s'         => \$user,
            'pass|p=s'         => \$pass,
			'slice=s'          => \$slice,
            'species|d=s'      => \$species,
			'schema_build|s=s' => \$schema_build,
		  );



my $dnadb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
												-dbname  => "${species}_core_${schema_build}",
												-host    => 'ens-staging',
												-port    => 3306,
												-user    =>  'ensro',
												-species => $species,
											   );

die ("You have not provided sufficient arguments to make a DB connection\n".
	 "-dbname  => ${species}_funcgen_${schema_build}, host   => $host, port   => $port, user   => $user, pass   => $pass") if ! $dnadb;


my $efg_db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
														  -dbname  => "${species}_funcgen_${schema_build}",
														  -host    => $host,
														  -port    => $port,
														  -user    => $user,
														  -pass    => $pass,
														  -species => $species,
														  -dnadb   => $dnadb,
														 );


die ("You have not provided sufficient arguments to make a DB connection\n".
	 "-dbname  => ${species}_funcgen_${schema_build}, host   => $host, port".
	 "  => $port, user   => $user, pass   => $pass") if ! $efg_db;

my $fset_adaptor = $efg_db->get_FeatureSetAdaptor();
my $dset_adaptor = $efg_db->get_DataSetAdaptor();
my $slice_adaptor = $efg_db->get_SliceAdaptor();
my @tl_slices = @{$slice_adaptor->fetch_all('toplevel')};

#reg features first
my @sets = @{$fset_adaptor->fetch_all_by_type('regulatory')};
print "Found ".scalar(@sets)." RegulatoryFeature sets:\n";

my (@states, @feats);

foreach my $set(@sets){
  $num_feats = 0;

  @states = sort(@{$set->get_all_states()});

  print "\n".$set->name.((@states) ? " with states:\t".(join("\t", @states))."\n" : "\n");

  #Now actually test features
  if(! $set->is_displayable()){
	print "Skipping feature counts for non-DISPLAYABLE set\n";
  }
  else{
	
	foreach my $slice(@tl_slices){
	  @feats = @{$set->get_Features_by_Slice($slice)};
	  $num_feats .= scalar(@feats);
	  
	  #we should count the types here?
	  if($slice->name =~ /chromosome/ || @feats){
		print 'Found '.scalar(@feats).' RegulatoryFeatures on slice '.$slice->name."\n";
	  }
	}
  }
  
  print "Total:\t$num_feats\n";
}

#now other data_sets
if(defined $slice){
  $slice = $slice_adaptor->fetch_by_name($slice);
}else{
  #default to chromosome 17;
  $slice = $slice_adaptor->fetch_by_region('chromosome', 17);
}

@sets = @{$dset_adaptor->fetch_all_by_feature_type('regulatory')};
my $non_disp_count = 0;

print "Checking DataSets\n";

foreach my $set(@sets){

  if($set->is_displayable()){
	print "\nFound DataSet ".$set->name." with states:\t".join("\t", sort(@{$set->get_all_states()}))."\n";

	my $fset = $set->get_displayable_product_FeatureSet();

	if(defined $fset){
	  print 'FeatureSet '.$fset->name." has states:\t".join("\t", sort(@{$fset->get_all_states()}))."\n"
	}else{
	  print "WARNING:\tNo DISPLAYABLE FeatureSet\n";
	}

	my @supporting_sets = @{$set->get_displayable_supporting_sets()};

	if(! @supporting_sets){
	  print "WARNING:\tNo DISPLAYABLE suppoorting sets\n";
	}else{
	  
	  foreach my $sset(@supporting_sets){
		
		if($sset->isa('Bio::EnsEMBL::Funcgen::ResultSet')){
		  @feats = @{$sset->get_ResultFeatures_by_Slice($slice)};
		}
		else{
		  @feats = @{$sset->get_Features_by_Slice($slice)};
		}
		$num_feats = scalar(@feats);

		#bah

		print ''.(($num_feats) ? '' : "WARNING:\t").'Supporting '.ref($sset).' '.$sset->name().
		  " has ${num_feats} for test slice ".$slice->name."\n";
	  }
	}
  }
  else{
	$non_disp_count++;
  }
}

print "Found $non_disp_count DataSets with no DISPLAYABLE status\n";
