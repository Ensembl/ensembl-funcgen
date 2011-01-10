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

=cut


use strict;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;


my $reg = "Bio::EnsEMBL::Registry";

my $port = $ENV{'EFG_PORT'};
my $user = $ENV{'EFG_READ_USER'};
my $host = $ENV{'EFG_HOST'};
my ($help, $pass, $species, $schema_build, $slice, $num_feats, $skip_counts);
my $usage;

GetOptions( 'help'             => \$help,
            'host|h=s'         => \$host,
            'port=i'           => \$port,
            'user|u=s'         => \$user,
            'pass|p=s'         => \$pass,
			'slice=s'          => \$slice,
            'species|d=s'      => \$species,
			'skip_counts'      => \$skip_counts,
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

	  if($skip_counts){
		print "Skipping counts for regulatory FeatureSet ".$set->name."\n";
	  }
	  else{
		
		@feats = @{$set->get_Features_by_Slice($slice)};
		$num_feats .= scalar(@feats);
		
		#we should count the types here?
		if($slice->name =~ /chromosome/ || @feats){
		  print 'Found '.scalar(@feats).' RegulatoryFeatures on slice '.$slice->name."\n";
		}
	  }
	}
  }
  
  print "Total:\t$num_feats\n" if ! $skip_counts;
}

#now other data_sets
if(defined $slice){
  $slice = $slice_adaptor->fetch_by_name($slice);
}else{
  #default to chromosome 17;
  $slice = $slice_adaptor->fetch_by_region('chromosome', 17);
}

@sets = @{$dset_adaptor->fetch_all()};
my $non_disp_count = 0;
my $warn;
print "Checking DataSets\n";

foreach my $set(@sets){

  if($set->is_displayable()){
	print "\nFound DataSet ".$set->name." with states:\t".join("\t", sort(@{$set->get_all_states()}))."\n";

	my $fset = $set->get_displayable_product_FeatureSet();

	if(defined $fset){
	  print 'FeatureSet '.$fset->name." has states:\t".join("\t", sort(@{$fset->get_all_states()}))."\n";
	  
	  if($skip_counts){
		print "Skipping features counts for product ".ref($fset).' '.$fset->name()."\n";
	  }
	  else{
		@feats = @{$fset->get_Features_by_Slice($slice)};
		$num_feats = scalar(@feats);
		$warn = ($num_feats) ? '' :  "WARNING:\t";
		
		print $warn.'Product '.ref($fset).' '.$fset->name().
		  " has ${num_feats} for test slice ".$slice->name."\n";
	  }
	}
	else{
	  print "WARNING:\tNo DISPLAYABLE FeatureSet\n";
	}

	my @supporting_sets = @{$set->get_displayable_supporting_sets()};

	if(! @supporting_sets){
	  print "WARNING:\tNo DISPLAYABLE supporting sets\n";
	}else{
	  
	  foreach my $sset(@supporting_sets){
		
		if($skip_counts){
		  print "Skipping features counts for supporting ".ref($sset).' '.$sset->name()."\n";
		}
		else{

		  if($sset->isa('Bio::EnsEMBL::Funcgen::ResultSet')){
			@feats = @{$sset->get_ResultFeatures_by_Slice($slice)};
		  }
		  else{
			@feats = @{$sset->get_Features_by_Slice($slice)};
		  }

		  $num_feats = scalar(@feats);
		  $warn = ($num_feats) ? '' :  "WARNING:\t";
		  
		  print $warn.'Supporting '.ref($sset).' '.$sset->name().
			" has ${num_feats} for test slice ".$slice->name."\n";
		}
	  }
	}
  }
  else{
	$non_disp_count++;
  }
}

print "Found $non_disp_count DataSets with no DISPLAYABLE status\n";
