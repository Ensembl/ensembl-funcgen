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

=head1 NAME

get_sets_for_regulatory_build.pl

=head1 SYNOPSIS

get_sets_for_regulatory_build.pl -s species -efgdbname efg_db_name  [-c cell_type -t type -r release $CONNECTION_PARAMS]

=head1 DESCRIPTION

Lists all feature_sets that should be used for the regulatory build

=head1 OPTIONS

=over

=item B<help>

Gives this help menu

=item B<-s>

Species for the dataset (e.g. homo_sapiens)

=item B<-c>

When specified, only display sets for given cell type (e.g. K562)

=item B<-t>

When specified, only displays sets for given type: 'focus' or 'non-focus'

=item B<-r>

When specified, only displays sets annotated for a given Ensembl release number

=item B<-efgdbhost>

Host for the EFG database (defaults to 'ens-genomics1')

=item B<-efgdbuser>

Read User for the EFG database (defaults to 'ensro')

=item B<-efgdbport>

Port for the EFG database (defaults to 3306)

=item B<-efgdbname>

Name of the EFG database (required)

=item B<-trackdbhost>

Host where the data tracking database is (defaults to 'ens-genomics1')

=item B<-trackdbuser>

Read User of the data tracking database (defaults to 'ensro')

=item B<-trackdbport>

Port of the host where the data tracking database is (defaults to 3306)

=item B<-trackdbname>

Name of the data tracking database (defaults to "efg_data_tracking")

=back

=cut

use warnings;
use strict;
use Bio::EnsEMBL::DBSQL::DBConnection;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Getopt::Long;
use Pod::Usage;
print "$0 @ARGV\n";

my ($species,$cell_type,$type,$release,$help);
my $trackdbhost = 'ens-genomics1';
my $trackdbport = 3306;
my $trackdbuser = 'ensro';
my $trackdbname = 'efg_data_tracking';
my $efgdbhost = 'ens-genomics1';
my $efgdbport = 3306;
my $efgdbuser = 'ensro';
my $efgdbname;

GetOptions (
	    'species|s=s'      => \$species,
	    'cell_type|c=s'    => \$cell_type,
	    'type|t=s'         => \$type,
	    'release|r=s'      => \$release,
	    'trackdbhost=s'    => \$trackdbhost,
	    'trackdbport=s'    => \$trackdbport,
	    'trackdbuser=s'    => \$trackdbuser,
	    'trackdbname=s'    => \$trackdbname,
	    'efgdbhost=s'      => \$efgdbhost,
	    'efgdbport=s'      => \$efgdbport,
	    'efgdbuser=s'      => \$efgdbuser,
	    'efgdbname=s'      => \$efgdbname,
	    "help|h"           => \$help,
	   )  or pod2usage( -exitval => 1 ); #Catch unknown opts

pod2usage(1) if ($help);
pod2usage(1) if (! $species );
pod2usage(1) if (! $efgdbname );
if($type && ($type ne 'focus') && ($type ne 'non_focus')){
  print "type must be 'focus' or 'non_focus'\n";
  exit 1;
}

my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(
						 -user     => $trackdbuser,
						 -dbname   => $trackdbname,
						 -host     => $trackdbhost,
						 );

my $efgdba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
  (
   -host    => $efgdbhost,
   -port    => $efgdbport,
   -user    => $efgdbuser,
   -dbname  => $efgdbname,
   -species => $species,
   -group   => 'funcgen',
  );
if(!$efgdba){ print "Could not connect to EFG DB\n"; return 1; }


my $ea = $efgdba->get_ExperimentAdaptor();
my $isa = $efgdba->get_InputSetAdaptor();
my $dsa = $efgdba->get_DataSetAdaptor();

#get all datasets from efg_data_tracking where release not NULL
my $sth = $dbc->prepare("SELECT `experiment_name`, cell_type, `feature_type`, `release_version`, efgdb_set_name  from `dataset` where `species`='$species' AND `release_version` is not null;");
$sth->execute();
while(my ($exp, $ct, $ft, $rel, $exp_name) = $sth->fetchrow_array()){

  #warn "$exp, $ct, $ft, $rel, $exp_name\n";
  #These vars are all badly named!


  if($release){ next if $rel != $release; }
  if($cell_type){ next if lc($cell_type) ne lc($ct); }
  
  #hack to get around efgdb_set_name bug
  #Need to fix this.
  my $exp_name = $ct."_".$ft."_".$exp;

  if(!$exp_name){ warn "Set for $species $ct $ft $exp does not seem to be in efgdb\n"; next; } 

  my $exp_obj = $ea->fetch_by_name($exp_name);
  if(!$exp_obj){ 
    warn "Experiment not found:\t$exp_name\n"; 
  } 
  else {
    my $found = 0;
    
	foreach my $is ( @{$isa->fetch_all_by_Experiment($exp_obj)} ){
    
	  foreach my $ds (@{$dsa->fetch_all_by_supporting_set($is)}){
		my $fs = $ds->product_FeatureSet();
	
		if($fs){
		  $found = 1;


		  #This only works with the previous imported focus sets
		  #Not new sets from the tracking DB!!


		  if($type){
			last if $fs->is_focus_set ? ($type eq 'non_focus') : ($type eq 'focus') ;
		  }

		  print $fs->name."\n";
		} 
		else {
		  warn "No Feature Set found\n";
		}
      }
    }
    if(!$found){ warn $exp_name." not found"; }
  }
}
$sth->finish();

