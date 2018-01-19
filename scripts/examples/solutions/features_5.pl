#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2018] EMBL-European Bioinformatics Institute

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

=cut

use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);


#Grab the eFG adaptor
my $fsa        = $registry->get_adaptor('Human', 'funcgen', 'featureset');
my $vista_fset = $fsa->fetch_by_name('VISTA enhancer set');
my @rf_fsets   = @{$fsa->fetch_all_by_type('regulatory')};

my $reg_feat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'regulatoryfeature');
my $ext_feat_adaptor = $registry->get_adaptor('Human', 'funcgen', 'externalfeature');


my @vista_enhancers = @{$ext_feat_adaptor->fetch_all_by_FeatureSets([$vista_fset])};


print "Got ".scalar(@vista_enhancers)." VISTA Enhancers\n";



# Do we need some summary methods from the MultiCell feature to tell us what the associated classifications are?
my %regfeat_ctype_ftypes;


foreach my $enhancer(@vista_enhancers){
  
  my @reg_feats = @{$reg_feat_adaptor->fetch_all_by_Slice_FeatureSets
					  (
					   $enhancer->feature_Slice,
					   \@rf_fsets
					  )};

  foreach my $reg_feat(@reg_feats){
	
	next if $reg_feat->cell_type->name eq 'MultiCell';

	$regfeat_ctype_ftypes{$reg_feat->cell_type->name}{$reg_feat->feature_type->name} ||=0;
	$regfeat_ctype_ftypes{$reg_feat->cell_type->name}{$reg_feat->feature_type->name} ++;
  }
}


foreach my $ctype(keys %regfeat_ctype_ftypes){

  print $ctype."\n";

  foreach my $ftype(keys %{$regfeat_ctype_ftypes{$ctype}}){
	print "\t".$ftype."\t".$regfeat_ctype_ftypes{$ctype}{$ftype}."\n";
  }
  
}


__END__


>perl features_5.pl
Got 212 VISTA Enhancers
K562
        Promoter Associated     3
        Unclassified    30
NHEK
        Promoter Associated     6
        Gene Associated 5
        Unclassified    29
GM06990 
        Promoter Associated     7
        Gene Associated 5
        Unclassified    12
IMR90   
        Promoter Associated     6
        Gene Associated 4
        Unclassified    28
HSMM
        Promoter Associated     7
        Gene Associated 2
        Unclassified    25
HepG2   
        Promoter Associated     3
        Gene Associated 5
        Unclassified    27
HeLa-S3 
        Promoter Associated     5
        Gene Associated 2
        Unclassified    23
NH-A
        Promoter Associated     6
        Gene Associated 3
        Unclassified    23
CD4
        Promoter Associated     7
        Gene Associated 1
        Non-Gene Associated     1
        Unclassified    2
HMEC
        Promoter Associated     8
        Gene Associated 3
        Unclassified    17
HUVEC   
        Promoter Associated     6
        Gene Associated 7
        Unclassified    26
H1ESC   
        Promoter Associated     9
        Gene Associated 3
        Unclassified    20
GM12878 
        Promoter Associated     4
        Gene Associated 5
        Unclassified    33
