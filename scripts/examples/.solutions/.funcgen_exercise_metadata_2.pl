#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# 2. Feature Types
# Get the list of all feature types in the Mouse eFG. For each one print its name, class and description
# Get the list of all Mouse feature types with class 'Transcription Factor'. How many transcription factors are there?

#Grab the eFG adaptor
my $fta = $registry->get_adaptor('Mouse', 'funcgen', 'featuretype');

my @mouse_feature_types = @{$fta->fetch_all()};
foreach my $feature_type (@mouse_feature_types){
	print $feature_type->name."\t".$feature_type->class."\t".$feature_type->description."\n";
}

my @tfs = @{$fta->fetch_all_by_class('Transcription Factor')};
print scalar(@tfs)." transcription factors in the mouse eFG\n";

