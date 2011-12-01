#!/usr/bin/perl
use strict;
use warnings;
use Bio::EnsEMBL::Registry;

my $registry = 'Bio::EnsEMBL::Registry';

$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous'
);

# 1. Cell Types
# Get the list of all cell types in the Human eFG. For each one print its name, gender and description
# Hint: not all cells have a determined gender

#Grab the eFG adaptor
my $cta = $registry->get_adaptor('Human', 'funcgen', 'celltype');

my @human_cell_types = @{$cta->fetch_all()};
foreach my $cell (@human_cell_types){
	my $gender = $cell->gender || "Undetermined gender"; 
	print $cell->name."\t".$gender."\t".$cell->description."\n";
}

