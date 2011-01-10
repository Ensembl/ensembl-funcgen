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

=head1 NAME

ensembl-efg convert_htilist_to_features.pl
  
=head1 SYNOPSIS

convert_hitlist_to_features.pl [options]

Options:

Mandatory


Optional


=head1 OPTIONS

=over 8

=item B<-name|n>

Mandatory:  Instance name for the data set, this is the directory where the native data files are located

=item B<-format|f>

Mandatory:  The format of the data files e.g. nimblegen

=over 8

=item B<-group|g>

Mandatory:  The name of the experimental group

=over 8

=item B<-data_root>

The root data dir containing native data and pipeline data, default = $ENV{'EFG_DATA'}

=over 8

=item B<-fasta>

Flag to turn on dumping of all probe_features in fasta format for the remapping pipeline

=item B<-norm>

Normalisation method, deafult is the Bioconductor vsn package which performs generalised log ratio transformations

=item B<-species|s>

Species name for the array.

=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> takes a input redundant probe name fasta file and generates an NR probe dbID fasta file.

=cut


#add @INC stuff here, or leave to .bashrc/.efg?

BEGIN{
  if (! defined $ENV{'EFG_DATA'}) {
	if (-f "~/src/ensembl-functgenomics/scripts/.efg") {
	  system (". ~/src/ensembl-functgenomics/scripts/.efg");
	} else {
	  die ("This script requires the .efg file available from ensembl-functgenomics\n".
		   "Please source it before running this script\n");
	}
  }
}
	

#use Bio::EnsEMBL::Root; #Only used for rearrange see pdocs
#Roll own Root object to handle debug levels, logging, dumps etc.

### MODULES ###
use Getopt::Long;
#use Carp;#For dev only? cluck not exported by default Remove this and implement in Helper
use Pod::Usage;
#POSIX? File stuff
use File::Path;
#use Bio::EnsEMBL::Registry;
use Bio::EnsEMBL::Utils::Exception qw( throw warning );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file run_system_cmd backup_file);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::FeatureSet;
use Bio::EnsEMBL::Funcgen::DataSet;
use Bio::EnsEMBL::Funcgen::AnnotatedFeature;
use Bio::EnsEMBL::Analysis;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use strict;

$| = 1;							#autoflush
my ($pass, $dbname, $help, $man, $ftname, $file, $species, $ctname);
my ($format, $set_name, $data_version, $clobber, $is_ucsc, $ctype, $ftype);
my $data_dir = $ENV{'EFG_DATA'};
my $user = $ENV{'EFG_WRITE_USER'};
my $host = $ENV{'EFG_HOST'};
my $port = $ENV{'EFG_PORT'};

#Definitely need some sort of Defs modules for each array?

$main::_debug_level = 0;
$main::_tee = 0;


warn "Hardcoded for FeatureSet type 'annotated'\n";

#Use some sort of DBDefs for now, but need  to integrate with Register, and have put SQL into (E)FGAdaptor?
#Use ArrayDefs.pm module for some of these, class, vendor, format?
#ArrayDefs would also contain paths to data and vendor specific parse methods?

GetOptions (
			"file|f=s"       => \$file,
			"pass|p=s"     => \$pass,
			"port=s"     => \$port,
			"host|h=s"     => \$host,
			"user|u=s"     => \$user,
			"dbname|d=s"   => \$dbname,
			"species=s"    => \$species,
			"help|?"       => \$help,
			"man|m"        => \$man,
		  	"feature_type|t=s" => \$ftname,
			"cell_type|c=s"    => \$ctname,
			"set_name|n=s"   => \$set_name,
			"data_version|s=s" => \$data_version,
			'clobber'          => \$clobber,
			"format=s"       => \$format,
			"is_ucsc=i"         => \$is_ucsc,#need true defined booelan rather than implicit
		   );




pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


#do mandatory params and checkign here
$is_ucsc = (! defined $is_ucsc && $format eq 'Wiggle') ? 1 : 0;





my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
											  -host => 'ensembldb.ensembl.org',
											  -user => 'anonymous',
											  -dbname => $species."_core_".$data_version,
											  -species => $species,
											 );



my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													   -dbname => $dbname,
													   -port   => $port,
													   -pass   => $pass,
													   -host   => $host,
													   -user   => $user,
													   -dnadb  => $cdb,
													 );


my $fset_a = $db->get_FeatureSetAdaptor();
my $dset_a = $db->get_DataSetAdaptor();
my $anal_a = $db->get_AnalysisAdaptor();
my $ft_adaptor = $db->get_FeatureTypeAdaptor();
my $ct_adaptor = $db->get_CellTypeAdaptor();
my $slice_a = $db->get_SliceAdaptor();
my $pfa = $db->get_AnnotatedFeatureAdaptor();

my @fsets = @{$fset_a->fetch_all_by_name($set_name)};
throw ("Found more than one FeatureSet with name:\t$set_name") if (scalar(@fsets) >1 );
my $fset = $fsets[0];

my %analyses = (
				Wiggle  => {(
							 logic_name => 'Wiggle',
							 description => 'Solexa clusters from wiggle data',
							 display_label   => 'Parzen',
						   )},

				HitList => {(
							 logic_name    => 'HitList',
							 description   => 'Custom hit list',
							 display_label => 'HitList',
							)},

			   );



throw("The file format you have specified is not accomodated:\t$format") if(! exists $analyses{$format});

my $anal = Bio::EnsEMBL::Analysis->new(
										-logic_name      => $analyses{$format}->{'logic_name'},
										-db              => 'NULL',
										-db_version      => 'NULL',
										-db_file         => 'NULL',
										-program         => 'NULL',
										-program_version => 'NULL',
										-program_file    => 'NULL',
										-gff_source      => 'NULL',
										-gff_feature     => 'NULL',
										-module          => 'NULL',
										-module_version  => 'NULL',
										-parameters      => 'NULL',
										-created         => 'NULL',
										-description     => $analyses{$format}->{'description'},
										-display_label   => $analyses{$format}->{'display_label'},
										-displayable     => 1,
									   );

$anal_a->store($anal);


if($ftname){
  warn "FeatureType name is $ftname";
  
  $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new(
													  -name => $ftname,
													 );
  
  ($ftype) = @{$ft_adaptor->store($ftype)};
}

if($ctname){
  
  warn "CellType name is $ctname";

  $ctype = Bio::EnsEMBL::Funcgen::CellType->new(
												   -name => $ctname,
												  );
  
  ($ctype) = @{$ct_adaptor->store($ctype)};
}





if($fset && ! $clobber){
  throw("Found pre-existing FeatureSet:\t$set_name\nUse -clobber to overwrite");
}elsif($fset && $clobber){
  my $cs_id = $db->get_FGCoordSystemAdaptor->fetch_by_name('chromosome')->dbID();
  my $sql = 'DELETE from predicted_feature where feature_set_id='.$fset->dbID().' and coord_system_id='.$cs_id;
  $db->dbc->do($sql) || throw('Failed to roll back predicted_features for feature_set_id'.$fset->dbID());
}else{#no fset
  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
												 -name         => $set_name,
												 -analysis     => $anal,
												 -feature_type => $ftype,
												 -cell_type    => $ctype,
												 -feature_class=> 'annotated',
												);

  ($fset) = @{$fset_a->store($fset)};
}

my @tmp = @{$dset_a->fetch_all_by_FeatureSet($fset)};

if(! @{$dset_a->fetch_all_by_FeatureSet($fset)}){
  my $dset = Bio::EnsEMBL::Funcgen::DataSet->new(
												 -feature_set => $fset,
												 -name        => $set_name,
												);

  $dset_a->store($dset);

}


my $infile = open_file($file);
my ($line, $chr, $start, $end, $score, $step, %slices);

my $cnt = 0;

#my $method = 'parse_'.$format;

#&$method;

&parse_Wiggle() if $format eq 'Wiggle';
&parse_HitList() if $format eq 'HitList';


close($infile);

print "Finished loading $cnt $ftname AnnotatedFeatures into FeatureSet $set_name\n";


sub parse_Wiggle{

  while ($line = <$infile>) {
	
	next if $line !~ /^[0-9f]/;

	chomp $line;
	
	if($line =~ /^f/){
	  
	  if(defined $score){
		$score = ($score * $step) / ($end - $start);
		&build_and_store_feature();
	  }
	  
	  if($line =~ /^fixedStep/){
		(undef, undef, $chr, undef, $start, undef, $step) = split/\s+|\=/o, $line;
		$chr =~ s/chr//;
		#warn 'slice is '.$slice_a->fetch_by_region('chromosome', $chr);
		
		$slices{$chr} ||= $slice_a->fetch_by_region('chromosome', $chr);
		
		#warn "Got $chr slice".$slices{$chr};
		$end = $start;
		
	  }else{
		throw("Found unexpected line in wiggle file:\n$line");
	  }
	}else{#got stepped score
	  $score += $line;
	  $end += $step;
	}
  }

  &build_and_store_feature();
}


sub parse_HitList{
  
  while ($line = <$infile>) {
	
	throw("Found unrecognised line:\n$line") if($line !~ /^chr/);
	chomp $line;

	($chr, $start, $end, $score) = split /\t/o, $line;

	$chr =~ s/chr//;
	$slices{$chr} ||= $slice_a->fetch_by_region('chromosome', $chr);


	&build_and_store_feature();

  }

}


sub build_and_store_feature{
  
  if($is_ucsc){
	$start ++;
	$end++;
  }


  my $pf = Bio::EnsEMBL::Funcgen::AnnotatedFeature->new(
														-feature_set => $fset,
														-strand      => 0,
														-start       => $start,
														-end         => $end,
														-score       => $score,
														-slice       => $slices{$chr},
													   );
  $pfa->store($pf);
  $cnt ++;

}
