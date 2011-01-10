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

use warnings;
use strict;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::ResultSet;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (open_file strip_param_args generate_slices_from_names);


#To do
# 1 Integrate into Exporter module
# 2 Genericise this to dump_features, use various format parsers for output
# 3 enable different set dumps, i.e. incorporate get_data.pl

$| =1;

my ($file, $ofile, $pass, $line, $fset_name, $anal_name);
my ($exp_name, $help, $pids, $man, @features, $chr, $not_status);
my ($dbhost, $port, $user, $dbname, $cdbname, $species, @slices, @skip_slices);
my ($dnadb_pass, $dnadb_user, $dnadb_name, $dnadb_port, $dnadb_host);
my $no_zip = 0;
my $out_dir = ".";
my $keep_colons = 0;
my @tmp_args = @ARGV;

GetOptions 
  (
   "feature_set=s"    => \$fset_name,
   "pass=s"           => \$pass,
   'user=s'           => \$user,
   "port=s"           => \$port,
   'species=s'        => \$species,
   "dbname=s"         => \$dbname,
   "dbhost=s"         => \$dbhost,
   
   "dnadb_pass=s"     => \$dnadb_pass,
   'dnadb_user=s'     => \$dnadb_user,
   "dnadb_port=s"     => \$dnadb_port,
   "dnadb_name=s"     => \$dnadb_name,
   "dnadb_host=s"     => \$dnadb_host,
   
   "outdir=s"         => \$out_dir,
   'slices=s{,}'      => \@slices,
   'skip_slices=s{,}' => \@skip_slices,
   'no_zip|z'         => \$no_zip,
   'keep_colons'      => \$keep_colons,
   "help|?"           => \$help,
   "man|m"            => \$man,
  ) 
  or pod2usage( -exitval => 1,
				-message => "Params are:\t@tmp_args"
			  );


pod2usage(0) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


### Set up adaptors and FeatureSet 

if(! $dbname ){
  throw("You must provide a funcgen -dbname paramter");
}
throw("Must define your funcgen dbhost -dbhost") if ! $dbhost;
#throw("Must supply an input file with -file") if ! $file;
#throw("Must supply a password for your fungen db") if ! $pass;

throw("Must pass a feature_set name via -feature_set, i.e. 'RegulatoryFeatures'") if ! $fset_name;

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
													  -host => $dbhost,
													  -dbname => $dbname,
													  -species => $species,
													  -user => $user,
													  -pass => $pass,
													  -port => $port,
													  -dnadb_pass => $dnadb_pass,
													  -dnadb_user => $dnadb_user,
													  -dnadb_port => $dnadb_port,
													  -dnadb_name => $dnadb_name,
													  -dnadb_host => $dnadb_host,
													 );


#Check DB connections
$db->dbc->db_handle;
$db->dnadb->dbc->db_handle;

#Check out_dir
#Use Helper!
system('mkdir -p '.$out_dir) if(! -d $out_dir);


#ChrX  . operon   XXXX YYYY  .  +  . ID=operon01;name=my_operon

my $slice_a = $db->get_SliceAdaptor();
my $fset_a = $db->get_FeatureSetAdaptor();
my $fset = $fset_a->fetch_by_name($fset_name);

if(! defined $fset){
  die("Could not fetch FeatureSet with name:\t$fset_name");
}

print "Dumping GFF features for FeatureSet:\t$fset_name\n";

my ($outline, @output);
my $fset_ftype = ucfirst($fset->feature_class).'Feature';



if(! @slices){
  print "No slices defined defaulting to current toplevel\n";
}
else{
  print "Restricting to slices:\t\t\t".join(', ', @slices)."\n";
}

if(@skip_slices){
  print "Skipping slices:\t\t\t".join(', ', @skip_slices)."\n";
}


@slices = @{&generate_slices_from_names($slice_a, \@slices, \@skip_slices, 1)};#toplevel flag


foreach my $slice(@slices){
  my $cnt = 0;

  print "\nDumping:\t\t\t\t".$slice->name."\n";

  my $seq_name = $slice->seq_region_name();

  if($slice->coord_system->name eq 'chromosome'){
	$seq_name = 'Chr'.$seq_name;
  }

  my $ofile_name = $fset_name.'.'.$seq_name.'.gff';
  
  if(! $keep_colons){
	#Remove colons from file patch which can cause problems
	#with scp
	$ofile_name =~ s/\:/_/go;
  }

  $ofile_name = $out_dir.'/'.$ofile_name;
  my $ofile = open_file($ofile_name, '>', 0775);

  foreach my $feature(@{$fset->get_Features_by_Slice($slice)}){
	$cnt++;
	
	#seqid source type start end score strand phase attrs

	#Let's push this onto an array here and only print OUT every 5000 lines to save I/O?
	#We are no handling strand here

	$outline = join("\t", ($seq_name, $dbname, $fset_ftype, $feature->start(), $feature->end(), '.', '.', '.', 'Name='.$feature->feature_type->name.';'));

	#associated_feature_types?
	#http://www.sequenceontology.org/gff3.shtml

	if($fset_ftype eq 'RegulatoryFeature'){
	  my @attrs;
		
	  foreach my $reg_attr(@{$feature->regulatory_attributes()}){
		
		if($reg_attr->isa('Bio::EnsEMBL::Funcgen::AnnotatedFeature')){

		  #Only need cell type here for MultiCell
		  my $attr_name = $reg_attr->feature_type->name;
		  $attr_name .= ':'.$reg_attr->cell_type->name if $fset->cell_type->name eq 'MultiCell';


		  #We should add MotifFeatures here so they have there AF context.
		  my @pwm_names;

		  foreach my $assoc_mf(@{$reg_attr->get_associated_MotifFeatures()}){
			#mf display_label or binding_matrix id/name here?
			push @pwm_names, $assoc_mf->binding_matrix->name;
		  }

		  if(@pwm_names){
			$attr_name .= '('.join(',', @pwm_names).')';
		  }

		  push @attrs, $attr_name;
		}
		elsif(! $reg_attr->isa('Bio::EnsEMBL::Funcgen::MotifFeature')){
		  #warn we have an unsupported ftype
		  warn "Found unsupported RegulatoryFeature attribute Feature type:\t".ref($reg_attr);
		}
	  }

	  $outline .= join('; ', (' ID='.$feature->stable_id(), 'bound_start='.$feature->bound_start, 
							  'bound_end='.$feature->bound_end, 'Note=Consists of following features: '.
							  join(',', @attrs)));
	}
	elsif($fset_ftype eq 'AnnotatedFeature'){
	  #Add peak summit here?

	  my @pwm_names;
	  
	  foreach my $assoc_mf(@{$feature->get_associated_MotifFeatures()}){
		#mf display_label or binding_matrix id/name here?
		push @pwm_names, $assoc_mf->binding_matrix->name;
	  }

	  if(@pwm_names){
		$outline .= 'Note=Contains the following PWMs:'.join(',', @pwm_names);
	  }
	}				
	elsif($fset_ftype eq 'ExternalFeature'){
	  #It may be more appropriate to have this as the Name for external_features?
	  $outline .= ' Alias='.$feature->display_label;
	}

	push @output, $outline;


	if(scalar(@output) == 1000){
	  print $ofile join("\n", @output)."\n";
	  @output = ();
	}
  }

  print $ofile join("\n", @output);
  @output = ();

  close($ofile);

  #Remove empty files and zip

  if(-z $ofile_name){
	print "No features found\n";
	unlink $ofile_name;
  }
  else{
	print "Features dumped:\t\t\t\t$cnt\n";

	if(! $no_zip){
	  system("gzip -f $ofile_name");
	}
  }
}
