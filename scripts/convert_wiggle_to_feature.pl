#!/software/bin/perl -w


####!/opt/local/bin/perl -w


=head1 NAME

ensembl-efg convert_wiggle_to_feature.pl
  
=head1 SYNOPSIS

make_nr_probe_fasta.pl [options]

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
use Bio::EnsEMBL::Funcgen::PredictedFeature;
use Bio::EnsEMBL::Analysis;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use strict;

$| = 1;							#autoflush
my ($pass, $dbname, $help, $man, $ftname, $file, $species, $set_name, $data_version, $clobber);
#my $reg = "Bio::EnsEMBL::Registry";
my $data_dir = $ENV{'EFG_DATA'};
my $user = "ensadmin";
my $host = 'localhost';
my $port = '3306';

#Definitely need some sort of Defs modules for each array?

$main::_debug_level = 0;
$main::_tee = 0;


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
			"set_name|n=s"   => \$set_name,
			"data_version|s=s" => \$data_version,
			'clobber'          => \$clobber,
		   );




pod2usage(1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;


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
my $slice_a = $db->get_SliceAdaptor();
my $pfa = $db->get_PredictedFeatureAdaptor();

my @fsets = @{$fset_a->fetch_all_by_name($set_name)};
throw ("Found more than one FeatureSet with name:\t$set_name") if (scalar(@fsets) >1 );
my $fset = $fsets[0];

my $wanal = Bio::EnsEMBL::Analysis->new(
										-logic_name      => 'Wiggle',
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
										-description     => 'Solexa clusters from wiggle data',
										-display_label   => 'Parzen',
										-displayable     => 1,
									   );

$anal_a->store($wanal);

warn "ftypename is $ftname";

my $ftype = Bio::EnsEMBL::Funcgen::FeatureType->new(
												  -name => $ftname,
												 );

($ftype) = @{$ft_adaptor->store($ftype)};

if($fset && ! $clobber){
  throw("Found pre-existing FeatureSet:\t$set_name\nUse -clobber to overwrite");
}elsif($fset && $clobber){
  my $sql = 'DELETE from predicted_feature where feature_set_id='.$fset->dbID();
  $db->dbc->do($sql) || throw('Failed to roll back predicted_features for feature_set_id'.$fset->dbID());
}else{#no fset
  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
												 -name         => $set_name,
												 -analysis     => $wanal,
												 -feature_type => $ftype,
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


my $wfile = open_file($file);
my ($line, $chr, $start, $end, $score, $step, $pf, %slices);

my $cnt = 0;

while ($line = <$wfile>) {

  next if $line !~ /^[0-9f]/;

  chomp $line;

  if($line =~ /^f/){

	if(defined $score){
	  $score = ($score * $step) / ($end - $start);

	  $pf = Bio::EnsEMBL::Funcgen::PredictedFeature->new(
														 -feature_set => $fset,
														 -strand      => 0,
														 -start       => ($start+1),
														 -end         => ($end+1),
														 -score       => $score,
														 -slice       => $slices{$chr},
														 #-display_label
														);
	  $pfa->store($pf);
	  $cnt ++;
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


$pf = Bio::EnsEMBL::Funcgen::PredictedFeature->new(
												   -feature_set => $fset,
												   -strand      => 0,
												   -start       => ($start+1),
												   -end         => ($end+1),
												   -score       => $score,
												   -slice       => $slices{$chr},
												  );
$pfa->store($pf);
$cnt ++;


close($wfile);

print "Finished loading $cnt $ftname PredictedFeatures into FeatureSet $set_name\n";
