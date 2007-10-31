#!/software/bin/perl
##!/usr/bin/perl

=head1 NAME

run_TileMap.pl - run TileMap on datasets in eFG database

=head1 SYNOPSIS

run_TileMap.pl -dbhost=host -dbport=port -dbuser=user -dbpass=XXXXXX \
	-dbname=homo_sapiens_funcgen_48_36j \
	-logic_name=TileMap -module=TileMap \
	-species homo_sapiens -data_version 46_36h \
	-input_name=21 -verbose

=head1 DESCRIPTION

This script runs TileMap on given eFG database result sets. 
To configure your analysis you need to set the following environment
variables (here bash syntax).

# General config

ANALYSIS_WORK_DIR='/tmp'
EXPERIMENT='ctcf_ren'
NORM_ANALYSIS='VSN_GLOG' # 'SANGER_PCR'
RESULT_SET_REGEXP='_IMPORT'
DATASET_NAME='TileMap'

export EXPERIMENT RESULT_SET_REGEXP NORM_ANALYSIS DATASET_NAME

# TilMap config
TILEMAP_DIR="$HOME/src/tilemap"
TM_LOGIC_NAME='TileMap'
TM_MODULE=$TM_LOGIC_NAME
TM_PROGRAM=$TM_LOGIC_NAME
TM_PROGRAM_FILE='tilemap'
TM_VERSION='2.0'
TM_PARAMETERS="$TILEMAP_DIR/efg/efg_runnable_tilemap_arg.txt"

export TM_LOGIC_NAME TM_MODULE TM_PROGRAM TM_PROGRAM_FILE TM_VERSION TM_PARAMETERS

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf, Ensembl Functional Genomics

=head1 CONTACT

Please post comments/questions to the Ensembl development list
<ensembl-dev@ebi.ac.uk>

=cut

use strict;
use warnings;
use Data::Dumper;
use Getopt::Long;

use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::Utils::Encode qw(get_encode_regions);
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::Utils::Exception qw(verbose throw warning info);
use Bio::EnsEMBL::Analysis::Tools::Logger qw(logger_verbosity logger_info);

$| = 1;

my $dbhost;
my $dbuser;
my $dbpass;
my $dbport;
my $dbname;
my $species;
my $data_version;
my $input_id;
my $input_name;
my $logic_name = $ENV{TM_LOGIC_NAME} || undef;
my $norm_analysis = $ENV{NORM_ANALYSIS} || undef; # 'VSN_GLOG' or 'SANGER_PCR';
my $experiment;
my $check  = 0;
my $output_dir;
my $write = 0;
my $help = 0; 
my $verbose = 0;
my $module;
my $analysis;
my $perl_path = 'Bio/EnsEMBL/Analysis/RunnableDB/Funcgen';
my $utils_verbosity = 'WARNING'; #how verbose do you want the 
my $logger_verbosity = 'OFF'; #how verbose do you want the 
#Bio::EnsEMBL::Utils::Exceptions module to be by default it is set to
#WARNING as this gives warning and throws but not deprecates or infos

my $assembly_version = 'NCBI36';

my @command_args = @ARGV;
&GetOptions( 
	'dbhost=s'      => \$dbhost,
	'dbname=s'      => \$dbname,
	'dbuser=s'      => \$dbuser,
	'dbpass=s'      => \$dbpass,
	'dbport=s'      => \$dbport,
	'species=s'     => \$species,
	'data_version=s'=> \$data_version,
	'input_id=s'    => \$input_id,
	'input_name=s'    => \$input_name,
	'logic_name|analysis=s'  => \$logic_name,
	'norm_analysis=s' => \$norm_analysis,
	'experiment=s' => \$experiment,
	'check'       => \$check,
	'write!' => \$write,
	'help!' => \$help,
	'verbose!' => \$verbose,
	'module=s'    => \$module,
	'runnabledb_path=s' => \$perl_path,
	'utils_verbosity=s' => \$utils_verbosity,
	'logger_verbosity=s' => \$logger_verbosity,
	) or ($help = 1);

$| = 1;

verbose($utils_verbosity);
logger_verbosity($logger_verbosity);

if($check ) {
  print STDERR "args: $dbhost : $dbuser : " . ($dbpass or '') . " : $dbname : $experiment : $input_id : $logic_name\n";
  exit 0;
}


throw("Must provide a dbhost, dbname, and dbuser!")
	if( !$dbhost || !$dbuser || !$dbname );

throw("Must provide a logic_name for analysis!")
	if( !$logic_name );

&usage(\@command_args) if($help);

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => 'ensembldb.ensembl.org',
     -port => 3306,
     -user => 'anonymous',
     -dbname => $species.'_core_'.$data_version,
     -species => $species,
     );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $dbhost,
     -user   => $dbuser,
     -dbname => $dbname,
     -pass   => $dbpass,
     -port   => $dbport,
     -dnadb  => $cdb,
     );

# get slice to work on
my $sa = $db->get_SliceAdaptor();
my $slice;

if ($input_name) {
	my @input_name = split(/[:,]/, $input_name);
    $slice = $sa->fetch_by_region('chromosome', @input_name);
	throw("Input name $input_name didn't return a slice. Must specify mandatory chromosome\n".
		  " name (start and end are optional), like '-input_name=22:1,20000'.") if (! $slice);
	$input_id = $slice->id;
} elsif (defined $ENV{LSB_JOBINDEX}) {
    warn("Performing whole genome analysis on toplevel slices using the farm (LSB_JOBINDEX: ".
         $ENV{LSB_JOBINDEX}.").\n");
    
    my $toplevel = $sa->fetch_all('toplevel');
    my @chr = sort (map $_->seq_region_name, @{$toplevel});
    #print Dumper @chr;
    
    my @slices;
    foreach my $chr (@chr) {
        
        next if ($chr =~ m/^NT_/);
        
        push @slices, $sa->fetch_by_region('chromosome', $chr);
    }
    #print Dumper @slices;

    $slice=$slices[$ENV{LSB_JOBINDEX}-1];      
    print Dumper ($ENV{LSB_JOBINDEX}, $slice->name);
        
	$input_id = $slice->id;
} else {

    throw("Must specify mandatory chromosome name (-input_name) or set\n ".
          "LSF environment variable LSB_JOBINDEX to perform whole\n ".
          "genome analysis on toplevel slices using the farm.\n");

}

throw("No input_id defined!")
    unless (defined $input_id);

print Dumper $input_id;

### setup analysis object
my $aa = $db->get_AnalysisAdaptor;
$analysis = $aa->fetch_by_logic_name($logic_name);

if(!$analysis){

    unless ($module) { 
        throw("The analysis with logic_name \"$logic_name\" is not stored".
              " in the database and the -module option wasn't used.\n". 
              "Either add the analysis to the db or use the -module flag".
              " to specify a module to use.\n".
              "(the analysis will be stored in the db if the features 
              are written)\n" ) ; 
        
    }
    
    print "Creating analysis object ".$logic_name." ".$module. "\n" if $verbose;
    
    print "This object will be stored in the database when the ".
        "features are written\n" if($verbose && $write);
    $analysis = Bio::EnsEMBL::Analysis->new
        (
         -logic_name => $logic_name,
         -module     => $module
         );
    
}

my ($runnable, $file);

if($analysis->module =~ "Bio::"){ 
  $runnable = $analysis->module; 
  ($file = $runnable) =~ s/::/\//g;
}else{
  $file = $perl_path."/".$analysis->module; 
  ($runnable = $file) =~ s/\//::/g;
}
eval{
  require "$file.pm";
};
if($@){
  throw("Couldn't require $file $@");
}
print STDERR "Creating runnable ".$file."\n" if($verbose);


$runnable =~ s/\//::/g;
print Dumper $runnable;
my $runobj = "$runnable"->new(-db            => $db,
                              -input_id      => $input_id,
                              -analysis      => $analysis,
                             );
print STDERR "Instantiated ".$runnable." runnabledb\n" if ($verbose);

$runobj->fetch_input;
print STDERR "Fetched input\n" if($verbose);

$runobj->run unless($runobj->input_is_void);
print STDERR "Run ".$runobj."\n" if($verbose);
print "Input is void not running\n" if($runobj->input_is_void);

if($write){
    $runobj->write_output;
    print STDERR "Written output\n" if($verbose);
} else {
    print STDERR "Output received\n" if($verbose);
    my @output = @{$runobj->output};
    #print Dumper @output;
}

sub usage{
  my ($command_args) = @_;
  print "Your commandline was :\n".
    "test_RunnableDB ".join("\t", @$command_args), "\n\n";
	exec('perldoc', $0);
	exit;
}


1;
