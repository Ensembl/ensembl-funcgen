#!/software/bin/perl

=head1 NAME

run_Nessie.pl - run Nessie (TilingHMM) on datasets in eFG database

=head1 SYNOPSIS

run_Nessie.pl -dbhost=dbhost -dbname=dbname -dbuser=dbuser -dbpass=password \
    -input_name=ENr333 -e H3ac-HeLa -result_set_analysis=SangerPCR \
    -logic_name Nessie_PCRarray

=head1 DESCRIPTION

This script will run Nessie on given eFG database result sets.

=head1 LICENCE

This code is distributed under an Apache style licence. Please see
http://www.ensembl.org/info/about/code_licence.html for details.

=head1 AUTHOR

Stefan Graf <graef@ebi.ac.uk>, Ensembl Functional Genomics

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
my $input_id;
my $input_name;
my $logic_name = 'Nessie';
my $experiment;
my $result_set_analysis;
#my $result_set_analysis = 'SangerPCR';
#my $result_set_analysis = 'VSN_GLOG';
my $check  = 0;
my $output_dir;
my $write = 0;
my $help = 0; 
my $verbose = 0;
my $module = 'Nessie';
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
             'input_id=s'    => \$input_id,
             'input_name=s'    => \$input_name,
             'logic_name|analysis:s'  => \$logic_name,
             'experiment:s'  => \$experiment,
             'result_set_analysis:s' => \$result_set_analysis,
             'check'       => \$check,
             'write!' => \$write,
             'help!' => \$help,
             'verbose!' => \$verbose,
             'module:s'    => \$module,
             'runnabledb_path:s' => \$perl_path,
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

if(!$dbhost || !$dbuser || !$dbname){
    throw("Must provide a dbhost, dbname, and dbuser!");
}

&usage(\@command_args) if($help);


use Bio::EnsEMBL::Registry;
Bio::EnsEMBL::Registry->load_registry_from_db
    (
     -host => 'ens-livemirror',
     -user => 'ensro',
#     -verbose => "1" 
     );
my $cdb = Bio::EnsEMBL::Registry->get_DBAdaptor('human', 'core');

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $dbhost,
     -user   => $dbuser,
     -dbname => $dbname,
     -pass   => $dbpass,
     -port   => $dbport,
     -dnadb  => $cdb,
     );

if (! defined $input_id) {
    warn("No input_id provided! Using input_name ".
         "to select an Encode region.");
    my $encode_regions = &get_encode_regions($cdb, $assembly_version);
    
    if (! defined $input_name) {
        warn("No input_name provided! Using LSB_JOBINDEX ".
             "to select an Encode region.");

        throw("LSF environment variable LSB_JOBINDEX not defined.") 
            if (! defined $ENV{LSB_JOBINDEX});
    
        my @encode_regions = sort keys %{$encode_regions};
        $input_name = $encode_regions[$ENV{LSB_JOBINDEX}-1];
    }
    $input_id = $encode_regions->{$input_name};
    print Dumper ($input_name, $input_id);
}
throw("No input_id defined!")
    unless (defined $input_id);

### check for experiment to be analysed
#print Dumper $experiment;
if (! defined $experiment ) {
    throw("Must provide an experiment name to be analysed.");
}

### check for result_set_analysis
#print Dumper $result_set_analysis;
if (! defined $result_set_analysis ) {
    throw("Must provide an result_set analysis name.");
}

### setup analysis object
my $aa = $db->get_AnalysisAdaptor;
$analysis = $aa->fetch_by_logic_name($logic_name);

if(!$analysis){
    
    unless ($module) { 
        throw("The analysis with logic_name \"$logic_name\" is not stored".
              " in the database and the -module option wasn't used.\n". 
              "Either add the analysis to the db or use the -module flag".
              " to specify a module to use.\n".
              "( the analysis will be stored in the db if the features 
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
#print Dumper $analysis;

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
my $runobj = "$runnable"->new(-db         => $db,
                              -input_id   => $input_id,
                              -analysis   => $analysis,
                              -experiment => $experiment,
                              -result_set_analysis => $result_set_analysis
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
    my @output = @{$runobj->output};
    print Dumper @output;
}

sub usage{
  my ($command_args) = @_;
  print "Your commandline was :\n".
    "test_RunnableDB ".join("\t", @$command_args), "\n\n";
	exec('perldoc', $0);
	exit;
}


1;
