#!/software/bin/perl
##!/usr/bin/env perl

=head1 NAME

run_Runnable.pl - run Funcgen runnables on datasets in eFG database

=head1 SYNOPSIS

run_Runnable.pl -host host -port port -user user -pass XXXXXX \
	-dbname database_name \
	-logic_name logic_name -module module_name \
	-species species -data_version version \
	-input_name chr_name

=head1 DESCRIPTION

This script runs a specified Runnable on given eFG database result sets. 
There is a second wrapper shell script called run_Runnable.sh that is 
supposed to simplify the use of this script. 

To configure your analysis you need to set a couple of environment 
variables as described below in bash syntax. These can be put together 
in config files by runnable that will automatically be sourced by the 
shell script, if the config file is named like ".efg_<runnable>_Runnable".

There is a set of variables that always need to be set for a certain 
analysis, described in the General config section:

    ### General config
    ANALYSIS_WORK_DIR='/tmp'
    EXPERIMENT='experiment_name'
    NORM_ANALYSIS='normalization_method'
    RESULT_SET_REGEXP='regular_expression'
    DATASET_NAME='data_set_name'
    
    export EXPERIMENT RESULT_SET_REGEXP NORM_ANALYSIS DATASET_NAME
    
    # Program specific config
    PROGRAM_DIR='path/to/program/directory'
    PATH="$PROGRAM_DIR:$PATH"
    LOGIC_NAME='logic_name'
    MODULE='module'
    PROGRAM='program_name'
    PROGRAM_FILE='program_filename'
    VERSION='version'
    
    export PROGRAM_DIR PATH LOGIC_NAME MODULE PROGRAM \
           PROGRAM_FILE VERSION 
    
Depending on the used runnable you need to specify further parameters to 
run the analysis. Below are listed the variables by runnable that are 
currently needed.

    ################################
    ### ACME specific parameters ###
    ################################
    # window size (usually 2-3 times the expected fragment size from the 
    # experiment and large enough to include about 10 probes, at least), and 
    WINDOW=1000
    # a threshold, which will be used to determine which probes are counted 
    # as positive in the chi-square test.
    THRESHOLD=0.95
    export WINDOW THRESHOLD
    
    ###########################
    ### Chipotle parameters ###
    ###########################
    PARAMETERS='-windowSize=>990, -stepSize=>198, -adjustPvalue=>BH, -alpha=>0.05'
    export PARAMETERS
    
	####################################
    ### Splitter specific parameters ###
	####################################
    #available Splitter parameters (for details see http://zlab.bu.edu/splitter)
    # -Input already sorted by genomic positions of probes
    #  'sorted=yes|no'
    # -Perform normalization between replicates
    #  'norm=none|qnorm|rnorm'
    # -Combine replicates
    #  'combine=mean|median'
    # -Signal cutoff for the combined intensities>=
    #  'cutoff=splitter|5pct|1pct|2sd|2.5sd'
    #   'splitterfrom='
    #   'splitterto='
    #   'splitterstep=10'
    #   'splitterratio=2'
    # -Determine how the probes are clustered
    #  Gap (Maxgap) <= n base pairs
    #   'maxgap=100'
    #  Clustering (Minrun) >= n probes
    #   'minrun=4'
    PARAMETERS="'norm=none' 'combine=median' 'cutoff=2sd' 'maxgap=100' 'minrun=4'"
    export PARAMETERS
    
	##################################
    ### TileMap specific paramters ###
	##################################
    # default config file used as template to write actual config file
    PARAMETERS="$PROGRAM_DIR/readme/sample2_tilemap_arg.txt"
    # Method to combine neighboring probes (0:HMM, 1:MA)
    METHOD=0
    # Posterior probability > val (default: 0.5)
    POSTPROB=0.5
    # Maximal gap allowed (1000: default)
    MAXGAP=250
    # Expected hybridization length (28: default)
    HYBLENGTH=28
    export PARAMETERS METHOD POSTPROB MAXGAP HYBLENGTH
    
Happy data analysis!

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

my ($host, $user, $pass, $port, $dbname, $species, $data_version);
my $input_id;
my $input_name;
my $input_chip;
my $logic_name = $ENV{LOGIC_NAME} || undef;
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
	'host=s'         => \$host,
	'dbname=s'       => \$dbname,
	'user=s'         => \$user,
	'pass=s'         => \$pass,
	'port=s'         => \$port,
	'species=s'      => \$species,
	'data_version=s' => \$data_version,
	'input_id=s'     => \$input_id,
	'input_name=s'   => \$input_name,
	'input_chip=s'   => \$input_chip,
	'logic_name|analysis=s'  => \$logic_name,
	'norm_analysis=s' => \$norm_analysis,
	'experiment=s'   => \$experiment,
	'check'          => \$check,
	'write!'         => \$write,
	'help!'          => \$help,
	'verbose!'       => \$verbose,
	'module=s'       => \$module,
	'runnabledb_path=s' => \$perl_path,
	'utils_verbosity=s' => \$utils_verbosity,
	'logger_verbosity=s' => \$logger_verbosity,
	) or ($help = 1);

$| = 1;

print "database: ", $dbname, "\n";

verbose($utils_verbosity);
logger_verbosity($logger_verbosity);

### check options ###
throw("Must specify mandatory database hostname (-host).\n") if (!$host);
throw("Must specify mandatory database username. (-user)\n") if (!$user);
throw("Must specify mandatory database name (-dbname).\n") if (!$dbname);
throw("Must specify mandatory database data version, like 47_36i (-data_version).\n") 
    if (!$data_version);
$ENV{DATA_VERSION}=$data_version;
throw("Must specify mandatory logic name for analysis (-logic_name).\n") 
	if (!$logic_name);

&usage(\@command_args) if($help);

my $cdb = Bio::EnsEMBL::DBSQL::DBAdaptor->new
    (
     -host => 'ensembldb.ensembl.org',
     -port => 3306,
     -user => 'anonymous',
     #-host => '127.0.0.1',
     #-port => 33064,
     #-user => 'ensro',
     -dbname => $species.'_core_'.$data_version,
     -species => $species,
     );

my $db = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new
    (
     -host   => $host,
     -user   => $user,
     -dbname => $dbname,
     -pass   => $pass,
     -port   => $port,
     -dnadb  => $cdb,
     );

# get slice to work on
my $sa = $db->get_SliceAdaptor();
my $slice;

if ($ENV{MODULE} eq 'MAT' && 
    (defined $input_chip || defined $ENV{LSB_JOBINDEX})) {

    $input_id = $input_chip || $ENV{LSB_JOBINDEX};

    warn("Performing analysis on all $ENV{NOCHIPS} chips on the farm (LSB_JOBINDEX: ".
         $ENV{LSB_JOBINDEX}.").\n") if (defined $ENV{LSB_JOBINDEX});

} elsif ($input_name) {

	my @input_name = split(/[:,]/, $input_name);
    $slice = $sa->fetch_by_region('chromosome', @input_name);
	throw("Input name $input_name didn't return a slice. Must specify mandatory chromosome\n".
		  " name (start and end are optional), like '-input_name=22:1,20000'.") if (! $slice);
	$input_id = $slice;

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
    #print Dumper ($ENV{LSB_JOBINDEX}, $slice->name);
    
    $input_id = $slice;

} else {
    
    throw("Must specify mandatory chromosome name (-input_name) or chip\n".
          "number (-input_chip) or set LSF environment variable LSB_JOBINDEX\n".
          "to perform whole genome / chip set analysis using the farm.\n");

}

throw("No input_id defined!")
    unless (defined $input_id);

print 'input_id: ', $input_id, "\n";

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
print 'runnable: ', $runnable, "\n";
my $runobj = "$runnable"->new(
                              -db            => $db,
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
	print scalar(@output), " features to annotate.\n";
	#print Dumper @output;
}

sub usage{
  my ($command_args) = @_;
  print "Your commandline was :\n".
    $0.join("\t", @$command_args), "\n\n";
	exec('perldoc', $0);
	exit;
}

1;
