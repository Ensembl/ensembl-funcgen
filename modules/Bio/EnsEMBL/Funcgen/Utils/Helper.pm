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

Bio::EnsEMBL::Funcgen::Utils::Helper
  
=head1 SYNOPSIS


 e.g. 


 my $object = Bio::EnsEMBL::Object->new
 (
     logging     => 1,
     log_file    => "/tmp/Misc.log",
     debug_level => 2,
     debug_file  => "/tmp/Misc.dbg",
 );

 $object->log("This is a log message.");
 $object->debug(1,"This is a debug message.");
 $object->system("rmdir /tmp/test");


 ----------------------------------------------------------------------------


=head1 OPTIONS

=over 8


=item B<-debug>

Turns on and defines the verbosity of debugging output, 1-3, default = 0 = off

=over 8

=item B<-log_file|l>

Defines the log file, default = "${instance}.log"

=item B<-help>

Print a brief help message and exits.

=item B<-man>

Prints the manual page and exits.

=back

=head1 DESCRIPTION

B<This program> performs several debugging and logging functions, aswell as providing several inheritable EFGUtils methods.

=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Utils::Helper;

use Bio::Root::Root;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw (throw stack_trace);
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw (get_date);
#use Devel::Timer;
use Carp;#? Can't use unless we can get it to redirect
use File::Basename;


use strict;
use vars qw(@ISA);
@ISA = qw(Bio::Root::Root);

#List of valid rollback levels
#To be used in conjunction with -full_delete
my @rollback_tables = ('data_set', 'feature_set', 'result_set', 'input_set', 'experiment', 'array', 'array_chip', 'experimental_chip');

#Some local filevars to avoid assigning to package typeglobs
my ($DBGFILE, $LOGFILE);

################################################################################

=head2 new

 Description : Constructor method to create a new object with passed or
               default attributes.

 Arg  [1]    : hash containing optional attributes :-
                 log_file    - name of log file (default = undef -> STDOUT)
                 debug_level - level of detail of debug message [1-3] (default = 0 = off)
                 debug_file  - name of debug file (default = undef -> STDERR)

 ReturnType  : Helper

 Example     : my $Helper = new Bio::EnsEMBL::Helper(
                                                      debug_level => 3,
                                                      debug_file  => "/tmp/efg.debug",
                                                      log_file    => "/tmp/efg.log",
                                                     );

 Exceptions  : throws exception if failed to open debug file
             : throws exception if failed to open log   file

=cut

################################################################################

#To do , change to rearrange

sub new{
    my ($caller, %args) = @_;

    my ($self, %attrdata, $argname);
    my $class = ref($caller) || $caller;

    #Create object from parent class
    $self = $class->SUPER::new(%args);

	#we need to mirror ensembl behaviour here
	#use rearrange and set default afterwards if not defined

    # objects private data and default values
	#Not all of these need to be in main

    %attrdata = (
				 _tee          => $main::_tee,
				 _debug_level  => $main::_debug_level,
				 _debug_file   => $main::_debug_file,
				 _log_file     => $main::_log_file,#default should be set in caller
				 _no_log       => $main::_no_log,#suppresses log file generation if log file not defined
				 _default_log_dir => $main::_default_log_dir,
		);

    # set each class attribute using passed value or default value
    foreach my $attrname (keys %attrdata){
	  ($argname = $attrname) =~ s/^_//; # remove leading underscore
	  $self->{$attrname} = (exists $args{$argname}) ? $args{$argname} : $attrdata{$attrname};
    }


	$self->{'_tee'} = 1 if $self->{'_no_log'};
	#should we undef log_file here too?
	#This currently only turns off default logging

	$self->{_default_log_dir} ||= $ENV{'HOME'}.'/logs';
	$self->{'_report'} = [];
	

    # DEBUG OUTPUT & STDERR

	#should default to lowest or highest debug level here!

    if(defined $self->{_debug_level} && $self->{_debug_level}){
        $main::_debug_level = $self->{_debug_level};
		
        if(defined $self->{_debug_file}){
			$main::_debug_file = $self->{_debug_file};
			  			  
            open($DBGFILE, '>>', $self->{_debug_file})
			  or throw("Failed to open debug file : $!");

			#open (DBGFILE, "<STDERR | tee -a ".$self->{_debug_file});#Mirrors STDERR to debug file
        }
        else{
            open($DBGFILE, '>&STDERR');
        }

        select $DBGFILE; $| = 1;  # make debug file unbuffered

        $self->debug(1,"Debugging started ".localtime()." on $0 at level ".$self->{_debug_level}." ...");
    }

	my $log_file =  $self->{_log_file};

	
	# LOG OUTPUT
	if (defined $self->{_log_file}){

	  #This causes print on unopened file as we try and log in the DESTROY
	  throw('You have specified mutually exclusive parameters log_file and no_log') if($self->{'_no_log'});
	  $main::_log_file = $self->{_log_file};
	  
	  #we need to implment tee here
	  if($self->{'_tee'}){
	    open($LOGFILE, ' | tee -a '.$log_file);
	  }
	  else{
		open($LOGFILE, '>>', $log_file)
		  or throw("Failed to open log file : $log_file\nError: $!");
	  }
	}
	else{
	  #Change this to get the name of the control script and append with PID.out
	  #This is to ensure that we always capture output
	  #We need to also log params
	  #We will have to call this from the child class.
	  

	  #Only do this if we don't have supress default logs set
	  #To avoid loads of loags during testing
	  if(! $self->{'_no_log'}){

		my @stack = stack_trace();
		my $top_level = $stack[$#stack];
		my (undef, $file) = @{$top_level};
		$file =~ s/.*\///;

		$self->run_system_cmd('mkdir '.$self->{_default_log_dir}) if(! -e $self->{_default_log_dir});
		$self->{'_log_file'} = $self->{_default_log_dir}.'/'.$file.'.'.$$.'.log';
		warn "No log file defined, defaulting to:\t".$self->{'_log_file'}."\n";

		#we should still tee here
		if($self->{'_tee'}){
		  open($LOGFILE, '| tee -a '.$self->{'_log_file'});
		}
		else{
		  open($LOGFILE, '>', $self->{'_log_file'})
			or throw('Failed to open log file : '.$self->{'_log_file'}."\nError: $!");
		}
		
	  }
	  else{
		#Have to include STD filehandles in operator
		open($LOGFILE, '>&STDOUT');
	  }
	}

	select $LOGFILE; $| = 1;  # make log file unbuffered
	$self->log("\n\nLogging started at ".localtime()."...");

    # RESET STDOUT TO DEFAULT
    select STDOUT; $| = 1; 

    $self->debug(2,"Helper class instance created.");

    return $self;
}


################################################################################

=head2 DESTROY

 Description : Called by gargbage collection to enable tidy up before object deleted

 ReturnType  : none

 Example     : none - should not be called directly

 Exceptions  : none

=cut

################################################################################

sub DESTROY{
    my ($self) = @_;


	$self->report;

    if($self->{_log_file}){
        $self->log("Logging complete ".localtime().".");
		$self->log('Virtual Memory '.`ps -p $$ -o vsz |tail -1`);
		$self->log('Resident Memory '.`ps -p $$ -o rss |tail -1`);
		
		
				   

		#       close LOGFILE;  # if inherited object then cannot close filehandle !!!
    }

    if($self->{_debug_level}){
        $self->debug(1,"Debugging complete ".localtime().".");
		#       close DBGFILE;  # if inherited object then cannot close filehandle !!!
    }

	if(defined $self->{'_timer'}){
		$self->{'_timer'}->report();
	}

	$self->debug(2,"Bio::EnsEMBL::Helper class instance destroyed.");

    return;
}




##Need generic method in here to get stack and line info
###Use Root.pm stack methods!
# and replace this with caller line method for logging
sub _get_stack{
  my ($self) = shift;
  

  #need to resolve this method with that in debug, pass log or debug arg for different format

  my @prog = (caller(2)) ? caller(2) : (caller(1)) ? caller(1) : (undef,"undef",0);

  return "[".localtime()." - ".basename($prog[1]).":$prog[2]]";
}


################################################################################

=head2 log

 Arg[0]      : string  - log message.
 Arg[1]      : boolean - memory usage, appends current process memory stats
 Description : Method to write messages to a previously set up log file.
 Return type : none
 Example     : $root->log("Processing file $filename ...", 1);
 Exceptions  : none

=cut

################################################################################

sub log{
  my ($self, $message, $mem, $date, $no_return) = @_;

  if($mem){
	$message.= " :: ".`ps -p $$ -o vsz |tail -1`;
	chomp $message;
	$message .= " KB";
  }
  
  if($date){
	my $time = localtime();
	chomp($time);
	$message .= ' - '.localtime();
  }

  $message .= "\n" if ! $no_return;

  print $LOGFILE "::\t$message";

  # Add to debug file if not printing to STDERR?
  # only if verbose?
  # this would double print everything to STDOUT if tee and debug has not redefined STDERR

  $self->debug(1,$message);
}

################################################################################


=head2 report

 Arg[0]      : optional string  - log message.
 Arg[1]      : optional boolean - memory usage, appends current process memory stats
 Description : Wrapper method for log, which also stores message for summary reporting
 Return type : none
 Example     : $root->report("WARNING: You have not done this or that and want it reported at the end of a script");
 Exceptions  : none

=cut

################################################################################

sub report{
  my ($self, $message, $mem) = @_;

  if(defined $message){

	$self->log($message, $mem);

	push @{$self->{'_report'}}, $message;
  }
  elsif(scalar(@{$self->{'_report'}})){
	print $LOGFILE "\n::\tSUMMARY REPORT\t::\n";
	print $LOGFILE join("\n", @{$self->{'_report'}})."\n";

	$self->{'_report'} = [];
  }

  return;
}






################################################################################

=head2 log_header

 Arg[0]      : string  - log message.
 Arg[1]      : boolean - memory usage, appends current process memory stats
 Description : Wrapper method to format a log as a header line
 Return type : none
 Example     : $root->log("Processing file $filename ...", 1);
 Exceptions  : none

=cut

################################################################################

sub log_header{
  my ($self, $message, $mem, $date) = @_;

  print $LOGFILE "\n\n";
  $self->log("::\t$message\t::\t::", $mem, $date);
  print $LOGFILE "\n";
}





################################################################################

=head2 debug

 Description : Method to write debug info to a previously set up debug file.
               Over-rides Root.pm on/off style debugging

 Args        : int: debug level and string: log message.

 ReturnType  : none

 Example     : $root->debug(2,"dir=$dir file=$file");

 Exceptions  : none

=cut

################################################################################

sub debug{
    my ($self,$level,$message) = @_;



    #Can we not detect whther message is a scalar, array or hash and Dump or print accordingly?

    my (@call,$cnt,$prog_name,$prog_line,$call_name,$call_line);

    $prog_name = $call_name = "undef";
    $prog_line = $call_line = $cnt = 0;

    # if debug on at the requested level then output the passed message
    if (defined $self->{_debug_level} && $level <= $self->{_debug_level}){

		######Replace this with Carp method?
        while (@call = caller($cnt++)){

            if ($cnt == 2){
                $call_name = basename($call[1]);
                $call_line = $call[2]
            }
            
            $prog_name = basename($call[1]);
            $prog_line = $call[2];
        }
           
		#This still attempts to print if file not opened
        print $DBGFILE "debug $message\t: [$$ - $prog_name:$prog_line  $call_name:$call_line]\n";

		#carp("carping $message");
    }
}


################################################################################

=head2 debug_hash

 Description : Method to write the contents of passed hash to debug output.

 Args        : int: debug level and hashref.

 ReturnType  : none

 Example     : $Helper->debug_hash(3,\%hash);

 Exceptions  : none

=cut

################################################################################

sub debug_hash{
    my ($self,$level,$hashref) = @_;
    
    my ($attr);
    
    # if debug on at the requested level then output the passed hash
    if (defined $self->{_debug_level} && $level <= $self->{_debug_level}){
		print $DBGFILE Data::Dumper::Dumper(\$hashref)."\n";
	}
}



################################################################################

=head2 run_system_cmd

 Description : Method to control the execution of the standard system() command

 ReturnType  : none

 Example     : $Helper->debug(2,"dir=$dir file=$file");

 Exceptions  : throws exception if system command returns none zero

=cut

################################################################################


#Move most of this to EFGUtils.pm
#Maintain wrapper here with throws, only warn in EFGUtils

sub run_system_cmd{
  my ($self, $command, $no_exit) = @_;

  my $redirect = '';

  $self->debug(3, "system($command)");
  
  # decide where the command line output should be redirected

  #This should account for redirects
  #This just sends everything to 1 no?

  if (defined $self->{_debug_level} && $self->{_debug_level} >= 3){

    if (defined $self->{_debug_file}){
      $redirect = " >>".$self->{_debug_file}." 2>&1";
    }
    else{
      $redirect = "";
    }
  }
  else{
    #$redirect = " > /dev/null 2>&1";
  }

  # execute the passed system command
  my $status = system("$command $redirect");
  my $exit_code = $status >> 8; 
 
  if ($status == -1) {	
	warn "Failed to execute: $!\n";
  }    
  elsif ($status & 127) {
	warn sprintf("Child died with signal %d, %s coredump\nError:\t$!",($status & 127),($status & 128) ? 'with' : 'without');
  }    
  elsif($status != 0) {	
	warn sprintf("Child exited with value %d\nError:\t$!\n", $exit_code); #get the true exit code
  }
 
  if ($exit_code != 0){
		  
    if (! $no_exit){
      throw("System command failed:\t$command\nExit code:\t$exit_code\n$!");
    }
    else{
      warn("System command returned non-zero exit code:\t$command\nExit code:\t$exit_code\n$!");
    }
  }
  
  #reverse boolean logic for perl...can't do this anymore due to tab2mage successful non-zero exit codes :/

  return $exit_code;
}


#add sys_get method ehre to handle system calls which retrieve data?
#i.e.backtick commands `find . -name *fasta`
#or use want or flag with above method?
#should open pipe instead to capture error?

sub get_data{
  my ($self, $data_type, $data_name) = @_;

  #This method is just to provide standard checking for specific get_data/config methods

  if(defined $data_name){
    throw("Defs data name $data_name for type '$data_type' does not exist\n") if (! exists $self->{"${data_type}"}{$data_name});
  }else{
    throw("Defs data type $data_type does not exist\n") if (! exists $self->{"${data_type}"});
  }
  
  return (defined $data_name) ? $self->{"${data_type}"}{$data_name} : $self->{"${data_type}"};
}


#sub Timer{
#	my ($self) = shift;

#	$self->{'_timer'} = new Devel::Timer()  if(! defined $self->{'_timer'});

#	return $self->{'_timer'};
	
#}


sub set_header_hash{
  my ($self, $header_ref, $fields) = @_;
	
  my %hpos;

  for my $x(0..$#{$header_ref}){
    $hpos{$header_ref->[$x]} = $x;
  }	


  if($fields){

    foreach my $field(@$fields){
	  
      if(! exists $hpos{$field}){
	throw("Header does not contain mandatory field:\t${field}");
      }
    }
  }
  
  return \%hpos;
}

#Move this to EFGUtils?

sub backup_file{
  my ($self, $file_path) = @_;

  throw("Must define a file path to backup") if(! $file_path);

   if (-f $file_path) {
	 $self->log("Backing up:\t$file_path");
    system ("mv ${file_path} ${file_path}.".`date '+%T'`);
  }

  return;

}

#This should move to Utils
#as it is a simple string manipulation

sub get_schema_and_build{
  my ($self, $dbname) = @_;
  my @dbname = split/_/, $dbname;
  return [$dbname[($#dbname -1)], $dbname[($#dbname )]];
}

=head2 get_regbuild_set_states

  Arg [1]    : Bio::EnsEMBL::DBAdaptor
  Example    : my ($dset_states, $rset_states, $fset_states) = $helper->get_regbuild_set_states($db);
  Description: Returns Array refs of appropriate states for sets use din the regulatory build
  Returntype : Array
  Exceptions : Warns if cannot find chromosome CoordSystem
  Caller     : HealthChecker & regulatory build code
  Status     : At risk

=cut


sub get_regbuild_set_states{
  my ($self, $db) = @_;
  
  my $cs_a = $db->get_CoordSystemAdaptor;

  #These states need to be mirrored in RegulatorySets.java

  my $chrom_cs = $cs_a->fetch_by_name('chromosome');
  my (@dset_states, @rset_states, @fset_states);

  if(! defined $chrom_cs){
	#This species most likely does not have a regbuild
	#really just need to get the 'highest' level here
	warn "Could not find Chromosome CoordSystem. ".$db->dbc->dbname.". most likely does not contain a RegulatoryBuild";
  }
  else{
	my $imp_cs_status = 'IMPORTED_'.$cs_a->fetch_by_name('chromosome')->version;
		
	#What about non-chromosome assemblies?
	#top level will not return version...why not?
	@dset_states = ('DISPLAYABLE');
	@rset_states = (@dset_states, 'DAS_DISPLAYABLE', $imp_cs_status);
	@fset_states = (@rset_states, 'MART_DISPLAYABLE');
  }

  return (\@dset_states, \@rset_states, \@fset_states);
}



=head2 define_and_validate_sets

  Arg [1]    : hash - set constructor parameters:
                            -dbadaptor    Bio::EnsEMBL::Funcgen::DBAdaptor
                            -name         Data/FeatureSet/ResultSet name to create
                            -feature_type Bio::EnsEMBL::Funcgen::FeatureType
                            -cell_type    Bio::EnsEMBL::Funcgen::CellType
                            -analysis     FeatureSet Bio::EnsEMBL::Analysis
                            -feature_class e.g. annotated or regulatory
                            -description  FeatureSet description
                            -recovery     Allows definition of extant sets so long as they match
                            -append       Boolean - Forces import on top of previously imported data
                            -rollback     Rolls back product feature set.
                            -supporting_sets Complete set of pre-stored supporting or input sets for this DataSet
                            -slices       ARRAYREF of Slices to rollback
  Example    : my $dset = $self->define_and_validate_Set(%params);
  Description: Checks whether set is already in DB based on set name, rolls back features
               if roll back flag set. Or creates new DataSet and Feature|ResultSet if not present.
  Returntype : Bio::EnsEMBL::Funcgen::DataSet
  Exceptions : Throws if DBAdaptor param not valid
  Caller     : Importers and Parsers
  Status     : At risk

=cut

sub define_and_validate_sets{
  my $self = shift;

  #change slice to slices to support multi slice import from InputSet::define_sets
  #Can't do full rollback in slice mode
  #This may not be safe in slice mode as we will then have mixed inputs/outputs

  my ($name, $anal, $ftype, $ctype, $type, $append, $db, $ssets, $description, $rollback, $recovery, $slices, $display_label) = rearrange(['NAME', 'ANALYSIS', 'FEATURE_TYPE', 'CELL_TYPE', 'FEATURE_CLASS', 'APPEND',
				 'DBADAPTOR', 'SUPPORTING_SETS', 'DESCRIPTION', 'ROLLBACK', 'RECOVERY', 'SLICES', 'DISPLAY_LABEL'], @_);


  #VALIDATE CONFIG HASH
  #$config_hash ||= {};#default so exists will work without testing
  #if(keys %{$config_hash}){
  #	#There is a module to handle config hashes somewhere!
  #	throw('config_hash not yet implemented for define_and_validate_sets');
  #my @known_config = ('full_delete');#We never want full delete here as this is a create method!
  #Can we set vars from has by refs like getopts?
  #map {
  #  throw("Found unsupported config hash parameter:\t$_") if ! grep(/^${_}$/, @known_config);
  #} keys %{$config_hash};
  #  }
  
  #define rollback level
  #extract this to _set_rollback_level($rollback_mode, $feature_class)
  my $rollback_level = 0;
  
  #These should be globally defined so all rollback methods can use them
  my %valid_rollback_modes = 
	(
	 product_features => 1,
	 #Just product features and FeatureSet status, what about DataSet status?
	 #full delete does nothing here?
	 
	 sets => 2,
	 #Includes product_features and
	 #deletes supporting_sets entries unless we specify append
	 #revoke all states on Feature/Data/InputSets
	 #Full delete removes Feature/Data/InputSet entries
	 #Never includes ResultSets!
	 
	 supporting_features => 3,
	 #Includes product_feature and sets
	 #Removes all states and supporting features 
	 #inc. ResultSet results/ResultFeatures
	 #Full_delete remove supporting set entries
	 #Otherwise just rollback states for affected sets
	);
  
  if($rollback){
	if(! exists $valid_rollback_modes{$rollback}){
	  #Default to some sensible values
	  $rollback = 'product_features';#default for FeatureSets

	  #Always want overwrite supporting sets if there is a difference
	  $rollback = 'sets' if ($type eq 'regulatory');
	  $rollback = 'supporting_sets' if ($type eq 'result');

	  warn ("You have not set a valid rollback mode(product_features|sets|supporting_features), defaulting to $rollback for feature class $type\n");
	}

	$rollback_level = $valid_rollback_modes{$rollback};
  }


  if($slices && (ref($slices) ne 'ARRAY')){
	throw('-slices param must be an ARRAYREF of Bio::EnsEMBL::Slice objects');
	#Rest of validation done in other methods
  }
  

 
  #But how are we going to resolve the append behaviour when we also want to validate the ssets?
  #Can't, so append also functions to enable addition in the absence of some or all previous data/esets?
  #No this is not true, we want to be able to fetch an extant set for import,
  #we just need to be aware of sset IMPORTED status?
  #This should be a recovery thing, allow fetch, but validate sets?
  

  #Check mandatory params
  if(! (ref($db) && $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'))){
	throw('Must provide a valid Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  }

  throw('Must provide a -name ') if(! defined $name);
  
  #Not necessarily, just do rollback then append?
  #But then we'd potentially have a supporting set associated which has had it's data removed from the feature set.
  #Generating sets for an ExpSet will always have append set
  #This could be valid for generically grabing/creating sets for adding new supporting sets e.g. reg build
  throw('-append and -rollback are mutually exclusive') if $rollback_level && $append;
  
  #This will never happen due to previous test? append will always fail?
  #warn('You are defining a pre-existing FeatureSet without rolling back'.
  #	   ' previous data, this could result in data duplication') if $append && ! $rollback_level;
  #Is this really possible, surely the supporting set will fail to store due to unique key?


  #Should we warn here about append && recovery?
  #Aren't these mutually exclusive?
  #Do we know if we have new data? append should override recovery, or just specifiy append
  #This will stop the import and highlight the issue to the user
  #We need to be able to run with both otherwise the import will not work


  throw('Must provide a -feature_class e.g. annotated, external, result or regulatory') if(! defined $type);
  #Check for annotated, external, regulatory etc here?
  #Should never be external as we don't have DataSets for external sets?
  
  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType',  $ftype);
  if (defined $ctype){
	$db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType', $ctype);
  }
  elsif($type ne 'regulatory'){
	throw('Only Data/FeatureSets with type \'regulatory\' can have an undefined CellType');
	#Coudl extend this to core set by name eq 'RegulatoryFeatures'?
  }

  $db->is_stored_and_valid('Bio::EnsEMBL::Analysis',  $anal);

  my $dset_adaptor = $db->get_DataSetAdaptor;
  my $fset_adaptor = $db->get_FeatureSetAdaptor;
  my $rset_adaptor = $db->get_ResultSetAdaptor;
  
  #DataSet centric definition to enable multiple DataSets
  #to be generated from the same supporting sets
  my $dset = $dset_adaptor->fetch_by_name($name);
  my ($fset, $rset, @input_sets);

  #Validate stored vs passed set data 

  if(defined $dset){
	$self->log('Found Stored DataSet '.$dset->name);

	if($type ne 'result'){#i.e. annotated
	  
	  #Does this account for regulatory?

	  $fset = $dset->product_FeatureSet;
	  #Here we have the possiblity that a feature_set with a different name may have 
	  #been associated with the DataSet

	  if(defined $fset){
		$self->log("Found associated product FeatureSet:\t".$fset->name);
		
		#if(! $clobber && 
		if($fset->name ne $name){
		  throw('Invalid product FeatureSet name ('.$fset->name.') for DataSet ('.$name.'). Rollback will overwrite the FeatureSet and mismatched name will be retained.');
		  #Need to clobber both or give explicit name for datasets or rename dataset???
		  #Force this throw for now, make this fix manual as we may end up automatically overwriting data
		}  
	  }

	  #This needs to be modified to support InputSets in ResultSets?
	  #Would never have mixed Input/ResultSets so no need
	  #Could potential need to do it for mixed Result/FeatureSets
	  #if we ever use an analysis which uses both set types
	  
	  #check supporting_sets here if defined
	  #We have the problem here of wanting to add ssets to a previously existing dset
	  #we may not know the original sset, or which of the ssets are new
	  #Hence there is a likelihood of a mismatch.

	  #Much of this is replicated in store_udpated sets
	  

	  if(defined $ssets){
		my @sorted_ssets = sort {$a->dbID <=> $b->dbID} @{$ssets};
		my @stored_ssets = sort {$a->dbID <=> $b->dbID} @{$dset->get_supporting_sets};
		my $mismatch = 0;
		
		$mismatch = 1 if(scalar(@sorted_ssets) != scalar(@stored_ssets));
		
		if(! $mismatch){
		  
		  for my $i(0..$#stored_ssets){
			
			if($stored_ssets[$i]->dbID != $sorted_ssets[$i]->dbID){
			  $mismatch=1;
			  last;
			}
		  }
		}
		
		
		

		if($mismatch){
		  #We're really print this names here which may hide the true cell/feature/anal type differences.
		  my $mismatch = 'There is a (name/type/analysis) mismatch between the supplied supporting_sets and the'.
			' supporting_sets in the DB for DataSet '.$dset->name."\n\nStored:\n"
			  .join(', ', (map { $_->name } @stored_ssets))."\n\nSupplied supporting_sets:\n"
				.join(', ', (map { $_->name } @sorted_ssets));
		  

		  if($append){
			warn($mismatch."\n\nAppending supporting set data to unvalidated supporting sets");
		  }
		  elsif($rollback_level > 1){#supporting set rollback
			warn($mismatch."\n\nReplacing previously stored supporting sets with newly defined sets\n");
			
			if($slices){
			  warn("WARNING:\tPerforming supporting_set rollback in slice mode. This may corrupt the supporting_set definition for other slices in this DataSet if they are not re-generated using the same supporting_sets\n");
			}

			#Remove supporting_set entries
			#This should be in a rollback_DataSet method
			#This has moved to DataSetAdaptor::store_update_sets
				
			#Reset supporting sets
			$dset->{'supporting_sets'} = undef;
			$dset->add_supporting_sets(\@sorted_ssets);
			#Move this to last block?
			#This will currently fail as it test for product_FeatureSet
			#How do we get around this? Remove IMPORTED status and only throw if fset has IMPORTED status?

			#warn "pre store sset ".@{$dset->get_supporting_sets};

			#($dset) = @{$dset_adaptor->store_updated_sets([$dset], $rollback_level)};
			#$dset->adaptor->store_regbuild_meta_strings($dset, $rollback_level) if $type eq 'regulatory';
		  }
		  else{
			throw($mismatch);
		  }
		}
	  }
	  else{
		warn("No supporting sets defined, skipping supporting set validation for definition of DataSet:\t".$name);
	  }
	}
	else{#result_features from InputSet
	  #Do we ever pass supporting sets here?
	  #Do we need to test vs stored_sets?


	  #There is the potential for more than one ResultSet to be associated with DataSet
	  #But as we are using the same name, this restricts the number wrt the cardinality
	  #of the name field. i.e. 1 name per analysis/cell_type/feature_type.
	  #This now works slightly differently to the rest of this method as we
	  #need to treat the ResultSet as we are currently treating the FeatureSet below.
	  
	  #However, the use case of this method is for one InputSet giving rise to one ResultSet
	  #Hence just throw if we find more than one or have a name mismatch???
	  my @stored_sets = @{$dset->get_supporting_sets};



	  #THis assumes we will always have supporting sets
	  #and is failing as we have removed this test in DataSet::new
	  #But where are we storing it without the supporting set?

	  if(scalar(@stored_sets) > 1){
		throw('define_and_validate_sets does not yet support DataSets with multiple supporting ResultSets for result_features');
	  }
	  elsif(! @stored_sets){
		throw("DataSet($name) does not have any stored supporting sets. These should have been defined when storing the DataSet");
		#Or should we handle this?
	  }

	  $rset = $stored_sets[0];

	  if($rset->set_type ne 'result'){
		throw("DataSet already contains a supporting set which is not a ResultSet:\t".$rset->set_type."\t".$stored_sets[0]->name);
	  }
	  elsif($ssets){
		#Do we ever pass supporting sets, test for completeness

		#Just test we have the same supplied ssets if it is defined
		if(scalar(@$ssets) != 1){
		  throw("ResultFeature data sets currently only support one supporting ResultSet.\nSupproting sets:\t".
				join(', ', (map { $_->name.'('.$_->set_type } @$ssets)));
		}
		elsif(! ($rset->dbID == $ssets->[0]->dbID) && 
			  ($ssets->[0]->set_type eq 'result')){
		  throw('Supplied supporting set('.$ssets->[0]->name.') does not match stored supporting set('.$rset->name.')');
		}
	  }
	  
	  @input_sets = @{$rset->get_InputSets};
	}
  }



  if($type eq 'result'){
	
	#Validate the defined InputSets
	if (scalar(@$ssets) > 1) {
	  throw("define_and_validate_sets does not yet support multiple InputSets for defining a ResultSet:\t".$name);
		
	}

	if ($ssets->[0]->set_type ne 'input') {
	  throw("To define a ResultSet($name) containing result_features, you must provide and InputSet as a supporting set\nArray based ResultSets(i.e. experimental_chip/channel) are not defined using this method, see specific Import Parsers.");
	}		
	  
	  
	#Try and grab the rset just in case it has been orphaned somehow
	if (! defined $rset) {
	  $rset = $rset_adaptor->fetch_all_by_name($name, $ftype, $ctype, $anal)->[0];
	  #Should only ever be one given all parts of unique key
	  @input_sets = @{$rset->get_InputSets} if $rset;
		
	}


	if (defined $rset) {		#Validate stored InputSets

	  if (scalar(@input_sets) != scalar(@$ssets)) {
		throw('Found mismatch between number of previously stored InputSets('.scalar(@input_sets).') and defined InputSets('.scalar(@$ssets).'). You must provide a complete list of InputSets to define your ResultSet.');
	  }
		
	  if ($input_sets[0]->dbID != $ssets->[0]->dbID) {
		throw('Found dbID mismatch between previously stored InputSet('.$input_sets[0]->name.') and define InputSet('.$ssets->[0]->name.')');
	  }

	  #rollback ResultSet/InputSet here?
	  if($rollback_level > 2){
		warn "rollback not yet fully implemented for Result/InputSets";
		
		#Does this need to be by slice?
		#What about states if we are running in parallel?
		
		if($slices){
		  map {$self->rollback_ResultSet($rset, $rollback, $_)} @$slices;
		}
		else{
		  $self->rollback_ResultSet($rset, $rollback);
		}

	  }

	} 
	else{#define ResultSet
	  ($rset) = @{$rset_adaptor->store(Bio::EnsEMBL::Funcgen::ResultSet->new
									   (
										-name => $name,
										-feature_type => $ftype,
										-cell_type => $ctype,
										-table_name => 'input_set',
										-table_id   => $ssets->[0]->dbID,
										-analysis   => $anal
									   )
									  )};
		
	}
  } 
  else{#annotated/regulatory/external i.e. FeatureSet

	#Try and grab the fset just in case it has been orphaned somehow
	if(! defined $fset){
	  $fset = $fset_adaptor->fetch_by_name($name);
	  
	  if(defined $fset){
		#Now we need to test whether it is attached to a dset
		#Will be incorrect dset if it is as we couldn't get it before
		#else we test the types and rollback
		$self->log("Found stored orphan FeatureSet:\t".$fset->name);
		
		my $stored_dset = $dset_adaptor->fetch_by_product_FeatureSet($fset);
		
		if(defined $stored_dset){
		  throw('Found FeatureSet('.$name.') associated with incorrect DataSet('.$stored_dset->name.
				").\nTry using another -name in the set parameters hash");
		  
		}
	  }
	}
	
	#Rollback or create FeatureSet
	if(defined $fset){
	  
	  if($rollback_level){
		#Don't check for IMPORTED here as we want to rollback anyway
		#Not forcing delete here as this may be used as a supporting set itself.

		$self->rollback_FeatureSet($fset, undef, $slices);
	  } 
	  elsif ($append || $recovery) {
		#This is only true if we have an sset mismatch
		
		#Do we need to revoke IMPORTED here too?
		#This behaves differently dependant on the supporting set.
		#InputSet status refers to loading in FeatureSet, where as ResultSet status refers to loading into result table

		#So we really want to revoke it
		#But this leaves us vulnerable to losing data if the import crashes after this point
		#because we have no way of assesing which is complete data and which is incomplete data
		#within a feature set.
		#This means we need a status on supporting_set, not InputSet or ResultSet
		#as this has to be in the context of a dataset.
		#Grrr, this means we need a SupportingSet class which simply wraps the InputSet/ResultSet
		#We also need a single dbID for the supporting_set table
		#Which means we will have to do some wierdity with the normal dbID implementation
		#i.e. Have supporting_set_id, so we can still access all the normal dbID method for the given Set class
		#This will have to be hardcoded into the state methods
		#Also will need to specify when we want to store as supporting_status or normal set status.

		#This is an awful lot to protect against vulnerability
		#Also as there easy way to track what features came from which supporting set
		#There isn't currently a viable way to rollback, hence will have to redo the whole set.

		#Maybe we can enforce this by procedure?
		#By simply not associating the supporting set until it has been loaded into the feature set?
		#This may cause even more tracking problems

		#Right then, simply warn and do not revoke feature_set IMPORTED to protect old data?
		#Parsers should identify supporting_sets(InputSets) which exist but do not have IMPORTED
		#status and fail, specifying -recover which will rollback_FeatureSet which will revoke the IMPORTED status

		#This can mean a failed import can leave a partially imported feature set with the IMPORTED status!!!

		#We just need to handle InputSets and ResultSets differently.
		#In parsers or here?
		#Probably best in the parsers as this is where the states are set.
	  

		#Should we throw here for ResultSet?
		#Force rollback of FeatureSet first or create new one?
		#And throw for InputSet?
		#This again comes back to whether we will ever have more than one file 
		#for a give InputSet, currently not.

		$self->log("WARNING\t::\tAdding data to a extant FeatureSet:\t".$fset->name);
	  } else {
		throw('Found extant FeatureSet '.$fset->name.'. Maybe you want to specify the rollback, append or recovery parameter or roll back the FeatureSet separately?');
	  }
	} else {
	  #create a new one
	  $self->log("Creating new FeatureSet:\t".$name);

	  $fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
													 -name          => $name,
													 -feature_type  => $ftype,
													 -cell_type     => $ctype,
													 -analysis      => $anal,
													 -feature_class => $type,
													 -description   => $description,
													 -display_label => $display_label,
													);
	  ($fset) = @{$fset_adaptor->store($fset)};
	}
  }

  #Create/Update the DataSet
  if(defined $dset){	
	#Could do these updates above?
	#But delayed to reduce redundancy

	if($type ne 'result'){

	  if(! defined $dset->product_FeatureSet){
		$self->log("Updating DataSet with new product FeatureSet:\t".$fset->name);
		$dset->product_FeatureSet($fset);
	  }
		
	  $dset = $dset_adaptor->store_updated_sets([$dset], $rollback_level)->[0];
	  #This cannot store the focus sets as we don't know which are which yet
	  #Only the script knows this
	  # $dset->adaptor->store_regbuild_meta_strings($dset, $rollback_level) if $type eq 'regulatory';
	}
	else{
	  #We may have the case where we have a DataSet(with a FeatureSet) but no ResultSet
	  #i.e. Load result_features after peak calls
	  #So update dset with ResultSet

	  if(! @{$dset->get_supporting_sets}){
		$self->log("Updating DataSet with new ResultSet:\t".$rset->name);
		$dset->add_supporting_sets([$rset]);
		$dset = $dset_adaptor->store_updated_sets([$dset], $rollback_level)->[0];
	  }
	}
  }
  else{
	$self->log("Creating new ${type}_feature DataSet:\t".$name);

	if($type ne 'result'){
	  ($dset) = @{$dset_adaptor->store(Bio::EnsEMBL::Funcgen::DataSet->new
									   (
										-name => $name,
										-feature_set => $fset,
										-supporting_sets => $ssets,
									   ))};
	  #$dset->adaptor->store_regbuild_meta_strings($dset, $rollback_level) if $type eq 'regulatory';
	}
	else{
	  warn "creating dataset $name with supporting set $rset";
	  ($dset) = @{$dset_adaptor->store(Bio::EnsEMBL::Funcgen::DataSet->new
									   (
										-name => $name,
										-supporting_sets => [$rset],
									   ))};
	}
  }

  return $dset;
}


#Rollback/load methods migrated from DBAdaptor
#Move to SetAdaptors, better located and will remove cyclical dependancy

=head2 rollback_FeatureSet

  Arg [0]    : Bio::EnsEMBL::Funcgen::FeatureSet
  Arg [1]    : optional - boolean force delete flag, if this FeatureSet is use as a support 
               for another DataSet.
  Arg [2]    : optional - arrayref of Bio::EnsEMBL::Slice objects to rollback
  Arg [3]    : optional - boolean flag to perform full rollback i.e. default will just remove feature
               specifying this with also delete the feature_set record
  Example    : $self->rollback_FeatureSet($fset);
  Description: Deletes all status and feature entries for this FeatureSet.
               Checks whether FeatureSet is a supporting set in any other DataSet.
  Returntype : none
  Exceptions : Throws if any deletes fails or if db method unavailable
  Caller     : Importers and Parsers
  Status     : At risk

=cut


sub rollback_FeatureSet{
  my ($self, $fset, $force_delete, $slices, $full_delete) = @_;

  #Remove force delete and just throw?
  #Currently only used in project_feature_set.
  #May want to keep an old RegBuild for mapping/comparison?
  #Coudl get around this by simply deleting the data_set? Unknown impact.
  #Move to config hash?
  #No need for rollback_level here as we always want to do the same thing
  
  my ($sql, $slice_name);
  my $slice_join = '';
  my $table = $fset->feature_class.'_feature';
  my $adaptor = $fset->adaptor || throw('FeatureSet must have an adaptor');
  my $db = $adaptor->db;
  #Cyclical dpendancy here, so not strictly necessary.
  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);


  $self->log_header('Rolling back '.$fset->feature_class." FeatureSet:\t".$fset->name);

  if($slices){

	if($full_delete){
	  throw("Cannot specify a full_delete for a Slice based rollback:\t".$fset->name);
	}

	
	if(! ref($slices) eq 'ARRAY'){
	  throw('Slices must be an ARRAYREF of Slice objects');
	}

   	map { throw("Must pass a valid Bio::EnsEMBL::Slice") if (! (ref($_) && $_->isa('Bio::EnsEMBL::Slice'))) } @$slices;
	$self->log("Restricting to slices:\n\t\t".join("\n\t\t", (map { $_->name } @$slices)) );
	#Allow subslice rollback only for one slice at a time
	my $subslice = (scalar(@$slices) == 1) ? 1 : 0;
	my @sr_ids;

	foreach my $slice(@$slices){
	  my $efg_sr_id = $fset->get_FeatureAdaptor->get_seq_region_id_by_Slice($slice);
	  
	  if(! $efg_sr_id){
		$self->log("Slice is not present in eFG DB:\t".$slice->name);
	  }else{

		if(! $subslice){#Test is not subslice
		  my $full_slice = $slice->adaptor->fetch_by_region(undef, $slice->seq_region_name);

		  if(($slice->start != 1) ||
			 ($full_slice->end != $slice->end)){
			throw("Can only rollback subslices one at a time:\nRollback slice:\t"
				  .$slice->name."\nFull slice:\t".$full_slice->name);
		  }
		}

		push @sr_ids, $efg_sr_id;
	  }
	}
	
	if(scalar(@sr_ids) == 1){
	  #Allow sub slice rollback
	  #add range here from meta coord?
	  $slice_join = " and f.seq_region_id=$sr_ids[0] and f.seq_region_start<=".$slices->[0]->end.' and f.seq_region_end>='.$slices->[0]->start;
	}
	else{
	  $slice_join = ' and f.seq_region_id in ('.join(', ', @sr_ids).')';
	}
  }
  


  #Check whether this is a supporting set for another data_set
  my @dsets = @{$db->get_DataSetAdaptor->fetch_all_by_supporting_set($fset)};

  if(@dsets){
	my $txt = $fset->name." is a supporting set of the following DataSets:\t".join(', ', (map {$_->name} @dsets));

	if($force_delete){
	  $self->log("WARNING:\t$txt\n");
	}
	else{
	  throw($txt."\nPlease resolve or specify the force_delete argument")
	}
  }

  #Remove states
  if(! $slices){
	$fset->adaptor->revoke_states($fset);
	
	#Revoke InputSet states here as this refers to whether
	#they are imported in the FeatureSet
	#Do this in FeatureSet->revoke_states?

	my $dset = $db->get_DataSetAdaptor->fetch_by_product_FeatureSet($fset);

	#Account for absent dset if we have an external_feature set

	if((! defined $dset) &&
	   $fset->feature_class ne 'external'){
	  warn "WARNING:\tFeatureSet ".$fset->name." does not have an associated DataSet. Rollback may be incomplete";
	}

	if($dset){
	  	
	  foreach my $sset(@{$dset->get_supporting_sets}){	  
		#Maybe skip this if we defined slice?
		
		#??? Do we want to do this?
		#This is dependant on the feature_class of the InputSet
		#result InputSets may have been imported as ResultFeatureCollections
		#So we want to leave those in place
		#annotated feature_class InputSets are directly imports, so the status of these refers
		#to the FeatureSet import status
		#Where is the imported status set for SWEmbl?
		
		if(($sset->feature_class eq 'annotated') &&
		   $sset->isa('Bio::EnsEMBL::Funcgen::InputSet')){
		  
		  $self->rollback_InputSet($sset) if $sset->isa('Bio::EnsEMBL::Funcgen::InputSet');
		  $self->rollback_InputSet($sset);#add full delete here?
		  #Do not want to rollback here for other type of sset
		}
	  }
	}
  }
  else{
	$self->log('Skipping '.$fset->name.' revoke_states for partial Slice rollback, maybe revoke IMPORTED? ');
  }

  #should add some log statements here?

  my $row_cnt;

  #Rollback reg attributes
  if($fset->feature_class eq 'regulatory'){
	$sql = "DELETE ra from regulatory_attribute ra, $table f where f.${table}_id=ra.${table}_id and f.feature_set_id=".$fset->dbID.$slice_join;
	$self->rollback_table($sql, 'regulatory_attribute', undef, $db);



	if($full_delete){
	  #Now delete meta entries
	  #This is messy as we use the following meta_key nomencalture
	  #which do not match the fset names
	  #regbuild.feature_set_ids_v5
	  #regbuild.feature_type_ids_v5
	  #regbuild.focus_feature_set_ids 
	  #regbuild.initial_release_date_v6
	  #regbuild.last_annotation_update_v6 
	  #regbuild.version NEED TO ADD THIS
	  #Also need to revise how these are generated by build_reg_feats.
	  #WHat about new cell_type level feature sets?
	  #How will we model these in the meta table?
	  
	  warn "Need to revise meta table entries before we add a delete here, remove manually for now for:\t".$fset->name;
	  
	  #We would only remove meta entries if we are performing a full rollback
	  my $version;
	  ($version = $fset->name) =~ s/.*_v([0-9]+)$/$1/;
	  $version = ($version eq  $fset->name) ? '' : "_v${version}";

	  #These are versionless meta_keys and apply to all sets
	  #handle these in reg build script
	  #'regbuild.initial_release_date',
	  #'regbuild.last_annotation_update'
	  #'regbuild.version'

	  foreach my  $mkey('regbuild.%s.feature_set_ids',
						'regbuild.%s.feature_type_ids', 
						'regbuild.%s.focus_feature_set_ids'){
		
		my $meta_key = sprintf($mkey, $fset->cell_type->name).$version;
		$sql = "DELETE from meta where meta_key='${meta_key}'";
		$self->rollback_table($sql, 'meta', undef, $db);
	  }
	}
  }


  #Need to remove object xrefs here
  #Do not remove xrefs as these may be used by something else!
  $sql = "DELETE ox from object_xref ox, $table f where ox.ensembl_object_type='".ucfirst($fset->feature_class)."Feature' and ox.ensembl_id=f.${table}_id and f.feature_set_id=".$fset->dbID.$slice_join;
  $self->rollback_table($sql, 'object_xref', 'object_xref_id', $db);
    

  #Remove associated_feature_type records
  #Do not remove actual feature_type records as they may be used by something else.

  $sql ="DELETE aft from associated_feature_type aft, $table f where f.feature_set_id=".$fset->dbID." and f.${table}_id=aft.table_id and aft.table_name='".$fset->feature_class."_feature'".$slice_join;
  $self->rollback_table($sql, 'associated_feature_type', undef, $db);



  #Remove features
  $sql = "DELETE f from $table f where f.feature_set_id=".$fset->dbID.$slice_join;
  $self->rollback_table($sql, $table, "${table}_id", $db);

  if($full_delete){ #Also delete feature/data_set records
	
	$sql = "DELETE from feature_set where feature_set_id=".$fset->dbID;
	$self->rollback_table($sql, 'feature_set', 'feature_set_id', $db);
	$self->log("Deleted feature_set entry for:\t".$fset->name);


	$sql = "DELETE from data_set where feature_set_id=".$fset->dbID;
	$self->rollback_table($sql, 'data_set', 'data_set_id', $db);
	$self->log("Deleted associated data_set entry for:\t".$fset->name);
  }

  return;
}


=head2 rollback_ResultSet

  Arg[1]     : Bio::EnsEMBL::Funcgen::ResultSet
  Arg[2]     : Boolean - optional flag to roll back array results
  Example    : $self->rollback_ResultSet($rset);
  Description: Deletes all status. chip_channel and result_set entries for this ResultSet.
               Will also rollback_results sets if rollback_results specified.  This will also
               update or delete associated ResultSets where appropriate.
  Returntype : Arrayref containing the ResultSet and associated DataSet which have not been rolled back
  Exceptions : Throws if ResultSet not valid
               Throws is result_rollback flag specified but associated product FeatureSet found.
  Caller     : General
  Status     : At risk

=cut

#Need to change slice to slices ref here
#Need to add full rollback, which will specify to remove all sets
#as well as results and 
#These params need clarifying as their nature changes between input_set and array rsets
#Don't we always want to rollback_results?
#force should only really be used to rollback InputSet ResultFeature sets
#i.e. Read collections which are not used as direct input for the linked product FeatureSet
#This should fail with array data associated with a product feature set

#Do we want to separate ResultFeature rollback from result rollback?
#Currently the array based collection rollback is done by hand
#Could be done via the ResultFeature Collector, but should probably use this method.


#rollback_results is only used in the MAGE parser to identify sets which have an 
#associated product fset.
#Can't really separate due to integrated functionality

sub rollback_ResultSet{
  my ($self, $rset, $rollback_results, $slice, $force, $full_delete) = @_;
  
  if(! (ref($rset) && $rset->can('adaptor') && defined $rset->adaptor)){
	throw('Must provide a valid stored Bio::EnsEMBL::ResultSet');
  }

  if($slice && $rset->table_name ne 'input_set'){
	throw('Can only rollback_ResultSet by Slice if the ResultSet contains InputSets');	
  }

  #We're still validating against itself??
  #And reciprocating part of the test :|
  my $sql;
  my $db = $rset->adaptor->db;#This needs to be tested
  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);
  $self->log("Rolling back ResultSet:\t".$rset->name);
  my $dset_adaptor = $self->db->get_DataSetAdaptor;
  my $rset_adaptor = $self->db->get_ResultSetAdaptor;
  my @skipped_sets;
  
  ### Check if this ResultSet is part of a DataSet with a product feature set

  foreach my $dset(@{$dset_adaptor->fetch_all_by_supporting_set($rset)}){
	
	if (defined $dset){
	  $self->log('Found linked DataSet('.$dset->name.") for ResultSet:\t".$rset->log_label);
	  
	  if(my $fset = $dset->product_FeatureSet){
		@skipped_sets = ($rset,$dset);

		#What impact does this have on result_rollback?
		#None as we never get there
		#But what if we have specified rollback results?
		#We should throw here as we can't perform the rollback
		
		if($rollback_results){

		  if($rset->table_name ne 'input_set' ||
			(! $force)){#is an input_set/reads collection
			#This will always throws for non-input_set ResultSets

			throw("Could not rollback supporting ResultSet and results for:\t".$rset->log_label.
				  "\nEither manually resolve the supporting/feature set relationship or set the 'force' flag.\n");
			#  ."Alternatively omit the rollback_results argument if you simply want to redefine the ResultSet without loading any new data");
			#This last bit is no longer true
			#Remove rollback_results?
		  }
		  else{
			@skipped_sets = ();
			$self->log("Forcing results rollback for InputSet based ResultSet:\t".$rset->log_label);
		  }
		}

		if(@skipped_sets){
		  $self->log('Skipping rollback. Found product FeatureSet('.$fset->name.") for supporting ResultSet:\t".$rset->log_label);
		}

	  }
	  elsif((! defined $slice) &&
			$full_delete){
		#Found rset in dset, but not yet processed so can remove safely.
		$self->unlink_ResultSet_DataSet($rset, $dset);
	  }
	}
  }
   

  #Now do similar for all associated ResultSets
  if(! @skipped_sets){

	
	#Rollback results if required
	if($rollback_results){

	  $self->log("Rolling back results for ResultSet:\t".$rset->log_label);
	  #Check result_set_input_ids are present in other result sets.
	  my @assoc_rsets = @{$rset_adaptor->fetch_all_linked_by_ResultSet($rset)};
	  my $feature_supporting = 0;
	  
	  foreach my $assoc_rset(@assoc_rsets){
		
		foreach my $dset(@{$dset_adaptor->fetch_all_by_supporting_set($assoc_rset)}){
		 
		  #Check for other product_FeatureSets
		  if(my $fset = $dset->product_FeatureSet){
			$feature_supporting++;
			$self->log('Found product FeatureSet('.$fset->name.
					   ") for associated supporting ResultSet:\t".$rset->log_label);

			if($rset->table_name ne 'input_set' ||
			   (! $force)){#is an input_set/reads collection
			  $feature_supporting++;
			}
		  }
		}					
	  }
		

	  if(! $feature_supporting){

		#RollBack result_feature table first
		$self->rollback_ResultFeatures($rset, $slice);

		#Now rollback other states
		$rset->adaptor->revoke_states($rset);


		#This also handles Echip status rollback
		if ($rset->table_name ne 'input_set'){
		  $self->log("Rolling back result table for ResultSet:\t".$rset->log_label);
		  $self->rollback_results($rset->result_set_input_ids);
		}

		$self->log('Removing result_set_input entries from associated ResultSets') if @assoc_rsets;
		
		if((! $slice) &&
		   $full_delete){

		  #Now remove result_set_input_ids from associated rsets.
		  foreach my $assoc_rset(@assoc_rsets){
			$sql = 'DELETE from result_set_input where result_set_id='.$assoc_rset->dbID.
			  ' and result_set_input_id in('.join', ', @{$assoc_rset->result_set_input_ids}.')';
			$db->dbc->do($sql);
			
			# we need to delete complete subsets from the result_set table.
			my $subset = 1;
			
			foreach my $cc_id(@{$assoc_rset->result_set_input_ids}){
			  
			  if(! grep { /$cc_id/ } @{$rset->result_set_input_ids}){
				$subset = 0;
				last;
			  }
			}
		  
			#$assoc_rset is complete subset of $rset so can delete
			#We know this does not have an assoicated product feature set
			#Only if it is not derived from an input_set
			if($subset){
			  $self->log("Deleting associated subset ResultSet:\t".$assoc_rset->log_label);
			  
			  #Delete status entries first
			  $assoc_rset->adaptor->revoke_states($assoc_rset);
			  
			  #All cc records will have already been deleted
			  $sql = 'DELETE from result_set where result_set_id='.$assoc_rset->dbID;
			  $db->dbc->do($sql);
			}
		  }
		}


		#Now warn about Echips in Experiments which may need removing.
		if($rset->table_name ne 'input_set'){
		  my %experiment_chips;
		  
		  foreach my $echip(@{$rset->get_ExperimentalChips}){
			$experiment_chips{$echip->experiment->name}{$echip->unique_id} = undef;
		  }
		  
		  foreach my $exp(keys %experiment_chips){
			$self->log("Experiment $exp has had ".scalar(values %{$experiment_chips{$exp}}).
					   " ExperimentalChips rolled back:\t".join('; ', values %{$experiment_chips{$exp}}).
					   ".\nTo fully remove these, use the rollback_experiment.pl (with -chip_ids) script");
		  }
		}
		else{
		  #Should only be one to rollback
		  foreach my $iset(@{$rset->get_InputSets}){
			$self->rollback_InputSet($iset);
		  }
		}
	  }
	  else{
		#$self->log("Skipping result rollback, found $feature_supporting associated supporting ResultSets for:\t".$rset->log_label);
		#warn("Skipping result rollback, found $feature_supporting associated supporting ResultSets for:\t".$rset->log_label);
		#do we need to return this info in skipped_rsets?
		#This is just to allow importer to know which ones 
		#weren't rolled back to avoid naming clashes.
		#so no.

		#But the results persist on the same chip_channel_ids
		#So not returning this rset may result in loading of more data
		#This should fail as status entries will not have been removed
		#Still we should throw here as we'll most likely want to manually resolve this
		#Besides this would be obfuscating the function

		throw("Could not rollback ResultSet and results, found $feature_supporting associated supporting ".
			  "ResultSets for:\t".$rset->log_label."\nManually resolve the supporting/feature set relationship or omit the ".
			 "rollback_results argument if you simply want to redefine the ResultSet without loading any new data");
	  }
	}
	else{
	  $self->log('Skipping results rollback');
	  
	  if($rset->name =~ /_IMPORT$/){
		throw("Rolling back an IMPORT set without rolling back the result can result in ophaning result records for a whole experiment.  Specify the result_rollback flag if you want to rollback the results for:\t".$rset->log_label);
	  }
	}
	
	#Delete chip_channel and result_set records
	#This should only be done with full delete
	if((! $slice) &&
	  $full_delete){
	  $sql = 'DELETE from result_set_input where result_set_id='.$rset->dbID;
	  $self->rollback_table($sql, 'result_set_input', 'result_set_input_id', $db);

	  $sql = 'DELETE from result_set where result_set_id='.$rset->dbID;
	  $db->dbc->do($sql);
	  $self->rollback_table($sql, 'result_set', 'result_set_id', $db);
	}
  }

  return \@skipped_sets;
}



sub unlink_ResultSet_DataSet{
  my ($self, $rset, $dset, $new_name) = @_;

  #validate set vars

  my $db = $rset->adaptor->db;

  $self->log("Removing supporting ResultSet from DataSet:\t".$dset->name."\tResultSet:".$rset->log_label);
  my $sql = 'DELETE from supporting_set where data_set_id='.$dset->dbID.
	' and type="result" and supporting_set_id='.$rset->dbID;
  
  warn "Removing ".$rset->log_label." as a supporting set to DataSet:\t".$dset->name.
	"\nThis may result in a DataSet with no supporting sets";
  $db->dbc->do($sql);

  if($new_name){
	#We risk overwriting any previously renamed result sets.
	#Should use datestamp?
	$sql = 'UPDATE result_set set name="OLD_'.$rset->name.'" where result_set_id='.$rset->dbID;
	$self->db->dbc->do($sql);

	if($dset->product_FeatureSet){
	  $self->log('Associated DataSet('.$dset->name.') has already been processed. It is not wise to replace a supporting set without first rolling back the FeatureSet, as there may be additional supporting data');
	  warn 'Associated DataSet('.$dset->name.') has already been processed. It is not wise to replace a supporting set without first rolling back the FeatureSet, as there may be additional supporting data';
	}
  }

  return;
}

=head2 rollback_InputSet

  Arg[1]     : Bio::EnsEMBL::Funcgen::InputSet
  Example    : $self->rollback_InputSet($eset);
  Description: Deletes all status entries for this InputSet and it's Subsets
  Returntype : none
  Exceptions : Throws if any deletes fails or if db method unavailable
  Caller     : Importers and Parsers
  Status     : At risk

=cut


sub rollback_InputSet{
  my ($self, $eset, $force_delete, $full_delete) = @_;


  #Need to implement force_delete!!!!!!!!!!!!!!!!!!!!!!
  #Need to check this is not used in a DataSet/ResultSet

  my $adaptor = $eset->adaptor || throw('InputSet must have an adaptor');
  my $db = $adaptor->db;
  

  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::InputSet', $eset);

  $self->log("Rolling back InputSet:\t".$eset->name);

  #SubSets
  foreach my $esset(@{$eset->get_InputSubsets}){
	$esset->adaptor->revoke_states($esset);
  }

  #InputSet
  $eset->adaptor->revoke_states($eset);

  return;
}
  

=head2 rollback_results

  Arg[1]     : Arrayref of chip_channel ids
  Example    : $self->rollback_results($rset->chip_channels_ids);
  Description: Deletes all result records for the given chip_channel ids.
               Also deletes all status records for associated experimental_chips or channels
  Returntype : None
  Exceptions : Throws if no chip_channel ids provided
  Caller     : General
  Status     : At risk

=cut

#changed implementation to take arrayref

sub rollback_results{
  my ($self, $cc_ids) = @_;

  my @cc_ids = @{$cc_ids};
  
  #Need to test for $self->db here?


  if(! scalar(@cc_ids) >0){
	throw('Must pass an array ref of result_set_input_ids to rollback');
  }
  
  #Rollback status entries
  #Cannot use revoke_states here?
  #We can if we retrieve the Chip or Channel first
  #Add to ResultSet adaptor
  my $sql = 'DELETE s from status s, result_set_input rsi WHERE rsi.result_set_input_id IN ('.join(',', @cc_ids).
	') AND rsi.table_id=s.table_id AND rsi.table_name=s.table_name';
  
  if(! $self->db->dbc->do($sql)){
	throw("Status rollback failed for result_set_input_ids:\t@cc_ids\n".$self->db->dbc->db_handle->errstr());
  }


  #Rollback result entries
  $sql = 'DELETE from result where result_set_input_id in ('.join(',', @cc_ids).');';
  $self->rollback_table($sql, 'result', 'result_id', $self->db);
   return;
}


=head2 rollback_ResultFeatures

  Arg[0]     : Bio::EnsEMBL::Funcgen::ResultSet
  Arg[1]     : Optional - Bio::EnsEMBL::Slice
  Arg[2]     : Optional - no_revoke Boolean. This is only used when generating new windows
               from a 0 window size which has been projected from a previous assembly.
  Example    : $self->rollback_result_features($rset);
  Description: Deletes all result_feature records for the given ResultSet.
               Also deletes 'RESULT_FEATURE_SET' status.
  Returntype : None
  Exceptions : Throws if ResultSet not provided
  Caller     : General
  Status     : At risk

=cut

#Need to implement this by slice?
#Then we would need to add a force flag for creation to override the status check

sub rollback_ResultFeatures{
  my ($self, $rset, $slice, $no_revoke) = @_;


  #what about?
  if(! (ref($rset) && $rset->can('adaptor') && defined $rset->adaptor)){
	throw('Must provide a valid stored Bio::EnsEMBL::ResultSet');
  }
  
  if(! $slice && $no_revoke){
	throw("Cannot rollback_ResultFeatures with no_reovke unless you specify a Slice");
  }
  #else warn if slice and no_revoke?

  my ($sql, $slice_name, $slice_constraint);

  if($slice){

	if(ref($slice) && $slice->isa('Bio::EnsEMBL::Slice')){

	  #Need to test for full slice here?
	  my $full_slice = $slice->adaptor->fetch_by_region(undef, $slice->seq_region_name);
	  $slice_name = "\t".$slice->name;
	  $slice_constraint = ' and seq_region_id='.$rset->adaptor->db->get_ResultFeatureAdaptor->get_seq_region_id_by_Slice($slice);
	  #my $slice_srid = $rset->adaptor->db->get_ResultFeatureAdaptor->get_seq_region_id_by_Slice($slice);
	  $slice_constraint = ' and seq_region_id='.$rset->adaptor->db->get_ResultFeatureAdaptor->get_seq_region_id_by_Slice($slice);

	  if(($slice->start != 1) ||
		 ($slice->end != $full_slice->end)){

		#Need to test whether we have non-0 wsize collections without the exact seq_region values
		$sql='SELECT window_size from result_feature where result_feature_id='.$rset->dbID.
		  ' and window_size!=0 and seq_region_start!='.$slice->start.' and seq_region_end!='.$slice->end.$slice_constraint;
		
		throw("rollback_ResultFeatures does not yet support non-full length Slices:\t".$slice_name);
	  }

	}
	else{
	  throw('slice argument must be a valid Bio::EnsEMBL::Slice');
	}
  }

  #We're still validating against itself??
  #And reciprocating part of the test :|
  my $db = $rset->adaptor->db;
  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);

  #Do this conditionally on whether it is a result_feature_set?
  #This may break if we have removed the status but not finished the rollback so no!
  $self->log("Rolling back result_feature table for ResultSet:\t".$rset->name.$slice_name);


  #Rollback status entry
  if($rset->has_status('RESULT_FEATURE_SET') && ! $no_revoke){
	$rset->adaptor->revoke_status('RESULT_FEATURE_SET', $rset);
  }


  #Cannot use revoke_states here?
  #We can if we retrieve the Chip or Channel first
  #Add to ResultSet adaptor
  $sql = 'DELETE from result_feature where result_set_id='.$rset->dbID.$slice_constraint;
  $self->rollback_table($sql, 'result_feature', 'result_feature_id', $db);

  return;
}



=head2 rollback_ArrayChips

  Arg[1]     : ARRAYREF: Bio::EnsEMBL::Funcgen::ArrayChip objects
  Example    : $self->rollback_ArrayChips([$achip1, $achip2]);
  Description: Deletes all Probes, ProbeSets, ProbeFeatures and 
               states associated with this ArrayChip
  Returntype : None
  Exceptions : Throws if ArrayChip not valid and stored
               Throws if ArrayChips are not of same class
  Caller     : General
  Status     : At risk

=cut

#This should be tied to a CS id!!!
#And analysis dependant?
#We may not want to delete alignment by different analyses?
#In practise the slice methods ignore analysis_id for this table
#So we currently never use this!
#So IMPORTED status should be tied to CS id and Analysis id?

sub rollback_ArrayChips{
  my ($self, $acs, $mode, $force, $keep_xrefs, $no_clean_up, $force_clean_up) = @_;
 
  #no_clean_up and force_clean_up allow analyze/optimize to be skipped until the last rollback
  #We could get around this by specifying all ArrayChips for all formats at the same time?
  #Need to implement in RollbackArrays
 
  $mode ||= 'probe';
  
  if($mode && ($mode ne 'probe' &&
			   $mode ne 'probe_feature' &&
			   $mode ne 'ProbeAlign' &&
			   $mode ne 'ProbeTranscriptAlign' &&
			   $mode ne 'probe2transcript')){
	throw("You have passed an invalid mode argument($mode), you must omit or specify either 'probe2transcript', 'probe', 'ProbeAlign, 'ProbeTranscriptAlign' or 'probe_feature' for all of the Align output");
  }
  
  if($force && ($force ne 'force')){
	throw("You have not specified a valid force argument($force), you must specify 'force' or omit");
  }

   if($keep_xrefs && ($keep_xrefs ne 'keep_xrefs')){
	throw("You have not specified a valid keep_xrefs argument($keep_xrefs), you must specify 'keep_xrefs' or omit");
  }


  if($keep_xrefs){

	if($mode eq 'probe' || $mode eq 'probe2transcript'){
	  throw("You cannot specify 'keep_xrefs' with mode $mode, you can only rollback features e.g. probe_feature, ProbeAlign or ProbeTranscriptAlign");
	}

	if($force){
	  throw("You cannot 'force' delete the probe2transcript xrefs and 'keep_xrefs' at the same time. Please specify just one.");
	}
  }




  my ($adaptor, $db, %classes);

  foreach my $ac(@$acs){
	$adaptor ||= $ac->adaptor || throw('ArrayChip must have an adaptor');
	$db      ||= $adaptor->db;
	$db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ArrayChip', $ac);

	if(! $ac->get_Array->class){
	  throw('The ArrayChip you are trying to rollback does not have a class attribute');
	}


	$classes{$ac->get_Array->class} = undef;

	#if($class && ($class ne $ac->get_Array->class)){
	#  throw('You can only rollback_ArrayChips for ArrayChips with the same class');
	#}
  }


  #This is always the case as we register the association before we set the Import status 
  #Hence the 2nd stage of the import fails as we have an associated ExperimentalChip
  #We need to make sure the ExperimentalChip and Channel have not been imported!!! 
  warn "NOTE: rollback_ArrayChips. Need to implement ExperimentlChip check, is the problem that ExperimentalChips are registered before ArrayChips imported?";  
  #Check for dependent ExperimentalChips
  #if(my @echips = @{$db->get_ExperimentalChipAdaptor->fetch_all_by_ArrayChip($ac)}){
#	my %exps;
#	my $txt = "Experiment\t\t\t\tExperimentalChip Unique IDs\n";
	
#	foreach my $ec(@echips){
#	  $exps{$ec->get_Experiment->name} ||= '';
	
#	  $exps{$ec->get_Experiment->name} .= "\t".$ec->unique_id;
#	}
	
#	map {$txt.= "\t".$_.":".$exps{$_}."\n"} keys %exps;
	
#	throw("Cannot rollback ArrayChip:\t".$ac->name.
#		  "\nFound Dependent Experimental Data:\n".$txt);
#  }
 
  
  my $ac_names = join(', ', (map { $_->name } @$acs));
  my $ac_ids   = join(', ', (map { $_->dbID } @$acs));


  $self->log("Rolling back ArrayChips $mode entries:\t$ac_names");
  my ($row_cnt, $probe_join, $sql);
  #$ac->adaptor->revoke_states($ac);#This need to be more specific to the type of rollback
  my $species = $db->species;
 
  if(!$species){
	throw('Cannot rollback probe2transcript level xrefs without specifying a species for the DBAdaptor');
  }
  #Will from registry? this return Homo sapiens?
  #Or homo_sapiens
  ($species = lc($species)) =~ s/ /_/;

  my $transc_edb_name = "${species}_core_Transcript";
  my $genome_edb_name = "${species}_core_Genome";

  #Maybe we want to rollback ProbeAlign and ProbeTranscriptAlign output separately so we 
  #can re-run just one part of the alignment step.

 
  #We want this Probe(Transcript)Align rollback available in the environment
  #So we can do it natively and before we get to the RunnableDB stage, 
  #where we would be trying multiple rollbacks in parallel
  #Wrapper script?
  #Or do we keep it simple here and maintain probe_feature wide rollback
  #And just the ProbeAlign/ProbeTranscriptAlign roll back in the environment?

    
  #We can restrict the probe deletes using the ac_id
  #We should test for other ac_ids using the same probe_id
  #Then fail unless we have specified force delete

  #These should be deleted for all other modes but only if force is set?
  #This may delete xrefs for other ArrayChips

  #The issues is if we need to specify force for one delete but don't want to delete something else?
  #force should only be used to delete upto and including the mode specified
  #no mode equates to probe mode
  #if no force then we fail if previous levels/modes have xrefs etc...


  #Let's grab the edb ids first and use them directly, this will avoid table locks on edb
  #and should also speed query up?


  if($mode eq 'probe2transcript' ||
	 $force){
	
	#Delete ProbeFeature UnmappedObjects	  
	$self->log("Deleting probe2transcript ProbeFeature UnmappedObjects");
	$sql = "DELETE uo FROM analysis a, unmapped_object uo, probe p, probe_feature pf, external_db e WHERE a.logic_name ='probe2transcript' AND a.analysis_id=uo.analysis_id AND p.probe_id=pf.probe_id and pf.probe_feature_id=uo.ensembl_id and uo.ensembl_object_type='ProbeFeature' and uo.external_db_id=e.external_db_id AND e.db_name ='${transc_edb_name}' AND p.array_chip_id IN($ac_ids)";
	$self->rollback_table($sql, 'unmapped_object', 'unmapped_object_id', $db, $no_clean_up);
	

	#Delete ProbeFeature Xrefs/DBEntries
	$self->log("Deleting probe2transcript ProbeFeature Xrefs");
	$sql = "DELETE ox FROM xref x, object_xref ox, probe p, probe_feature pf, external_db e WHERE x.external_db_id=e.external_db_id AND e.db_name ='${transc_edb_name}' AND x.xref_id=ox.xref_id AND ox.ensembl_object_type='ProbeFeature' AND ox.ensembl_id=pf.probe_feature_id AND pf.probe_id=p.probe_id AND ox.linkage_annotation!='ProbeTranscriptAlign' AND p.array_chip_id IN($ac_ids)";
	$self->rollback_table($sql, 'object_xref', 'object_xref_id', $db, $no_clean_up);


	#Probe/Set specific entries
	for my $xref_object('Probe', 'ProbeSet'){
	  $probe_join = ($xref_object eq 'ProbeSet') ? 'p.probe_set_id' : 'p.probe_id';
	  
	  #Delete Probe/Set UnmappedObjects

	  $self->log("Deleting probe2transcript $xref_object UnmappedObjects");
	 
	  $sql = "DELETE uo FROM analysis a, unmapped_object uo, probe p, external_db e WHERE a.logic_name='probe2transcript' AND a.analysis_id=uo.analysis_id AND uo.ensembl_object_type='${xref_object}' AND $probe_join=uo.ensembl_id AND uo.external_db_id=e.external_db_id AND e.db_name='${transc_edb_name}' AND p.array_chip_id IN($ac_ids)";
	  #.' and edb.db_release="'.$schema_build.'"'; 
	  $self->rollback_table($sql, 'unmapped_object', 'unmapped_object_id', $db, $no_clean_up);

	  #Delete Probe/Set Xrefs/DBEntries
	  $sql = "DELETE ox FROM xref x, object_xref ox, external_db e, probe p WHERE x.xref_id=ox.xref_id AND e.external_db_id=x.external_db_id AND e.db_name ='${transc_edb_name}' AND ox.ensembl_object_type='${xref_object}' AND ox.ensembl_id=${probe_join} AND p.array_chip_id IN($ac_ids)";
	  $self->log("Deleting probe2transcript $xref_object xref records");
	  $self->rollback_table($sql, 'object_xref', 'object_xref_id', $db, $no_clean_up);
	}
  }
  elsif(! $keep_xrefs){#Need to check for existing xrefs if not force
	#we don't know whether this is on probe or probeset level
	#This is a little hacky as there's not way we can guarantee this xref will be from probe2transcript
	#until we get the analysis_id moved from identity_xref to xref
	#We are also using the Probe/Set Xrefs as a proxy for all other Xrefs and UnmappedObjects
	#Do we need to set a status here? Would have problem rolling back the states of associated ArrayChips
	
	for my $xref_object('Probe', 'ProbeSet'){
	  
	  $probe_join = ($xref_object eq 'ProbeSet') ? 'p.probe_set_id' : 'p.probe_id';
	  
	  $row_cnt = $db->dbc->db_handle->selectrow_array("SELECT COUNT(*) FROM xref x, object_xref ox, external_db e, probe p WHERE x.xref_id=ox.xref_id AND e.external_db_id=x.external_db_id AND e.db_name ='${transc_edb_name}' and ox.ensembl_object_type='${xref_object}' and ox.ensembl_id=${probe_join} AND p.array_chip_id IN($ac_ids)");
	  
	  if($row_cnt){
		throw("Cannot rollback ArrayChips($ac_names), found $row_cnt $xref_object Xrefs. Pass 'force' argument or 'probe2transcript' mode to delete");
	  }
	  else{
		#$self->log("Found $row_cnt $xref_object Xrefs");
	  }
	}
  }
	

  #ProbeFeatures inc ProbeTranscriptAlign xrefs

  if($mode ne 'probe2transcript'){
	
	if(($mode eq 'probe' && $force) ||
	   $mode eq 'probe_feature'  ||
	   $mode eq 'ProbeAlign' ||
	   $mode eq 'ProbeTranscriptAlign'){


	  #Should really revoke some state here but we only have IMPORTED
   
	  #ProbeTranscriptAlign Xref/DBEntries
	  
	  #my (@anal_ids) = @{$db->get_AnalysisAdaptor->generic_fetch("a.module='ProbeAlign'")};
	  #Grrrr! AnalysisAdaptor is not a standard BaseAdaptor implementation
	  #my @anal_ids = @{$db->dbc->db_handle->selectall_arrayref('select analysis_id from analysis where module like "%ProbeAlign"')};
	  #@anal_ids = map {$_= "@$_"} @anal_ids;
	
	  if($mode ne 'ProbeAlign'){
		my $lnames = join(', ', (map { "'${_}_ProbeTranscriptAlign'" } keys(%classes)));

		$sql = "DELETE ox from object_xref ox, xref x, probe p, probe_feature pf, external_db e WHERE ox.ensembl_object_type='ProbeFeature' AND ox.linkage_annotation='ProbeTranscriptAlign' AND ox.xref_id=x.xref_id AND e.external_db_id=x.external_db_id and e.db_name='${transc_edb_name}' AND ox.ensembl_id=pf.probe_feature_id AND pf.probe_id=p.probe_id AND p.array_chip_id IN($ac_ids)";
		$self->log("Deleting ProbeFeature Xref/DBEntry records for:\t$lnames");
		$self->rollback_table($sql, 'object_xref', 'object_xref_id', $db, $no_clean_up);


		#Can't include uo.type='ProbeTranscriptAlign' in these deletes yet as uo.type is enum'd to xref or probe2transcript
		#will have to join to analysis and do a like "%ProbeTranscriptAlign" on the the logic name?
		#or/and ur.summary_description='Promiscuous probe'?

		$sql = "DELETE uo from unmapped_object uo, probe p, external_db e, analysis a WHERE uo.ensembl_object_type='Probe' AND uo.analysis_id=a.analysis_id AND a.logic_name in (${lnames}) AND e.external_db_id=uo.external_db_id and e.db_name='${transc_edb_name}' AND uo.ensembl_id=p.probe_id AND p.array_chip_id IN($ac_ids)";
	
		$self->log("Deleting UnmappedObjects for:\t${lnames}");
		$self->rollback_table($sql, 'unmapped_object', 'unmapped_object_id', $db, $no_clean_up);


		#Now the actual ProbeFeatures		
		$sql = "DELETE pf from probe_feature pf, probe p, analysis a WHERE a.logic_name in(${lnames}) AND a.analysis_id=pf.analysis_id AND pf.probe_id=p.probe_id AND p.array_chip_id IN($ac_ids)";
		$self->log("Deleting ProbeFeatures for:\t${lnames}");
		$self->rollback_table($sql, 'probe_feature', 'probe_feature_id', $db, $no_clean_up);
	  }

	  if($mode ne 'ProbeTranscriptAlign'){
		my $lnames = join(', ', (map { "'${_}_ProbeAlign'" } keys(%classes)));

		$sql = "DELETE uo from unmapped_object uo, probe p, external_db e, analysis a WHERE uo.ensembl_object_type='Probe' AND uo.analysis_id=a.analysis_id AND a.logic_name=(${lnames}) AND e.external_db_id=uo.external_db_id and e.db_name='${genome_edb_name}' AND uo.ensembl_id=p.probe_id AND p.array_chip_id IN($ac_ids)";
		$self->log("Deleting UnmappedObjects for:\t${lnames}");
		$self->rollback_table($sql, 'unmapped_object', 'unmapped_object_id', $db, $no_clean_up);
		

		$sql = "DELETE pf from probe_feature pf, probe p, analysis a WHERE a.logic_name in(${lnames}) AND a.analysis_id=pf.analysis_id AND pf.probe_id=p.probe_id AND p.array_chip_id IN($ac_ids)";
		$self->log("Deleting ProbeFeatures for:\t${lnames}");
		$self->rollback_table($sql, 'probe_feature', 'probe_feature_id', $db, $no_clean_up);
	  }
	}
	else{
	  #Need to count to see if we can carry on with a unforced probe rollback?
	  #Do we need this level of control here
	  #Can't we assume that if you want probe you also want probe_feature?
	  #Leave for safety, at least until we get the dependant ExperimetnalChip test sorted
	  #What about if we only want to delete one array from an associated set?
	  #This would delete all the features from the rest?
	  
	  $sql = "select count(*) from object_xref ox, xref x, probe p, external_db e WHERE ox.ensembl_object_type='ProbeFeature' AND ox.linkage_annotation='ProbeTranscriptAlign' AND ox.xref_id=x.xref_id AND e.external_db_id=x.external_db_id and e.db_name='${transc_edb_name}' AND ox.ensembl_id=p.probe_id AND p.array_chip_id IN($ac_ids)";
	  $row_cnt =  $db->dbc->db_handle->selectrow_array($sql);
	  
	  if($row_cnt){
		throw("Cannot rollback ArrayChips($ac_names), found $row_cnt ProbeFeatures. Pass 'force' argument or 'probe_feature' mode to delete");
	  }
	   else{
		 $self->log("Found $row_cnt ProbeFeatures");
	  }
	}
	
	if($mode eq 'probe'){
	  #Don't need to rollback on a CS as we have no dependant EChips?
	  #Is this true?  Should we enforce a 3rd CoordSystem argument, 'all' string we delete all?
	  
	  foreach my $ac(@$acs){
		$ac->adaptor->revoke_states($ac);#Do we need to change this to revoke specific states?
		#Current states are only IMPORTED, so not just yet, but we could change this for safety?
	  }
	  
	  #ProbeSets
	  $sql = "DELETE ps from probe p, probe_set ps where p.array_chip_id IN($ac_ids) and p.probe_set_id=ps.probe_set_id";
	  $self->rollback_table($sql, 'probe_set', 'probe_set_id', $db, $no_clean_up);

	  #Probes
	  $sql = "DELETE from probe where array_chip_id IN($ac_ids)";  
	  $self->rollback_table($sql, 'probe', 'probe_id', $db, $no_clean_up);	  
	}
  }
  
  $self->log("Finished $mode roll back for ArrayChip:\t$ac_names");
  return;
}


#This will just fail silently if the reset value
#Is less than the true autoinc value
#i.e. if there are parallel inserts going on
#So we can never assume that the $new_auto_inc will be used


sub rollback_table{
  my ($self, $sql, $table, $id_field, $db, $no_clean_up, $force_clean_up) = @_;

  my $row_cnt;
  eval { $row_cnt = $db->dbc->do($sql) };
  
  if($@){
  	throw("Failed to rollback table $table using sql:\t$sql\n$@");
  }

  $row_cnt = 0 if $row_cnt eq '0E0';
  $self->log("Deleted $row_cnt $table records");
  
  if($force_clean_up ||
	 ($row_cnt && ! $no_clean_up)){
	$self->refresh_table($table, $id_field, $db);
  }

  return;
}

#Now separated so that we can do this once at the end of a rollback of many Sets

sub refresh_table{
  my ($self, $table, $id_field, $db) = @_;

  #This only works if the new calue is available
  #i.e. do not need lock for this to be safe
  $self->reset_table_autoinc($table, $id_field, $db) if $id_field;
  
  $self->log("Optimizing and Analyzing $table");	
  
  $db->dbc->do("optimize table $table");#defrag data, sorts indices, updates table stats
  $db->dbc->do("analyze  table $table");#analyses key distribution
  
  return;
}



sub reset_table_autoinc{
  #Is this called elsewhere or can we merge with
  my($self, $table_name, $autoinc_field, $db) = @_;

  if(! ($table_name && $autoinc_field && $db)){
	throw('You must pass a table_name and an autoinc_field to reset the autoinc value');
  }

  if(! (ref($db) && $db->isa('Bio::EnsEMBL::DBSQL::DBAdaptor'))){
	throw('Must pass a valid Bio::EnsEMBL::DBSQL::DBAdaptor');
  }

  #my $sql = "show table status where name='$table_name'";
  #my ($autoinc) = ${$db->dbc->db_handle->selectrow_array($sql)}[11];
  #11 is the field in the show table status table
  #We cannot select just the Auto_increment, so this will fail if the table format changes

  #Why do we need autoinc here?

  my $sql = "select $autoinc_field from $table_name order by $autoinc_field desc limit 1";
  my ($current_auto_inc) = $db->dbc->db_handle->selectrow_array($sql);
  my $new_autoinc = ($current_auto_inc) ? ($current_auto_inc + 1) : 1;  
  $sql = "ALTER TABLE $table_name AUTO_INCREMENT=$new_autoinc";
  $db->dbc->do($sql);
  return;
}

#$qry=mysql_query("show table status where name='newblog'") or die (mysql_error());
#$row=mysql_fetch_array($qry);
#$newtid=$row[10];
#Reset autoincrement

#ALTER TABLE tablename AUTO_INCREMENT = 1




=head2 get_core_display_name_by_stable_id

  Args [1]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Args [2]   : stable ID from core DB.
  Args [3]   : stable feature type e.g. gene, transcript, translation
  Example    : $self->validate_and_store_feature_types;
  Description: Builds a cache of stable ID to display names.
  Returntype : string - display name
  Exceptions : Throws is type is not valid.
  Caller     : General
  Status     : At risk

=cut

# --------------------------------------------------------------------------------
# Build a cache of ensembl stable ID -> display_name
# Return hashref keyed on {$type}{$stable_id}
#Need to update cache if we're doing more than one 'type' at a time
# as it will never get loaded for the new type!

sub get_core_display_name_by_stable_id{
  my ($self, $cdb, $stable_id, $type) = @_;

  $type = lc($type);

  if($type !~ /(gene|transcript|translation)/){
	throw("Cannot get display_name for stable_id $stable_id with type $type");
  }
  
  if(! exists $self->{'display_name_cache'}->{$stable_id}){
	($self->{'display_name_cache'}->{$stable_id}) = $cdb->dbc->db_handle->selectrow_array("SELECT x.display_label FROM ${type}_stable_id s, $type t, xref x where t.display_xref_id=x.xref_id and s.${type}_id=t.${type}_id and s.stable_id='${stable_id}'");
  }

  return $self->{'display_name_cache'}->{$stable_id};
}


=head2 get_core_stable_id_by_display_name

  Args [1]   : Bio::EnsEMBL::DBSQL::DBAdaptor
  Args [2]   : display name (e.g. from core DB or GNC name)
  Example    : 
  Description: Builds a cache of stable ID to display names.
  Returntype : string - gene stable ID
  Exceptions : None
  Caller     : General
  Status     : At risk

=cut

# --------------------------------------------------------------------------------
# Build a cache of ensembl stable ID -> display_name
# Return hashref keyed on {$type}{$stable_id}
#Need to update cache if we're doing more than one 'type' at a time
# as it will never get loaded for the new type!

sub get_core_stable_id_by_display_name{
  my ($self, $cdb, $display_name) = @_;

  #if($type !~ /(gene|transcript|translation)/){
#	throw("Cannot get display_name for stable_id $stable_id with type $type");
#  }
  
  if(! exists $self->{'stable_id_cache'}->{$display_name}){
	($self->{'stable_id_cache'}->{$display_name}) = $cdb->dbc->db_handle->selectrow_array("SELECT s.stable_id FROM gene_stable_id s, gene g, xref x where g.display_xref_id=x.xref_id and s.gene_id=g.gene_id and x.display_label='${display_name}'");
  }

  return $self->{'stable_id_cache'}->{$display_name};
}



#Can we do this for several common sets of params
#And therefore cut down the code that we need to write for every new script?

#sub generate_efgdb_from_params{
#  my ($self, $argv) = @_;
#}


1;

