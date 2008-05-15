
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

=head1 NOTES


=head1 AUTHOR(S)

Nathan Johnson, njohnson@ebi.ac.uk


=cut

################################################################################

package Bio::EnsEMBL::Funcgen::Utils::Helper;
#put in Utils?
use Bio::Root::Root;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw (throw stack_trace);
#use Devel::Timer;
use Carp;#? Can't use unless we can get it to redirect
use File::Basename;

use strict;
use vars qw(@ISA);
@ISA = qw(Bio::Root::Root);

################################################################################

=head2 new

 Description : Constructor method to create a new object with passed or
               default attributes.

 Arg  [1]    : hash containing optional attributes :-
                 log_file    - name of log file (default = undef -> STDOUT)
                 debug_level - level of detail of debug message [1-3] (default = 0 = off)
                 debug_file  - name of debug file (default = undef -> STDERR)

 ReturnType  : Helper

 Example     : my $Helper = Bio::EnsEMBL::Helper->new(
                                                      debug_level => 3,
                                                      debug_file  => "/tmp/efg.debug",
                                                      log_file    => "/tmp/efg.log",
                                                     );

 Exceptions  : throws exception if failed to open debug file
             : throws exception if failed to open log   file

=cut

################################################################################

sub new{
    my ($caller, %args) = @_;

    my ($self,%attrdata,$attrname,$argname);
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
    foreach $attrname (keys %attrdata){
        ($argname = $attrname) =~ s/^_//; # remove leading underscore
        $self->{$attrname} = (exists $args{$argname}) ? $args{$argname} : $attrdata{$attrname};
    }

	$self->{'_tee'} = 1 if $self->{'_no_log'};
	#should we undef log_file here too?
	#This currently only turns off default logging

	$self->{_default_log_dir} ||= $ENV{'HOME'}.'/logs';


    # DEBUG OUTPUT & STDERR
    if(defined $self->{_debug_level} && $self->{_debug_level}){
        $main::_debug_level = $self->{_debug_level};
		
        if(defined $self->{_debug_file}){
			$main::_debug_file = $self->{_debug_file};
			  			  
            open(DBGFILE,">>".$self->{_debug_file})
			  or throw("Failed to open debug file : $!");

			#open (DBGFILE, "<STDERR | tee -a ".$self->{_debug_file});#Mirrors STDERR to debug file
        }
        else{
            open(DBGFILE,">&STDERR");
        }

        select DBGFILE; $| = 1;  # make debug file unbuffered

        $self->debug(1,"Debugging started ".localtime()." on $0 at level ".$self->{_debug_level}." ...");
    }

	# LOG OUTPUT
	if (defined $self->{_log_file}){
	  $main::_log_file = $self->{_log_file};
		
	  my $log_file = '>>'.$self->{'_log_file'};

	  #we need to implment tee here
	  if($self->{'_tee'}){
	    #we're not resetting $main::_tee here, we only use it once.
	    $log_file = '| tee -a '.$self->{_log_file};
	  }

	  open(LOGFILE, $log_file)
	    or throw("Failed to open log file : $log_file\nError: $!");
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
		my $log_file = '>>'.$self->{'_log_file'};
		warn "No log file defined, defaulting to:\t".$self->{'_log_file'}."\n";

		#we should still tee here
		if($self->{'_tee'}){
		  #we're not resetting $main::_tee here, we only use it once.
		  $log_file = '| tee -a '.$log_file;
		}

		open(LOGFILE,  $log_file)
		  or throw("Failed to open log file : $log_file\nError: $!");

	  }
	  else{
		open(LOGFILE,">&STDOUT");
	  }
	}

	select LOGFILE; $| = 1;  # make log file unbuffered

	$self->log("\n\nLogging started at ".localtime()."...");

    # RESET STDOUT TO DEFAULT
    select STDOUT; $| = 1; 

    $self->debug(2,"Helper class instance created.");

    return ($self);
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

    if($self->{_log_file}){
        $self->log("Logging complete ".localtime().".");

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
  my ($self, $message, $mem) = @_;

  if($mem){
	$message.= " :: ".`ps -p $$ -o vsz |tail -1`;
	chomp $message;
	$message .= " KB";
  }

  print LOGFILE "::\t$message\n";

  # Add to debug file if not printing to STDERR?
  # only if verbose?
  # this would double print everything to STDOUT if tee and debug has not redefined STDERR

  $self->debug(1,$message);
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
  my ($self, $message, $mem) = @_;

  $self->log("::\t$message\t::\t::\n", $mem);
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
        print DBGFILE "debug $message\t: [$$ - $prog_name:$prog_line  $call_name:$call_line]\n";

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
		print DBGFILE Data::Dumper::Dumper(\$hashref)."\n";
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
  return [$dbname[($#dbname -1)], $dbname[($#dbname -1)]];
}



=head2 define_and_validate_sets

  Arg [1]    : hash - set constructor parameters:
                            -dbadaptor    Bio::EnsEMBL::Funcgen::DBAdaptor
                            -name         Data/FeatureSet name to create
                            -feature_type Bio::EnsEMBL::Funcgen::FeatureType
                            -cell_type    Bio::EnsEMBL::Funcgen::CellType
                            -analysis     FeatureSet Bio::EnsEMBL::Analysis
                            -type         e.g. annotated or regulatory
                            -description  FeatureSet description
                            -rollback     Boolean - Forces rollback of previously imported data
                            -append       Boolean - Forces import on top of previously imported dataç
  Example    : my $fset = $self->define_and_validate_Set(%params);
  Description: Checks whether set is already in DB based on set name, rolls back features 
               if roll back flag set. Or creates new DB if not present.
  Returntype : Bio::EnsEMBL::Funcgen::DataSet
  Exceptions : Throws if DBAdaptor param not valid
  Caller     : Importers and Parsers
  Status     : At risk

=cut

sub define_and_validate_sets{
  my $self = shift;

  my ($name, $anal, $ftype, $ctype, $type, $append, $rollback, $db, $ssets, $description)
    = rearrange(['NAME', 'ANALYSIS', 'FEATURE_TYPE', 'CELL_TYPE', 'TYPE', 'APPEND', 'ROLLBACK', 'DBADAPTOR', 'SUPPORTING_SETS', 'DESCRIPTION'], @_);

  #We need an append flag to allow addition of Features to a pre-existing feature set
  #We should implement rearrange here, will this capture any ill-defined parameters
  #add db, rollback and append to params


  #Check mandatory params
  if(! (ref($db) && $db->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'))){
	throw('Must provide a valid Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor');
  }

  throw('Must provide a -name ') if(! defined $name);
  
  #Not necessarily, just do rollback then append?
  #Generating sets for an ExpSet will always have append set
  #This could be valid for generically grabing/creating sets for adding new supporting sets e.g. reg build
  #throw('-append and -rollback are mutually exclusive') if $rollback && $append;
  
  #warn only for append?
  warn('You are defining a pre-existing FeatureSet without rolling back'.
	   ' previous data, this could result in data duplication') if $append;


  throw('Must provide a -type e.g. annotated, external or regulatory') if(! defined $type);
  #Check for annotated, external, regulatory etc here?
  #Should never be external as we don't have DataSets for external sets?
  
  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureType',  $ftype);
  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::CellType',  $ctype);
  $db->is_stored_and_valid('Bio::EnsEMBL::Analysis',  $anal);




  my $dset_adaptor = $db->get_DataSetAdaptor;
  my $fset_adaptor = $db->get_FeatureSetAdaptor;
  
  my $dset = $dset_adaptor->fetch_by_name($name);


  my ($fset);

  if(defined $dset){
	$self->log('Found Stored DataSet '.$dset->name);

	$fset = $dset->product_FeatureSet;
	#Here we have the possiblity that a feature_set with a different name may have been associated with the DataSet


	if(defined $fset){
	  $self->log("Found associated product FeatureSet:\t".$fset->name);
	  
	  #if(! $clobber && 
	  if($fset->name ne $name){
		throw('Invalid product FeatureSet name ('.$fset->name.') for DataSet ('.$name.'). Rollback will overwrite the FeatureSet and mismatched name will be retained.');
		#Need to clobber both or give explicit name for datasets or rename dataset???
		#Force this throw for now, make this fix manual as we may end up automatically overwriting data
	  }  
	}

	#check supporting_sets here if defined
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
		throw('There is a (name/type/anal) mismatch between the supplied supporting_sets and the'.
			  'supporting_sets in the DB for DataSet '.$dset->name."\nStored:\t"
			  .join(', ', (map $_->name, @stored_ssets))."\nSupplied supporting_sets:\t"
			  .join(', ', (map $_->name, @sorted_ssets)));
	  }
	}
	else{
	  warn("Skipping validating of supporting sets for Data/FeatureSet definition:\t".$name);
	}
  }


  if(! defined $fset){
	#Try and grab it anyway just in case it has been orphaned somehow
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

	if($rollback){
	  $self->rollback_FeatureSet($fset);
	  #Not forcing delete here as this may be used as a supporting set itself.
	}else{
	  throw('Found pre-existing FeatureSet '.$fset->name.'. Maybe you want to specify the rollback flag?');
	}
  }
  else{
	#create a new one
	$self->log("Creating new FeatureSet:\t".$name);

	$fset = Bio::EnsEMBL::Funcgen::FeatureSet->new(
												   -name         => $name,
												   -feature_type => $ftype,
												   -cell_type    => $ctype,
												   -analysis     => $anal,
												   -type         => $type,
												   -description  => $description,
												  );
	($fset) = @{$fset_adaptor->store($fset)};
  }


  #Create/Update the DataSet
  if(defined $dset){
	
	if(! defined $dset->product_FeatureSet){
	  $self->log("Updating DataSet with new product FeatureSet:\t".$fset->name);
	  ($dset) = @{$dset_adaptor->store_updated_sets($dset->product_FeatureSet($fset))};
	}
  }
  else{
	$self->log("Creating new DataSet:\t".$name);
	$dset = Bio::EnsEMBL::Funcgen::DataSet->new(
												-name => $name,
												-feature_set => $fset,
												-supporting_sets => $ssets,
											   );
	($dset) = @{$dset_adaptor->store($dset)};
  }
  
  return $dset;
}


#Rollback/load methods migrated from DBAdaptor

#Do we need to add a rolling back status?
#Set before and remove afterwards?

#These assume the parent class has a db attr
#do we need a $self->can(db) test here

=head2 rollback_FeatureSet

  Arg [1]    : Bio::EnsEMBL::Funcgen::FeatureSet
  Arg [2]    : boolean - Force delete flag
  Example    : $self->rollback_FeatureSet($fset);
  Description: Deletes all status and feature entries for this FeatureSet.
               Checks whether FeatureSet is a supporting set in any other DataSet.
  Returntype : none
  Exceptions : Throws if any deletes fails or if db method unavailable
  Caller     : Importers and Parsers
  Status     : At risk

=cut


sub rollback_FeatureSet{
  my ($self, $fset, $force_delete) = @_;

  #Need to test before we do adaptor call? Cyclical dependency here :|
  #We need to implement this method locally
  #This is because we don't force Helper to have a DB attribute
  #Maybe we should?
  #We're always going to have a DB so why not?
  #Because we might want to use the Helper to Log before we can create the DB?

  my $adaptor = $fset->adaptor || throw('FeatureSet must have an adaptor');
  my $db = $adaptor->db;
  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::FeatureSet', $fset);

  $self->log('Rolling back '.$fset->type." FeatureSet:\t".$fset->name);

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
  $fset->adaptor->revoke_states($fset);

  #Remove feature
  my $sql = 'DELETE from '.$fset->type.'_feature where feature_set_id='.$fset->dbID;  

  if(! $db->dbc->do($sql)){
	throw('Failed to rollback '.$fset->type."_features for FeatureSet:\t".$fset->name.' (dbID:'.$fset->dbID.')');
  }

  return;
}


=head2 rollback_ResultSet

  Arg[1]     : Bio::EnsEMBL::Funcgen::ResultSet
  Arg[2]     : Boolean - optional flag to roll back IMPORT set results
  Example    : $self->rollback_ResultSet($rset);
  Description: Deletes all status. chip_channel and result_set entries for this ResultSet.
               Will also rollback_results sets if rollback_results specified.  This will also
               update or delete associated ResultSets where appropriate.
               If an associated
               I
  Returntype : Arrayref 
  Exceptions : Throws if ResultSet not valid
               Throws is result_rollback flag specified but associated product FeatureSet found.
  Caller     : General
  Status     : At risk

=cut



sub rollback_ResultSet{
  my ($self, $rset, $rollback_results) = @_;
  
  #what about?
  if(! (ref($rset) && $rset->can('adaptor') && defined $rset->adaptor)){
	throw('Must provide a valid stored Bio::EnsEMBL::ResultSet');
  }
  
  #We're still validating against itself??
  #And reciprocating part of the test :|
  my $sql;
  my $db = $rset->adaptor->db;
  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ResultSet', $rset);
  $self->log("Rolling back ResultSet:\t".$rset->name);
  
  #Is this ResultSet used in a DataSet?
  my $dset_adaptor = $self->db->get_DataSetAdaptor;
  my $rset_adaptor = $self->db->get_ResultSetAdaptor;
  my @skipped_sets;
  
  foreach my $dset(@{$dset_adaptor->fetch_all_by_supporting_set($rset)}){
	
	if (defined $dset){
	  $self->log('Found linked DataSet('.$dset->name.") for ResultSet:\t".$rset->log_label);
	  
	  if(my $fset = $dset->product_FeatureSet){
		$self->log('Skipping rollback. Found product FeatureSet('.$fset->name.") for supporting ResultSet:\t".$rset->log_label);
		warn('Add more info on logs here on which script to use to edit the DataSet');
		
		@skipped_sets = ($rset,$dset);

		#What impact does this have on result_rollback?
		#None as we never get there
		#But what if we have specified rollback results?
		#We should throw here as we can't perform the rollback
		
		if($rollback_results){
		  throw("Could not rollback supporting ResultSet and results for:\t".$rset->log_label.
				"\nManually resolve the supporting/feature set relationship or omit the ".
				"rollback_results argument if you simply want to redefine the ResultSet without loading any new data");
		}
	  }
	  else{
		#Found rset in dset, but not yet processed so can remove safely.
		$self->log("Removing supporting ResultSet from DataSet:\t".$dset->name."\tResultSet:".$rset->log_label);
		$sql = 'DELETE from supporting_set where data_set_id='.$dset->dbID.
		  ' and type="result" and supporting_set_id='.$rset->dbID;
		$db->dbc->do($sql);
	  }
	}
  }
  

  #Now do similar for all associated ResultSets
  if(! @skipped_sets){

	#Rollback status first

	#We need to move this as we could still fail after here, but before we've deleted anything else


	#Use revoked states here!
	$rset->adaptor->revoke_states($rset);
	
	#Rollback results if required
	if($rollback_results){
	  $self->log("Rolling back result for ResultSet:\t".$rset->log_label);
	  #First we need to check whether these cc_ids are present in other result sets.
	  #Get all associated data_sets
	  #checking for other product_FeatureSets for given cc_ids
	  #If we are rolling back results then we delete the cc_ids from 
	  #associated result_sets and dompletely delete any ResultSets 
	  #which are a subset of this one
	  #Warn if we delete an off target set
	  #Warn if we don't rollback results due associated supporting ResultSet
	  #which has been used to produce a product FeatureSet.
		
	  my @assoc_rsets = @{$rset_adaptor->fetch_all_linked_by_ResultSet($rset)};
	  my $feature_supporting = 0;
	  
	  foreach my $assoc_rset(@assoc_rsets){
		
		foreach my $dset(@{$dset_adaptor->fetch_all_by_supporting_set($assoc_rset)}){
		  
		  if(my $fset = $dset->product_FeatureSet){
			$feature_supporting++;
			$self->log('Found product FeatureSet('.$fset->name.
					   ") for associated supporting ResultSet:\t".$rset->log_label);
		  }
		}					
	  }
		

	  if(! $feature_supporting){
		#This also handles Echip status rollback
		$self->rollback_results($rset->chip_channel_ids);
		$self->log('Removing chip_channel entries from associated ResultSets');
		
		#Now remove cc_ids from associated rsets.
		foreach my $assoc_rset(@assoc_rsets){
		  $sql = 'DELETE from chip_channel where result_set_id='.$assoc_rset->dbID.
			' and chip_channel_id in('.join', ', @{$assoc_rset->chip_channel_ids}.')';
		  $db->dbc->do($sql);
		  
		  # we need to delete complete subsets from the result_set table.
		  my $subset = 1;
		  
		  foreach my $cc_id(@{$assoc_rset->chip_channel_ids}){
			
			if(! grep/$cc_id/, @{$rset->chip_channel_ids}){
			  $subset = 0;
			  last;
			}
		  }
			
		  #Found complete subset so can delete
		  if($subset){
			$self->log("Deleting associated subset ResultSet:\t".$assoc_rset->log_label);
			
			#Delete status entries first
			$assoc_rset->adaptor->revoke_states($assoc_rset);
			
			#All cc records will have already been deleted
			$sql = 'DELETE from result_set where result_set_id='.$assoc_rset->dbID;
			$db->dbc->do($sql);
		  }
		}


		#Now warn about Echips in Experiments which may need removing.
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
	$sql = 'DELETE from chip_channel where result_set_id='.$rset->dbID;
	$db->dbc->do($sql);
	$sql = 'DELETE from result_set where result_set_id='.$rset->dbID;
	$db->dbc->do($sql);
  }

  return \@skipped_sets;
}





=head2 rollback_ExperimentalSet

  Arg[1]     : Bio::EnsEMBL::Funcgen::ExperimentalSet
  Example    : $self->rollback_ExperimentalSet($eset);
  Description: Deletes all status entries for this ExperimentalSet and it's ExperimentalSubSets
  Returntype : none
  Exceptions : Throws if any deletes fails or if db method unavailable
  Caller     : Importers and Parsers
  Status     : At risk

=cut


sub rollback_ExperimentalSet{
  my ($self, $eset, $force_delete) = @_;


  #Need to implement force_delete!!!!!!!!!!!!!!!!!!!!!!

  my $adaptor = $eset->adaptor || throw('ExperimentalSet must have an adaptor');
  my $db = $adaptor->db;
  

  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ExperimentalSet', $eset);

  $self->log("Rolling back ExperimentSet:\t".$eset->name);

  #ExperimentalSubSets
  foreach my $esset(@{$eset->get_subsets}){
	$esset->adaptor->revoke_states($esset);
  }

  #ExperimentalSet
  $eset->adaptor->revoke_states($eset);

  
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
  

  if(! scalar(@cc_ids) >0){
	throw('Must pass an array ref of chip_channel ids to rollback');
  }
  
  #Rollback status entries
  #Cannot use revoke_states here?
  #We can if we retrieve the Chip or Channel first
  #Add to ResultSet adaptor
  my $sql = 'DELETE s from status s, chip_channel cc WHERE cc.chip_channel_id IN ('.join(',', @cc_ids).
	') AND cc.table_id=s.table_id AND cc.table_name=s.table_name';
  
  if(! $self->db->dbc->do($sql)){
	throw("Status rollback failed for chip_channel_ids:\t@cc_ids\n".$self->db->dbc->db_handle->errstr());
  }


  #Rollback result entries
  $sql = 'DELETE from result where chip_channel_id in ('.join(',', @cc_ids).');';
  if(! $self->db->dbc->do($sql)){
	throw("Results rollback failed for chip_channel_ids:\t@cc_ids\n".$self->db->dbc->db_handle->errstr());
  }

  return;
}


=head2 rollback_ArrayChip

  Arg[1]     : Bio::EnsEMBL::Funcgen::ArrayChip
  Example    : $self->rollback_ArrayChip($achip);
  Description: Deletes all Probes, ProbeSets, ProbeFeatures and 
               states associated with this ArrayChip
  Returntype : None
  Exceptions : Throws if ArrayChip not valid and stored
  Caller     : General
  Status     : At risk

=cut


sub rollback_ArrayChip{
  my ($self, $ac) = @_;

  my $adaptor = $ac->adaptor || throw('ArrayChip must have an adaptor');
  my $db = $adaptor->db;
  $db->is_stored_and_valid('Bio::EnsEMBL::Funcgen::ArrayChip', $ac);


  #This is always the case as we register the association before we set the Import status 
  #Hence the 2nd stage of the import fails as we have an associated ExperimentalChip
  #We need to make sure the ExperimentalChip and Channel have not been imported!!!

  warn "Need to implement EC/Channel import check before rolling back ArrayChip";

  

  #Check for dependent ExperimentalChips
#  if(my @echips = @{$db->get_ExperimentalChipAdaptor->fetch_all_by_ArrayChip($ac)}){
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
  


  $self->log("Rolling back ArrayChip:\t".$ac->name);
  $ac->adaptor->revoke_states($ac);

  
  #Need to rollback feature here too!
  #Don't need to rollback on a CS as we have no dependant EChips?
  #Is this true?  Should we enforce a 3rd CoordSystem argument, 'all' string we delete all?


  my $sql = 'DELETE ps from probe p, probe_set ps where p.array_chip_id='.$ac->dbID().' and p.probe_set_id=ps.probe_set_id';
  if(! $db->dbc->do($sql)){
	throw("ProbeSet rollback failed for ArrayChip:\t".$ac->name()."\n".$self->dbc->db_handle->errstr());
  }

  $sql = 'DELETE pf from probe p, probe_feature pf where p.array_chip_id='.$ac->dbID().' and p.probe_id=pf.probe_id';
  if(! $db->dbc->do($sql)){
	throw("ProbeFeature rollback failed for ArrayChip:\t".$ac->name()."\n".$self->dbc->db_handle->errstr());
  }


  $sql = 'DELETE from probe where array_chip_id='.$ac->dbID();
  
  if(! $db->dbc->do($sql)){
	throw("Probe rollback failed for ArrayChip:\t".$ac->name()."\n".$self->dbc->db_handle->errstr());
  }
     
  return;
}




1;

