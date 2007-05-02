
=head1 NAME

Bio::EnsEMBL::Funcgen::Helper
  
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

=item B<-instance|i>

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

package Bio::EnsEMBL::Funcgen::Helper;
#put in Utils?
use Bio::Root::Root;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Exception qw (throw);
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


    # objects private data and default values
    %attrdata = (
		 _tee          => $main::_tee,
		 _debug_level  => $main::_debug_level,
		 _debug_file   => $main::_debug_file,
		 _log_file     => $main::_log_file,#default should be set in caller
		);

    # set each class attribute using passed value or default value
    foreach $attrname (keys %attrdata){
        ($argname = $attrname) =~ s/^_//; # remove leading underscore
        $self->{$attrname} = (exists $args{$argname}) ? $args{$argname} : $attrdata{$attrname};
    }

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
		
	  my $log_file = ">>".$self->{'_log_file'};

	  #we need to implment tee here
	  if($self->{'_tee'}){
	    #we're not resetting $main::_tee here, we only use it once.
	    $log_file = "| tee -a ".$self->{_log_file};
	  }

	  open(LOGFILE, $log_file)
	    or throw("Failed to open log file : $log_file\nError: $!");
	}
	else{
	  open(LOGFILE,">&STDOUT");
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

  if ($self->{_debug_level} >= 3){

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

  #This method is just to provide standard checking for specific get_data methods
	

  


  if(defined $data_name){
    throw("Defs data name $data_name for type '$data_type' does not exist\n") if (! exists $self->{"${data_type}"}{$data_name});
  }else{
    throw("Defs data type $data_type does not exist\n") if (! exists $self->{"${data_type}"});
  }
  
  return (defined $data_name) ? $self->{"${data_type}"}{$data_name} : $self->{"${data_type}"};
}


sub Timer{
	my ($self) = shift;

	$self->{'_timer'} = new Devel::Timer()  if(! defined $self->{'_timer'});

	return $self->{'_timer'};
	
}


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



1;

