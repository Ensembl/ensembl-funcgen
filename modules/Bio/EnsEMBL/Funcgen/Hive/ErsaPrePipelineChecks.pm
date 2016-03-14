package Bio::EnsEMBL::Funcgen::Hive::ErsaPrePipelineChecks;

use strict;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
    my $self = shift;

    use Bio::EnsEMBL::Utils::Logger;
    my $logger = Bio::EnsEMBL::Utils::Logger->new();    
    my $error_msg;
    
    # - Has PERL5LIB been exported
    # - Has JAVA_HOME been exported
    # Are certain directories writeable? Like the one with the chromosome file in it
    
    my @exported_variables = qw( PERL5LIB JAVA_HOME R_LIBS CLASSPATH );
    my $cmd;
    
    ENVIRONMENT_VARIABLE:
    foreach my $exported_variable (@exported_variables) {
    
      $cmd = qq(set | grep $exported_variable);
      my $set_command = `$cmd`;
      if (! $set_command) {
	$error_msg .= "$exported_variable has not been set!\n";
	next ENVIRONMENT_VARIABLE;
      }

      #
      # If "foo" has been exported, the command
      #     
      # export | grep foo
      #
      # will yield this:
      #
      # declare -x foo="bar"
      # 
      $cmd = qq(export | grep $exported_variable);
      my $declare_command = `$cmd`;
      
      if (! $declare_command) {
	$error_msg .= "$exported_variable has not been exported!\n";
      }
      
    }

    system("which bwa > /dev/null");    
    if ($?) {
      $error_msg .= "Can't find bwa in path.\n";
    }
    
    system("which R > /dev/null");    
    if ($?) {
      $error_msg .= "Can't find R in path.\n";
    }
    
    my $picard_output = `java picard.cmdline.PicardCommandLine`;
    if ($picard_output !~ /^USAGE: PicardCommandLine/) {
      $error_msg .= qq(Can't run picard! "java picard.cmdline.PicardCommandLine".\n);
    }
    
    system("which Rscript > /dev/null");    
    if ($?) {
      $error_msg .= "Can't find Rscript in path.\n";
    } else {
    
      # Should be something like 
      # "R scripting front-end version 3.2.2 (2015-08-14)"
      #
      my $version_string = `Rscript --version`;
      
      $version_string_found = $version_string =~ /R scripting front-end version (.+?) /;
      
      if ($version_string_found) {
      
	my $version = $1;
	
	# If version string starts with 3.2, we are probably ok.
	#
	my $version_ok =
	     $version =~ /^3\.2/
	  || $version =~ /^3\.3/
	  || $version =~ /^3\.4/
	;
	
	if (!$version_ok) {
	  $error_msg .= "Rscript on PATH has version $version, this may cause trouble.\n";
	}
	
      } else {
	$error_msg .= "Can't find version from Rscript in string: '$version'\n";
      }
    }
    

    
    if ($error_msg) {
      die($error_msg);
    }
    $logger->info("Pre pipeline checks have completed successfully. You can now run the rest of the pipeline.\n");
}

1;
