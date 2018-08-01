package Bio::EnsEMBL::Funcgen::Utils::GoodUtils;

use warnings;
use strict;
use Carp;

use base qw( Exporter );
use vars qw( @EXPORT_OK );

@EXPORT_OK = qw(
  run_cmd
  create_production_name
  create_species_assembly_path
);

sub create_production_name {
    my $name = shift;
    
    my $max_length = 70;
    
    $name =~ s/[^a-zA-Z0-9_]/_/g;
    my $shortened = substr( $name, 0, $max_length );

    return $shortened;
}

=head2 create_species_assembly_path

  Bio::EnsEMBL::Registry->load_all($registry);
  my $species_assembly_path = create_species_assembly_path($species);

=cut
sub create_species_assembly_path {

  my $species = shift;
  
    my $coordsystem_adaptor = Bio::EnsEMBL::Registry->get_adaptor($species, 'core', 'coordsystem');
    my $default_chromosome_coordsystem = $coordsystem_adaptor->fetch_by_name('chromosome');
    my $default_assembly = $default_chromosome_coordsystem->version;

    #die ("\n\n------------------> $species, $default_assembly\n\n");
    
    return $species . '/' . $default_assembly;
}

sub run_cmd {

    my $cmd  = shift;
    my $test_for_success = shift;

    my $param_ok = (!defined $test_for_success) || (ref $test_for_success eq 'ARRAY');

    confess("Parameter error! If test_for_success is set, it must be an array of hashes!") 
      unless ($param_ok);
    
    my $stdout = `$cmd`;

    my $execution_failed = $? == -1;
    confess("Could not execute command:\n$cmd\n")
      if ($execution_failed);

    my $program_died = $? & 127;
    confess(
        sprintf (
            "Child died with signal %d, %s coredump\n",
            ($? & 127), ($? & 128) ? 'with' : 'without'
        )
    ) if ($program_died);

    my $exit_value = $? >> 8;
    my $program_completed_successfully = $exit_value == 0;
    confess("exited with value $exit_value")
        if (!$program_completed_successfully);

    if ($test_for_success) {

      foreach my $current_test_for_success (@$test_for_success) {

        confess('Type error') unless(ref $current_test_for_success eq 'HASH');

        use Hash::Util qw( lock_hash );
        lock_hash(%$current_test_for_success);

        my $current_test = $current_test_for_success->{test};
        confess('Test must be a sub!') unless (ref $current_test eq 'CODE');

        my $test_succeeded = $current_test->();

        confess(
            "The following command failed:\n"
            . "\n" . $cmd . "\n\n"
            . "Reason: " . $current_test_for_success->{fail_msg} . "\n"
        ) unless($test_succeeded);
      }
    }
    return $stdout;
}

1;

