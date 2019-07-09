#!/usr/bin/env perl

use strict;
use Data::Dumper;

=head1 find_species_for_probemapping

=head2 Usage

  find_species_for_probemapping.pl $(st1 details script)

=cut

use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::ApiVersion;

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

use Bio::EnsEMBL::Utils::CliHelper;
my $cli_helper = Bio::EnsEMBL::Utils::CliHelper->new();

my $optsd = $cli_helper->get_dba_opts();
my $opts  = $cli_helper->process_args($optsd, \&usage);

$opts->{dbname} ||= '';

$cli_helper->load_registry_for_opts( $opts );

my $species_on_staging = Bio::EnsEMBL::Registry->get_all_species;
my $species_hash = find_species_that_need_probemapping_update();

if (scalar keys %{$species_hash}) {
  print Dumper($species_hash);
  $logger->info("Done.\n");
}
else {
  my $api_version = Bio::EnsEMBL::ApiVersion->software_version();

  my $msg = 'No species found. This is unusual. Please make sure that you are '
      . 'using the latest \'release\' branch of the main ensembl repository. '
      . 'This script is using API version '
      . $api_version . "\n";
  $logger->warning($msg);
}

$logger->finish_log;

sub find_species_that_need_probemapping_update {

  my $species_hash = {};

  SPECIES:
  foreach my $species (sort @$species_on_staging) {
    
    my $funcgen_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
    
    if (! $funcgen_dba) {
      next SPECIES;
    }

    my $core_dba = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'core');
    my $meta_adaptor = $core_dba->get_MetaContainer;
    my $gene_build_version_from_core = $meta_adaptor->single_value_by_key('genebuild.last_geneset_update', 0);
    
    my $probemapping_adaptor = $funcgen_dba->get_ProbemappingAdaptor;
    my $probemappings = $probemapping_adaptor->fetch_all;
    my $gene_build_version = $probemappings->[0]->gene_build_version;
    
    if ($gene_build_version_from_core ne $gene_build_version) {
      
      $species_hash->{$species} = {
        gene_build_version_core    => $gene_build_version_from_core,
        gene_build_version_funcgen => $gene_build_version,
        species                    => $species,
        dbname                     => $funcgen_dba->dbc->dbname,
      };
    }
  }
  return $species_hash;
}

sub usage {
  print <<EOF; exit(0);

Usage: perl $0 <options>

  --host    Database host to connect to

  --port    Database port to connect to

  --user    Database username

  --pass    Password for user

  --release EnsEMBL release version

  --help    This message

 Example: perl $0 --host=mysql-ens-sta-1 --port=4519 --user=ensro

          or

          perl $0 \$(st1 details script)

EOF
}

