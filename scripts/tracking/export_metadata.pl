#!/usr/bin/env perl

# Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
# Copyright [2016-2018] EMBL-European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

=head1

./scripts/tracking/export_metadata.pl -c ./scripts/tracking/register.conf -o registry_script.csv

=cut

use strict;
use warnings;
use autodie;
use feature qw(say);

use Cwd 'abs_path';
use Config::Tiny;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::OntologyXref;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(uniq);

main();

sub main {
    my ( $csv, $config );

    # ----------------------------
    # read command line parameters
    # ----------------------------
    GetOptions(
        'o=s'  => \$csv,
        'c=s'  => \$config,
    );
    # -----------------
    # initialize logger
    # -----------------
    my $logger = Bio::EnsEMBL::Utils::Logger->new();
    $logger->init_log();

    # -------------------------------------
    # check that config and csv files exist
    # -------------------------------------
    if ( !-e $config ) {
        $logger->error(
            'Config file ' . abs_path($config) . ' doesn\'t exist!',
            0, 1 );
    }

    $logger->info( 'CSV file used: ' . abs_path($csv) . "\n",       0, 1 );
    $logger->info( 'Config file used: ' . abs_path($config) . "\n", 0, 1 );

    # --------------------
    # read the config file
    # --------------------
    my $cfg = Config::Tiny->read($config);

    # ------------------------------------------------------------
    # connect to funcgen tracking db, fetch all necessary adaptors
    # ------------------------------------------------------------
    $logger->info( 'Connecting to ' . $cfg->{efg_db}->{dbname} . '... ',
        0, 0 );
    my $dba = fetch_DBAdaptor($cfg);
    
    my $dbh = $dba->dbc->db_handle;
    my $sql = &get_big_sql_statement;
    my $column_headers = &get_column_headers;
    
    my $sth = $dbh->prepare($sql);
    $sth->execute;
    
    open my $csv_fh, '>', $csv;
    
    my $separator = "\t";
    
    my $header_line = join $separator, @$column_headers;
    
    $csv_fh->print($header_line);
    $csv_fh->print("\n");

    while (my $ary_ref = $sth->fetchrow_arrayref) {
      my $csv_line = join $separator, map {
        # Replace undefined values with empty strings.
        if (! defined $_) { '' } else { $_ }
      } @$ary_ref;
      $csv_fh->print($csv_line);
      $csv_fh->print("\n");
    }
    $logger->finish_log( "done\n", 0, 1 );

    return 1;
}

sub fetch_DBAdaptor {
    my ($cfg) = @_;
    my %adaptors;

    # Tracking DB hidden from user, hence no get_TrackingAdaptor method.
    # TrackingAdaptor->new() does not YET accept DBAdaptor object
    my $tracking_adaptor = Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor->new(
        -user    => $cfg->{efg_db}->{user},
        -pass    => $cfg->{efg_db}->{pass},
        -host    => $cfg->{efg_db}->{host},
        -port    => $cfg->{efg_db}->{port},
        -dbname  => $cfg->{efg_db}->{dbname},
        -species => $cfg->{general}->{species},
    );

    my $dba = $tracking_adaptor->db();
    return $dba;
}

sub get_column_headers {
  my @column_headers = qw(
    accession
    epigenome
    feature_type
    biological_replicate
    technical_replicate
    gender
    md5_checksum
    local_url
    analysis
    experimental_group
    ontology_xrefs
    xrefs
    epigenome_description
    control_id
    download_url
    info
  );
  return \@column_headers;
}

sub get_big_sql_statement {

  return <<SQL
select 
      CASE 
          WHEN experiment.is_control=1 THEN experiment.name
          ELSE input_subset.name
      END,
      epigenome.name, 
      feature_type.name, 
      input_subset.biological_replicate, 
      input_subset.technical_replicate,
      CASE 
          WHEN epigenome.gender is null THEN "unknown"
          ELSE epigenome.gender
      END,
      md5sum,
      local_url,
      analysis.logic_name,
      experimental_group.name,
      concat(linkage_annotation, '-', group_concat(xref.dbprimary_acc)),
      '-',
      '-',
      control.name,
      input_subset_tracking.download_url,
      input_subset_tracking.notes
    from 
      experiment 
      join input_subset using (experiment_id, epigenome_id) 
      join input_subset_tracking using (input_subset_id)
      join analysis on (input_subset.analysis_id=analysis.analysis_id)
      join epigenome using (epigenome_id) 
      join feature_type on (feature_type.feature_type_id=experiment.feature_type_id)
      join object_xref on (epigenome_id=object_xref.ensembl_id and object_xref.ensembl_object_type="epigenome")
      join xref using (xref_id)
      join external_db on (xref.external_db_id=external_db.external_db_id and external_db.db_name="EFO")
      join experimental_group using (experimental_group_id)
      left join experiment control on (experiment.control_id=control.experiment_id)
    group by 
      experiment.name, 
      epigenome.name, 
      feature_type.name, 
      input_subset.biological_replicate, 
      input_subset.technical_replicate,
      epigenome.gender,
      md5sum,
      local_url,
      analysis.logic_name
SQL
;
}

