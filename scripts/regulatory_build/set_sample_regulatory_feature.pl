#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2025] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <ensembl-dev@ebi.ac.uk>.

  Questions may also be sent to the Ensembl help desk at
  <helpdesk@ensembl.org>.

=head1 NAME

  A script to set the sample regulatory feature. The sql is based on patch_87_88_c.sql.

  set_sample_regulatory_feature.pl \
      --registry /homes/mnuhn/work_dir_regbuild_testrun/lib/ensembl-funcgen/registry.with_previous_version.human_regbuild_testdb19.pm \
      --species homo_sapiens \

=cut

use strict;
use Getopt::Long;
use Bio::EnsEMBL::Registry;
use Data::Dumper;
use Bio::EnsEMBL::Utils::Logger;

my %options;
GetOptions (
    \%options,
    "species|s=s",
    "registry|r=s",
 );

my $species  = $options{'species'};
my $registry = $options{'registry'};

Bio::EnsEMBL::Registry->load_all($registry);

my $logger = Bio::EnsEMBL::Utils::Logger->new();
$logger->init_log;

my $funcgen_adaptor          = Bio::EnsEMBL::Registry->get_DBAdaptor($species, 'funcgen');
my $db_connection = $funcgen_adaptor->dbc;

my $sth;

$logger->info("Step 1\n");
$sth = $db_connection->prepare(qq(
  drop table if exists `temp_sample_id_max_activities`;
));
$sth->execute;

$logger->info("Step 2\n");
$sth = $db_connection->prepare(qq(
  drop table if exists `temp_sample_id`;
));
$sth->execute;

$logger->info("Step 3\n");
$sth = $db_connection->prepare(qq(
  create table temp_sample_id_max_activities as 
  select 
    regulatory_build_id, max(num_activities) as max_num_activities 
  from (
    select 
      regulatory_build_id, regulatory_feature_id, count(distinct activity) num_activities 
    from 
      regulatory_build 
      join regulatory_feature using (regulatory_build_id) 
      join regulatory_activity using (regulatory_feature_id) 
    group by 
      regulatory_build_id, regulatory_feature_id 
    order by 
      num_activities desc
  ) a group by regulatory_build_id;
));
$sth->execute;

$logger->info("Step 4\n");
$sth = $db_connection->prepare(qq(
  create table temp_sample_id as 
  select 
    temp_sample_id_max_activities.regulatory_build_id, min(regulatory_feature_id) as sample_regulatory_feature_id 
  from 
    temp_sample_id_max_activities join (
      select 
        regulatory_build_id, regulatory_feature_id, count(distinct activity) num_activities 
      from 
        regulatory_build 
        join regulatory_feature using (regulatory_build_id) 
        join regulatory_activity using (regulatory_feature_id) 
      group by 
        regulatory_build_id, regulatory_feature_id 
      order by 
        num_activities desc
  ) a on (
    a.regulatory_build_id=temp_sample_id_max_activities.regulatory_build_id 
    and temp_sample_id_max_activities.max_num_activities=a.num_activities
  )
  group by 
    temp_sample_id_max_activities.regulatory_build_id
  ;
));
$sth->execute;

$logger->info("Step 5\n");
$sth = $db_connection->prepare(qq(
  update 
    regulatory_build, temp_sample_id 
  set 
    regulatory_build.sample_regulatory_feature_id=temp_sample_id.sample_regulatory_feature_id 
  where 
    regulatory_build.regulatory_build_id=temp_sample_id.regulatory_build_id
  ;
));
$sth->execute;

$logger->info("Step 6\n");
$sth = $db_connection->prepare(qq(
  drop table `temp_sample_id_max_activities`;
));
$sth->execute;

$logger->info("Step 7\n");
$sth = $db_connection->prepare(qq(
  drop table `temp_sample_id`;
));
$sth->execute;

$logger->finish_log;
