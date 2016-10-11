#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016] EMBL-European Bioinformatics Institute

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
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=head1 NAME

add_hive_url_to_meta.pl - Sets/checks the hive url in the meta table

=head1 SYNOPSIS

add_hive_url_to_meta.pl -host <string> -user <string> -pass <string> -dbname <string> -url <string> -pipeline <String> [ -port <int> -help -man ]

=head1 PARAMETERS

  Mandatory:
    -host     <string>  DB host
    -user     <string>  DB user name
    -pass     <string>  DB password
    -dbname   <string>  DB name
    -url      <string>  The url of your hive DB e.g. 
                            mysql://DB_USER:DB_PASS@DB_HOST:DB_PORT/DB_NAME
    -pipeline <string>  Name of the pipeline (i.e. $ENV_NAME from your derivative pipeline.env)

  Optional:
    -port     <int>     DB port (default = 3306)  
    -species  <string>  Only required where species is not already set in meta or for multi species DBs
    -help
    -man

=head1 DESCRIPTION

This script checks whether a hive url meta tabe entry has been set for the 
given pipeline. If it is absent it sets it, if it does not if dies, if it is 
absent, it sets the value.

The script is used during initialisation of the pipeline environment. The aim 
is to ensure that it is not possible to run two versions of the same pipeline
on the same output DB.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils        qw( add_DB_url_to_meta  );
use Bio::EnsEMBL::Funcgen::Utils::DBAdaptorHelper qw( get_DB_options_config 
                                                      create_DBAdaptor_from_options );

my ($species, $url, $pname, $help, $man);
my $db_opts = get_DB_options_config(['funcgen']);

my @tmp_args = @ARGV;

GetOptions (
            %$db_opts,
            'url=s'      => \$url,
            'help'             => sub { pod2usage(-exitval => 0); }, #do we need verbose here?
            #removed ~? frm here as we don't want to exit with 0 for ?
            'man|m'            => sub { pod2usage(-exitval => 0, -verbose => 2); },
           ) or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 
           
#Drop this script in favour of a 1 liner in the env?

#We just want to create a core DBAdaptor here so we avoid defaulting to ensembldb for the dnadb
#and we don't want to force specification of the dnadb params
add_DB_url_to_meta('hive', $url, create_DBAdaptor_from_options($db_opts, 'funcgen', 'pass'));

1;
