=head1 LICENSE

  Copyright (c) 1999-2013 The European Bioinformatics Institute and
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

Bio::EnsEMBL::Funcgen::Utils::DBAHelper

=head1 DESCRIPTION

This module provides a set of useful methods for handling Ensembl DBAdaptors.

=head1 SYNOPSIS

use Bio::EnsEMBL::Funcgen::Utils::DBAHelpers qw(list your required methods here);

=cut


# These methods were removed from EFGUtils do to a cyclical dependancy
# caused by the DBAdaptor using EFGUtils

###############################################################################

package Bio::EnsEMBL::Funcgen::Utils::DBAdaptorHelper;

require Exporter;
@ISA = qw(Exporter);
@EXPORT_OK = qw(               
                create_DBAdaptor_from_params
                mysql_args_from_DBAdaptor_params
                process_DB_options
                process_funcgen_DB_options
               );

#After global @ declarations to avoid warnings
use strict;
use warnings;

use Getopt::Long                   qw( GetOptionsFromArray );
use Bio::EnsEMBL::Utils::Exception qw( throw               );
use Bio::EnsEMBL::Utils::Scalar    qw( assert_ref          );
use Data::Dumper                   qw( Dumper              );
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor; #This was causing cyclical requirements as it requires EFGUtils
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;

#Keys here map to distinct DBAdaptors
my %db_type_options = ('funcgen' => ['user=s', 'pass=s', 'port=i', 'host=s', 'dbname=s'],   
                       'core'    => ['dnadb_user=s', 'dnadb_pass=s', 'dnadb_port=i', 'dnadb_host=s', 'dnadb_name=s'],
                       'hive'    => ['pdb_user=s', 'pdb_pass=s', 'pdb_port=i', 'pdb_host=s', 'pdb_name=s']);

#Add in some aliases
$db_type_options{dna}  = $db_type_options{core};
#These follwing two don't map to distinct DBAdaptors
#Hence will need manually will need manually re-assigning
#before use with create_DBAdaptors_from_params
$db_type_options{db}   = $db_type_options{funcgen};
$db_type_options{pipeline} = $db_type_options{hive}; 



=head2 process_DB_options

  Arg [1]    : Arrayref - DB types to handle e.g. funcgen|db, core|dna, hive|pipeline
  Arg [2]    : Hashref  - Options hash after GetOptions processing
  Arg [3]    : Boolean (optional) - Flag to set whether the DB type parameters are optional or not
  Example    : my $db_params = get_DBAdaptor_params_from_ARGV(\@ARGV, ['funcgen', 'core']);
  Description: Returns either an option hash for use with GetOptions, or processes the resulting 
               GetOptions processes hash to return DB constructor parameters
  Returntype : Hashref - Either DB options keys with param ref values (when no DB option hash is passed)
                         or
                         DB type keys, values are Hashrefs of DBAdaptor constructor params
  Exceptions : Throws if arguments are not valid.
               Throws if DB type is not valid.
               Throws if optional boolean is not specified and mandatory params
               for each given DB type are not found.
  Caller     : General, scripts.
  Status     : At risk

=cut
  


sub process_DB_options{
  my ($db_types, $db_opts, $optional, $maintain_options) = @_;  
  
  if(! defined $db_types){
    $db_types = ['funcgen', 'pipeline', 'core'];
  }
  elsif( (ref($db_types) ne 'ARRAY') ||
          (scalar(@$db_types) < 1) ){
    throw("DB types arg must be an arrayref of valid cmdline DB types:\t". 
      join("\t", keys(%db_type_options)) );            
  }
     
  my %return_hash = ();
  
  if($db_opts){
    
    assert_ref($db_opts, 'HASH', 'db_opts');
    my $param;
    
    foreach my $db_type(@$db_types){
      
      if(! exists $db_type_options{$db_type}){
        throw("$db_type is not valid. Valid DB types are:\t". 
          join("\t", keys(%db_type_options)) ); 
      }
      
      
      foreach my $opt(@{$db_type_options{$db_type}}){
        my $value = ${$db_opts->{$opt}};
             
        if(defined $value){
          #will already 'exist' from the first call of process_DB_options
          $param = $opt;
          
          if(! $maintain_options){
            ($param = $opt) =~ s/.*_//;
          }
          
          $param =~ s/\=.*//;    
          #handle dbname oddity
          $param = 'dbname' if $param eq 'name';
          $return_hash{$db_type}{'-'.$param} = $value;
        }
      }
     
      if(! $optional){
        
        if(! ($return_hash{$db_type}{'-user'} && 
              $return_hash{$db_type}{'-host'} && 
              $return_hash{$db_type}{'-dbname'}) ){
          throw("Could not find mandatory (user|host|dbname) $db_type DB param options:\n".Dumper($db_opts));
        }
      }
    }
  }
  else{
  
    foreach my $db_type(@$db_types){
    
      if(! exists $db_type_options{$db_type}){
        throw("$db_type is not valid. Valid DB types are:\t". 
          join("\t", keys(%db_type_options)) ); 
      }
    
      map {my $param; $return_hash{$_} = \$param } @{$db_type_options{$db_type}};  
    }
  }
  
  return \%return_hash;
}

#This simply allows dnadb_params to be optional

sub process_funcgen_DB_options {
  my $db_opts = $_[0];  
  
  return {%{&process_DB_options(['funcgen'], $db_opts)},
          %{&process_DB_options(['core'],    $db_opts, 1, 1)}};
}



=head2 create_DBAdaptors_from_params

  Arg [1]    : Hashref - DB contructor params
  Arg [2]    : String  - DB type
  Example    : my $db = create_DBAdaptors_from_params($funcgen_params, 'funcgen');
  DESCRIPTION: Creates and tests the relevant DBAdaptors from the passed parameters.
  Returntype : Bio::EnsEMBL::DBSQL::DBAdaptor (Or Funcgen or Hive DBAdaptor)
  Exceptions : Throws if arguments not valid 
               Throws if DB type is not supported
  Caller     : General, scripts
  Status     : At risk

=cut
    
sub create_DBAdaptor_from_params {
  my ($db_params, $db_type) = @_;
  
  assert_ref($db_params, 'HASH', 'db_params');
  my %dba_modules = (
                     funcgen => 'Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor',
                     core    => 'Bio::EnsEMBL::DBSQL::DBAdaptor',
                     hive    => 'Bio::EnsEMBL::Hive::DBSQL::DBAdaptor',
                    );    
  
             
  if( scalar(keys %$db_params) < 1 ){
    throw('You must pass a hashref DB parameters.');
  }  
  
  if(! (defined $db_type && 
        exists $dba_modules{$ db_type}) ){
    throw("DB type argument is invalid, valid DB types are:\t".join("\t", (keys %dba_modules)) );        
  }
  
  my $dba = $dba_modules{$db_type}->new(%$db_params);
    
  #Test connection
  $dba->dbc->db_handle;    
  return $dba; 
}    
   

sub mysql_args_from_DBAdaptor_params {
  my $db_params = $_[0];    
  assert_ref($db_params, 'HASH', 'db_params');
  
  my @db_mysql_params = (['-host',  '-h', 1],
                         ['-user',  '-u', 1],
                         ['-pass',  '-p'],
                         ['-port',  '-P'],
                         ['-dbname', ' ', 1]);
  
  
  my $mysql_params= '';
  
  foreach my $param_info(@db_mysql_params){
    my ($db_param, $mysql_param, $mandatory) = @$param_info;
    
    my $param = '';
    
    if(! exists ${$db_params}{$db_param}) {
      if($mandatory) {
        throw($db_param.' is a mandatory DB param');  
      }
    }
    else{
      $param .= ' '.$mysql_param.$db_params->{$db_param};  
    }
    
    $mysql_params .= $param;
  }  
  
  return $mysql_params;
}
   
sub parse_DB_url {
  my $url = $_[0];
  
  my ($user, $pass, $host, $dbname);

#  if($url =~ /(.*):\/\/(.*):(.*)@(.*):(.*)\/(.)/){
   
#   if($1 ne 'mysql://'){ #dbtype
#     pod2usage(-exitval => 1, -message => 'This script currently only supports MySQL urls');
#   }
  
#   if(! ($2 && $4 && $6)){
#     die('The -url  must contain at least user, host and name elements');
#   }

#   ($user, $host, $dbname) = ($2, $4, $6);
#   ($pass, $port) = ($3, $5);
#}
#else{
#  pod2usage(-exitval => 1, -message => "Unrecognized hive -url format:\t$url"); 
#} 

}

1;
