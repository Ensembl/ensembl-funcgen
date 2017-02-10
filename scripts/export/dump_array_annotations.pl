#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2017] EMBL-European Bioinformatics Institute

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

=head1 NAME dump_array_annotations.pl

=head1 SYNOPSIS

dump_array_annotations.pl -host <string> -user <string> -pass <string> -dbname <string> \
                          (-arrays|class) <string> [<string> ...]  (-xrefs|features) \
                          [ OPTIONAL PARAMTERS ]

=head1 PARAMETERS

  Mandatory:
    -user     Funcgen DB user
    -pass     Funcgen DB password
    -port     Funcgen DB port
    -host     Funcgen DB host
    -dbname   Funcgen DB name
    
    -arrays   List of array names to dump
     OR
    -class    Class of arrays to dump
  
  
    -xrefs    Dump tab delimited file of probe/set transcript xref annotations
     OR
    -features Dumper bed file of genomic (including gapped cdna) probe alignments

  Optional:
    -dnadb_user Core DB user
    -dnadb_pass Core DB
    -dnadb_port Core DB
    -dnadb_host Core DB
    -dnadb_name Core DB
    -outdir     Output directory (default=./)
    -merged     Flag to merge the dump output into non-redundant rows wrt probe
                features/xrefs (probe name and/or array fields will be comma
                separated list of values)
    -prefix     Dump prefix to be used with -merged. This will be used in place 
                of the array names(s) in the output file name, and also the bed 
                track line name field.
    -list       Lists all the array names and classes in the given DB.
    -help       Prints some help
    -man        Prints the full man page

=head1 DESCRIPTION

Dumps either the probe features or the transcript annotations for the specified 
arrays. For speed, this is done using a direct SQL approach, as opposed to the 
dump_features.pl script which uses the API and is hence more robust but much slower.

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Funcgen::Utils::DBAdaptorHelper qw(get_DB_options_config
                                                     create_Funcgen_DBAdaptor_from_options
                                                     get_MYSQL_args_from_options);
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils  qw(run_system_cmd
                                               open_file);
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;

my @tmp_args = @ARGV;
my (@arrays, $class, $dump_features, $list,
    $outdir, $dump_xrefs, $merged, $prefix);
my $db_opts = get_DB_options_config(['funcgen', 'core']);

GetOptions (
            #Mandatory
            %{$db_opts},
            'arrays=s{,}' => \@arrays,
            'class=s'     => \$class,
            'outdir=s'    => \$outdir, 
            
            'features'   => \$dump_features,
            'xrefs'      => \$dump_xrefs,
                 
            #Optional
            'prefix=s' => \$prefix,
            'merged'        => \$merged,
            #'species=s'     => \$species,
            'list'          => \$list,
            
            'help'             => sub { pod2usage(-exitval => 0); }, #do we need verbose here?
            #removed ~? frm here as we don't want to exit with 0 for ?
            
            'man|m'            => sub { pod2usage(-exitval => 0, -verbose => 2); },
           ) 
           or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 

# VALIDATE PARAMETERS

my $db         = create_Funcgen_DBAdaptor_from_options($db_opts);
my $mysql_args = get_MYSQL_args_from_options($db_opts, 'funcgen');

if(! $outdir){
  $outdir = '.';
}
elsif(! -d $outdir){
    pod2usage(-exitval => 1, -message => "-outdir does not exist or is not a directory:\t$outdir");
}           
           
if(($dump_features && $dump_xrefs) || 
  ((! $dump_features) && (! $dump_xrefs)) ){
  pod2usage(-exitval => 1, -message => "Please specify -features or -xrefs to dump");
}          
     
if($prefix && ! $merged){
  warn 'Your have specified a -prefix but not -merged, -prefix will be ignored'; 
}     
           

my $array_adaptor = $db->get_ArrayAdaptor;
 
my (%class_arrays, $edb_name, $schema_build);

foreach my $array(@{$array_adaptor->fetch_all}){
  $class_arrays{$array->class} ||= [];  
  push @{$class_arrays{$array->class}}, $array;
} 

if($list){
  
  #do this here to ensure we group by class
  foreach my $c(keys %class_arrays){ 
    foreach my $a(@{$class_arrays{$c}}){
      print $c."\t".$a->name."\n";
    }
  }
  exit;  
}

if($dump_xrefs){
  $edb_name     = $db->species.'_core_Transcript'; 
}
else{
  $schema_build = $db->_get_schema_build($db);
}  
  


if(! @arrays){
  
  if(defined $class){
    @arrays = map { $_->name } @{$class_arrays{$class}};
    
    warn "No arrays defined, defaulting to all for $class:\n\t".
      join("\n\t", (@arrays)."\n");
  }
  else{
    pod2usage(-exitval => 1, 
              -message =>'You must provide either a list of -arrays or an array -class to dump');  
  }
}
else{
  
  if(! defined $class){
    #Just find the class of the first array
    my $aname = $arrays[0];
    
    foreach my $c(keys %class_arrays){
      if(grep { /^$aname$/ } (map { $_->name } @{$class_arrays{$c}}) ){
         $class = $c;
         last;
      }
    }
  
    if(! defined $class){
      pod2usage(-exitval => 1, 
                -message =>"$aname is not a valid array name. Please use -list to see all the available arrays");
    } 
  } 
   
   
  foreach my $aname(@arrays){
 
      
    if(! grep { /^$aname$/ } (map { $_->name } @{$class_arrays{$class}}) ){
      pod2usage(-exitval => 1, 
                -message => "$aname is not a valid array name. Valid $class array names:\n\t".
                            join("\n\t", (map { $_->name } @{$class_arrays{$class}}) ).
                            "\nPlease select array names from one only array class.\n"); 
    }
  }
}
  
my $redirect = '>';
#Could do a lot of this SQL generation just once
 
foreach my $array(@arrays){
  my $table_sql = 'array a, array_chip ac, probe p';
  my ($constraint_sql, $select_sql, $outfile);
  
  if($merged){
    $prefix ||= $class.'.'.join('_', @arrays);
    $constraint_sql = 'WHERE a.name in("'.join('", "', @arrays).'")';
  }
  else{
    $prefix = $class.'.'.$array;
    $constraint_sql = 'WHERE a.name ="'.$array.'"';
  }  
    
  if($class =~ /AFFY/){
    $table_sql      .= ', probe_set ps';
    $constraint_sql .= ' AND ps.probe_set_id=p.probe_set_id'; 
  }
  
  $constraint_sql .=  ' AND a.array_id=ac.array_id AND ac.array_chip_id=p.array_chip_id';
 
 
  if($dump_features){
    # DEFINE OUTPUT & PRINT TRACK LINE
    $outfile = $outdir.'/'.$prefix.'.bed';
    my $of = open_file($outfile, '>');
    print $of  "track name=$prefix description=\"Ensembl ".$class.':'.$prefix." mappings\" useScore=0\n";
    close($of);
    $redirect = '>>';
    
    # DEFINE SQL
    #Adding in schema_build clause here to avoid the product wrt nr seq_region entries   
    $constraint_sql .= " AND pf.seq_region_id=sr.seq_region_id AND sr.schema_build=\"${schema_build}\" AND p.probe_id=pf.probe_id GROUP by pf.probe_feature_id";
    $table_sql      .= ', probe_feature pf, seq_region sr';   
    my $name_sql     = 'p.name';
    
    if($class =~ /AFFY/){
      $name_sql = 'group_concat(a.name, ":", ps.name, ":", p.name)'; #nice 
    }
    
    
    $select_sql = 'SELECT sr.name, pf.seq_region_start, pf.seq_region_end, '.
    
    #$name_sql.', pf.mismatches, (CASE pf.seq_region_strand WHEN "1" THEN "+" WHEN "-1" THEN "-" ELSE ".")';   
    #CASE is only for stored_procedures!
    #triple replace for now until, replace this replace with something better?
    $name_sql.', pf.mismatches, replace(replace(replace(pf.seq_region_strand, "-1", "-"), "1", "+"), "0", ".")';   
  
    #currently calc score calc on the fly 
    #todo write function/stroed procedure to parse the cigarline and calculate the %ID
 
  }
  else{ #dump_xrefs
    # DEFINE OUTPUT 
    $outfile = $outdir.'/'.$prefix.'.xrefs.txt';
    
    # DEFINE SQL
    $table_sql  .= ', object_xref ox, xref x, external_db edb';
    my $id_field;    
    
    if($class =~ /AFFY/){
      $id_field      = 'ps.name,';
      $constraint_sql .= ' AND ps.probe_set_id=ox.ensembl_id AND ox.ensembl_object_type="ProbeSet"'; 
    }
    else{
      $id_field      = 'p.name,';
      $constraint_sql .= ' AND p.probe_id=ox.ensembl_id AND ox.ensembl_object_type="Probe"';   
    }
    
    $constraint_sql .= ' AND ox.xref_id=x.xref_id AND x.external_db_id = edb.external_db_id '.
      "AND edb.db_name=\"$edb_name\" GROUP BY $id_field x.xref_id";
      
    #Can't restrict by single db_release as we may have mixed releases if we have imported a new array
    #This should be fine as the array mapping pipline always removes xrefs before importing new ones
    
    $select_sql .= "SELECT $id_field group_concat(a.name), x.dbprimary_acc, x.display_label, ox.linkage_annotation ";
  }
 
  # BUILD THE SQL CMDLINE
  my $sql = "mysql --skip-column-names --quick -e'${select_sql} FROM ${table_sql} ${constraint_sql}' ".
    $mysql_args;
  
  #warn "$sql $redirect $outfile";
  
  # DO THE DUMP
  run_system_cmd("$sql $redirect $outfile");
  last if $merged;
}
  
   
 
