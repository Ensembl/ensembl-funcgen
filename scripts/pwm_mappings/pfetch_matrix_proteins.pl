#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2013] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

 pfetch_jaspar_proteins.pl
  
=head1 SYNOPSIS

 pfetch_jaspar_proteins.pl [options]

=head1 OPTIONS

 DB Connection:
  --host            DB host
  --dbname          DB name (default = JASPAR_v5_0)
  --user            DB user name
  --pass            DB password (optional)
  --port            DB port (optional)
  --out_dir         Output directory (optional)
  --collections     List of Jaspar collections query, default is CORE and PBM
  --man
  --help

=head1 DESCRIPTION

B<This program>

=cut

use warnings;
use strict;

#use Bio::SeqIO;
#use Bio::DB::SwissProt;

#my $sp = new Bio::DB::SwissProt;

use Getopt::Long;
use Pod::Usage;
use DBI;

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_backtick_cmd 
                                               run_system_cmd
                                               open_file );
my ($pass, $user, $host, @cols);
my $dbname  = 'JASPAR_v5_0';
my $driver  = 'mysql';
my $port    = 3306;
my $out_dir = '.';

#Probably need to add more here?
my @tmp_args = @ARGV;

GetOptions 
  (
   'dbname=s'  => \$dbname,
   'pass=s'    => \$pass,
   'port=s'    => \$port,
   'host=s'    => \$host,
   'user=s'    => \$user,
   #"ssh"             => \$ssh,
   'out_dir=s' => \$out_dir,
   'collections=s{,}'  => \@cols,
   #'list_collections' => \$list_cols,
   
   #Other params
   'man'       => sub { pod2usage(-exitval => 0, -verbose => 2); },
   'help|?'    => sub { pod2usage(-exitval => 0, -verbose => 1, -message => "Params are:\t@tmp_args"); }
  )
  or pod2usage( -exitval => 1,
                -message => "Params are:\t@tmp_args"
              );

print "pfetch_jaspar_proteins.pl @tmp_args\n";


#? will catch unkown options?
if (@ARGV){
  pod2usage( -exitval =>1,
			 -message => "You have specified unknown options/args:\t@ARGV");
}

die('Mandatory --user param not specified') if ! $user;
die('Mandatory --host param not specified') if ! $host;
die("-out_dir does not exist:\t$out_dir")   if ! -d $out_dir; 

my $dbh;
my $dsn    = sprintf( "DBI:%s:%s:host=%s;port=%s",
                      $driver, $dbname, $host, $port);
eval {
  $dbh = DBI->connect( $dsn, $user, $pass, { 'RaiseError' => 1 } );
};

my $error = $@;

if ( !$dbh || $error || !$dbh->ping() ) {
  die("Could not connect to database $dbname as user $user using [$dsn] as a locator:\n".$error);
}


if(@cols){ #Validate vs DB

  die();
}
else{
  @cols = qw( CORE PBM );
  #CORE isn't strictly a collection
}




#What about species resrictions?

my $sql = 'SELECT m.BASE_ID as id, m.NAME name, mp.ACC uniprot_acc, '.
  'm.COLLECTION as collection, t.SPECIES as lspecies '.
  'FROM MATRIX m, MATRIX_PROTEIN mp, MATRIX_SPECIES ms, TAX t '.
  'WHERE m.ID=mp.ID and m.ID=ms.ID and ms.TAX_ID=t.TAX_ID '.
  'and m.COLLECTION in("'.join('", "', @cols).'")';

#warn $sql;
  
my $sth = $dbh->prepare($sql);
$sth->execute;

my ($last_id, @fails);
my $skipped = 0;
my $wrote   = 0;

while(my $hashref = $sth->fetchrow_hashref){
  
  if(! $last_id){
    $last_id = $hashref->{id};  
  }
  elsif($last_id eq $hashref->{id}){
    #The data suggests that this is actually a many to 1 relationship
    #so catch just in case
    warn('Found join product, likely due to >1 version or the many to many matrix - species '.
      'relationship actually being used '.$last_id."\n");
    
    #Actually these are more likely caused by a new version in the DB
    #These may have different MATRIX_ANNOTATION
    #but should have the same species
    #VERSION should really be in the MATRIX_ANNOTATION table. 
    $skipped ++;
    next; 
  }
  
  
  #pfetch may return a header with a .1 suffix, and maybe others
  #> pfetch Q6LBK7
  #>Q6LBK7.1 Q6LBK7_HUMAN C-myc protein (Fragment)
  #QIPELENNEKAPKVVILKKATAYILSVQAEEQKLISEEDLLRKRREQLKHKLEQLRNSCA
  
  
  #UP accs are cat'd in the field, sigh     
  my @up_accs = split(', ', $hashref->{uniprot_acc});
    
  foreach my $j_acc(@up_accs){
    #Assuming that these are all uniprot accessions 
    #pfetch will bring anything back from UniProtKB (inc current data) and EMBL
    my $pf_fasta = run_backtick_cmd('pfetch '.$hashref->{uniprot_acc});
  
    if(! defined $pf_fasta){
      push @fails, $hashref->{id};
    }
    else{
      my ($pf_header, @seq) = split("\n", $pf_fasta);      
      (my $pf_acc = $pf_header) =~ s/>([^ ])+ /$1/;
      #Keep pf accession for provenance
      
      $pf_fasta = '>'.join("\t", ($hashref->{id}, $j_acc, $pf_acc, $hashref->{name}, 
                                  $hashref->{collection}, $hashref->{lspecies}))."\n".join("\n", @seq);
     
      my $fh = open_file($out_dir.'/'.$hashref->{id}.'_'.$j_acc.'.fa', '>');
      print $fh $pf_fasta;
      $fh->close;
      $wrote++;
    }
  }
}
  
if(@fails){
  print STDERR 'Failed to retrieve '.scalar(@fails)." records for following protein accessions:\n".
    join("\n", @fails)."\n";  
} 
 
print "Wrote $wrote protein fasta files\n";
print "Skipped $skipped Jaspar records due to > 1 version being found\n"; 


#Check how I did that for Kaz work
#blastall -p blastp -d pepdb/Homo_sapiens.GRCh37.58.pep.all.fa -M BLOSUM80 -m 8 -b 1 -i example.fas
