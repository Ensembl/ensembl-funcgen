#!/usr/bin/env perl

=head1 LICENSE


  Copyright (c) 1999-2011 The European Bioinformatics Institute and
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

  dump_genes.pl

=head1 SYNOPSIS
 
  dump_genes.pl

=head1 DESCRIPTION

dump_genes.pl dumps all the genes in a database.
#It\'s a stripped down version of gene2flat.???

=cut

use strict;

use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::SeqIO;
use Getopt::Long;

my $dbname    = '';
my $dbhost    = '';
my $dbuser    = 'ensro';
my $dbport    = 3306;
my $dbpass;#    = undef;

my $dnadbname    = '';
my $dnadbhost    = '';
my $dnadbuser    = 'ensro';
my $dnadbport    = 3306;
my $dnadbpass    = undef;

my ($species, $multi_species, $dnadbmulti_species);
my $species_id=1;
my $dnadbspecies_id = 1;

my $cdna;
my $stable_id = 0;
my $db_id = 0;
my $file;
my $slicename;
my $verbose;
my $ev;
my $biotype;
my $help;
my @stable_id_list;
my $stable_id_file;

GetOptions(
		   'dbhost=s'    => \$dbhost,
		   'dbname=s'    => \$dbname,
		   'dbuser=s'    => \$dbuser,
		   'dbpass=s'    => \$dbpass,
		   'dbport=s'    => \$dbport,
		   'dnadbhost=s' => \$dnadbhost,
		   'dnadbname=s' => \$dnadbname,
		   'dnadbuser=s' => \$dnadbuser,
		   'dnadbport=s' => \$dnadbport,
		   'dnadbpass=s' => \$dnadbpass,
		   'dnadbspecies_id=i', => \$dnadbspecies_id,
		   'dnadbmulti_species!' => \$dnadbmulti_species,
		   'species=s'   => \$species,
		   'multi_species' => \$multi_species,
			 'species_id=i' => \$species_id,
		   'stable_id!'  => \$stable_id,
		   'db_id!'      => \$db_id,
		   'file=s'      => \$file,
       'slicename=s' => \$slicename,
       'verbose'     => \$verbose,
       'cdna!'       => \$cdna,
       'evidence!'   => \$ev,
		   'biotype=s'   => \$biotype,
		   'stable_id_file=s' => \$stable_id_file,
		   'help!'       => \$help,
		  ) or die ("Couldn't get options");

die ("'dbname=s'    => $dbname,
'dbhost=s'   => $dbname,
'dbuser=s'    => $dbuser,
'dbpass=s'    => $dbpass,
'dbport=s'    => $dbport,
'dnadbhost=s' => $dnadbhost,
'dnadbname=s' => $dnadbname,
'dnadbuser=s' => $dnadbuser,
'dnadbport=s' => $dnadbport,
'dnadbpass=s' => $dnadbpass,
'dnadbspecies_id=i' => $dnadbspecies_id,
'dnadbmulti_species!' => $dnadbmulti_species,
'species=s    => $species,
'species_id=i' => $species_id,
'multi_species!' => $multi_species,
'stable_id!'  => $stable_id,
'db_id!'      => $db_id,
'file=s'      => $file,
'slicename=s' => $slicename,
'verbose'     => $verbose,
'cdna!'       => $cdna,
'evidence!'   => $ev, use the transcript supporting feature description as the display_id
'biotype=s'   => $biotype,
")    unless( $dbhost &&  $dbname);

if ($stable_id and $db_id){
  $verbose and print STDERR "Entry ids will be db_id.stable_id\n";
} elsif ($stable_id) {
  $verbose and print STDERR "Entry ids will be translation stable_ids\n";
} else {
  $verbose and print STDERR "Entry ids will be translation dbIDs\n";
}

my $db;

if ($dnadbname) {
  if (not $dnadbhost or not $dnadbuser) {
    die "Fine. Your DNA is not in '$dbname' but in '$dnadbname'. But you must give a user and host for it\n";
  }
 
  my $dnadb = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                                 '-host'   => $dnadbhost,
                                                 '-user'   => $dnadbuser,
                                                 '-dbname' => $dnadbname,
                                                 '-pass'   => $dnadbpass,
                                                 '-port'   => $dnadbport,
												 '-species' => $species,
												 '-multispecies_db' => $dnadbmulti_species,
												 '-species_id' => $dnadbspecies_id
                                              );

  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
										   '-host'   => $dbhost,
										   '-user'   => $dbuser,
										   '-dbname' => $dbname,
										   '-pass'   => $dbpass,
										   '-port'   => $dbport,
										   '-dnadb' => $dnadb,
										   '-species' => $species,
										   '-multispecies_db' => $multi_species,
											 '-species_id' => $species_id
                                              );
} else {
  $db = new Bio::EnsEMBL::DBSQL::DBAdaptor(
                                              '-host'   => $dbhost,
                                              '-user'   => $dbuser,
                                              '-dbname' => $dbname,
                                              '-pass'   => $dbpass,
                                              '-port'   => $dbport,
										      '-group'  => 'core',
										   '-species' => $species,
										   '-multispecies_db' => $multi_species,
											 '-species_id' => $species_id
                                              );
}


#Test DBConnections here before we generate file.
$db->dbc->db_handle;
$db->dnadb->dbc->db_handle;



my $fh;
if($file){
  $verbose and print STDERR "Going to write peptides to '$file'\n";
  open (FH, '>'.$file) or die "couldn't open file ".$file." $!";
  $fh = \*FH;
} else{
  $verbose and print STDERR "Going to write peptides to stdout\n";
  $fh = \*STDOUT;
}

if($stable_id_file){
  open (FH, "$stable_id_file") or die "couldn't open file ".$stable_id_file." $!";
  while(<FH>){
    chomp;
    push @stable_id_list, $_;
  }
}


my $seqio = Bio::SeqIO->new('-format' => 'Fasta' , -fh => $fh ) ;
my $gene_adaptor = $db->get_GeneAdaptor();
my @genes;

if (defined $slicename) {
  my $slice = $db->get_SliceAdaptor->fetch_by_name($slicename);
  @genes = @{$slice->get_all_Genes};
}
if (scalar(@stable_id_list)){
  foreach my $id (@stable_id_list){
    my $gene = $gene_adaptor->fetch_by_stable_id($id);
    push  @genes, $gene if $gene
  }
}

unless (scalar(@genes) > 0 ){
  print STDERR "No genes found so far fetching all...\n";

  #my $gene_ids = $gene_adaptor->list_dbIDs();
  #@genes = @{$gene_adaptor->fetch_all_by_dbID_list($gene_ids)};

   @genes = @{$gene_adaptor->fetch_all};


}





foreach my $gene (@genes) {
  if ($biotype){
    next unless $gene->biotype eq $biotype;
  }
  my $gene_id = $gene->dbID();

  foreach my $trans ( @{$gene->get_all_Transcripts}) {
   my $identifier;
   my $tseq;
   
   if ($cdna){
   $identifier = $trans->stable_id . " ". $gene->description ;
   $identifier .= $gene->seq_region_name.":".$gene->start.":".$gene->end.":".$gene->strand;
     if ($ev){
       my $feat = $trans->get_all_supporting_features->[0];
       $identifier .= " ".$feat->hseqname." ";
       open (ID,"mfetch -f des ".$feat->hseqname." |") or die ("Cannot do the mfetch thing\n");
       while (<ID>){
	 chomp;
	 if ($_ =~ /DE\s+(.+)/){
	   $identifier .= " $1";
	 }
       }
       close ID;
     }
   $tseq = $trans->seq;
   } else {
     next if (!$trans->translation);
     if($db_id){
       $identifier = $trans->translation->dbID;
     }
     if($stable_id){
       if(!$db_id){
	 $identifier = $trans->stable_id ." " . $gene->stable_id;
       } else {
	 $identifier .= ".".$trans->stable_id;
       }
     }
     if ($ev){
       my $feat = $trans->get_all_supporting_features->[0];
       if ($feat){
	 open (ID,"mfetch -f des ".$feat->hseqname." |") or die ("Cannot do the mfetch thing\n");
	 $identifier .= " ".$feat->hseqname." ";
	 while (<ID>){
	   chomp;
	   if ($_ =~ /DE\s+(.+)/){
	     $identifier .= " $1";
	   }
	 }
	 close ID;
       }
     }
     $tseq = $trans->translate();
     if ( $tseq->seq =~ /\*/ ) {
       print STDERR "Translation of $identifier has stop codons ",
	 "- Skipping! (in ",$trans->slice->name(),")\n";
       next;
     }
   }
   
   $tseq->display_id($identifier);
#   $tseq->desc("Translation id $identifier gene $gene_id");
   $seqio->write_seq($tseq);
 }
}
close($fh);
