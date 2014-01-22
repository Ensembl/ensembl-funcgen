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
  --jdb_user    Jaspar DB user name
  --jdb_host    Jaspar DB host
  --jdb_name    Jaspar DB name (default = JASPAR_v5_0)
  --jdb_pass    Jaspar DB password (optional)
  --jdb_port    Jaspar DB port (optional)
  --out_dir     Output directory (optional)
  --user        DB user name (optional)
  --host        DB host (optional)
  --dbname      DB name (optional)
  --pass        DB password (optional)
  --port        DB port (optional)
  --collections List of Jaspar collections to query, default is CORE and PBM
  --man
  --help

=head1 DESCRIPTION

B<This program>

=cut



#ISSUES/CAVEATS

# 1 Theoretically have to validate xrefs after they are run by core
#   for the genes/translations we are interested in here, it is very 
#   unlikely they will change. Would need to check stable_id, uniprot accs
#   and optionally FeatureType name vs Gene external names
#   For release 75 xrefs only finished just before funcgen handoffs. 
#   Which doesn't give much time to rerun and patch in the new motifs. 
#   If we ever use the motifs for some other analysis, then this will cause 
#   further delays which will probably result in delaying the hand over

# 2 Full FeatureType > PFM > Gene display_id/name validation is not required
#   This may mean we may match a gene which does not have the feature type name as 
#   an eternal name. Need to monitor how many this affects



#TODO

# 1 Review counts

# 2 Make sure we update binding_matrix, as it may have change between Jaspar releases
#   do we store this version info anywhere?  Will actual matrix have changed if we enounter
#   different versions. Probably!

# 3 Attempt blast match for ftypes which have no species specific pfms
#   But have a ftype matched pfm on a different species.
#   This will require these we can identify a unique protein based on the ftype name

# 4 Add mode which will allow associations between fuzzy match pairs

# 5 Load new analysis based on jdb name. Load binding matrices. Leave old ones, and remove them afterwards based on the analysis_id.
#   Add version to description. Currently just set to 'Jaspar Matrix', which is redundant as this described by the analysis

# 6 DONE add support for Jaspar DB version? Add to analysis description (is already in analysis logic name).

# 7 Add support for restricting to a single FeatureType or PFM

# 8 Consider Ftype gene xrefs which are pre-existing. These used to be located in a flat file.
#   Actually, there are associations which can be made using Jsapar wich can't be done directly via the gene
#   Should we sub out the Gene assoication bit, so this can be reused?

# 9 Handle non species hits? Or just add species to the where clause and let the blast step handle that

use warnings;
use strict;

#use Bio::SeqIO;
#use Bio::DB::SwissProt;

#my $sp = new Bio::DB::SwissProt;

use Getopt::Long;
use Pod::Usage;
use DBI qw( :sql_types );

use Bio::EnsEMBL::Utils::Scalar                   qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument                 qw( rearrange );
use Bio::EnsEMBL::Analysis;
use Bio::EnsEMBL::DBEntry;

use Bio::EnsEMBL::Funcgen::Utils::Helper;
use Bio::EnsEMBL::Funcgen::Utils::DBAdaptorHelper qw( get_DB_options_config
                                                      create_Funcgen_DBAdaptor_from_options ); 
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils        qw( run_backtick_cmd 
                                                      run_system_cmd
                                                      open_file
                                                      add_external_db );

my @tmp_args = @ARGV;


### Some adventures in GetOptions processing ###

#my ($pass, $user, $host, @cols);
#my $dbname  = 'JASPAR_v5_0';
#my $driver  = 'mysql';
#my $port    = 3306;
#my $out_dir = '.';

#As we never use this vars now, do we even need to declare vars which dont have defaults?
#can we not simply specify an empty ref? [] {} or \undef

#Can we also add a preceding - to the options config such that we can use these will rearrange?
#or will this handle striping the absent -?

#my ($odir, @cols);

#Don't really need to keep these separate, but don eas an example
#my $db_opts    = get_DB_options_config(['funcgen', 'dna', 'jdb'], 1);
#my $other_opts = {'-out_dir=s' => \$odir,
#                  '-collections=s{,}'  => \@cols};#\[]};

#Dang this cause read only modification error!!!!
#We will never used these vars!
 
#GetOptions
# (%{$db_opts},
#  %{$other_opts},
#  man  => sub { pod2usage(-exitval => 0, -verbose => 2); },
#  help => sub { pod2usage(-exitval => 0, 
#                          -verbose => 1,
#                          -message => "Params are:\t@tmp_args"); }
# ) or #Catch unkown options
# pod2usage( -exitval => 1,
#             -message => "Params are:\t@tmp_args");

#Don't do this for $db_opts as process_DB_options does this
#We need to nest these options, such that we can pass the original options spec
#i.e. main rearrange doesn't have to handle them
#in fact it wouldn't anyway, it would just ignore them
#do we even need to keep these separate
#The only advantage is so we don't have to hadnle the individual vars in rearrange
#and we don't have to pass @_ to the create_Funcgen_DBAdaptor method

#foreach my $key(keys %$other_opts){
#  (my $new_key = $key) =~ s/\=.*$//o;
#  $other_opts->{$new_key}      = delete $other_opts->{$key};
#}

#$other_opts->{'-db_options'} = $db_opts;

#${,} does not work and s@ only pushes multiple -option specs onto same array
#Multiple names should only be specified in the @opts_config not the %opts                                   
my %opts = ('jdb_name' => 'JASPAR_v5_0',
            'jdb_version' => '5.0',
            'man'      => sub { pod2usage(-exitval => 0, -verbose => 2); },
            'help'     => sub { pod2usage(-exitval => 0, 
                                          -verbose => 1,
                                          -message => "Params are:\t@tmp_args"); });
my @db_opts  = keys %{get_DB_options_config(['funcgen', 'dna', 'jdb'], 1)}; #allow custom flag                                        
my @opts_config = (@db_opts, 
                   qw( out_dir=s  collection=s@ pfetch_pfms
                       create_edb help ));

GetOptions(\%opts, @opts_config) or #Catch unkown options
  pod2usage( -exitval => 1,
             -message => "Params are:\t@tmp_args");



#Delete the help and man subs as we don't want to pass those around
#Getopt::Long qw(:config auto_help) would do this for us but the help 
#is not verbose enough as it only prints the synopsis
#Not much point as it doesn't handle man, but do we need man?
#or could we rely on perldoc here?
delete $opts{help};
delete $opts{man};

#Rename collections array, so it makes sense if we ever require in this code
$opts{collections} = delete $opts{collection};

foreach my $key(keys %opts){
  $opts{"-${key}"} = delete $opts{$key};
}

#Gash! Neither way works the we need here so either I stripe the =.* spec off the keys
#or I have to add the - prefix GRRR!

#Probably better the other way, as it will be less fiddling and will allow us
#to update specific opts hashes to localize groups of options in predefined hashes 
#Also wont have to delete the man and help entries, and we will be able to use =s{,}
#No no no, we still have to deref the opts values!
#Better this way


if (@ARGV){ #Catch trailing args
  pod2usage( -exitval =>1,
             -message => "You have specified unknown options/args:\t@ARGV");
}


print "pfetch_matrix_proteins.pl @tmp_args\n";


#Sub everythign so we can re-use it in other scripts/modules if required

&main(\%opts);

sub main{
  my $opts = shift;    
  my ($out_dir, $colls, $create_edb, $jdb_version, $pfetch_pfms) = rearrange
    (['OUT_DIR', 'COLLECTIONS', 'JDB_VERSION','CREATE_EDB', 'PFETCH_PFMS'], %$opts );
  
  #assert_ref($colls, 'ARRAY');
  
  if(! $out_dir){
    $out_dir = '.';  
  }   
  elsif(! -d $out_dir){  
    die("-out_dir does not exist:\t$out_dir");
  }


  die('Mandatory --jdb_user param not specified') if ! defined $opts->{'-jdb_user'};
  die('Mandatory --jdb_host param not specified') if ! defined $opts->{'-jdb_host'};
  die('Mandatory --jdb_name param not specified') if ! defined $opts->{'-jdb_name'};
  

  #Currently get undefs here
  #add support in process_DB_options?

  my $jdbh;
  my $dsn    = sprintf( "DBI:%s:%s:host=%s;port=%s",
                        'mysql', $opts->{'-jdb_name'} , $opts->{'-jdb_host'}, 
                        $opts->{-jdb_port} || 3306);
                        
  eval {
    $jdbh = DBI->connect( $dsn, $opts->{'-jdb_user'}, $opts->{'-jdb_pass'}, { 'RaiseError' => 1 } );
  };

  my $error = $@;

  if ( ! $jdbh || $error || ! $jdbh->ping ) {
    die('Could not connect to database '.$opts->{'-jdb_name'}.' as user '.
      $opts->{'-jdb_user'}." using [$dsn] as a locator:\n".$error);
  }


  if($colls){ #Validate vs DB
    assert_ref($colls, 'ARRAY', 'Collections');

    die('collection validation not yet implemented');
  }
  else{
    $colls = [qw( CORE PBM )];
    #CORE isn't strictly a 'collection', although it is stored as such in the Jaspar DB
  }

  my $efg_db  = create_Funcgen_DBAdaptor_from_options($opts, 'optional', 1);#Conditional validate dnadb flag
  #Don't need to test this

  my $helper = Bio::EnsEMBL::Funcgen::Utils::Helper->new(no_log => 1);
  my $db_species = $efg_db->species;
  my $species = $db_species;
  
  if($species !~ /_/){ #assume we have failed to get a latin name
    die("$species latin species name, need to add a -species option");    
  }
  else{ #Convert to TAX format
    ($species = ucfirst(lc($species)) ) =~ s/_/ /;  
  }
  
  my $analysis   = $efg_db->get_AnalysisAdaptor->fetch_by_logic_name($opts->{'-jdb_name'});
  my $db_version = $efg_db->_get_schema_build($efg_db->dnadb);
  
  if(! defined $analysis){
    #Using Jaspar as the logic name here will prevent being able to load two versions along 
    #side each other which could complicate things
    
    $analysis = Bio::EnsEMBL::Analysis->new
     (-logic_name  => $opts->{'-jdb_name'},
      -db_version  => $jdb_version,
      -description => 'Position Frequency Matrices from the '.
        "<a href='http://jaspar.genereg.net/'>Jaspar Database</a> ($jdb_version)",
      -display_label => 'Jaspar'
        ); 
    $efg_db->get_AnalysisAdaptor->store($analysis); 
  }
  else{   #rollback!
    my $delete_sql = 'DELETE from object_xref where analysis_id='.$analysis->dbID;  
    #warn "Rolling back with:\t$delete_sql";
    $efg_db->dbc->db_handle->do($delete_sql);  
    
    #Do we need finer grain control of rollback here with different analyses?
    #We might just want to rerun the blast step or this step?
    #is blast step dependant on the output of this step
    #is yes then we need to do everythign at the same time
    
    #todo optionally rollback binding_matrices and motif_features?
  }
    
  if($create_edb){
    add_external_db($efg_db, $db_species.'_core_Gene', $db_version, 'EnsemblGene');
    #todo add translation here too?      
  }
  
  
  #Turn this into a two subs (PFM first, and Gene first)
  #Each will take the counts hashes generated by the other
  #Final output should be details of the matches for use with final
  #blast step. This will match all known ftype transcript/translations (based on external_name match)
  #to all the matrix proteins
  
  #So we really want main to generate the DBAdaptor
  #and then call the other subs
  #conditional on the input opts.

    
    
  #Let's make this optionally do everything i.e. do the name/xref match, load the matrix entry, store/clean/update the xref
  #then output a file of TF names that it can't resolve.

  #Only do those ftypes we have data for?    
  my $fset_a = $efg_db->get_FeatureSetAdaptor;
  my @fsets = grep { ($_->feature_type->class eq 'Transcription Factor') ||
                      ($_->feature_type->class eq 'Transcription Factor Complex') } 
                      @{$fset_a->fetch_all_by_feature_class('annotated')};
  #Make ftypes non-redundant    
  my %seen = ();
  my @ftypes = grep { ! $seen{$_->name}++ } (map { $_->feature_type } @fsets);
 
  my $dnadb      = $efg_db->dnadb; 
  my $gene_a     = $dnadb->get_GeneAdaptor;
  my $dbentry_a  = $efg_db->get_DBEntryAdaptor;
  my $bm_adaptor = $efg_db->get_BindingMatrixAdaptor;
  
  #These are the minimum fields(and aliases) required to support preprocess_jaspar_rows
  #So could add more, but can't remove any  
  #todo change this to use COLLATE for lc match
  my $sql = 'SELECT m.BASE_ID as id, m.VERSION as version, t.SPECIES as species, mp.ACC as accs, m.NAME as name '.
    'FROM MATRIX m, MATRIX_PROTEIN mp, MATRIX_SPECIES ms, TAX t '.
    'WHERE m.ID=mp.ID and m.ID=mp.ID and m.ID=ms.ID and ms.TAX_ID = t.TAX_ID and LOWER(m.NAME) like ?';
  my $sth  = $jdbh->prepare($sql); 

  #Can't do an IN statement as the accs are cat'd within the field!!!
  my $sql2 = 'SELECT m.BASE_ID as id, m.VERSION as version, t.SPECIES as species, mp.ACC as accs, m.NAME as name '.
    'FROM MATRIX m, MATRIX_PROTEIN mp, MATRIX_SPECIES ms, TAX t '.
    'WHERE m.ID=mp.ID and m.ID=mp.ID and m.ID=ms.ID and ms.TAX_ID = t.TAX_ID and mp.ACC like concat("%", ? "%")';
  my $sth2 = $jdbh->prepare($sql2);

  #Old method was to get the gene based on the ftype name, then blast that against 
  #the pfm protein (or vice versa) to get the best blast hit.
  #We will do this for those TFs which we do not have any hits for
  #but we need to set a threshold?
  
  #Currently doing this the other way around
  #using TF names to bring back matrices directly
  #then pulling back proteins based on the uniprot accs
  #need to validate they have the ft as an external name
  #but can probably take them regardless if there isn't a better hit
  #Will any of these be identity xrefs? No, where do these annotations come from?
  
  #external_dbs
  #ftypes tend to be Uniprot_gn
  #Jaspar uniprot accessions tend to be Uniprot/SWISSPROT
  #and can be mapped to many translations
  #Jaspar names tend to be Uniprot_gn
    
  my (%ftype_matches, %no_matches, %counts);
  #todo put these counters in the hash? or remove hash in favour of individual counts
  my $ftype_to_pfm_complex    = 0;
  my $gene_to_pfm_complex     = 0;
  my $ftype_to_gene_only_full = 0;  
  my $ftype_to_gene_full      = 0;
  my $ftype_to_gene_no_full   = 0;
  my $ftype_to_gene           = 0;
  my $ftype_to_pwm            = 0;
  my $ftype_pep_fasta = $out_dir.'/'.$db_species.'_TranscriptionFactor_FeatureTypes.pep.fasta';
  my $fasta_fh = open_file($ftype_pep_fasta, '>');
  
    foreach my $ftype(@ftypes){
      my $ftype_name = $ftype->name;
      
      ####  FeatureType name > PFM name > Gene accession/external_name ###
      $sth->bind_param(1, '%'.lc($ftype_name).'%', SQL_VARCHAR);
      $sth->execute;
      #Filters out older versions and non species specific rows
      #Also handles fuzzy matching
      my ($pfm_matches, $non_species_rows) = &preprocess_jaspar_rows($sth, $species, $ftype_name);
   
      if(! scalar(@$pfm_matches)){
        $no_matches{$ftype_name}{ftype_to_pfm} = $ftype;    
      }
      else{
        $ftype_matches{$ftype_name}{feature_type} = $ftype;
        
        #todo Do we still want to load pfm matches that don't validate
        #even id we get others which do?
        #consider in blast first? 
              
            
        foreach my $row_ref(@$pfm_matches){
          my $found_gene = 0;
          my ($id, $species, $accs, $name) = 
            ($row_ref->{id}, $row_ref->{species}, $row_ref->{accs}, $row_ref->{name});
       
          if($name =~ /::/){
            $ftype_to_pfm_complex++;  
          }

          foreach my $i(0..$#{$accs}){  
            my $matched_display_label;
            my $seen_multi;
            my @genes = @{$gene_a->fetch_all_by_external_name($accs->[$i])};
      
            if(scalar (@genes) >1){
              $counts{multi_acc_to_gene_hits}++;
            }
            elsif(! scalar(@genes) 
              $ftype_matches{$ftype_name}{ftype_pfm_name} ||= [$id, []]; 
              push @{$ftype_matches{$ftype_name}{ftype_pfm_name}}, $accs->[$i]; #No gene matches
              #This is very unlikely, as we are filtering out the non-species specific stuff                                             
            }  
          
            if($name =~/::/){
              print "Found ".scalar(@genes)." genes for TF complex hit $name with acc ".$accs->[$i]."\n";  
            }
          
            foreach my $gene(@genes){
              $found_gene = 1;
              my @xrefs = grep  {lc($_->display_id) eq lc($ftype_name)} @{$gene->get_all_DBLinks};
              
              #Match the ftype name to the gene external names
              if(@xrefs){     
                $counts{ftype_pwm_gene_and_name}++;    #This is a redundant count   
                print "PFM (name) FeatureType: $ftype_name > PFM $id $name > ".$gene->display_id.' '.$accs->[$i]."\n"; 
                
                if($matched_display_label){
                  
                  if(! $seen_multi){
                    $counts{multi_gene_to_display_label_hits}++; 
                    #This is likely due to paralogs/gene families
                    $seen_multi = 1;
                  }
                }
                
                $matched_display_label = 1;
                $ftype_matches{$ftype_name}{pwm_gene_name} ||= [];
                push @{$ftype_matches{$ftype_name}{pwm_gene_name}}, [$gene, $id, $name];
                
                #This is quicker than using the API
                #my $display_name = $helper->get_core_display_name_by_stable_id($dnadb, $gene->stable_id, 'gene');
                #my $dbentry = Bio::EnsEMBL::DBEntry->new
                # (-dbname                 => $db_species.'_core_Gene',
                #  -release                => $db_version,
                #  -status                 => 'KNOWNXREF',
                #  #-display_label_linkable => 1,
                #  #-db_display_name        => $self->db->dnadb->dbc->dbname,
                #  -db_display_name        => 'EnsemblGene',
                #  -type                   => 'MISC',
                #  -primary_id             => $gene->stable_id,
                #  -display_id             => $display_name,
                #  -info_type              => 'MISC',
                #  -info_text              => 'GENE',
                #  -linkage_annotation     => 'FeatureType > Jaspar PFM > Gene (name match)',
                #  -analysis               => $analysis,
                # );
                #DBI->trace(2);
                #$dbentry_a->store($dbentry, $ftype->dbID, 'FeatureType');#1 is ignore release flag 
                #DBI->trace(0);     
                
                    #Seems like non of these are on the Translation level
              #We could load the Transcript xrefs instead, which are almost as specific
              #We shoudl really identify which of the objects have the acc hit
              #and root the rest of the search off that
              #too much
          
             #foreach my $xref(@$xrefs){
              #  #next if $xref->ensembl_object_type eq 'Gene';
              #  #Here is appears we can have various Gene and Transcript xrefs
              #  #No easy way to resolve the relationships between them and Translations
              #  #to make them nr, and also ignore transcripts if we have a more specific translation
              #  
              #  if($xref->ensembl_object_type eq 'Translation'){
              #    print 'Found Tranlsation xref!!! '.$xref->ensembl_id;
              #  }
              #}
                    
              }
              else{
                print "PFM (acc): FeatureType $ftype_name > PFM $id $name > ".$gene->display_id.' '.$accs->[$i]."\n"; 
                $ftype_matches{$ftype_name}{pwm_gene_accession} ||= [];
                push @{$ftype_matches{$ftype_name}{pwm_gene_accession}}, [$gene, $id, $name];
              }
            }
          }
           
           
          if(! $found_gene){
            #cchange this to cache in no_matches?
            #so we know which ones we need to consider for the blast?
            $no_matches{$ftype_name}{pfm_to_gene} = $ftype;          
            $ftype_to_pfm++;    
          }
        }
      
       
      #LOAD THE MATRIX FIRST OR LAST We may not have any hits
      
      #if(@display_label_matches ||
      #   @gene_acc_matches){
      
      #  my $matrix =        
      #     
      #     
      #}
                      
      if(exists $ftype_matches{$ftype_name}{pwm_gene_name} ){ 
        
        if(! exists $ftype_matches{$ftype_name}{pwm_gene_accesssion}){
          $ftype_to_gene_only_full ++;  
        }  
          
        $ftype_to_gene_full ++;  
      }
          
      if(exists $ftype_matches{$ftype_name}{pwm_gene_accesssion}){
              
        if(! exists $ftype_matches{$ftype_name}{pwm_gene_name}){
          $ftype_to_gene_no_full ++;  
        }  
          
        $ftype_to_gene++;  
        
        #foreach my $gene(@gene_acc_matches){
        #  #This is quicker than using the API
        #  my $display_name = $helper->get_core_display_name_by_stable_id($dnadb, $gene->stable_id, 'gene');
  
        #  my $dbentry = Bio::EnsEMBL::DBEntry->new
        #   (-dbname                 => $db_species.'_core_Gene',
        #    -release                => $db_version,
        #    -status                 => 'KNOWNXREF',
        #    #-display_label_linkable => 1,
        #    #-db_display_name        => $self->db->dnadb->dbc->dbname,
        #    -db_display_name        => 'EnsemblGene',
        #    -type                   => 'MISC',
        #    -primary_id             => $gene->stable_id,
        #    -display_id             => $display_name,
        #    -info_type              => 'MISC',
        #    -info_text              => 'GENE',
        #    -linkage_annotation     => 'FeatureType > Jaspar PFM > Gene (accession match)',
        #    -analysis               => $analysis,
        #   );
    
          #DBI->trace(2);
        #  $dbentry_a->store($dbentry, $ftype->dbID, 'FeatureType');#1 is ignore release flag 
        #}
          
          #This is potentially a multi acc > gene hit
          #Based only on the uniprot acc i.e. no display_id - ftype hits
          
        #Load and try an ID in reverse too?
        #This may identify a pfm which does not match the ftype name
        #we should list all of these such that we can manually eye ball them
        #most likely just slight variations on name or complex
        
        #or should we just do evertyhing in reverse too!
        #so long as we make the ftype to pfm association
        #we don't care how many xrefs there are
        
        #We need to highlight  
        #links which are only made via the ftype name to the pfm
        #with no further validation?           
      }
   

      ####  FeatureType name > Gene external_name > PFM accession/name ###
      
      #we want to avoid duplicating the fully validated hits
      #$counts{any_matches}{$ftype_name} = 'full';
  
      #Use different hashes for reverse match so we now what came from where
  
  
  
 
      my @genes = @{$gene_a->fetch_all_by_external_name($ftype_name)};
      
    
      if(scalar(@genes > 1)){
        $counts{'multi_ftype_to_gene_hits'} ++;  
      }
      elsif(scalar(@genes)){
        $ftype_gene_matches{$ftype_name}{feature_type} = $ftype;
      }
      else{
        $no_matches{$ftype_name}{ftype_to_gene} = $ftype;
      }
     
      
      foreach my $gene(@genes){
        my $seen_pwm = 0;
        my %seen     = ();
        my @dbprimary_accs =  map {$_->primary_id } 
                                (grep { ! $seen{$_->primary_id}++ } @{$gene->get_all_DBLinks});
        #dbprimary accs maybe redundant as they can be associated to 
        #gene|transcript|translation from the same external_db
        
        
        foreach my $xref_acc(@dbprimary_accs){
          $sth2->bind_param(1, $xref_acc, SQL_VARCHAR);
          $sth2->execute;
          my ($rows, $non_species_rows) = &preprocess_jaspar_rows($sth2, $species, undef, $xref_acc);
        
          foreach my $row_ref(@$rows){
            my ($id, $species, $name) = 
              ($row_ref->{id}, $row_ref->{species}, $row_ref->{name});
            $seen_pwm = 1;
           
            #Need to account for fuzzy TF complex matching here
            #as these will have already been caught about
              
              
                
            if($name =~ /(^|::)${ftype_name}(::|$)/i){
              #should have already seen these matches
              #Count, so we can snity check
              $counts{ftype_gene_pwm_and_name}++;  
            }
            else{ #These are the new ones we're interested in!
              #$counts{ftype_gene_pwm_no_name}++; 
              
              if($name =~ /::/){
                $gene_to_pfm_complex++;  
              }  
                           
              $ftype_gene_matches{$ftype_name}{pwm_accession} ||= [];
              
              push @{$ftype_gene_matches{$ftype_name}{pwm_accession}}, [$gene, $id, $name]; 
              #Need to remove these from the no match cache
              print "GENE (acc): FeatureType $ftype_name > ".$gene->display_id." ($xref_acc) > PFM $id ${name}\n";
              
            }
          }
          
          
          #Now dump the pep fasta for the transcript which matches the ftype name
          my @xrefs = grep {lc($_->display_id) eq lc($ftype_name)} @{$gene->get_all_DBLinks};
          
          my $seen_xrefs;
          
          foreach my $xref(@xrefs){
            
            if($xref->ensembl_object_type eq 'Transcript'){
              
              if(! exists $seen_xrefs{'Transcript:'.$xref->ensembl_id}){
               my $pep = $gene->get_Transcript->get_Tranlsation;
          
               print $fasta_fh '>'.$pep->stable_id.'('.$ftype_name.")\n".
                $pep->seq."\n";
               $seen_xrefs{'Transcript:'.$xref->ensembl_id} = 1; 
            }
            
            if($xref->ensembl_object_type eq 'Translation'){
              die("Found unexpected Translation xref for FeatureType name $ftype_name");
            }
          }     
        }
        
        #Now cache those with no acc/name hits
        if(! $seen_pwm){
          #or just count here?
          $no_matches{$ftype_name}{gene_to_pwm} = $ftype;     
        }
      }
      
      #Now we handle no hits here
      #If we truly have no hits for this ftype
      #then we need to add the genes returned from the 2nd bit
      #to a cache for dumping
      #we also need to add the ftype name to the headers
      #there can be >1 due to close paralogs
      
      #Then we can blast these against the whole of the jaspar proteins
      #taking the top hit above a threshold
      #optionally across species!
      
      #todo Count complex hits     
     
      #Let's do the blast step separately?
      #Which can optionally take a list of gene stable_ids?
      #Actually, these stable_ids need to come from the same level of as the 
      #DBEntry which matched the ftype name
      #this is likely Uniport_gn and therefore gene level
      
      #Let's post pone this for now
      #and do the all vs all
      
      #We still need to get the gene names in the fasta
      #so we can match these to the feature_types?
      #So should inputs be a list of ftype names
      #or a list of stable_ids
      #probably take both?
      #and sub out code which IDs genes/DBEntries based on ftype name
      
      
       
    }
  
    #Are we interested in the Gene only matches for xrefs?

  
  
    #write info to file for TFs with no hits or poor hits
    #This can then be used by the pfetch script to do the relevant blast
 
    #these are nr!
    #my $ftype_to_gene_only_full = 0;  
    #my $ftype_to_gene_full      = 0;
    #my $ftype_to_gene_no_full   = 0;
    #my $ftype_to_gene           = 0;
    
  
    #break this down into full and non-display_id matches
    print "FeatureType > PFM > Gene matches:\n";
    print "Total redundant fully validated matches, FeatureType > PFM > Gene name:\t".$counts{ftype_pwm_gene_and_name}."\n"; 
    print "Total non-redundant wrt FeatureType matches (all/only name match/only no name match):\t\t".
      $ftype_to_gene_full.'/'.$ftype_to_gene_only_full.'/'.$ftype_to_gene_no_full."\n";
    #.scalar(keys %{$counts{any_matches}})."\n";
    
    
    print "Found no FeatureType - PFM matches:\t\t\t\t".scalar(keys %no_pfm_matches)."\n";
    print "Found only Feature - PFM matches (i.e. no PFM-Gene acc match):\t".scalar(keys %no_gene_matches)."\n";
    
    print 'Both of the following counts ere probably the result of different genes giving'.
      " rise to the same protein product through identical splicing:\n"; 
    print "Multiple Gene hits for a single acc:\t\t\t\t".$counts{multi_acc_to_gene_hits}."\n";
    print "Including a matching display label:\t\t\t\t".$counts{multi_gene_to_display_label_hits}."\n";
    print "Total FeatureType to PFM TF complexes matches:\t\t\t$gene_to_pfm_complex\n"; 
      
    
    print "\nFeatureType > Gene > PFM matches:\n";
    print "Ignoring fully validated matches(nr), FeatureType > Gene > PFM name:\t".$counts{ftype_gene_pwm_and_name}."\n";
  
    #push @{$gene_to_pfms{$ftype_name}{$gene->display_id}}, "$id $name"; 
  
    #will there be redundant gene hit in here between ftypes?
  
 
    print "Distinct FeatureTypes with PFM matches (not name validated):\t\t".scalar(keys %gene_to_pfms)."\n";
    
    my $total_gene_matches = 0;
    my $total_pfm_matches  = 0;
    
    foreach my $href(values %gene_to_pfms){
      
      foreach my $gid(keys %$href){
        $total_gene_matches++;
        $total_pfm_matches += scalar(@{$href->{$gid}});
      }
    }   
      
    
    print "Total Gene matches with PFM matches (not name validated):\t\t".$total_gene_matches."\n";
    print "Total PWM matches (not name validated):\t\t\t\t\t".$total_pfm_matches."\n";
    print "Total Gene to PFM TF complexes matches:\t\t\t\t\t$gene_to_pfm_complex\n"; 
       
    
  
    #print "All PFM matches (not distinct wrt ftype or gene, not full name match):\t".$counts{ftype_gene_pwm_no_name}."\n";
    
     
    
    #Not getting the converage of the previous approach...yet!
    
    #Before reverse search, multi acc and TF complex support
    
    #Found FeatureType - Gene matches (inc non-display_id):  29
    #Found no FeatureType - PFM matches:     61
    #Found only Feature - PFM matches (i.e. no PFM-Gene acc match):  1
    #Both of the following counts ere probably the result of different genes giving rise to the same protein product through identical splicing:
    #Multiple Gene hits for a single acc:    1
    #Including a matching display label:     1
      
    #Now with gene first search? and fuzzy matching and cat'd acc handling
    # but no complex::processing 
    
    #FeatureType > PFM > Gene matches:
    #Found FeatureType - Gene matches (inc non-display_id):  36
    #Found no FeatureType - PFM matches:     54
    #Found only Feature - PFM matches (i.e. no PFM-Gene acc match):  1 
    #Both of the following counts ere probably the result of different genes giving rise to the same protein product through identical splicing:
    #Multiple Gene hits for a single acc:    5
    #Including a matching display label:     5
    
    #FeatureType > Gene > PFM matches:
    #Full (should already have these in the above numbers):  85
    #All PFM matches (not distinct wrt ftype or gene):       23
    
      
 
}
 
 
sub preprocess_jaspar_rows {
  my ($executed_sth, $species, $name, $acc) = @_;
  
 
  if($name && $acc){
    die('Name and Accession arguments are mutually exclusive');  
  }
 
  #todo validate params
  my %rows;
  my @non_species_rows;
  
  #todo probably want to process non species rows properly too? 
  #'SELECT m.BASE_ID as id, t.SPECIES as species, mp.ACC as acc, m.NAME as name '.
  
  #Count complex '::' hits in caller not here!
    
  while(my $href = $executed_sth->fetchrow_hashref){
    
    if($href->{species} ne $species){
      push @non_species_rows, $href;
    }
    else{
      
      if($name){ #Allow TF complex matches
        #next if $href->{name} !~ /(^|::)${name}(::|$)/; #erroneous match
        if ($href->{name} !~ /(^|::)${name}(::|$)/i){
          warn "Skipping fuzzy name match:\t\t$name (Ensembl) vs ".$href->{name}." (Jaspar)\n";
        }
        elsif($href->{name} =~ /::/){
          print "Matched complex:\t$name (Ensembl) vs ".$href->{name}." (Jaspar)\n";
            
        }
      }
      
      if( (! exists $rows{$href->{id}}) ||
            $rows{$href->{id}}{version} < $href->{version}) {
              
        #Sanity check it isn't the same version?!
        #Anything could happen
      
        if((exists $rows{$href->{id}}) &&
            $rows{$href->{id}}{version} == $href->{version}) {
          die('Found two '.$href->{id}." with the same version:\t".$href->{version});        
        }        
     
        my @accs = split(/,\s+/, $href->{accs});
      
        if($acc){
        
          if(grep /^${acc}$/, @accs){
            $href->{accs} = [$acc];
          }
          else{  
            warn "Erroneous fuzzy accession match:\t$acc (Ensembl) vs @accs (Jaspar)\n";
            next; } #This was an erroneous match
        }
        else{
          $href->{accs} = \@accs;    
        }
     
        $rows{$href->{id}} = $href;
      }
    }         
  }
  
  my @rows = values(%rows);
  return (\@rows, \@non_species_rows);
}
 
 
sub pfetch_jaspar_pfm_fastas{
  my $jdbh    = shift;
  my $out_dir = shift; 
  my $colls   = shift;
  #my $pfms   = shift;
 
  assert_ref($jdbh, 'DBI', 'Jaspar DB Handle');
 
  if(! $out_dir){
    $out_dir = '.';  
  }   
  elsif(! -d $out_dir){  
    die("Output directory does not exist:\t$out_dir");
  }
 
  #Now need to change this based on what we have already done
  #probably need to sub some of it out


  #We probably need to handle a funcgen DB here too
  #use DBAdaptor::create_Funcgen_DB_from_options
  #This will be an issue if we only want the core DB
  #but we defo want the funcgenDB too.
  #do we need to make core DB mandatory?
  #default core DB will be unsafe
  #shoudl probably just remove that auto load stuff
  #or at least make is an opt in mode?
  
  #What about species restrictions?
  
  my $sql = 'SELECT m.BASE_ID as id, m.NAME name, mp.ACC uniprot_acc, '.
    'm.COLLECTION as collection, t.SPECIES as lspecies '.
    'FROM MATRIX m, MATRIX_PROTEIN mp, MATRIX_SPECIES ms, TAX t '.
    'WHERE m.ID=mp.ID and m.ID=ms.ID and ms.TAX_ID=t.TAX_ID '.
    'and m.COLLECTION in("'.join('", "', @$colls).'")';
  
  #warn $sql;
    
  my $sth = $jdbh->prepare($sql);
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
  
}

1;