#!/usr/bin/env perl

=head1 LICENSE

Copyright [1999-2014] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute

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

load_jaspar_matrices.pl
  
=head1 SYNOPSIS

load_jaspar_matrices.pl [options]

=head1 OPTIONS

 DB Connection:
  --user               DB user name 
  --host               DB host 
  --dbname             DB name 
  --pass               DB password 
  --port               DB port (optional)  --jdb_user           Jaspar DB user name
  --dnadb_user         Core DB user
  --dnadb_host         Core DB host
  --dnadb_name         Core DB name
  --dnadb_pass         Core DB password (optional)
  --dnadb_port         Core DB port (optional)
  --jdb_user           Jaspar DB user
  --jdb_host           Jaspar DB host
  --jdb_name           Jaspar DB name (default = JASPAR_v5_0). This will be used as the analysis logic_name
  --jdb_pass           Jaspar DB password (optional)
  --jdb_port           Jaspar DB port (optional)
  
 Inputs/Outputs: 
  --out_dir            Output directory (optional)
  --pep_fasta          Location of the relevant species whole proteome fasta DB (i.e. use formatdb)
  --pfm_file           Location of the file containing all the Position Frequency Matrices
 
 Other:
  --collections        List of Jaspar collections to query, default is CORE and PBM
  --dump_jaspar_fasta  Dumps fasta files for the Jaspar accessions
  --skip_blast  	     Skips blast step if blast results file already exists (mainly for testing)
  --man
  --help

=head1 DESCRIPTION

B<This program>


For a Jaspar analysis which has been previosuly imported, binding matrix entries and associated xrefs are 
deleted before being re-imported.

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

# 10 Use helper to log/debug, and add helper options

# 11 Update to use Uniprot rest web service rather than flaky (and proprietary) pfetch
#    dbfetch service which has some latency, so probably better to update to use that directly using LWP

# 12 Consider ftype_to_pfm name only matches if we can't even get a blast result. This may be because 
#    the feature type name is not present as a gene external name 

# 13 Does this handle rollback of previously loaded xrefs/associations? Yes it does.

# 14 This script is dependant on core xrefs being in place. As this is the first
#    step in the motif pipeline, it effectively holds up the alignments.
#    Supporting a post-import tidyup, is probably more trouble than it's worth.
#    It is probably only needed with a new assembly, or a significantly revised genebuild
#    A validate mode maybe useful, which would do the associations, and check whether they
#    match what is in the DB already. So decoupling the association from the DB import would be nice



use warnings;
use strict;

#use Bio::SeqIO;
#use Bio::DB::SwissProt; #This uses dbfetch not uniprot rest service (which is always up)

#my $sp = new Bio::DB::SwissProt;

use Data::Dumper qw( Dumper );

use Getopt::Long;
use Pod::Usage;
use DBI qw( :sql_types );

use Bio::EnsEMBL::Utils::Scalar   qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument qw( rearrange );
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
my %opts = (
  'jdb_name'    => 'JASPAR_v5_0',
  'jdb_version' => '5.0',
  'man'         => sub { pod2usage( -exitval => 0, -verbose => 2 ); },
  'help'        => sub {
    pod2usage( -exitval => 0,
               -verbose => 1,
               -message => "Params are:\t@tmp_args" );
  } );
my @db_opts = keys %{ get_DB_options_config( [ 'funcgen', 'dna', 'jdb' ], 1 ) };    #allow custom flag
my @opts_config = (
  @db_opts, qw( out_dir=s  collection=s@ dump_jaspar_fasta
    pfm_file=s pep_fasta=s   skip_blast        help ) );

GetOptions( \%opts, @opts_config ) or    #Catch unknown options
  pod2usage( -exitval => 1, -message => "Params are:\t@tmp_args" );

#Delete the help and man subs as we don't want to pass those around
#Getopt::Long qw(:config auto_help) would do this for us but the help
#is not verbose enough as it only prints the synopsis
#Not much point as it doesn't handle man, but do we need man?
#or could we rely on perldoc here?
delete $opts{help};
delete $opts{man};

#Rename collections array, so it makes sense if we ever require in this code
$opts{collections} = delete $opts{collection};

foreach my $key ( keys %opts ) {
  $opts{"-${key}"} = delete $opts{$key};
}

#Gash! Neither way works the we need here so either I stripe the =.* spec off the keys
#or I have to add the - prefix GRRR!

#Probably better the other way, as it will be less fiddling and will allow us
#to update specific opts hashes to localize groups of options in predefined hashes
#Also wont have to delete the man and help entries, and we will be able to use =s{,}
#No no no, we still have to deref the opts values!
#Better this way

if (@ARGV) {    #Catch trailing args
  pod2usage( -exitval => 1,
             -message => "You have specified unknown options/args:\t@ARGV" );
}

print "load_jaspar_matrices.pl @tmp_args\n";

#Sub everythign so we can re-use it in other scripts/modules if required

&main( \%opts );

sub main {
  my $opts = shift;
  my ($out_dir, $colls, $pfm_file, $jdb_version, 
      $pep_fasta, $dump_pfm_fasta, $skip_blast) = rearrange
      (['OUT_DIR', 'COLLECTIONS', 'PFM_FILE', 'JDB_VERSION',
        'PEP_FASTA', 'DUMP_JASPAR_FASTA', 'SKIP_BLAST'], %$opts );

  #assert_ref($colls, 'ARRAY');

  if ( !$out_dir ) {
    $out_dir = '.';
  }
  elsif ( !-d $out_dir ) {
    die("-out_dir does not exist:\t$out_dir");
  }

  if ( !( $pfm_file && -f $pfm_file ) ) {
    die( '-pfm_file needs to be a valid position frequency matrix file,' .
         " containing all the required PFMs:\t$pfm_file" );
  }

  die('Mandatory --jdb_user param not specified')
    if !defined $opts->{'-jdb_user'};
  die('Mandatory --jdb_host param not specified')
    if !defined $opts->{'-jdb_host'};
  die('Mandatory --jdb_name param not specified')
    if !defined $opts->{'-jdb_name'};
    
  if((! defined $pep_fasta) ||
      ! -f $pep_fasta){
    die("-pep_fasta is not defined or does not exists:\t$pep_fasta");      
  }  

  #Currently get undefs here
  #add support in process_DB_options?

  my $jdbh;
  my $dsn = sprintf( "DBI:%s:%s:host=%s;port=%s",
                     'mysql', $opts->{'-jdb_name'}, $opts->{'-jdb_host'},
                     $opts->{-jdb_port} || 3306 );

  eval {
    $jdbh = DBI->connect( $dsn, $opts->{'-jdb_user'}, $opts->{'-jdb_pass'},
                          { 'RaiseError' => 1 } );
  };

  my $error = $@;

  if ( !$jdbh || $error || !$jdbh->ping ) {
    die( 'Could not connect to database ' . $opts->{'-jdb_name'} . ' as user ' .
         $opts->{'-jdb_user'} . " using [$dsn] as a locator:\n" . $error );
  }

  if ($colls) {    #Validate vs DB
    assert_ref( $colls, 'ARRAY', 'Collections' );

    die('collection validation not yet implemented');
  }
  else {
    $colls = [qw( CORE PBM )];

#CORE isn't strictly a 'collection', although it is stored as such in the Jaspar DB
  }

  my $efg_db = create_Funcgen_DBAdaptor_from_options( $opts, 'pass', 1 ); 
  #Require as password here, but write user may not require one? 

  my $helper     = Bio::EnsEMBL::Funcgen::Utils::Helper->new( no_log => 1 );
  my $db_species = $efg_db->species;
  my $species    = $db_species;

  if ( $species !~ /_/ ) {    #assume we have failed to get a latin name
    die("$species latin species name, need to add a -species option");
  }
  else {                      #Convert to TAX format
    ( $species = ucfirst( lc($species) ) ) =~ s/_/ /;
  }

  my $analysis =
    $efg_db->get_AnalysisAdaptor->fetch_by_logic_name( $opts->{'-jdb_name'} );
  my $db_version = $efg_db->_get_schema_build( $efg_db->dnadb );

  if ( ! defined $analysis ) {

#Using Jaspar as the logic name here will prevent being able to load two versions along
#side each other which could complicate things

    $analysis =
      Bio::EnsEMBL::Analysis->new(
      -logic_name  => $opts->{'-jdb_name'},
      -db_version  => $jdb_version,
      -description => 'Position Frequency Matrices from the ' .
"<a href='http://jaspar.genereg.net/'>Jaspar Database</a> ($jdb_version)",
      -display_label => 'Jaspar' );
    $efg_db->get_AnalysisAdaptor->store($analysis);
  }
  else {    #rollback!
    my $delete_sql =
      'DELETE from object_xref where analysis_id=' . $analysis->dbID;

    #warn "Rolling back with:\t$delete_sql";
    $efg_db->dbc->db_handle->do($delete_sql);

    #Do we need finer grain control of rollback here with different analyses?
    #We might just want to rerun the blast step or this step?
    #is blast step dependant on the output of this step
    #is yes then we need to do everythign at the same time

    $delete_sql =
      'DELETE from binding_matrix where analysis_id=' . $analysis->dbID;
    $efg_db->dbc->db_handle->do($delete_sql);

#todo rollback motif_features and associated_motif_feature and regulatory_attribute ?
#add this rollback function to the Helper based on the analysis.
  }


  #Just do this anyway. It will warn if it already exists

  add_external_db( $efg_db, $db_species . '_core_Gene',
                     $db_version, 'EnsemblGene' );

  if ($dump_pfm_fasta) {
    &dump_jaspar_fasta( $jdbh, $out_dir, $colls );
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
  my @fsets  = grep {
    ( $_->feature_type->class eq 'Transcription Factor' ) ||
      ( $_->feature_type->class eq 'Transcription Factor Complex' )
  } @{ $fset_a->fetch_all_by_feature_class('annotated') };

  #Make ftypes non-redundant
  my %ftypes;

  #my @ftypes = grep { ! $seen{$_->name}++ } (map { $_->feature_type } @fsets);
  
  #warn "hardcoding for just 10 ftypes";
  
  map { $ftypes{ $_->feature_type->name } = $_->feature_type } @fsets; #[0-9];

  my $dnadb      = $efg_db->dnadb;
  my $gene_a     = $dnadb->get_GeneAdaptor;
  my $ftype_a    = $efg_db->get_FeatureTypeAdaptor;
  my $dbentry_a  = $efg_db->get_DBEntryAdaptor;
  my $bm_adaptor = $efg_db->get_BindingMatrixAdaptor;

#These are the minimum fields(and aliases) required to support preprocess_jaspar_rows
#So could add more, but can't remove any
#todo change this to use COLLATE for lc match
  my $sql =
'SELECT m.BASE_ID as id, m.VERSION as version, t.SPECIES as species, mp.ACC as accs, m.NAME as name '
    . 'FROM MATRIX m, MATRIX_PROTEIN mp, MATRIX_SPECIES ms, TAX t '
    . 'WHERE m.ID=mp.ID and m.ID=ms.ID and ms.TAX_ID = t.TAX_ID and LOWER(m.NAME) like ?';
  my $sth = $jdbh->prepare($sql);

  #Can't do an IN statement as the accs are cat'd within the field!!!
  my $sql2 =
'SELECT m.BASE_ID as id, m.VERSION as version, t.SPECIES as species, mp.ACC as accs, m.NAME as name '
    . 'FROM MATRIX m, MATRIX_PROTEIN mp, MATRIX_SPECIES ms, TAX t '
    . 'WHERE m.ID=mp.ID and m.ID=ms.ID and ms.TAX_ID = t.TAX_ID and mp.ACC like concat("%", ? "%")';
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

  my ( %ftype_matches, %no_matches, %peps_to_gene_ftype_names, %counts, %pfm_hits);

#todo put these counters in the hash? or remove hash in favour of individual counts
  my $ftype_to_pfm_complex    = 0;
  my $gene_to_pfm_complex     = 0;
  my $ftype_to_gene_only_full = 0;
  my $ftype_to_gene_full      = 0;
  my $ftype_to_gene_no_full   = 0;
  my $ftype_to_gene           = 0;
  my $ftype_to_pfm            = 0;
  my $multi_pep_dump          = 0;


  #warn "Hardcoding 10 ftypes for testing!";
  #@ftypes = @ftypes[0..9];

  #could probably change this to use keys
  #and then only access the %ftype hash when absolutely necessary

  foreach my $ftype ( values %ftypes ) {
    my $ftype_name = $ftype->name;
    print "Finding associations for:\t$ftype_name\n";
    ####  FeatureType name > PFM name > Gene accession/external_name ###
    $sth->bind_param( 1, '%' . lc($ftype_name) . '%', SQL_VARCHAR );
    $sth->execute;

    #Filters out older versions and non species specific rows
    #Also handles fuzzy matching
    my ( $pfm_matches, $non_species_rows ) =
      &preprocess_jaspar_rows( $sth, $species, $ftype_name );

    if ( ! scalar(@$pfm_matches) ) {
      $no_matches{$ftype_name}{ftype_to_pfm} = $ftype;
    }
    else {
      #We never actually use this
      #$ftype_matches{$ftype_name}{feature_type} = $ftype;

      foreach my $row_ref (@$pfm_matches) {
        my $found_gene = 0;
        my ( $id, $species, $accs, $name, $version ) = 
          ($row_ref->{id},   $row_ref->{species}, $row_ref->{accs},
           $row_ref->{name}, $row_ref->{version} );

        my $pfm_key = $id . '.' . $version . ' ' . $name;

        if ( $name =~ /::/ ) {
          $ftype_to_pfm_complex++;
        }

        foreach my $i ( 0 .. $#{$accs} ) {
          my $matched_display_label;
          my $seen_multi;
          my @genes = @{ $gene_a->fetch_all_by_external_name( $accs->[$i] ) };

          if ( scalar(@genes) > 1 ) {
            $counts{multi_acc_to_gene_hits}++;
          }
          elsif ( !scalar(@genes) ) {
            $ftype_matches{$ftype_name}{ftype_pfm_name} ||= [ $id, [] ];
            push @{ $ftype_matches{$ftype_name}{ftype_pfm_name} },
              $accs->[$i];    #No gene matches
             #This is very unlikely, as we are filtering out the non-species specific stuff
          }

          if ( $name =~ /::/ ) {
            print "Found ".scalar(@genes)." genes for TF complex hit $name with acc ".
              $accs->[$i]."\n";
          }

          foreach my $gene (@genes) {
            $found_gene = 1;
            my @xrefs = grep { lc( $_->display_id ) eq lc($ftype_name) }
              @{ $gene->get_all_DBLinks };

            #Match the ftype name to the gene external names
            if (@xrefs) {
              $counts{ftype_pwm_gene_and_name}++;    #This is a redundant count
              print "PFM (name) FeatureType: $ftype_name > PFM $id $name > " .
                $gene->display_id . ' ' . $accs->[$i] . "\n";

              if ($matched_display_label) {

                if ( !$seen_multi ) {
                  $counts{multi_gene_to_display_label_hits}++;

                  #This is likely due to paralogs/gene families
                  $seen_multi = 1;
                }
              }

              $matched_display_label = 1;
              $pfm_hits{$pfm_key}{ $accs->[$i] }{$ftype_name}{pwm_gene_name} ||= [];
              push @{ $pfm_hits{$pfm_key}{ $accs->[$i] }{$ftype_name}{pwm_gene_name} }, $gene;

            } ## end if (@xrefs)
            else {
              print "PFM (acc): FeatureType $ftype_name > PFM $id $name > " .
                $gene->display_id . ' ' . $accs->[$i] . "\n";

              $pfm_hits{$pfm_key}{ $accs->[$i] }{$ftype_name}{pwm_gene_accession} ||= [];
              push @{ $pfm_hits{$pfm_key}{ $accs->[$i] }{$ftype_name}{pwm_gene_accession} }, $gene;
            }
          } ## end foreach my $gene (@genes)
        } ## end foreach my $i ( 0 .. $#{$accs...})

        if ( ! $found_gene ) {
          $no_matches{$ftype_name}{pfm_to_gene} = $ftype;
          $ftype_to_pfm++;
        }
      } ## end foreach my $row_ref (@$pfm_matches)
    } ## end else [ if ( !scalar(@$pfm_matches...))]

    ### DO SOME COUNTS 
    
    if ( $ftype_matches{$ftype_name} ) { #Avoids auto-vivification of hash value

      if ( exists $ftype_matches{$ftype_name}{pwm_gene_name} ) {

        if ( !exists $ftype_matches{$ftype_name}{pwm_gene_accesssion} ) {
          $ftype_to_gene_only_full++;
        }

        $ftype_to_gene_full++;
      }

      if ( exists $ftype_matches{$ftype_name}{pwm_gene_accesssion} ) {

        if ( !exists $ftype_matches{$ftype_name}{pwm_gene_name} ) {
          $ftype_to_gene_no_full++;
        }

        $ftype_to_gene++;

      } ## end if ( exists $ftype_matches...)
    } ## end if ( $ftype_matches{$ftype_name...})


    ####  FeatureType name > Gene external_name > PFM accession/name ###
    my @genes = @{ $gene_a->fetch_all_by_external_name($ftype_name) };

    if ( scalar( @genes > 1 ) ) {
      $counts{'multi_ftype_to_gene_hits'}++;
    }
    elsif ( scalar(@genes) ) {

      foreach my $gene (@genes) {

        foreach my $pep ( map { $_->translation } @{ $gene->get_all_Transcripts } ){

          if ( defined $pep ) {
            $peps_to_gene_ftype_names{ $pep->stable_id } = [ $gene, $ftype_name ];
          }
        }
      }
    }
    else {
      $no_matches{$ftype_name}{ftype_to_gene} = $ftype;
    }

    my @validated_genes = ();

    foreach my $gene (@genes) {
      my $seen_pwm       = 0;
      my %seen           = ();
      my @dbprimary_accs = map { $_->primary_id }
        ( grep { !$seen{ $_->primary_id }++ } @{ $gene->get_all_DBLinks } );

      #dbprimary accs maybe redundant as they can be associated to
      #gene|transcript|translation from the same external_db

      foreach my $xref_acc (@dbprimary_accs) {
        $sth2->bind_param( 1, $xref_acc, SQL_VARCHAR );
        $sth2->execute;
        my ( $rows, $non_species_rows ) =
          &preprocess_jaspar_rows( $sth2, $species, undef, $xref_acc );

        foreach my $row_ref (@$rows) {
          my ( $id, $species, $name, $version ) = 
           ($row_ref->{id},   $row_ref->{species},
            $row_ref->{name}, $row_ref->{version});
          
          my $pfm_key = $id.'.'.$version.' '.$name;  
          $seen_pwm = 1;

          if ( $name =~ /(^|::)${ftype_name}(::|$)/i )
          {    #Ignore previously seen fully validated
            $counts{ftype_gene_pwm_and_name}++;    #Count for sanity check

            #These are being pushed for each acc
            #not for each ftype
            #
            print "Gene (name) FeatureType: $ftype_name > Gene $xref_acc " .
              $gene->display_id . " > PFM $id $name\n";
            push @validated_genes, $gene;

            #              $gene_first_validated = 1;
          }
          else {    #These are the new ones we're interested in!
            $counts{ftype_gene_pwm_no_name}++;

            if ( $name =~ /::/ ) {
              $gene_to_pfm_complex++;
            }

#Need to handle versions here
#We don't yet know how many this is affecting
#$ftype_matches{$ftype_name}{gene_pwm_accession} ||= [];
#push @{$ftype_matches{$ftype_name}{gene_pwm_accession}}, [$gene, $id, $name, $version, $xref_acc];
#Need to remove these from the no match cache? Or handle this later

            $pfm_hits{$pfm_key}{$xref_acc}{$ftype_name}{gene_pwm_accession} ||=[];
            push @{$pfm_hits{$pfm_key}{$xref_acc}{$ftype_name}{gene_pwm_accession}}, $gene;#do we even need the name here?
              
              print "GENE (acc): FeatureType $ftype_name > ".$gene->display_id." ($xref_acc) > PFM $id ${name}\n";
            }
          }
                 
          #Now dump the pep fasta for the transcript which matches the ftype name
          my @xrefs = grep {lc($_->display_id) eq lc($ftype_name)} @{$gene->get_all_DBLinks};
          my %seen_xrefs;
          
     
            #Now cache those with no acc/name hits
          if ( !$seen_pwm ) {
            #or just count here?
            $no_matches{$ftype_name}{gene_to_pfm} = $ftype;
          }
        }
      }

      #Validate we get full validation in both directions

      #This appears to be a failure to take the highest version
      #This is because for the Gene first we are preprocessing rows on separate
      #accs, so it never sees the the other version
      #where as for the PWM first approach we are preprocessing the rows based
      #on the name, which bring back both rows.

     #We need to handle version in the main cache for the Gene first approach
     #but for acc only matches
      warn "TODO !!! Can we convert match validation to use pfm_matches or use booleans !!!";

    #if(exists $ftype_matches{$ftype_name}){
    #  if(@validated_genes &&
    #     ! exists $ftype_matches{$ftype_name}{pwm_gene_name}){
    #    die("Found Gene first validated match, but not PFM first validated match:\t".$ftype_name);
    #  }
    #  elsif(exists $ftype_matches{$ftype_name}{pwm_gene_name} &&
    #        ! @validated_genes){
    #    die("Found PWM first validated match, but not Gene first validated match:\t".$ftype_name);
    #  }
    #  elsif(@validated_genes && #implicit the pwm_gene_name exists
    #    (scalar(@validated_genes) != scalar(@{$ftype_matches{$ftype_name}{pwm_gene_name}}))){
    #
    #    #This currently happens for only for E2F1 in human
    #
    #   warn("Found mismatch between the number of validate PWM/Gene first matches for $ftype_name:\t".
    #    scalar(@{$ftype_matches{$ftype_name}{pwm_gene_name}}).' vs '.scalar(@validated_genes) );
    #  }
    #}
   }
      
          #Not getting the converage of the previous approach...yet!

          #Before reverse search, multi acc and TF complex support

#Found FeatureType - Gene matches (inc non-display_id):  29
#Found no FeatureType - PFM matches:     61
#Found only Feature - PFM matches (i.e. no PFM-Gene acc match):  1
#Both of the following counts ere probably the result of different genes giving rise to the same protein product through identical splicing:
#Multiple Gene hits for a single acc:    1
#Including a matching display label:     1

          #Now with gene first search? and fuzzy matching and catd acc handling
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

          #Now after fuzzy name match bug fix

          #should probably do this before the loop
          #so we can do the blast in line

   (my $blast_results = $pep_fasta) =~ s/.*\///g;
   $blast_results =~ s/\.fa(asta)//;
   #this subing does not work, fasta remains in string 
   $blast_results = $out_dir."/Jaspar_vs_${blast_results}.blastp.txt";
   
   warn "HARDCODING PEP FAST TO $out_dir/Jaspar_vs_homo_sapiens_core_76_38.pep.blastp.txt";
   $blast_results = $out_dir.'/Jaspar_vs_homo_sapiens_core_76_38.pep.blastp.txt'; 

   if($skip_blast){
      
     if(! -f $blast_results){
       die("Cannot -skip_blast as blast output does not exist:\t".$blast_results);  
     }
       
     print "SKIPPING BLAST STEP!\nUsing existing blast results:\t$blast_results\n";
   }
   else{
     warn "TODO need to log blast output";
     &blast_jaspar_matrices( $pep_fasta, $out_dir, $blast_results );
   }



   my $blast_fh = open_file($blast_results);
   my ( $pfm_id, $line, $ens_sid, $pfm_acc, $ftype, %blast_hits, $name );

   warn "ARE WE handling versions here?";

   #This is currently returning many many PFM hits for each gene
   #many of varying quality
   #previous data had multi species hits of >90% id
   #and good coverage
   #clearly some went down to 70%
   #but the are lots of 100% hits
   #This is a by product of only blasting against selected proteins and not all of them!!!!


  #Do some no_match preprocessing before we consider the blast hits
  foreach my $ftype(keys %no_matches){
    
    if(scalar keys(%{$no_matches{$ftype}}) != 2){# we have at least 1 match
      delete $no_matches{$ftype};
    }  
  }
  

   my $used_blast_hits = 0;

   


   while ( ( $line = $blast_fh->getline ) && defined $line ) {
     ( $pfm_id, $ens_sid ) = split( /\t/, $line );
     ( $pfm_id, $pfm_acc, $name ) = split(/__/, $pfm_id );
     my $pfm_key = $pfm_id.' '.$name;
          
     if ( exists $peps_to_gene_ftype_names{$ens_sid} ) {
       my ( $gene, $ftype ) = @{ $peps_to_gene_ftype_names{$ens_sid} };
       #Only set this if we don't already have a good xref match
       #already contains version from dumps 

       if(! exists $pfm_hits{$pfm_key}{$pfm_acc}){
         #This will autovivify the $pfm_id key, bu thtat fine as we are setting it in here anyway
         print "Found no xref matches for PFM $pfm_key($pfm_acc), defaulting to best blast hit:\t $ens_sid(".$gene->stable_id.") $ftype\n"; 
         $pfm_hits{$pfm_key}{$pfm_acc}{$ftype}{blast} = [$gene];
         $used_blast_hits++;
    
         if(exists $no_matches{$ftype}){
           delete $no_matches{$ftype};
         }
         
       }
       else{
         warn "TODO: Log if we have a blast hit which differs or matches the existing hit(s)";
       }
     }
     else{
       print "Best PFM $pfm_key($pfm_acc) blast hit $ens_sid is not a recognisable TranscriptionFactor FeatureType\n";
     }     
  }
  
  #break this down into full and non-display_id matches
  print "FeatureType > PFM > Gene matches:\n";
  print "Total redundant fully validated matches, FeatureType > PFM > Gene name:\t".$counts{ftype_pwm_gene_and_name}."\n"; 
  print "Total non-redundant wrt FeatureType matches (all/only name match/only no name match):\t\t".
    $ftype_to_gene_full.'/'.$ftype_to_gene_only_full.'/'.$ftype_to_gene_no_full."\n";
  #.scalar(keys %{$counts{any_matches}})."\n";
  
  
  print "Found no FeatureType - PFM matches:\t\t\t\t".scalar(keys %{$no_matches{ftype_to_pfm}})."\n";
  print "Found only Feature - PFM matches (i.e. no PFM-Gene acc match):\t".scalar(keys %{$no_matches{pfm_to_gene}})."\n";
  
  print 'Both of the following counts ere probably the result of different genes giving'.
    " rise to the same protein product through identical splicing:\n"; 
  print "Multiple Gene hits for a single acc:\t\t\t\t".$counts{multi_acc_to_gene_hits}."\n";
  print "Including a matching display label:\t\t\t\t".$counts{multi_gene_to_display_label_hits}."\n";
  print "Total FeatureType to PFM TF complexes matches:\t\t\t$gene_to_pfm_complex\n"; 
      
  my $total_gene_matches = 0;
  my $total_pfm_matches  = 0;
  my $distinct_ftype_pfm = 0;
  warn "TODO fix counts";   
  #foreach my $href(values %ftype_matches){
    
  #  if(exists $href->{gene_pwm_accession}){
  #    $distinct_ftype_pfm++;
    
      #really need to count those which do not currently have any pwm_gene matches
      #These will be the new ones!
      #Theoretically these shoudl awlays be previously unseen
      #as they are not fully validated
      #and if there is both a gene_pwm_accession and a pwm_gene_accession
      #path, then this should have been captured as a pwm_gene_name hit
  #    my $last_gene;
      
  #    foreach my $gene_pwm_hit(@{$href->{gene_pwm_accession}}){
  #      $last_gene ||= $gene_pwm_hit->[0];
        
  #      if($last_gene->stable_id ne $gene_pwm_hit->[0]->stable_id){
  #        $total_gene_matches++;
  #      }
     
  #    $total_pfm_matches ++;
  #    }
  #  }  
  #}
 
  print "\nFeatureType > Gene > PFM matches:\n";
  print "Ignoring fully validated matches(nr), FeatureType > Gene > PFM name:\t".$counts{ftype_gene_pwm_and_name}."\n";
  #print "Distinct FeatureTypes with PFM matches (not name validated):\t\t${distinct_ftype_pfm}\n";
  #print "Total Gene matches with PFM matches (not name validated):\t\t".$total_gene_matches."\n";
  #print "Total PWM matches (not name validated):\t\t\t\t\t".$total_pfm_matches."\n";
  print "Total Gene to PFM TF complexes matches:\t\t\t\t\t$gene_to_pfm_complex\n"; 
  print "All PFM matches (not distinct wrt ftype or gene, not full name match):\t".
    $counts{ftype_gene_pwm_no_name}."\n";
  print "Multi pep dumps:\t$multi_pep_dump\n";   
  print "Total PFM blast hits used:\t$used_blast_hits\n";
  print "Total distinct PFMs with any hit ".scalar(keys(%pfm_hits))."\n";
  

  
  ### CACHE PFMS ###
  my $fh = open_file($pfm_file);
  my (%pfms, $id);
  
  while(($line = $fh->getline) && defined $line){
    
    if($line =~ /^>([^ ]+) /){
      $id = $1;
     # warn "caching pfm freqs for $id";
      $pfms{$id} = ''; 
    }
    else{
      $pfms{$id} .= $line;    
    }
  }
  
  #warn "cached ".scalar(keys %pfms)." pfms";
  #if(! exists $pfms{'MA0492.1'}){
  #  die('BLARRT');  
  #}
  
  
 
  #my ($version, $acc);
  #The multi acc entries all appear to be TF complexes
  #So long as there aren't >1 ftype hit per acc do we care?
  #Different accs may have the same ftype hit
  #We should probably warn about this for now
  #but it probably won't happen
  #should we consider gene here too?
  
  #In the cse of >1 ftype hit/acc, then take the one with the best hit?
  #or throw?

  my @hit_types = qw(pwm_gene_name pwm_gene_accession gene_pwm_accession blast);
  my %bm_descs =
   (pwm_gene_name      => 'Name/Accession association',
    pwm_gene_accession => 'Accession association',
    gene_pwm_accession => 'Accession association',
    blast              => 'Blast association');

  my @colls_skipped;
  my $bms_loaded;

  foreach my $pfm_key(keys %pfm_hits){
    #This is now a lot more simple as we have already keyed on accession!
    
    my ($pfm_id_version, $name) = split(/ /, $pfm_key);
    my ($pfm_id, $version)      = split(/\./, $pfm_id_version);
    my $acc_hits                = $pfm_hits{$pfm_key};
    my (%seen_ftypes, $bm_ftype, %validation_types, $store_assoc_ftypes, @assoc_ftypes);
    my @accs = keys %{$acc_hits};
    
    if(! scalar(@accs)){
      die("Found pfm_hit entry for $pfm_key with no acc hits");  
    }
    elsif(scalar(@accs) > 1){
      
      #Deal with TF complexes here!
      #i.e. we need to store an ftype based on the name
      #associate that with the binding_matrix
      #then store associated ftypes to the bm based on what we find below.
      
      #These can then be used to pull back the relevant bms when doing the mapping for a given feature set 
      
      
      if($name !~ /::/){
        die("Multiple accession are not yet supported for PFM which are not Transcription Factor Complexes:\t$name");  
      }
   
      $store_assoc_ftypes = 1;
      $bm_ftype           = $ftype_a->fetch_by_name($name);   
      
        
      if(! defined $bm_ftype){
        #define and store here
       
        print "Storing new Transcription Factor Complex FeatureType:\t$name\n";
        
        $bm_ftype = Bio::EnsEMBL::Funcgen::FeatureType->new
                     (-name         => $name,
                      -class        => 'Transcription Factor Complex',
                      -description  => '$name binding',
                      -so_accession => 'SO:0000235',
                      -so_name      => 'TF_binding_site');
        $ftype_a->store($bm_ftype);
        
      }  
      #validate it doesn't already exist?
      
      $ftypes{$bm_ftype->name} = $bm_ftype;
      $bm_ftype = $bm_ftype->name;
    }
    

    
    foreach my $acc(@accs){
      my @ftypes = keys %{$acc_hits->{$acc}};
      #We need to restrict this to 1 so we can assoicate the bm with just 1 ftype        
      
      if(scalar(@ftypes) > 1){
        die("Found pfm_hit acc hit for $pfm_key $acc with more than 1 FeatureType:\t @ftypes");
        #If this happens we could arbitrarily take the best match type  
      }
      elsif(! @ftypes){
        die("Found empty acc hit for $pfm_key $acc");
      }

      my $ftype = $ftypes[0];
      $bm_ftype ||= $ftype;

      if(exists $seen_ftypes{$ftype}){
        die( "Already seen match for FeatureType $ftype:\t$pfm_key $acc"); 
        #CHANGE THIS TO A WARN AND STORE DUPLICATES IF IT IS A PROBLEM
        #as it maybe possible that we map to the same ftype via different genes/accs
        #if we are dealing with a gene family/paralog
        #we may also see an ftype mapping to > 1 gene for 1 acc 
        #for the same reason
         
      }
      
      if($store_assoc_ftypes){
        push @assoc_ftypes, $ftypes{$ftype};        
      }
      
      my ($gene, %seen_genes, $last_ht);
      
          
      for my $ht(@hit_types){
        
        if(exists $acc_hits->{$acc}->{$ftype}->{$ht}){
          my @genes = @{$acc_hits->{$acc}->{$ftype}->{$ht}}; 
         
          if(! @genes){
            die("Failed to find any gene hits for feature_type hit for $pfm_key $acc $ftype hit_type $ht");
          }
         
          foreach my $gene(@genes){
         
          if(! defined $gene){
            die("Feature type gene hit for $pfm_key $acc $ftype hit_type $ht is not defined");
          }
          elsif(exists $seen_genes{$gene->stable_id}){
            warn "Skipping $ht hit for $pfm_key $acc $ftype as we have already seen a better hit";
            next; #$ht 
          }  
          
          ### STORE THE XREF! ###
          #This is quicker than using the API
          #especially if we see it more than once
          #Could just use display_xref here instead
          my $display_name = $helper->get_core_display_name_by_stable_id($efg_db->dnadb, $gene->stable_id, 'gene');
          $last_ht         = $ht;
          my $linkage_txt;
         
         
          if($ht eq 'pwm_gene_name'){
            $linkage_txt = 'FeatureType($name) > PFM($acc) > Gene($name) : $pfm_key';   
          }
          elsif($ht eq 'pwm_gene_accession'){
            $linkage_txt = 'FeatureType($name) > PFM($acc) > Gene : $pfm_key';   
          }
          elsif($ht eq 'gene_pwm_accession'){
            $linkage_txt = 'FeatureType($name) > Gene($acc) > PFM : $pfm_key';   
          }
          elsif($ht eq 'blast'){
            $linkage_txt = 'External name association to support $pfm_key $acc $name PFM blast'; 
          }
          else{
            die("$ht hit type is currently not supported");  
          }
     
          my $dbentry = Bio::EnsEMBL::DBEntry->new
           (-dbname                 => $db_species.'_core_Gene',
            -release                => $db_version,
            -status                 => 'KNOWNXREF',
            -display_label_linkable => 1,
            -db_display_name        => $efg_db->dnadb->dbc->dbname,
            -db_display_name        => 'EnsemblGene',
            -type                   => 'MISC',
            -primary_id             => $gene->stable_id,
            -display_id             => $display_name,
            -info_type              => 'MISC',
            -info_text              => 'GENE',
            -linkage_annotation     => $linkage_txt,
            -analysis               => $analysis);
    
          #DBI->trace(2);
          $dbentry_a->store($dbentry, $ftypes{$ftype}->dbID, 'FeatureType');#1 is ignore release flag 
        }
        }
      }      
    
      if(! $last_ht){# keys(%seen_genes)){
        die("Failed to find expected hit_type for feature_type hit for $pfm_key $acc $ftype:\t".
          join(keys(%{$acc_hits->{$acc}->{$ftype}})));
      }
      
      $validation_types{$last_ht} = 1;
    }

    #We may have hit types from multiple gene hits on one ftype 
    #How can we easily represent that in the description
    #xref validated or blast validated match
    #We really need another field for this
    #can we slip this in as another last minute patch
    #we need to handle interspecies blast hits
    #put this in the description for now
    #or last minute patch?

    ### STORE THE BINDING_MATRIX ###
    #validate we have a bm_ftype?
    my $vtype;
    
    #Get best validation method 
    
    foreach my $ht(@hit_types){
    
      if(exists $validation_types{$ht}){
        $vtype = $bm_descs{$ht};
        last;  
      }
    }    
      
    $bm_ftype = $ftypes{$bm_ftype}; 
       
    if(@assoc_ftypes){
      
      $bm_ftype->adaptor->store_associated_feature_types($bm_ftype, \@assoc_ftypes, 1);#rollback flag  
    }   
       
    if(! exists $pfms{$pfm_id_version}){
      warn("Could not acquire frequecies for PFM $pfm_id_version. This is probably a collection not contained in the Jaspar CORE set");  
      #These will need fetching from the ARCHIVE!!!
      push @colls_skipped, $pfm_id_version;
      next;
    }    
       
    my $bm = Bio::EnsEMBL::Funcgen::BindingMatrix->new
     (-name         => $pfm_id.'.'.$version,
      -analysis     => $analysis,
      -feature_type => $bm_ftype,
      -description  => $vtype,
      #description is not currently used in the web interface
      #this qualitative info should probably be loaded as another xref?
      #although this may be useful in the interface as a confidence level?
      #we should have an explicit version
      -frequencies  => $pfms{$pfm_id_version});
   $bm_adaptor->store($bm);
   
   
   $bms_loaded++; 

     #might need this if we have to pick a best hit across  different ftypes
     #   foreach my $sid(keys %{$pfm_hits{$pfm_id}}){
     #   
     #     for my $i(0..$#hit_types){
     # 
     #       if(exists $pfm_hits{$pfm_id}{$sid}{$hit_types[$i]}){
     #         
     #         if(! defined $last_hit_type ||
     #            ($last_hit_type > $i)){
     #           $hit_info = $pfm_hits{$pfm_id}{$sid}{$hit_types[$i]};
     #           $hit_id = $sid;
     #           $last_hit_type = $i;
     #           last;
     #         }
     #       }
     #     }    
  }  
 
  #we probably want to print those ftype which match the pfm name and nothing else
  #it's unlikely that these won't have blast hits
  #but we may be able to see some obvious matches which aren't handled by this pipeline
  
  print "Failed to identify any matches for FeatureTypes:\n\t".join("\n\t", keys %no_matches)."\n";
  print "Skipped ".scalar(@colls_skipped)." non CORE PFMs:\n@colls_skipped\n"; 
  print "Loaded $bms_loaded Binding Matrices\n";
  return;
}
 

#todo dirs and filenames in here probably not the best
#but rushing now
#could probably boot strap this a bit, but using the DB
#to check that all the files are there?


#if we want to do this in line
#we really need to dump individual gene peps
#and have just one joint matrix fasta
#Doing it in line will probably be quite a bit slower?

 
sub blast_jaspar_matrices{
  my $pep_fasta  = shift;
  my $matrix_dir = shift;
  my $blast_out  = shift;
  my $matrix_ids = shift;
  
  if(-f $blast_out){
    print "Removing existing blast output file:\t$blast_out\n";
    unlink $blast_out || die("Failed to remove pre-existing blast output file:\t$blast_out\n$!\n");  
  }
  
  
  opendir(DIR, $matrix_dir) or die("Cannot open matrix fasta dir:\t$matrix_dir\n$!");
  my @files = grep {$_ =~ /\.fa$/} readdir DIR;
  closedir DIR;
  
  
  if(! @files){
    die("Failed to identify any fasta files (.fa) in the input directory:\t".$matrix_dir);
  }
  else{
    print 'Blasting '.scalar(@files)." Jaspar PFM proteins against $pep_fasta\n";  
  }
  
  foreach my $pfm_fasta(@files){
    
    next if $pep_fasta =~ /${pfm_fasta}$/;
    #This should happen as the suffixes are different
    
    
    #this does not currently validate that the fasta files in this dir 
    
    #-b 1 top hit
    #-m 8 Tabular format
    
    #do we even need to write to file?
    
    #it appear get_gene_from_protein.pl was used for this
    #but it didn't handle coverage mismatches or % ID
    
    #change this to gapped alignment and we need to ensure coverage of query
    
    my $cmd = "blastall -p blastp -d $pep_fasta -M BLOSUM80 -m 8 -b 1 -i $matrix_dir/$pfm_fasta >> $blast_out";
    
    if(run_system_cmd($cmd, 1) != 0){ #No exit flag
      #This can be caused by a gene accession being used instead of a uniprot accession
      #This means pfetch returns the EMBL (cdna/dna) seq instead of the UniProtKB prot seq
      #It seesm like all the yeast accs are gene names rather than protein accessions
      #e.g. YHL020C
      #SwissProt is broken/old
      #Need to update to use unirpot rest services
      #warn "FAILED TO BLAST $pfm_fasta\n";
      
      #we are also getting the following error for each blast
      
      #but i can't reproduce them on the cmdline
      #Error:  No such file or directory
  
      #and also these?
      #The invlaid query seqs are due to pfetch returning EMBL c/dna records for genes accessions
      #need to implement uniprot rest API for this
      #[blastall] WARNING: Sequence number 1 had length 0
      #[blastall] WARNING: Sequence number 1 had length 0
      #[blastall] WARNING: MA0434.1: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
      #[blastall] WARNING: MA0434.1: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
      #Child process exited with value 1
      #Error:  No such file or directory
      #System command returned non-zero exit code:     blastall -p blastp -d /nfs/users/nfs_n/nj1/scratch_prd/fasta/homo_sapiens/homo_sapiens_core_75_37.pep.fa -M BLOSUM80 -m 8 -b 1 -i .//MA0434.1_YPR013C.fa >> .//Jaspar_vs_homo_sapiens_core_75_37.pep.blastp.24811.txt at /nfs/users/nfs_n/nj1/src/ensembl-funcgen/modules/Bio/EnsEMBL/Funcgen/Utils/EFGUtils.pm line 1433.
      
      #[blastall] WARNING: Sequence number 1 had length 0
      #[blastall] WARNING: MA0365.1: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
      #[blastall] WARNING: MA0365.1: Could not calculate ungapped Karlin-Altschul parameters due to an invalid query sequence or its translation. Please verify the query sequence(s) and/or filtering options
      #Child process exited with value 1
      #Error:  No such file or directory
      #System command returned non-zero exit code:     blastall -p blastp -d /nfs/users/nfs_n/nj1/scratch_prd/fasta/homo_sapiens/homo_sapiens_core_75_37.pep.fa -M BLOSUM80 -m 8 -b 1 -i .//MA0365.1_YLR176C.fa >> .//Jaspar_vs_homo_sapiens_core_75_37.pep.blastp.24811.txt at /nfs/users/nfs_n/nj1/src/ensembl-funcgen/modules/Bio/EnsEMBL/Funcgen/Utils/EFGUtils.pm line 1433.
      
      
    }  
  }
  
  return $blast_out; 
} 


#todo count/warn skipped versions
 
sub preprocess_jaspar_rows {
  my ($executed_sth, $species, $name, $acc) = @_;
  #todo validate params
  my %rows;
  my @non_species_rows;
  my $all_species = 0;
  
  if($name && $acc){
    die('Name and Accession arguments are mutually exclusive');  
  }
 
  if(defined $species &&
     ($species =~ /^all$/i)){
    undef $species;
    $all_species = 1;    
  }
  #todo probably want to process non species rows properly too? 
  #Can't drop species from this, as we would get merging of accs across species
  
  while(my $href = $executed_sth->fetchrow_hashref){
    
    if(defined $species && 
       ($href->{species} ne $species)){
      push @non_species_rows, $href;
    }
    else{
      
      #Matrices with different species are given different versions
      #and different internal IDs, but the same BASE_ID
      #This allows different MATRIC_PROTEIN/ANNOTATION records
      #per species
      #But the MATRIX_SPECIES table suggestes that we also get > 1 species per internal ID
      #I suspect this table is/should never be used and is effectively redundant
      
      if((! defined $species) && 
         (! $all_species)){ #Only relevant if we haven't already filtered     
      
        if(exists $rows{$href->{id}} && 
           ($href->{species} ne $rows{$href->{id}}{species})){
          die('Found muliple associated species for PFM '.$href->{id}.
            ":\t".$href->{species}.' '.$rows{$href->{id}}{species});
          #This is not yet supported in dump_jaspar_pfm_fastas   
          
          #So it seems different version can be on teh same species or different
          #species
          #We need to be able to ignore this if we want to dump 
          #them in a species agnostic way
              
        }
      }
      
      if($name){ #Validate the query name 
 
        if ($href->{name} !~ /(^|::)${name}(::|$)/i){  #Allow TF complex matches
          warn "Skipping fuzzy name match:\t\t$name (Ensembl) vs ".$href->{name}." (Jaspar)\n";
          next;
        }
        elsif($href->{name} =~ /::/){
          print "Matched TF complex:\t$name (Ensembl) vs ".$href->{name}." (Jaspar)\n";
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
      
        if($acc){ #Validate the query acc
        
          if(grep /^${acc}$/, @accs){
            $href->{accs} = [$acc];
          }
          else{  
            warn "Skipping fuzzy accession match:\t$acc (Ensembl) vs @accs (Jaspar)\n";
            next; #$href
          }
        }
        else{
          $href->{accs} = \@accs;    
        }
        $rows{$href->{id}} = $href;
      }
    }         
  }
  
  return ([values %rows], \@non_species_rows);
}
 
 
sub dump_jaspar_fasta{
  my $jdbh    = shift;
  my $out_dir = shift; 
  my $colls   = shift;
  #my $pfms   = shift;
 
  assert_ref($jdbh, 'DBI::db', 'Jaspar DB Handle');
 
 #todo remove any pre-existing jaspar files?
 #these will be picked up by the blast step
 #although their IDs should not existing in the DB
 #else they would be over-written here
 
  if(! $out_dir){
    $out_dir = '.';  
  }   
  elsif(! -d $out_dir){  
    die("Output directory does not exist:\t$out_dir");
  }

  print "Dumping Jaspar PFM proteins for collections:\t".join(' ', @$colls)."\n";

  my $sql = 'SELECT m.BASE_ID as id, m.VERSION as version, t.SPECIES as species, mp.ACC as accs, m.NAME as name '.
   'FROM MATRIX m, MATRIX_PROTEIN mp, MATRIX_SPECIES ms, TAX t '.
   'WHERE m.ID=mp.ID and m.ID=ms.ID and ms.TAX_ID = t.TAX_ID and m.COLLECTION in("'.join('", "', @$colls).'")';
  my $sth = $jdbh->prepare($sql);
  $sth->execute;
  my ($rows) = preprocess_jaspar_rows($sth, 'all');#No species restriction!  
  my ($last_id, @fails);
  my $skipped = 0;
  my $wrote   = 0;
  
  foreach my $row(@$rows){  
    #pfetch may return a header with a .1 suffix, and maybe others
    #> pfetch Q6LBK7
    #>Q6LBK7.1 Q6LBK7_HUMAN C-myc protein (Fragment)
    #QIPELENNEKAPKVVILKKATAYILSVQAEEQKLISEEDLLRKRREQLKHKLEQLRNSCA
    
    foreach my $j_acc(@{$row->{accs}}){
      #Assuming that these are all uniprot accessions 
      #pfetch will bring anything back from UniProtKB (inc current data) and EMBL
      
      #actually pfetch will alternately return a translation or a cdna seq (or no match)
      #as it alternates between DBs if there are entries in both.
      #and it's hard to tell which will be returned first so have to try twice
      #Could use magical undocumented -d pep argument 
      #-d pep is not uniprot?
      #So some translations must be in EMBL, and -d pep does not bring these back
      #This also requires a version number to be specified, which we do not have.
         
      #Need to use Uniprot webservice here!
    
      my ($pf_fasta, @seq);
      
      for(0..1){
        $pf_fasta = run_backtick_cmd('pfetch '.$j_acc);
        (undef, @seq) = split("\n", $pf_fasta);   
       
        if((! defined $seq[0]) ||
           ($seq[0] =~ /^[agct]+$/i)){
          @seq = ();
        } 
        else{
          last;  
        }
      }
         
      #Should  probably merge fastas from complexes into the same fasta
      #and remove acc from file name.
    
      if(! @seq){
        push @fails, $row->{id};
      }
      else{
        #print 'Writing fasta for '.$row->{id}.' '.$j_acc."\n";
     
        #We really only need the id and the accession in the header
        $pf_fasta = '>'.$row->{id}.'.'.$row->{version}.'__'.$j_acc.'__'.$row->{name}."\n".join("\n", @seq);
       
        #can we do this from one query file, instead of many?
        #currently this is overwriting fasta for matrices which have > 1 acc (complexes?)
            
        my $fh = open_file($out_dir.'/'.$row->{id}.'.'.$row->{version}.'_'.$j_acc.'.fa', '>');
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
   
  print "Wrote $wrote Jaspar PFM protein fasta files\n";
  print "Skipped $skipped Jaspar records due to > 1 version being found\n"; 

  return;
}

1;
