#!/usr/bin/env perl

# Copyright [1999-2018] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
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

use strict;
#use warnings::unused;
use autodie;
use feature qw(say);

use Cwd 'abs_path';
use Config::Tiny;
use Bio::EnsEMBL::Utils::Logger;
use Bio::EnsEMBL::Utils::Exception qw(throw);
#use Bio::EnsEMBL::OntologyXref;
use Bio::EnsEMBL::DBEntry;
use Bio::EnsEMBL::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::DBSQL::DBEntryAdaptor;
use Getopt::Long;
use File::Basename;
use List::MoreUtils qw(uniq);
use Bio::EnsEMBL::Funcgen::ReadFile;
use Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration;
use Bio::EnsEMBL::Funcgen::FeatureType;
use DBI;
use File::Basename;

# use DateTime;

#TODO Dry run implementation
#TODO Rollback function
#TODO POD for every subroutine
#TODO Logger
#TODO PerlCritic 
#TODO check for external db availability in verify_basic_objects
#TODO healthchecks, ie.invalid br/tr values, invalid local/download url
#TODO use state variable for $control_db_ids

main();

sub main {
    my ( $csv, $config, $help, $dry );

    # ----------------------------
    # read command line parameters
    # ----------------------------
    GetOptions(
        'i=s'  => \$csv,
        'c=s'  => \$config,
        'h'    => \$help,
        'help' => \$help,
        'n'    => \$dry,      #not implemented
    );

    # ------------------------------------------------------
    # display usage and exit if anything critical is missing
    # ------------------------------------------------------
    if ( $help || !$csv || !$config ) {
        usage();
        exit;
    }

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
    if ( !-e $csv ) {
        $logger->error(
            'Input csv file ' . abs_path($csv) . ' doesn\'t exist!',
            0, 1 );
    }

    $logger->info( 'CSV file used: ' . abs_path($csv) . "\n",       0, 1 );
    $logger->info( 'Config file used: ' . abs_path($config) . "\n", 0, 1 );

    # --------------------
    # read the config file
    # --------------------
    my $cfg = Config::Tiny->read($config);

    # ----------------
    # get current date
    # ----------------
    my ( $sec, $min, $hour, $mday, $mon, $year, $wday, $yday, $isdst )
        = localtime(time);
    $cfg->{date} = ( 1900 + $year ) . '-' . ++$mon . '-' . $mday;

    # -------------------
    # empty data 
    # -------------------
    #empty_data($cfg);
    


    # -------------------
    # read input csv file
    # -------------------
    $logger->info( 'Reading metadata from ' . abs_path($csv) . '... ', 0, 0 );
    my ( $control_data, $signal_data, $control_sets )
        = get_data_from_csv( abs_path($csv), $logger );

	$logger->info( "done\n", 0, 1 );

    my $control_entries = keys %{$control_data};
    my $signal_entries  = keys %{$signal_data};
	
	
    $logger->info( $control_entries . " control entries found\n", 1, 1 );
    $logger->info( $signal_entries . " signal entries found\n",   1, 1 );
	

	### first test
	#exit;


    # ------------------------------------------------------------
    # connect to funcgen tracking db, fetch all necessary adaptors
    # ------------------------------------------------------------
    $logger->info( 'Connecting to ' . $cfg->{efg_db}->{dbname} . '... ',
        0, 0 );
    my $adaptors = fetch_adaptors($cfg);
    $logger->info( "done\n", 0, 1 );

	# second test
	#exit;
    # -----------------------------------------------
    # verify that fundamental objects exist in the db
    # -----------------------------------------------
    verify_basic_objects_in_db( $control_data, $signal_data, $adaptors,
        $logger );

 	#Thirth test
	#exit;

    # ------------
    # registration
    # ------------
	
    my %epi_Ihec;
    
	
	#register control files 
	my $control_experiments_ids = register_ctl($logger, $control_data, $adaptors, $cfg, \%epi_Ihec);

	#fourth control
	#exit;
	#register signal files
    for my $entry ( values %{$signal_data} ) {
        register_signal ( $logger, $entry, $adaptors, $cfg, $control_experiments_ids, \%epi_Ihec);
    }




	print "FIN\n";
    return 1;
}


sub empty_data {
    my $cfg = shift;

    my $connection = DBI->connect("DBI:mysql:".$cfg->{efg_db}->{dbname}.":".$cfg->{efg_db}->{host}.":".$cfg->{efg_db}->{port}, $cfg->{efg_db}->{user}, $cfg->{efg_db}->{pass});
    
    my $sql = "DELETE FROM experiment";
    my $sth = $connection->prepare($sql);
    $sth->execute();

    $sql = "DELETE FROM read_file_experimental_configuration";
    $sth = $connection->prepare($sql);
    $sth->execute();

    $sql = "DELETE FROM read_file";
    $sth = $connection->prepare($sql);
    $sth->execute();

    $connection->disconnect();
    return 1;

}


sub get_data_from_csv {
    my ( $csv, $logger ) = @_;
    my ( %control_data, %signal_data );
	
    open my $csv_fh, '<', $csv;

    while ( readline $csv_fh ) {
        chomp;

        if (/^epi_accession/i) {
            next;    # ignore input file header
        }




        my ($epi_accession,         $accession,             $exp_accession,			$epigenome_name,    	
            $feature_type_name,     $br,                    $new_br,				$tr,					
            $new_tr,     			$gender,	            $md5,                   $local_url,         	
            $analysis_name,         $exp_group_name,        $assay_xrefs,			$ontology_xref_accs, 	
            $xref_accs, 			$epigenome_description,	$controlled_by, 		$paired,				
            $paired_end_tag,		$read_length,			$multiple,				$paired_with,           
            $download_url,          $info,                  $derived
        ) = split /\t/;

		
		#Columns $exp_accession, $br and $tr are ignored



        my $entry = {};
        $entry->{epi_accession}         = $epi_accession;
        $entry->{accession}             = $accession;
        $entry->{ex_accession}          = $exp_accession;
        $entry->{epigenome_name}        = create_epigenome_production_name( $epigenome_name );
        $entry->{feature_type_name}     = $feature_type_name;
        $entry->{br}                    = $new_br;
        $entry->{tr}                    = $new_tr;
		$entry->{gender}                = $gender;
		$entry->{md5}                   = $md5;
        $entry->{local_url}             = $local_url;
        $entry->{analysis_name}         = $analysis_name;
        $entry->{exp_group_name}        = $exp_group_name;
        $entry->{assay_xrefs}        	= $assay_xrefs;
		$entry->{ontology_xref_accs}    = $ontology_xref_accs;
        $entry->{xref_accs}             = $xref_accs;
        $entry->{epigenome_description} = $epigenome_description;
        $entry->{controlled_by}         = $controlled_by;
        $entry->{paired}         		= $paired;
		$entry->{paired_end_tag}        = $paired_end_tag;
		$entry->{read_length}         	= $read_length;
		$entry->{multiple}         		= $multiple;
		$entry->{download_url}          = $download_url;
        $entry->{info}                  = $info;
        $entry->{derived}               = $derived;
		

        if (uc $entry->{exp_group_name} eq 'ROADMAP'){
            $entry->{exp_group_name} = 'Roadmap Epigenomics';
        }
        


		if (uc $entry->{paired} eq 'YES'){
			$entry->{paired} = 1;
		}else{
			$entry->{paired} = undef;
		}

        

        if ( $control_data{$accession} || $signal_data{$accession} ) {
            $logger->error( 'Accession ' . $accession . ' is NOT unique!',
                0, 0 );
        }
		

		
        if ( $feature_type_name eq 'WCE' ) {
            $entry->{is_control} = 1;
            $control_data{$accession} = $entry;
        }
        else {
            $entry->{is_control} = 0;
            $signal_data{$accession} = $entry;
        }

        $entry->{derived} =~ s/^\s+|\s+$//g; #remove espaces
        $entry->{paired_end_tag} =~ s/^\s+|\s+$//g; #remove espaces

        if ($entry->{derived} ne '-' && $entry->{derived} ne ''){
            my($filename, $dirs, $suffix) = fileparse($entry->{local_url});
            my $file_new_name = $entry->{derived};
            if ( $entry->{paired_end_tag} ne '-' && $entry->{paired_end_tag} ne '' ){
                $file_new_name .= '_'.$entry->{paired_end_tag};
            }
            $file_new_name .= '.fastq.gz';
            $entry->{local_url} = $dirs.$file_new_name;
        }



		$entry = verify_entry_metadata( $entry, $logger );

    }

    close $csv_fh;
	
    return ( \%control_data, \%signal_data );
}

sub verify_entry_metadata {
    my ( $entry, $logger ) = @_;

    my @mandatory = qw(
        accession         epigenome_name
        feature_type_name br
        tr                md5
        local_url         analysis_name
        exp_group_name    multiple
    );
	
	

    for my $man (@mandatory) {
        if ( !$entry->{$man} ) {
            $logger->error(
                'There is no ' . $man . ' value for ' . $entry->{accession},
                0, 0 );
        }
    }

    my @optional = qw(
        gender        			assay_xrefs
		ontology_xref_accs		xref_accs     
		epigenome_description	controlled_by 
		download_url	        info
		read_length
    );

    for my $opt (@optional) {
        if ( $entry->{$opt} eq '-' ) {
            $entry->{$opt} = undef;
        }
    }

	#if paired is true, paired_end_tag must be informed
	if ($entry->{paired}){
		if (!$entry->{paired_end_tag} || $entry->{paired_end_tag} eq '-'){
            $logger->error(
                'There is no paired_end_tag value for ' . $entry->{accession},
                0, 0 );
		}
	}else{
		$entry->{paired_end_tag} = undef;
	}



    return $entry;
}

# avoid using a register.conf file, use the existing mechanism
#
# use registry instead
# list and expose TrackingAdaptor in DBAdaptor
sub fetch_adaptors {
    my ($cfg) = @_;
    my %adaptors;

    my $dba = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new(
        -user    => $cfg->{efg_db}->{user},
        -pass    => $cfg->{efg_db}->{pass},
        -host    => $cfg->{efg_db}->{host},
        -port    => $cfg->{efg_db}->{port},
        -dbname  => $cfg->{efg_db}->{dbname},
        -species => $cfg->{general}->{species},
    );


    my $dbCoreAdaptor = Bio::EnsEMBL::DBSQL::DBAdaptor->new(
        -user    => $cfg->{dna_db}->{user},
        -host    => $cfg->{dna_db}->{host},
        -port    => $cfg->{dna_db}->{port},
        -dbname  => $cfg->{dna_db}->{dbname}
    );

	$adaptors{epigenome}    = $dba->get_adaptor('Epigenome');
	$adaptors{feature_type}    = $dba->get_adaptor('FeatureType');
	$adaptors{analysis}    = $dba->get_adaptor('Analysis');
	$adaptors{exp_group}    = $dba->get_adaptor('ExperimentalGroup');
	$adaptors{experiment}    = $dba->get_adaptor('Experiment');
	$adaptors{read_file_experimental_configuration} = $dba->get_adaptor('ReadFileExperimentalConfiguration');
	$adaptors{read_file}    = $dba->get_adaptor('ReadFile');
		
    $adaptors{db}       = $dba;
    $adaptors{db_entry} = Bio::EnsEMBL::DBSQL::DBEntryAdaptor->new($dba);

    $adaptors{Gene} = $dbCoreAdaptor->get_adaptor("Gene");

    return \%adaptors;
}

sub verify_basic_objects_in_db {
    my ( $control_data, $signal_data, $adaptors, $logger ) = @_;
    my $abort = 0;
    my %to_register;

    for my $entry ( values %{$control_data}, values %{$signal_data} ) {


        my @lstAnalysis = split (',', $entry->{analysis_name});

        my $analysis;

        foreach my $analysis_name (@lstAnalysis){
            $analysis_name =~ s/\s//g; #remove blank spaces
            
            my $analysis = $adaptors->{analysis}->fetch_by_logic_name( $analysis_name );
            if ( !defined $analysis ) {
                $to_register{Analysis} //= [];
                push @{ $to_register{Analysis} }, $entry->{analysis_name};
                $abort = 1;

            }else{
                my $feature_type = fetch_feature_type($entry, $adaptors, $analysis);
                if (!defined $feature_type){
                    $to_register{Feature_Type} //= [];
                    push @{ $to_register{Feature_Type} }, $entry->{feature_type_name};
                    #$abort = 1;
                }
            }
        }


			
        my $exp_group = $adaptors->{exp_group}->fetch_by_name( $entry->{exp_group_name} );

        if ( !defined $exp_group ) {
            $to_register{Experimental_Group} //= [];
            push @{ $to_register{Experimental_Group} },
                $entry->{exp_group_name};
            $abort = 1;

            
        }else{
        	$entry->{exp_group_name} = $exp_group->name;
        }
    }
	
	##########################################
	## extra checks must be implemented here ##
	##########################################
	
	
	
    if ($abort) {
        for my $object ( keys %to_register ) {
            for my $missing ( uniq @{ $to_register{$object} } ) {
                $logger->warning(
                    'Register ' . $object . ': ' . $missing . "\n",
                    0, 1 );
            }
        }

        $logger->error( 'Aborting registration' . "\n", 0, 1 );
    }


    return 1;
}

sub get_epigenome_unique_name{
    my $epigenome_name = shift;
    my $adaptors = shift;

    my $epi_counter = 1; 
    my $epi_name;
    do{
        $epi_name = "new_".$epigenome_name.'_'.$epi_counter;
        $epi_counter++;

    }while($adaptors->{epigenome}->fetch_by_name( $epi_name ));   

    return $epi_name;
}

sub register_ctl {
	my ( $logger, $ctl_data, $adaptors, $cfg, $epi_name_ihec) = @_;
	
	my $registered_file = {};
	my $control_experiments_ids = {};

    my $external_exp = {};
    my $file_exp = {};
	my %control_data = %{$ctl_data};
    my $experiment;
    #my %epi_name_ihec = %{$epi_ihec};
   
    foreach my $ctl (keys %control_data){
        
        my $entry = $control_data{$ctl};

        my @lstAnalysis = split (',', $entry->{analysis_name});

        my $epi_new_name;

        foreach my $analysis_name (@lstAnalysis){

            $analysis_name =~ s/\s//g; #remove blank spaces
            my $analysis = $adaptors->{analysis}->fetch_by_logic_name( $analysis_name );

            #my $epigenome = $adaptors->{epigenome}->fetch_by_name( $entry->{epigenome_name} );
            my $epigenome;
            my $old_epigenome_name = $entry->{epigenome_name}; 
            $old_epigenome_name =~ s/\:/_/g;
            $old_epigenome_name =~ s/\+//g;
            $old_epigenome_name =~ s/\(//g;
            $old_epigenome_name =~ s/\)//g;
            $old_epigenome_name =~ s/\-/_/g;
            $old_epigenome_name =~ s/\./_/g;
            $old_epigenome_name =~ s/\//_/g;
            $old_epigenome_name =~ s/ /_/g;

            my $ihec_epiname = $entry->{xref_accs}."_".$old_epigenome_name;

            if (!(exists ($epi_name_ihec->{$ihec_epiname})) && !($epi_new_name)){
               
                $epi_new_name = get_epigenome_unique_name($entry->{epigenome_name}, $adaptors );
               
                $entry->{epigenome_name}=$epi_new_name;

                $epigenome = store_epigenome( $entry, $adaptors );

                $epi_name_ihec->{$ihec_epiname} = $epi_new_name;


            }else{
                if ($epi_new_name){
                    $epigenome = $adaptors->{epigenome}->fetch_by_name( $epi_new_name );
                }else{
                    $epigenome = $adaptors->{epigenome}->fetch_by_name( $epi_name_ihec->{$ihec_epiname} );
                }
                
            }

            my $subFix='';
            if($analysis_name eq 'ChIP-Seq_control_hist'){
                $subFix = '_His';
            }else{
                if ($analysis_name eq 'ChIP-Seq_control_TF'){
                    $subFix = '_TF';
                }
            }

            my $ex_accession = $entry->{ex_accession}.$subFix;

            if (!$external_exp->{$ex_accession}){
                my $exp_group = $adaptors->{exp_group}->fetch_by_name( $entry->{exp_group_name} );
                my $feature_type = fetch_feature_type( $entry, $adaptors, $analysis );
                my $experiment_name = create_control_experiment_name( $entry, $cfg, $epigenome, $adaptors );
                $experiment = store_experiment($entry, $experiment_name, $adaptors, $epigenome, $feature_type, $exp_group, '');
                $external_exp->{$ex_accession} = $experiment->dbID();


            }else{
                $experiment = $adaptors->{experiment}->fetch_by_dbID( $external_exp->{$ex_accession} );
            }



            store_db_xref( $entry, $adaptors, $epigenome );
            store_db_xref_encode_epigenome($entry, $adaptors, $epigenome);
            my $reg_file = store_read_file_only( $logger, $entry, $analysis, $adaptors );
            store_read_file_experimental_configuration ($reg_file, $entry, $adaptors, $experiment);
            store_db_xref_read_file($entry, $adaptors, $reg_file);
            my $accession = $entry->{accession}.$subFix;
            $file_exp->{$accession}=$external_exp->{$ex_accession};


        }

    }

    $logger->info( "Ctl Successful Registration\n", 1, 1 );

    return $file_exp;

}

sub register_signal {
    my ( $logger, $entry, $adaptors, $cfg, $control_experiments_ids, $epi_name_ihec) = @_;
    #my %epi_name_ihec = %{$epi_ihec};

    $logger->info( 'Registering ' . $entry->{accession} . "\n", 0, 1 );

    #my $analysis = $entry->{analysis_name};
	my $analysis = $adaptors->{analysis}->fetch_by_logic_name( $entry->{analysis_name} );

    #my $epigenome = $adaptors->{epigenome}->fetch_by_name( $entry->{epigenome_name} );
    my $epigenome;
    my $old_epigenome_name = $entry->{epigenome_name}; 
    $old_epigenome_name =~ s/\:/_/g;
    $old_epigenome_name =~ s/\+//g;
    $old_epigenome_name =~ s/\(//g;
    $old_epigenome_name =~ s/\)//g;
    $old_epigenome_name =~ s/\-/_/g;
    $old_epigenome_name =~ s/\./_/g;
    $old_epigenome_name =~ s/\//_/g;
    $old_epigenome_name =~ s/ /_/g;

    my $ihec_epiname = $entry->{xref_accs}."_".$old_epigenome_name;

    

    if (!(exists ($epi_name_ihec->{$ihec_epiname}))){
        my $epi_new_name = get_epigenome_unique_name($entry->{epigenome_name}, $adaptors );
        $entry->{epigenome_name}=$epi_new_name;
        $epigenome = store_epigenome( $entry, $adaptors );
        $epi_name_ihec->{$ihec_epiname} = $entry->{epigenome_name};
    }else{

        $epigenome = $adaptors->{epigenome}->fetch_by_name( $epi_name_ihec->{$ihec_epiname} );
    }





    #if ( $entry->{ontology_xref_accs} ) {
    #    store_ontology_xref( $entry, $adaptors, $epigenome );
    #}

    if ( $entry->{xref_accs} ) {
        store_db_xref( $entry, $adaptors, $epigenome );
    }



    my $feature_type = fetch_feature_type( $entry, $adaptors, $analysis );
    
    if (! defined $feature_type) {
        

        die("Couldn't fetch feature type with name " . $entry->{feature_type_name} . "!");
    }

	
    my $exp_group
        = $adaptors->{exp_group}->fetch_by_name( $entry->{exp_group_name} );


    my $experiment_name = create_signal_experiment_name( $entry, $cfg, $epigenome, $adaptors );
	
   
    my $experiment = $adaptors->{experiment}->fetch_by_name($experiment_name);
	

    # don't use control_db_ids, fetch directly from db for every signal file
    my $subFix;
    if ( !$experiment ) {

		my $controlled = undef;
		if ($entry->{controlled_by}){
			if ($entry->{controlled_by} ne '-'){
                my @lstctl_by = split (',', $entry->{controlled_by});
                my $ctl_by = $lstctl_by[0];
                $ctl_by =~ s/\s//g; #remove blank spaces
                
                my $ftClass = $feature_type->class;
                if (uc $ftClass eq 'HISTONE'){
                    $subFix = '_His';
                }else{
                   #if (uc $ftClass eq 'TRANSCRIPTION FACTOR'){
                        $subFix = '_TF';
                    #}
                }

               

                my $controlledSub = $ctl_by.$subFix;
 

				my $controlled_id = $control_experiments_ids->{$controlledSub};

                #print "Controlled_id: ".$controlled_id."\n";
				$controlled = $adaptors->{experiment}->fetch_by_dbID($controlled_id);
			}
		}
        $experiment = store_experiment($entry, $experiment_name, $adaptors, $epigenome, $feature_type,  $exp_group, $controlled );

    }


    my $reg_file = store_read_file_only( $logger, $entry, $analysis, $adaptors );
    store_read_file_experimental_configuration ($reg_file, $entry, $adaptors, $experiment);
    store_db_xref_read_file($entry, $adaptors, $reg_file);
	

    $logger->info( "Successful Registration\n", 1, 1 );

    return 1;
}

sub store_epigenome {
    my ( $entry, $adaptors ) = @_;

    my $production_name
        = create_epigenome_production_name( $entry->{epigenome_name} );



    my @epicomp = split (';', $entry->{epigenome_description});;
    my $epidesc = @epicomp[0];
    #my $epilabel = @epicomp [1];
    #my $epilabel = $epidesc;
    $epidesc =~ s/^\s+//;
    $epidesc =~ s/\s+$//;
    #$epilabel =~ s/^\s+//;
    #$epilabel =~ s/\s+$//;

    my $epilabel = $entry->{epigenome_name};

     my $epigenome = Bio::EnsEMBL::Funcgen::Epigenome->new(
        -name            => $entry->{epigenome_name},
        -short_name      => $epilabel,
        -description     => $epidesc,
        -production_name => $production_name,
        -gender          => $entry->{gender}
    );

    $adaptors->{epigenome}->store($epigenome);

    return $epigenome;
}

sub fetch_feature_type {
    my ( $entry, $adaptors, $analysis ) = @_;

    my $ft_name = $entry->{feature_type_name};

    my $feature_type;

    if ( $entry->{is_control} ) {
        $feature_type = $adaptors->{feature_type}->fetch_by_name( $ft_name, 'DNA', $analysis );
        #print "Analysis control: ".$analysis."\n";
        #print "Feature control: ".$ft_name."\n";
    }
    else {
        $feature_type = $adaptors->{feature_type}->fetch_by_name($ft_name);
            #print "Feature: ".$ft_name."\n";
    }

    if (!$feature_type){
        #my $is_TF =  check_TF($adaptors, $ft_name);

        #if ($is_TF == 1){
            my $class = 'Transcription Factor';
            my $description = $ft_name. ' TF binding';
            my $so_accession = 'SO:0000235';
            my $so_name = 'TF_binding_site';
            $feature_type = create_feature_type($adaptors, $ft_name, $class, $description, $so_accession, $so_name);
        #}
        
    }

    return $feature_type;
}

sub check_TF {
    my $adaptors = shift;
    my $ft_name = shift;

    my $is_TF = 0;
    
    my $lstgene = $adaptors->{Gene}->fetch_all_by_external_name($ft_name);

    LBL_GENE: {
        foreach my $gene (@{$lstgene}){
            my $lstDbEntries = $gene->get_all_DBLinks('GO');
            foreach my $dbentry (@{$lstDbEntries}){
                if (index(uc $dbentry->description(), 'DNA BINDING TRANSCRIPTION FACTOR ACTIVITY') != -1) {
                    $is_TF = 1;
                    last LBL_GENE;
                } 
            }
    
        }
    }

    return $is_TF;
}

sub create_feature_type {
    my ($adaptors, $ft_name, $class, $description, $so_accession, $so_name) = @_;

    my $feature = Bio::EnsEMBL::Funcgen::FeatureType->new(
        -name               => $ft_name,
        -class              => $class,
        -description        => $description,
        -so_name            => $so_name,
        -so_accession       => $so_accession
        );

    $adaptors->{feature_type}->store($feature);

    $feature= $adaptors->{feature_type}->fetch_by_name($ft_name);


    return $feature;
}

sub create_control_experiment_name {
    my ( $entry, $cfg, $epigenome, $adaptors ) = @_;

    my $experiment_name;
    my $number = 1;

    do {
        $experiment_name
            = $epigenome->{production_name} . '_'
            . $entry->{feature_type_name} . '_'
            . $entry->{analysis_name} . '_no'
            . $number . '_'
            . $entry->{exp_group_name};

        if ( $cfg->{general}->{release} ) {
            $experiment_name .= $cfg->{general}->{release};
        }

        $experiment_name =~ s/\s//g;
        $number++;
    } while ( $adaptors->{experiment}->fetch_by_name($experiment_name) );

    return $experiment_name;
}

sub create_signal_experiment_name {
    my ( $entry, $cfg, $epigenome, $adaptors ) = @_;

    my $experiment_name;

    $experiment_name
        = $epigenome->{production_name} . '_'
        . $entry->{feature_type_name} . '_'
        . $entry->{analysis_name} . '_'
        . $entry->{exp_group_name};

    if ( $cfg->{general}->{release} ) {
        $experiment_name .= $cfg->{general}->{release};
    }

    $experiment_name =~ s/\s//g;

    return $experiment_name;
}

sub store_experiment {
    my ($entry, $experiment_name, $adaptors, $epigenome, $feature_type, $exp_group, $control_experiment ) = @_;


    #my $control_experiment;
    #if ( !$entry->{is_control} && $entry->{controlled_by} ) {

    #    my $control_db_id = $control_db_ids->{ $entry->{controlled_by} };

    #    $control_experiment
    #        = $adaptors->{experiment}->fetch_by_dbID($control_db_id);
    #}

	if ($entry->{is_control}){
		$control_experiment=undef;
	}

	
    my $experiment = Bio::EnsEMBL::Funcgen::Experiment->new(
        -NAME               => $experiment_name,
        -EPIGENOME          => $epigenome,
        -FEATURE_TYPE       => $feature_type,
        -EXPERIMENTAL_GROUP => $exp_group,
        -IS_CONTROL         => $entry->{is_control},
        -CONTROL            => $control_experiment,
    );




    $adaptors->{experiment}->store($experiment);

    return $experiment;
}


sub store_read_file_only {
	my ( $logger, $entry, $analysis, $adaptors )= @_;

    my $read_file = $adaptors->{read_file}->fetch_by_name( $entry->{accession} );

    if ($read_file) {
        $logger->warning(
            'A read file entry for accession '
                . $entry->{accession}
                . ' already exists in DB! ',
            0, 1
        );
    }
	
	if (!$entry->{info}){
		$entry->{info} = undef;
	}
	
	my $is_paired = undef;
	if ($entry->{paired}){
		$is_paired=1;
	}


    my $read_file_name = $entry->{accession};
    if ($entry->{derived} ne '-' && $entry->{derived} ne ''){
        $read_file_name = $entry->{derived};
        if ($is_paired == 1){
            $read_file_name .='_'.$entry->{paired_end_tag};
        }
    }

    my $read_file = Bio::EnsEMBL::Funcgen::ReadFile->new(
        -name           => $read_file_name,
        -analysis       => $analysis,
        -is_paired_end  => $is_paired,
        -file_size      => undef,
        -read_length    => $entry->{read_length},
        -md5sum         => $entry->{md5},
        -file           => $entry->{local_url},
        -notes          => undef,
    );
	$adaptors->{read_file}->store($read_file);
	
	return $read_file;
}

sub store_read_file_experimental_configuration {
	 my ( $read_file, $entry, $adaptors, $experiment ) = @_;
	 
     my $paired_end_tag = $entry->{paired_end_tag};
	 
     if (!$paired_end_tag || $paired_end_tag eq '-') {
       $paired_end_tag = undef;
     }
	 
     
	 my $read_file_experimental_configuration = Bio::EnsEMBL::Funcgen::ReadFileExperimentalConfiguration->new(
     -read_file_id => $read_file->dbID,
     -experiment_id         => $experiment->dbID,
     -biological_replicate  => $entry->{br},
     -technical_replicate   => $entry->{tr},
     -paired_end_tag        => $paired_end_tag,
     -multiple              => $entry->{multiple},
   );
    
   
   
   $adaptors->{read_file_experimental_configuration}->store($read_file_experimental_configuration);

 
   return 1;
}



=pod
sub store_ontology_xref {
    my ( $entry, $adaptors, $epigenome_obj ) = @_;

    my @ontology_accessions = split /;/, $entry->{ontology_xref_accs};
    my %valid_linkage_annotations
        = ( 'SAMPLE' => 1, 'TISSUE' => 1, 'CONDITION' => 1 );

    for my $ontology_accession (@ontology_accessions) {

        my ( $linkage_annotation, $primary_id ) = split /-/,
            $ontology_accession;
        my ($dbname) = split /:/, $primary_id;

        if ( !$valid_linkage_annotations{$linkage_annotation} ) {
            throw
                'Invalid linkage_annotation, please use \'SAMPLE\', \'TISSUE\' or \'CONDITION\'';
        }

        my $ontology_xref = Bio::EnsEMBL::OntologyXref->new(
            -primary_id         => $primary_id,
            -dbname             => $dbname,
            -linkage_annotation => $linkage_annotation,
        );

        my $ignore_release = 1;

        my $epigenome_id = $epigenome_obj->dbID();

        $adaptors->{db_entry}
            ->store( $ontology_xref, $epigenome_id, 'epigenome',
            $ignore_release );
    }

    return 1;
}
=cut

sub store_db_xref {
    my ( $entry, $adaptors, $epigenome_obj ) = @_;

    my @xref_accessions = split /;/, $entry->{xref_accs};

    for my $xref_acc (@xref_accessions) {
        my ( $dbname, $primary_id ) = split /-/, $xref_acc;
        $dbname = "EpiRR";

        my $xref = Bio::EnsEMBL::DBEntry->new(
            -primary_id => $primary_id,
            -dbname     => $dbname,
            -display_id => undef, # if not set to undef, an empty string will be stored in the display_label column
            -info_text  => undef # if not set to undef, an empty string will be stored in the info_text column
        );

       #my $ignore_release = 1;
       my $ignore_release = 0;

        my $epigenome_id = $epigenome_obj->dbID();



        $adaptors->{db_entry}->store( $xref, $epigenome_id, 'epigenome', $ignore_release );

        

    }

    return 1;
}

sub store_db_xref_read_file {
    my ( $entry, $adaptors, $read_file_obj ) = @_;

    my $dbname = "ENCODE";

    my $xref = Bio::EnsEMBL::DBEntry->new(
        -primary_id => $entry->{accession},
        -dbname     => $dbname,
        -display_id => undef, # if not set to undef, an empty string will be stored in the display_label column
        -info_text  => undef # if not set to undef, an empty string will be stored in the info_text column
    );

    #my $ignore_release = 1;
    my $ignore_release = 0;

    my $read_file_id = $read_file_obj->dbID();

    $adaptors->{db_entry}->store( $xref, $read_file_id, 'ReadFile', $ignore_release );

    return 1;
}

sub store_db_xref_encode_epigenome {
    my ( $entry, $adaptors, $epigenome_obj ) = @_;

    my $dbname = "ENCODE";

    my $xref = Bio::EnsEMBL::DBEntry->new(
        -primary_id => $entry->{epi_accession},
        -dbname     => $dbname,
        -display_id => undef, # if not set to undef, an empty string will be stored in the display_label column
        -info_text  => undef # if not set to undef, an empty string will be stored in the info_text column
    );

    #my $ignore_release = 1;
    my $ignore_release = 0;

    my $epigenome_id = $epigenome_obj->dbID();

    $adaptors->{db_entry}->store( $xref, $epigenome_id, 'epigenome', $ignore_release );

    return 1;
}

sub create_epigenome_production_name {
    my ($epigenome_name) = shift;

    $epigenome_name =~ s/\:/_/g;
    $epigenome_name =~ s/\+//g;
    $epigenome_name =~ s/\(//g;
    $epigenome_name =~ s/\)//g;
    $epigenome_name =~ s/\-/_/g;
    $epigenome_name =~ s/\./_/g;
    $epigenome_name =~ s/\//_/g;
    $epigenome_name =~ s/ /_/g;

    return $epigenome_name;
}

sub usage {
    my $usage = << 'END_USAGE';

Usage: register_metadata.pl -i <input_file> -c <config_file>

Options:
-i input_file:  this is the tab delimited text file that contains the metadata
-c config_file: this is the configuration file that contains the database connection details
-n dry_run:     run the script without storing any data into the database, NOT YET IMPLEMENTED
-h help:        shows this help message
 
END_USAGE

    say $usage;

    return 1;
}
