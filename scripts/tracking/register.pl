#!/usr/bin/perl

use strict;
use warnings;
use Data::Dumper;
use Carp;

use JSON qw(decode_json);
use JSON::Path;
use HTTP::Tiny;

use Pod::Usage;
use Getopt::Long;
use Config::Tiny;
use feature qw(say);

use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( dump_data get_date);
use Bio::EnsEMBL::Utils::Exception qw(throw);
use Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor;
use Bio::EnsEMBL::Funcgen::CellType;
use Bio::EnsEMBL::Funcgen::FeatureType;
use Bio::EnsEMBL::Funcgen::Experiment;
use Bio::EnsEMBL::Funcgen::ExperimentalGroup;
use Bio::EnsEMBL::Funcgen::InputSubset;

use Bio::EnsEMBL::Utils::SqlHelper;


use constant CONFIG  => 'register.config.ini';
select((select(STDOUT), $|=1)[0]);

my $http = HTTP::Tiny->new();
my $server = 'https://www.encodeproject.org';
my $global_headers = { 'Content-Type' => 'application/json' };
my $last_request_time = Time::HiRes::time();
my $request_count = 0;

main();

sub main {
  my $cfg = Config::Tiny->new;
     $cfg = Config::Tiny->read(CONFIG);

  _connect_to_trackingDB($cfg);
  _get_trackingDB_adaptors($cfg);

  my $headers = {};
  my $ids = [];
     $ids = ['ENCSR000CEH', 'ENCSR000CDB','ENCSR000CBK'];
     # $ids = parse_flatfile($cfg->{file}->{input});
  say "Found " .scalar(@{$ids}) . ' Experiments to register';



  my $skip = {
#    ENCSR000ADN => 1,
#    ENCSR000COE => 1,
#    ENCSR000CNE => 1,
#    ENCSR000CAW => 1,
#    ENCSR000CNG => 1,
#    ENCSR000CNI => 1,
  };

  # $objects contains Ensembl API objects
  my $objects = {};
  $objects->{exp_group} = $cfg->{tr_a}->{eg}->fetch_by_name('ENCODE');
###############################################################################
# Each ID is one ENCODE experiment. An ENCODE experiment is different to our
# (Regulation) definition of an experiment, as it can contain different
# biological replicates.
# Each experiment has the same ExperimentalGroup (ENCODE).
#
# Each file in each experiment should have the same Analysis, CellType and
# FeatureType.
#
#

  for my $id(@{$ids}){
    next if (exists $skip->{$id});
    # next unless ($id eq 'ENCSR000ESH');
    say $id;

    my $type = 'experiments';
    my $url = "$cfg->{url}->{base}/$type/$id/?format=json";

    $headers->{'Content-Type'} = 'application/json';
    my $json = fetch_and_decode_json($url, $headers);

    # $data contains data parsed from JSON or data derived from JSON
    my $data = {};
    $data->{experiment_id} = $id;
    $data->{experiment_url} = $url;
    parse_ENCODE_json($cfg, $json, $data);


    compare_biosample_set_in_root($data);

    for my $file (@{$data->{files} }){
      fetch_analyis($cfg, $data, $objects);

      create_cell_type_name($cfg, $data, $file);
      create_experiment_name($data, $file);
      my $exp_name =  $data->{experiment_name};

      fetch_cell_type($cfg, $data, $objects);
      fetch_feature_type($cfg, $data, $objects, $json);
      fetch_experiment($cfg, $data, $objects, $exp_name);
      fetch_input_subset($cfg, $data, $objects, $file);

    }




    # say $data->{life_stage};
    # say dump_data($data,1,1);

  }
}
sub print_objects {
  my ($obj) = @_;

  foreach my $key (sort keys %{$obj}){
    say "$key: " . ref($obj->{$key});
  }
}

###############################################################################
#                              fetch_analyis
###############################################################################
sub fetch_analyis {
  my ($cfg, $data, $objects) = @_;

  my $anal = $cfg->{tr_a}->{an}->fetch_by_logic_name($data->{analysis});

  if(! defined $anal){
    throw "Register Analysis: '" . $data->{analysis}."'";
  }

  $objects->{analysis} = $anal;

  return;
}

###############################################################################
#                           fetch_cell_feature_type
###############################################################################
sub fetch_cell_type {
  my ($cfg, $data, $objects) = @_;

  my $ct_name = $data->{cell_type_name};
  my $cell_type = $cfg->{tr_a}->{ct}->fetch_by_name($ct_name);

  if(!defined $cell_type) {
    $cell_type = Bio::EnsEMBL::Funcgen::CellType->new(
     -name          => $data->{cell_type_name},
     -display_label => $data->{cell_type_name},
     -description   => $data->{ct_desc},
     -gender        => $data->{sex},
     -tissue        => $data->{cell_type},
     );
    $cfg->{tr_a}->{ct}->store($cell_type);
  }

  $objects->{cell_type} = $cell_type;

  return;
}


###############################################################################
#                           fetch_feature_type
###############################################################################
sub fetch_feature_type {
  my ($cfg, $data, $objects, $json) = @_;

  my $ft_name = $data->{feature_type};
  my $feature_type = $cfg->{tr_a}->{ft}->fetch_by_name($ft_name);

  if(!defined $feature_type) {
    if($json->{target}->{investigated_as}->[0] eq 'transcription factor') {
      $feature_type = Bio::EnsEMBL::Funcgen::FeatureType->new(
        -name         => $ft_name,
        -class        => 'Transcription Factor',
        -description  => $ft_name . ' Transcription Factor Binding',
        -so_accession => 'SO:0000235',
        -so_name      => 'TF_binding_site',
        );
      $cfg->{tr_a}->{ft}->store($feature_type);
    }
    else {
      my $id = $data->{experiment_id};
      my $url = $data->{experiment_url};
      throw "Register FeatureType: '$ft_name' URL: $id $url";
    }
  }

  $objects->{feature_type} = $feature_type;

  return;
}



###############################################################################
#                          create_and_store_input_subset
###############################################################################
sub fetch_input_subset {
  my ($cfg, $data, $objects, $file) = @_;

  my $iss = $cfg->{tr_a}->{iss}->fetch_by_name($file->{name});

  if(!defined $iss){
    my $is_control = ($data->{feature_type} eq 'WCE') ? 1 : 0;

    $iss = Bio::EnsEMBL::Funcgen::InputSubset->new (
         -analysis      => $objects->{analysis},
         -cell_type     => $objects->{cell_type},
         -experiment    => $objects->{experiment},
         -feature_type  => $objects->{feature_type},
         -is_control    => $is_control,
         -name          => $file->{name},
         -replicate     => $file->{replicate},
      );
    $cfg->{tr_a}->{iss}->store($iss);

    my $tr_info->{info} = {
         availability_date => $data->{availability_date},
         download_url      => $file->{download_url},
         download_date     => undef,
         local_url         => undef,
         md5sum            => $file->{md5sum},
         notes             => 'Experiment:'.$data->{experiment_id},
       };
      $cfg->{tr_a}->{tr}->store_tracking_info($iss, $tr_info);

   }

  push(@{$objects->{input_subsets}}, $iss);

  return;
}

###############################################################################
#                          create_and_store_experiment
###############################################################################
# Could cache experiment, but again, readability vs efficiency

sub fetch_experiment {
  my ($cfg, $data, $objects, $exp_name) = @_;

  my $experiment = $cfg->{tr_a}->{ex}->fetch_by_name($exp_name);
  if(! defined $experiment) {

    my $exp_name = $data->{experiment_name};
    $experiment = Bio::EnsEMBL::Funcgen::Experiment->new (
          -NAME                => $data->{experiment_name},
          -CELL_TYPE           => $objects->{cell_type},
          -FEATURE_TYPE        => $objects->{feature_type},
          -EXPERIMENTAL_GROUP  => $objects->{exp_group},
          # -DATE                => get_date('mysql'),
          -DESCRIPTION         => $data->{description},
       );

    $cfg->{tr_a}->{ex}->store($experiment);

    my $tr_info->{info} = {
     # experiment_id     => $experiment->dbID,
     notes             => 'Experiment:'.$data->{experiment_id},
    };
   $cfg->{tr_a}->{tr}->store_tracking_info($experiment, $tr_info);

  }
  $objects->{experiment} = $experiment;
}


###############################################################################
#                          create_cell_type_name
###############################################################################

sub create_cell_type_name {
  my ($cfg, $data, $file) = @_;

  my $name;
  $name .= $data->{cell_type}.':';
  $name .= $cfg->{life_stage}->{$data->{life_stage}};

  if( (defined $data->{age}) and (defined $data->{age_units})){
    $name .= $data->{age};
    $name .= $cfg->{age_units}->{$data->{age_units}};
  }

  $data->{cell_type_name}  = $name;
  return $name;

}

###############################################################################
#                          create_experiment_name
###############################################################################

sub create_experiment_name {
  my ($data, $file) = @_;
  # say dump_data($data,1,1);die;

  my $name;
  $name .= '_'         . $data->{feature_type};
  $name .= '_ENCODE_'  . $data->{lab};
  $name .= '_BR'       . $file->{bio_replicate};
  $data->{experiment_name} = $data->{cell_type_name} . $name;
  return $name;

}

###############################################################################
#                          compare_biosample_set_in_root
###############################################################################
# Some data which can be  unique in each ENCODE file must be unique for an
# Regulation experiment. This method compares elements we need to be identical
# across files. It also compares itself to itself, but as this is not time
# critical, readability has been choosen over efficiency

sub compare_biosample_set_in_root {
  my ($data) = @_;

  my $file_1 = $data->{files}->[0];
  my $acc    = $data->{accession};
  my $uuid_1 = $file_1->{uuid};

  for my $file_2 (@{$data->{files}}){
    my $uuid_2 = $file_2->{uuid};

    my @fields = qw(age age_units description life_stage sex);
    for my $f(@fields){
      if(defined $file_1->{$f}){
        if(! defined $file_2){
          croak "Missing '$f'. $acc: Files: $uuid_1 $uuid_2";
        }
        if($file_1->{$f} ne $file_2->{$f}){
          croak "Difference '$f'. $acc: Files: $uuid_1 $uuid_2";
        }
      }
    }
  }

  $data->{age}            = $file_1->{age};
  $data->{age_units}      = $file_1->{age_units};
  $data->{bio_replicate}  = $file_1->{bio_replicate};
  $data->{cell_type_desc} = $file_1->{description};
  $data->{life_stage}     = $file_1->{life_stage};
  $data->{sex}            = _get_sex($file_1->{sex});


}

###############################################################################
#                                 _get_sex
###############################################################################

sub _get_sex {
 my ($sex) = @_;

 my $sex_result;

  if($sex eq 'unknown'){
    $sex_result = undef;
  }
  elsif($sex =~ m!hermaphrodite|female|male|mixed!o){
   $sex_result  = $sex;
  }
  else{
    throw("Unkown type of sex: $sex");
  }
  return $sex_result;
}
###############################################################################
#                           assign_feature_type
###############################################################################
sub assign_feature_type {
  my ($json, $data, $cfg) = @_;

  my $ftype = undef;
  if(defined $json->{target}->{label} ){
    $ftype = $json->{target}->{label};
  }
  elsif(defined $json->{assay_term_name}){
    $ftype = $json->{assay_term_name};
  }
  else{
    throw("Implement ". $data->{experiment_url});
  }

  if(exists $cfg->{feature_type}->{$ftype}){
    $ftype = $cfg->{feature_type}->{$ftype};
  }

  return $ftype;
}


 sub assing_analysis {
  my ($json, $data) = @_;

  my $anal;
  if (defined $json->{assay_term_name}) {
    $anal = $json->{assay_term_name};
  }
  else {
    throw "No analysis found";
  }
  return $anal;
 }
###############################################################################
#                           fetch_and_decode_json
###############################################################################

sub fetch_and_decode_json {
  my ($url, $headers) = @_;

  if(!exists $headers->{'Content-Type'} ){
   $headers->{'Content-Type'} = 'application/json';
  }
  my $response = $http->get($url, {headers => $headers});

  _test_json_response($response);
  $response = JSON::XS::decode_json ($response->{content});

  return $response;
}

###############################################################################
#                           parse_ENCODE_json
###############################################################################

sub parse_ENCODE_json {
  my ($cfg, $json, $data) = @_;

  if($json->{status} ne 'released'){
    throw "Implement";
  }

  $data->{accession}          = $json->{accession};
  $data->{availability_date}  = $json->{date_released};
  $data->{cell_type}          = $json->{biosample_term_name};
  $data->{description}        = $json->{description};
  $data->{lab}                = $json->{lab}->{institute_label};
  # test if they are all defined, avoid autovivification
  test_root_elemets($cfg, $data);
  $data->{analysis}           = assing_analysis($json, $data);
  $data->{feature_type}       = assign_feature_type($json, $data, $cfg);

  $data->{files}              = parse_ENCODE_json_files($cfg, $json, $data);


  return $data;
}



###############################################################################
#                           test_root_elemets
###############################################################################
sub test_root_elemets {
  my ($cfg, $data) = @_;

  foreach my $key (sort keys %{$data}){
    if( (!defined $data->{$key}) or (length ($data->{$key} ) == 0) ) {
      my $json_key = $cfg->{json_root}->{$key};
      my $acc = $data->{accession};
      my $url = $cfg->{url}->{file};
      throw("$acc: JSON element $json_key not defined. $url");
    }
  }
}


###############################################################################
#                           parse_ENCODE_json_files
###############################################################################
sub parse_ENCODE_json_files {
  my ($cfg, $json, $data) = @_;

  if(scalar(@{$json->{files}}) < 1 ){
    throw 'Not an ARRAYREF';
  }
  if(! defined $json->{files}){
    throw('$json->{files} not defined');
  }

  my @files;

  for my $f (@{$json->{files}}){
    # say dump_data($file,1,1);
    next unless($f->{output_category} eq 'raw data');
    my $tmp = {};
    $tmp->{name}          = $f->{accession};
    $tmp->{md5sum}        = $f->{md5sum};
    $tmp->{download_url}  = $cfg->{url}->{base} . $f->{href};
    $tmp->{bio_replicate} = $f->{replicate}->{biological_replicate_number};
    $tmp->{replicate}     = $f->{replicate}->{technical_replicate_number};
    $tmp->{uuid}          = $f->{replicate}->{uuid};

    if(! defined $f->{replicate}->{uuid}){
      my $acc = $f->{accession};
      my $url = $cfg->{url}->{file};
      throw ("Missing UUID in file $acc $url");
    }

    parse_ENCODE_replicates($cfg, $tmp, $json, $data);
    push(@files, $tmp);
  }
  return \@files;
}


###############################################################################
#                            parse_ENCODE_replicates
###############################################################################
# iterates through all root.replicates. Stores data from Replicates matching
# the file currently processed (same uuid),

sub parse_ENCODE_replicates {
  my ($cfg, $tmp, $json, $data) = @_;

  if(!defined $json->{replicates}){
    croak 'json.replicates missing';
  }
  my $found = 0;
  for my $r(@{$json->{replicates}}){
    my $bs =  $r->{library}->{biosample};
    next unless ($r->{uuid} eq $tmp->{uuid});
    $found = 1;
    $tmp->{age}           = (defined $bs->{age}) ? $bs->{age} : undef;
    $tmp->{age_units}     = (defined $bs->{age_units})?$bs->{age_units}:undef;
    $tmp->{biosample_acc} = $bs->{accession};
    $tmp->{description}   = $bs->{description};
    $tmp->{life_stage}    = $bs->{life_stage};
    $tmp->{sex}           = $bs->{model_organism_sex};

  }
  if($found == 0){
    my $uuid = $tmp->{uuid};
    my $url = $data->{experiment_url};
    croak "Could not find replicate $uuid [$url]";
  }
  return;
}

###############################################################################
#                            assign_cell_type_name
###############################################################################

sub assign_cell_type_name {
  my ($json) = @_;

  my $name = undef;
  if(defined $json->{biosample_term_name}){
    $name = $json->{biosample_term_name}
  }
  elsif(defined $json->{replicates}) {
    say $json->{'@id'};die;
    for my $r(@{$json->{replicates}}) {
      my $tmp = $r->{library}->{biosample}->{biosample_term_name};
      if( (defined $name) && ($name ne $tmp) ){
        throw($json->{accession}. " $name vs $tmp");
      }
      $name = $tmp;
    }

  }

  if(!defined $name) {
    throw($json);
  }
  return $name;

}


###############################################################################
#                               parse_flatfile
###############################################################################

sub parse_flatfile {
  my ($filename) = @_;

  my @ids;
  open(my $fh,'<',$filename) or die "Can not open/access '$filename'\n$!";
    while(my $line = <$fh>){
      $line =~ /^(ENC\S+)\s/;
      push(@ids, $1);
    }
  close($fh);
  return \@ids;
}


###############################################################################
#                            _test_json_response
###############################################################################
sub _test_json_response {
  my ($response) = @_;

  my $status = $response->{status};

  if(!$response->{success}) {
    # Quickly check for rate limit exceeded & Retry-After (
    #  lowercase due to our client)
    if($status == 429 && exists $response->{headers}->{'retry-after'}) {
      my $retry = $response->{headers}->{'retry-after'};
      Time::HiRes::sleep($retry);
      # After sleeping see that we re-request
      # return perform_rest_action($endpoint, $parameters, $headers);
    }
    else {
      my ($status, $reason) = ($response->{status}, $response->{reason});
      die "Failed for endpoint! Status code: ${status}. Reason: ${reason}\n";
    }
  }

}




###############################################################################
#                            _get_trackingDB_adaptors
###############################################################################
sub _get_trackingDB_adaptors {
  my ($cfg) = @_;

# Tracking DB hidden from user, hence no get_TrackingAdaptor method.
# TrackingAdaptor->new() does not YET accept DBAdaptor object

  $cfg->{tr_a}->{tr} =
    Bio::EnsEMBL::Funcgen::DBSQL::TrackingAdaptor->new (
        -user       => $cfg->{efg_db}->{user},
        -pass       => $cfg->{efg_db}->{pass},
        -host       => $cfg->{efg_db}->{host},
        -port       => $cfg->{efg_db}->{port},
        -dbname     => $cfg->{efg_db}->{dbname},
        -species    => $cfg->{generic}->{species},
        -dnadb_user => $cfg->{dna_db}->{user},
        -dnadb_pass => $cfg->{dna_db}->{pass},
        -dnadb_host => $cfg->{dna_db}->{host},
        -dnadb_port => $cfg->{dna_db}->{port},
        -dnadb_name => $cfg->{dna_db}->{dbname},
        );

  my $db_a = $cfg->{tr_a}->{tr}->db;

  $cfg->{tr_a}->{ct} = $db_a->get_CellTypeAdaptor();
  $cfg->{tr_a}->{ft} = $db_a->get_FeatureTypeAdaptor();
  $cfg->{tr_a}->{an} = $db_a->get_AnalysisAdaptor();

  $cfg->{tr_a}->{eg} = $db_a->get_ExperimentalGroupAdaptor();

  $cfg->{tr_a}->{ex} = $db_a->get_ExperimentAdaptor();
  $cfg->{tr_a}->{iss} = $db_a->get_InputSubsetAdaptor();

  $cfg->{tr_a}->{rs} = $db_a->get_ResultSetAdaptor();
  $cfg->{tr_a}->{rf} = $db_a->get_RegulatoryFeatureAdaptor();

  $cfg->{tr_a}->{fs} = $db_a->get_FeatureSetAdaptor();
  $cfg->{tr_a}->{ds} = $db_a->get_DataSetAdaptor();
  $cfg->{tr_a}->{af} = $db_a->get_AnnotatedFeatureAdaptor();
}

###############################################################################
#                            _connect_to_trackingDB
###############################################################################
sub _connect_to_trackingDB {
  my ($cfg) = @_;

  # say dump_data($cfg->{efg_db},1,1);
  my $db_a = Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor->new (
    -user       => $cfg->{efg_db}->{user},
    -pass       => $cfg->{efg_db}->{pass},
    -host       => $cfg->{efg_db}->{host},
    -port       => $cfg->{efg_db}->{port},
    -dbname     => $cfg->{efg_db}->{dbname},
    -dnadb_name => $cfg->{dna_db}->{dbname},
    );
  $db_a->dbc->do("SET sql_mode='traditional'");
  say "\nConnected to trDB: " . $cfg->{efg_db}->{dbname}  ."\n";

  return($cfg->{dba_tracking} = $db_a);
}


