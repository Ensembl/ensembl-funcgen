use strict;
use warnings;
use feature qw(say);
use Carp;
use JSON;
use REST::Client;
use Data::Dumper;
use URI::Escape;
use Epigenome;
use Experiment;
use File;
use JSON::MaybeXS;

main();

# IHEC:IHECRE00001827.1
sub main {
  my $self = bless({}, __PACKAGE__);

  $self->init_REST();
  $self->{tsv} = 
  '/homes/juettema/src/ensembl-funcgen/scripts/registration/ReferenceEpigenome_Report_2018_10_18.tsv';
  $self->set_assay_type();
  $self->parse_encode_tsv();
  die;
}

sub init_REST {
  my ($self) = @_;
  my $client = REST::Client->new({
    host => 'https://www.encodeproject.org',
    });
  $client->addHeader('Accept', 'application/json');
  $client->addHeader('Content-Type', 'application/json');

  $self->{client} = $client;
}

sub set_assay_type {
  my ($self) = @_;

  $self->{assays}->{'ChIP-seq'}   = 1;
  $self->{assays}->{'DNase-seq'}  = 1;
  $self->{assays}->{'ATAC-seq'}   = 1;
}

# 0:Accession 1:Biosample term  2:Assay Type  3:Description 4:Project 5:Status  
# 6:Species 7:External identifiers  8:Assay term id   9:Biosample type  10:Assay synonyms  
# 11:Biosample term id  12:Biosample synonyms
sub parse_encode_tsv {
  my ($self) = @_;

  my $file = $self->{tsv};
  open(my $fh, '<:encoding(UTF-8)', $file) or croak "Could not open $file: $!";
    my @lines = <$fh>;
  close($fh);
  #date and link
  shift @lines;
  #Column headers, unfortunately with spaces
  shift @lines;
  foreach my $line(@lines) {
    chomp($line);
    my $epigenome = $self->epigenome_from_tsv($line);   
    next if(!defined $epigenome);
    say $epigenome->ihec;
    # say $epigenome->print_all_attributes;
    my $epigenome_json = $self->fetch_epigenome_json($epigenome);

    $self->iterate_datasets($epigenome_json);
    # my $files = $self->create_file_objects($epigenome_json);

    die;
  }
}

sub epigenome_from_tsv {
  my ($self, $tsv) = @_;
  my @f = split("\t", $tsv);
  return undef if($f[7] !~ /^IHEC:/);
  my $epigenome = new Epigenome(
      accession           => $f[0],
      biosample_term      => $f[1],
      assay_type          => $f[2],
      description         => $f[3],
      project             => lc($f[4]),
      status              => $f[5],
      species             => $f[6], 
      ihec                => $f[7],
      assay_term_id       => $f[8],
      biosample_type      => $f[9],
      assay_synonyms      => $f[10],
      biosample_term_id   => $f[11],
      biosample_synonyms  => $f[12],
    );
  return($epigenome);
}

# 
sub fetch_epigenome_json {
  my ($self, $epigenome) = @_;

  my $acc = $epigenome->{accession};
  my $url = "reference-epigenomes/$acc/?format=json";
  my $r = from_json($self->{client}->GET($url)->responseContent());
  return($r);
 }

sub iterate_datasets {
  my ($self, $epigenome_json) = @_;


  my @datasets = @{$epigenome_json->{related_datasets}};
  foreach my $ds (@datasets) {
    next unless (defined $self->{assays}->{$ds->{assay_term_name}});
    my $exp_endpoint = $ds->{'@id'};

    $self->create_experiment($exp_endpoint);
    # my @files = @{$ds->{files}};
    # foreach my $f (@files){
    #   next unless ($f->{file_type} eq 'fastq');
    # }
    # my $experiment = new Experiment(
    #   accession       => $rd->{accession},
    #   assay_term_name => $rd->{assay_term_name},
    #   status          => $rd->{status},
    # );

  }
}

sub create_experiment {
  my ($self, $exp_endpoint) = @_;

  my $url = $exp_endpoint . '?format=json';
  my $r = from_json($self->{client}->GET($url)->responseContent());
say $r->{target}->{label};
say $r->{target}->{name};
}

# sub create_file_objects {
#   my ($self, $epigenome_json) = @_;

#   my @files = @{$epigenome_json->{files}};

#   my @objects;
#   foreach my $f (@files) {
#     next unless($f->{file_type} eq 'fastq');
#     my $file = new File(
#       accession   => $f->{accession},
#       file_type   => $f->{file_type},
#       md5sum      => $f->{md5sum},
#       read_length => $f->{},
#       read_count  => $f->{},
#       run_type    => $f->{},
#       aliases     => $f->{},
#       );
#   }
# }

# sub add_experiment {
#   my ($self, $epigenome) = @_;

#   my $ihec = uri_escape($epigenome->ihec());

#   my $url = "search/?type=Experiment&related_series.dbxrefs=$ihec&format=json";
#   say $url;
#   $self->{client}->GET($url);
#   my $r = JSON->new->utf8->decode($self->{client}->responseContent());
#   my $experiments = $r->{'@graph'};
#   foreach my $e (@{$experiments}) {
#     next unless ( defined $self->{assays}->{ $e->{assay_title} } );

#     my $files = $e->{files};
#     foreach my $f(@{$files}){
#       next unless($f->{file_format} eq 'fastq');
#       my $file = $self->create_file($f);
#       say Dumper($file);
# die;
#     }

#     say Dumper($e);die;

#   }
# }

sub populate_experiment {

}

sub fetch_file {
  my ($self, $encsr) = @_;

  my $url = "search/?type=File&file_format=fastq&dataset=/experiments/$encsr/&format=json";

}

sub create_file {
  my ($self,$json) = @_;

  my $file = new File(
    accession   => $json->{accession},
    file_type   => $json->{file_type},
    md5sum      => $json->{md5sum},
    read_length => $json->{read_length},
    read_count  => $json->{read_count},
    run_type    => $json->{run_type},
    );
  
  if(defined $json->{aliases} and scalar( @{$json->{aliases}} ) >0){
    foreach my $a ( @{ $json->{aliases} } ) {
      $file->file_push($a);
    }
  }
  return($file);
}
