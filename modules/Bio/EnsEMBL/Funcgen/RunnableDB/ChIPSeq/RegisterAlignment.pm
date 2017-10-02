package Bio::EnsEMBL::Funcgen::RunnableDB::ChIPSeq::RegisterAlignment;

use strict;
use base 'Bio::EnsEMBL::Hive::Process';
use Data::Dumper;

sub run {

  my $self = shift;
  
  my $species = $self->param_required('species');
  my $plan    = $self->param_required('plan');
  
  my $align_plan = $plan
    ->{remove_duplicates}
    ->{align}
  ;
  
  my $alignment_name = $align_plan->{name};
  my $read_names     = $align_plan->{read_files};
  my $bam_file       = $align_plan->{bam_file}->{stored};
  
  $self->register_alignment(
    {
      bam_file_path  => $bam_file,
      read_names     => $read_names,
      species        => $species,
      alignment_name => $alignment_name,
    }
  );

  my $remove_duplicates_plan = $plan
    ->{remove_duplicates}
  ;
  
  my $alignment_name = $remove_duplicates_plan->{name};
  my $read_names     = $read_names;
  my $bam_file       = $remove_duplicates_plan->{bam_file}->{stored};
  
  $self->register_alignment(
    {
      bam_file_path  => $bam_file,
      read_names     => $read_names,
      species        => $species,
      alignment_name => $alignment_name,
    }
  );

  return;
}

sub register_alignment {

  my $self  = shift;
  my $param = shift;
  
  my $bam_file_path  = $param->{bam_file_path};
  my $read_names     = $param->{read_names};
  my $species        = $param->{species};
  my $alignment_name = $param->{alignment_name};

  my $read_file_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'ReadFile'
      );
  
  my @read_file_ids;
  foreach my $read_name (@$read_names) {
    my $read_file = $read_file_adaptor->fetch_by_name($read_name);
    push @read_file_ids, $read_file->dbID;
  }
  
  use Bio::EnsEMBL::Funcgen::DataFile;
  my $data_file = Bio::EnsEMBL::Funcgen::DataFile
    ->new(
      -table_id     => 0,
      -table_name   => 'alignment',
      -path         => $bam_file_path,
      -file_type    => 'BAM',
      -md5sum       => undef,
    );

  my $data_file_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'DataFile'
    );
  $data_file_adaptor->store($data_file);

  my $alignment_adaptor = Bio::EnsEMBL::Registry
    ->get_adaptor(
        $species, 
        'funcgen', 
        'Alignment'
    );
  
  my $alignment = Bio::EnsEMBL::Funcgen::Alignment->new(
    -name          => $alignment_name,
    -analysis_id   => 1,
    -read_file_ids => \@read_file_ids,
    -bam_file_id   => $data_file->dbID,
  );
  
  $alignment_adaptor->store($alignment);

  return;
}

1;
