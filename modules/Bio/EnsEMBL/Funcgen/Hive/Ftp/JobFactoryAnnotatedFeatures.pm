package Bio::EnsEMBL::Funcgen::Hive::Ftp::JobFactoryAnnotatedFeatures;

use strict;
use Data::Dumper;
use Bio::EnsEMBL::Registry;
use base ('Bio::EnsEMBL::Hive::Process');

sub run {
  my $self = shift;

  my $temp_dir = $self->param('temp_dir');
  my $species  = $self->param('species');

  use File::Path qw(make_path remove_tree);
  make_path($temp_dir);

  my $batch_size = 10_000_000;
  my $sql = <<SQL
    select 
      floor(annotated_feature_id / ?) as batch_number, 
      min(annotated_feature_id) as min_id, 
      max(annotated_feature_id) as max_id 
    from 
      annotated_feature 
    group by 
      batch_number 
    order by 
      batch_number asc;
SQL
;

  my $coord_system_adaptor = Bio::EnsEMBL::Registry->get_adaptor( $species, 'Core', 'CoordSystem' );
  if (!$coord_system_adaptor) {
    die("Can't get coord system adaptor! Please configure your registry accordingly.")
  }
  my ($cs) = @{$coord_system_adaptor->fetch_all()};
  my $assembly = $cs->version();
  if (!$assembly) {
    die("Can't work out assembly for $species!")
  }
#   use Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator;
#   my $ref_build_file_locator = Bio::EnsEMBL::Funcgen::Hive::RefBuildFileLocator->new;
  
  my $funcgen_adaptor = Bio::EnsEMBL::Registry->get_DBAdaptor( $species, 'Funcgen' );
  my $dbc = $funcgen_adaptor->dbc;
  
  use Bio::EnsEMBL::Utils::SqlHelper;
  my $helper = Bio::EnsEMBL::Utils::SqlHelper->new(
    -DB_CONNECTION => $funcgen_adaptor->dbc
  );
  
  my @all_annotated_feature_batches;
  
  $helper->execute_no_return(
    -SQL      => $sql,
    -PARAMS => [ $batch_size ],
    -CALLBACK => sub {
      my @row = @{ shift @_ };

      my $batch_number = $row[0];
      my $min_id       = $row[1];
      my $max_id       = $row[2];

      my $current_batch = {
	batch_number => $batch_number,
	min_id       => $min_id,
	max_id       => $max_id,
      };
      push @all_annotated_feature_batches, $current_batch;
    }
  );

  my $sql_fetch_feature_sets = <<SQL
    select 
      feature_set.feature_set_id as feature_set_id, 
      feature_type.name          as feature_type_name, 
      epigenome.production_name  as epigenome_production_name,
      epigenome.gender           as epigenome_gender
    from 
      feature_set 
      join feature_type using (feature_type_id) 
      join epigenome using (epigenome_id) 
    where 
      feature_set.type="annotated"
      and feature_set_id in (select distinct feature_set_id from annotated_feature)
SQL
;
  my $sth_fetch_feature_sets = $dbc->prepare($sql_fetch_feature_sets);
  $sth_fetch_feature_sets->execute();
  my $x = $sth_fetch_feature_sets->fetchall_hashref('feature_set_id');

  # @feature_set_list looks like this:
  # 
  #   {
  #     'feature_type_name' => 'ELF1',
  #     'feature_set_id' => '38',
  #     'epigenome_production_name' => 'MEL'
  #   },
  #   {
  #     'feature_type_name' => 'H3K4me1',
  #     'feature_set_id' => '4',
  #     'epigenome_production_name' => 'thymus_a8w'
  #   },
  #   {
  #     'feature_type_name' => 'H3K9ac',
  #     'feature_set_id' => '34',
  #     'epigenome_production_name' => 'heart_a8w'
  #   },
  # 
  my @feature_set_list = values %$x;
#   print Dumper(\@feature_set_list);
  use Hash::Util qw( lock_hash unlock_hash );
  
  foreach my $current_feature_set (@feature_set_list) {
  
    lock_hash( %$current_feature_set );
    
    my $gff_files_from_batches = $temp_dir . "/gff_files_from_batches.".$current_feature_set->{epigenome_production_name}.'.'.$current_feature_set->{feature_type_name}.".txt";
    open GFF_OUT, ">$gff_files_from_batches";

    my $bed_files_from_batches = $temp_dir . "/bed_files_from_batches.".$current_feature_set->{epigenome_production_name}.'.'.$current_feature_set->{feature_type_name}.".txt";
    open BED_OUT, ">$bed_files_from_batches";
  
    foreach my $current_annotated_feature_batch (@all_annotated_feature_batches) {
    
      lock_hash( %$current_annotated_feature_batch );
  
      my @directory_components;
      push @directory_components, $temp_dir;
      push @directory_components, reverse split '', $current_annotated_feature_batch->{batch_number};
      push @directory_components, $species;
      push @directory_components, $current_feature_set->{epigenome_production_name};
      push @directory_components, $current_feature_set->{feature_type_name};

      my $directory = join '/', @directory_components;
      my $gff_file      = join '_', (
	'annotated_features',
  #       $feature_set_id
	$current_annotated_feature_batch->{min_id},
	$current_annotated_feature_batch->{max_id},
      );
      $gff_file .= '.gff';
      my $bed_file      = join '_', (
	'annotated_features',
  #       $feature_set_id
	$current_annotated_feature_batch->{min_id},
	$current_annotated_feature_batch->{max_id},
      );
      $bed_file .= '.bed';
      print GFF_OUT "$directory/$gff_file\n";
      print BED_OUT "$directory/$bed_file\n";
      
      unlock_hash(%$current_annotated_feature_batch);
      
      $current_annotated_feature_batch->{temporary_directory} = $directory;
      $current_annotated_feature_batch->{partial_gff_file}  = $gff_file;
      $current_annotated_feature_batch->{partial_bed_file}  = $bed_file;
      $current_annotated_feature_batch->{feature_type_name}          = $current_feature_set->{feature_type_name};
      $current_annotated_feature_batch->{feature_set_id}             = $current_feature_set->{feature_set_id};
      $current_annotated_feature_batch->{epigenome_production_name}  = $current_feature_set->{epigenome_production_name};
      
      lock_hash(%$current_annotated_feature_batch);
      
      $self->dataflow_output_id($current_annotated_feature_batch, 2);
    }
    
    close(GFF_OUT);
    close(BED_OUT);
    
#     my $chromosome_sizes_file = $ref_build_file_locator->locate({
#       species          => $species,
#       assembly         => $assembly,
#       epigenome_gender => $current_feature_set->{epigenome_gender},
#       file_type        => 'chromosome_lengths_by_species_assembly',
#     });
#     my $chromosome_sizes_file_ucsc = "$temp_dir/$species.$assembly";

    $self->dataflow_output_id({
      gff_files_from_batches     => $gff_files_from_batches,
      bed_files_from_batches     => $bed_files_from_batches,
      merged_gff                 => $temp_dir . '/AnnotatedFeatures.'.$current_feature_set->{epigenome_production_name}.'.'.$current_feature_set->{feature_type_name}.'.gff',
      merged_bed                 => $temp_dir . '/AnnotatedFeatures.'.$current_feature_set->{epigenome_production_name}.'.'.$current_feature_set->{feature_type_name}.'.bed',
      converted_big_bed          => $temp_dir . '/AnnotatedFeatures.'.$current_feature_set->{epigenome_production_name}.'.'.$current_feature_set->{feature_type_name}.'.bb',
      species                    => $species,
      feature_type_name          => $current_feature_set->{feature_type_name},
      feature_set_id             => $current_feature_set->{feature_set_id},
      epigenome_production_name  => $current_feature_set->{epigenome_production_name},
      epigenome_gender           => $current_feature_set->{epigenome_gender},
      assembly                   => $assembly,
#       chromosome_sizes_file      => $chromosome_sizes_file,
#       chromosome_sizes_file_ucsc => $chromosome_sizes_file_ucsc
    }, 1);
  }
}

1;
