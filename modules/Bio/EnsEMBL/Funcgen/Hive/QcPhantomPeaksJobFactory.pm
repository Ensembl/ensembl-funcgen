
=pod 

=head1 NAME

Bio::EnsEMBL::Funcgen::Hive::QcPhantomPeaksJobFactory

=head1 DESCRIPTION

default_stuff='out_db => {"-dbname" => "mn1_faang2_tracking_homo_sapiens_funcgen_81_38","-host" => "ens-genomics1","-pass" => "ensembl","-port" => 3306,"-user" => "ensadmin"}, work_root_dir => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/debug", data_root_dir => "/lustre/scratch109/ensembl/funcgen/mn1/ersa/faang/", pipeline_name => "blah", use_tracking_db => 1, dnadb => {"-dnadb_host" => "ens-livemirror","-dnadb_name" => "homo_sapiens_core_82_38","-dnadb_pass" => "","-dnadb_port" => 3306,"-dnadb_user" => "ensro"}'

standaloneJob.pl Bio::EnsEMBL::Funcgen::Hive::QcPhantomPeaksJobFactory -input_id "{ $default_stuff, input_subset_id => 1234, }"

export R_LIBS=/software/ensembl/funcgen/R-modules

        {   -logic_name => 'QcRunPhantomPeaks',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LSF',
            -parameters => { 
		  cmd => 
		    qq( Rscript /software/ensembl/funcgen/spp_package/run_spp.R )
		  . qq(    -c=#bam_file# )
		  . qq(    -savp -out=#phantom_peak_out_file# )
            },
            -flow_into => { 
	      '1' => [ 'QCLoadPhantomPeaksToDB' ],
            },
        },
        
        {   -logic_name => 'QCLoadPhantomPeaksToDB',
            -module     => 'Bio::EnsEMBL::Hive::RunnableDB::SystemCmd',
            -meadow_type=> 'LOCAL',
	    -parameters => {
                'cmd' =>
		    qq( load_phantom_peak_file.pl )
		  . qq(    --result_set_id #result_set_id# )
		  . qq(    --result_file #phantom_peak_out_file# )
		  . qq(    --user #tracking_db_user# --pass #tracking_db_pass# --host #tracking_db_host# --dbname #tracking_db_name# )
            },
            
=cut

package Bio::EnsEMBL::Funcgen::Hive::QcPhantomPeaksJobFactory;

use warnings;
use strict;

use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

sub run {
  my $self = shift;
  my $input_subset_id = $self->param('input_subset_id');

  my $rset;
  
  my $set_type = $self->param_required('set_type');
  if ($set_type eq 'ResultSet') {
    $rset = $self->fetch_Set_input('ResultSet'); 
  } else {
    my $fset     = $self->fetch_Set_input('FeatureSet');
    my $analysis = $fset->analysis;
    $rset        = $self->ResultSet; 
  }

  my $input_id = $self->create_input_id($rset);
  
  use Data::Dumper;
  $Data::Dumper::Maxdepth = 0;
  print Dumper($input_id);

  $self->dataflow_output_id($input_id, 2);
  return;
}

sub create_input_id {

  my $self = shift;
  
  my $result_set = shift;
  my $result_set_id = $result_set->dbID;

  my $align_prefix   = $self->get_alignment_path_prefix_by_ResultSet($result_set, undef, 1);#validate aligned flag 
  my $control_prefix = $self->get_alignment_path_prefix_by_ResultSet($result_set, 1, 1);#and control flag 
  
  my $signal_bam_file  = $align_prefix   . '.bam';
  my $control_bam_file = $control_prefix . '.bam';
  
  if (! -e $signal_bam_file) {
    die("$signal_bam_file doesn't exist!");
  }
  if (! -e $control_bam_file) {
    die("$control_bam_file doesn't exist!");
  }
  
  my $work_dir = $self->param_required('work_root_dir');
  my $temp_dir = "$work_dir/temp/Qc/PhantomPeaks/$result_set_id";
  my $out_db = $self->param('out_db');  
  
  use File::Path qw( make_path );
  make_path($temp_dir);
  
  use File::Basename;
  (my $signal_bam_file_base_name,  my $signal_bam_directory)  = fileparse($signal_bam_file);
  (my $control_bam_file_base_name, my $control_bam_directory) = fileparse($control_bam_file);
  
  die unless($signal_bam_directory eq $control_bam_directory);

  my $input_id_common = [
      # Directory into which the bam files will be copied
      tempdir            => $temp_dir,      
      result_set_id      => $result_set_id,
      
      # Connection details for the db to which the results will be written
      tracking_db_user   => $out_db->dbc->user,
      tracking_db_pass   => $out_db->dbc->password,
      tracking_db_host   => $out_db->dbc->host,
      tracking_db_name   => $out_db->dbc->dbname,
      tracking_db_port   => $out_db->dbc->port,
  ];
  
  my $signal_phantom_peaks_file  = "$temp_dir/${signal_bam_file_base_name}.signal_phantom_peaks.txt";
  my $control_phantom_peaks_file = "$temp_dir/${control_bam_file_base_name}.control_phantom_peaks.txt";
  
  my $input_id = [
    {
      @$input_id_common,
      phantom_peak_out_file  => $signal_phantom_peaks_file,
      bam_file => $signal_bam_file,
      is_control => undef,
    },
    {
      @$input_id_common,
      phantom_peak_out_file  => $control_phantom_peaks_file,
      bam_file => $control_bam_file,
      is_control => 1,
    }
  ];
  return $input_id;
}

1;


