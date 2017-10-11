=pod

=head1 NAME

Bio::EnsEMBL::Funcgen::Sequencing::SeqTools

=head1 DESCRIPTION


=cut

package Bio::EnsEMBL::Funcgen::Sequencing::SeqTools;

use warnings;
use strict;

use Net::FTP;
use feature qw(say);
use DBI     qw(:sql_types);

use Bio::EnsEMBL::Funcgen::Experiment;

use File::Basename                         qw( dirname basename );
use Bio::EnsEMBL::Utils::Scalar            qw( assert_ref );
use Bio::EnsEMBL::Utils::Argument          qw( rearrange );
use Bio::EnsEMBL::Utils::Exception         qw( throw );
use Bio::EnsEMBL::Funcgen::Utils::EFGUtils qw( run_system_cmd
                                               run_backtick_cmd
                                               validate_package_path
                                               which_path );

use base qw( Exporter );
use vars qw( @EXPORT );

@EXPORT = qw(
  _init_peak_caller
  _run_peak_caller
  merge_bams
  remove_duplicates_from_bam
  post_process_IDR
  pre_process_IDR
  run_IDR
  run_peak_caller
  );

sub merge_bams {

  my $param = shift;
  
  my $bams    = $param->{input_bams};
  my $outfile = $param->{output_bam};
  my $debug   = $param->{debug};

  assert_ref($bams, 'ARRAY', 'bam files');
  
  if(! scalar(@$bams)){
    throw('Must provide an arrayref of bam files to merge');
  }

  my $cmd;
  my $skip_merge = 0;
  if(scalar(@$bams) == 1){
    #samtools merge cannot handle a single input!
    #Instead it throws a seemingly completely unrelated error message:
    #Note: Samtools' merge does not reconstruct the @RG dictionary in the header. Users
    #  must provide the correct header with -h, or uses Picard which properly maintains
    #  the header dictionary in merging.

    #Rather than having the caller have to handle this, let's just do the expected thing here
    #and warn.
    $skip_merge = 1;
    warn 'Only 1 bam file has been specified, merge will be skipped, '.
      "otherwise file will be processed accordingly\n";
  }
  if(! $skip_merge) {
    # -f option forces overwriting the outfile when it already exists. This
    # is useful, if a failed job is being rerun.
    #
    use Data::Dumper;
    print Dumper($bams);
    
    $cmd = "samtools merge -f $outfile ".join(' ', @$bams);
  } else{ 
    $cmd = 'cp '.$bams->[0].' '.$outfile;
  }
  warn "Merging with:\n$cmd\n" if $debug;
  run_system_cmd($cmd);  
  warn "Finished merge to $outfile" if $debug;
  return;
}

sub remove_duplicates_from_bam {

  my $param = shift;
  
  my $input_bam  = $param->{input_bam};
  my $output_bam = $param->{output_bam};
  my $debug      = $param->{debug};

#   my $metrics_file = "${output_bam}.merged_duplication_removal_metrics.tab";
  
  # Remove unmapped reads
  #
  # Must set -b or picard will complain: "Error parsing text SAM file. Empty 
  # sequence dictionary."
  #
  my $cmd_removeUnmappedReads = qq(samtools view -F 4 -b -o - -);
  
  # Picard must be in the classpath before running this module, e.g. like this:
  # export CLASSPATH=/software/ensembl/funcgen/picard.jar
  #
  # The output is always a sam file, even if bam was specified, hence 
  # SamFormatConverter is run on the output.
  #
#   my $cmd_MarkDuplicates = qq(java picard.cmdline.PicardCommandLine MarkDuplicates ) 
#   . qq( REMOVE_DUPLICATES=true ) 
#   . qq( VALIDATION_STRINGENCY=LENIENT ) 
#   . qq( ASSUME_SORTED=true     ) 
#   . qq( INPUT=$input_bam       ) 
#   . qq( COMPRESSION_LEVEL=0    ) 
#   . qq( OUTPUT=/dev/stdout ) 
#   # Prevent any messages from going into the pipe:
#   . qq( QUIET=true ) 
#   . qq( METRICS_FILE=$metrics_file );

  my $cmd_MarkDuplicates = qq(samtools rmdup $input_bam -);
  
  #my $cmd = qq(bash -o pipefail -c "$cmd_removeUnmappedReads | $cmd_MarkDuplicates | $cmd_SamFormatConverter");
  my $cmd = qq(bash -o pipefail -c "$cmd_MarkDuplicates | $cmd_removeUnmappedReads > $output_bam");
  warn "Running\n$cmd\n";
  run_system_cmd($cmd);

#   unlink($metrics_file);
  return;
}

#These are split into _init and _run methods to allow initialisation before running the
#peak caller proper. This fits well with the fetch_input/run jobs procedure of the hive.

#Can't do this easily, as we don't know exactly what we want to pass through to the
#PeakCaller. Isn't this already done with run_aligner?
#Actually, we have acces to get_files_by_format here
#so we could optionally defined the signal and control files here given the right input
#Probably to separate out get_alignment_files_by_ResultSet_formats
#If this is possible, then we can add FeatureSet support here and move all the Analysis
#processing here from RunPeaks::fetch_input

sub _init_peak_caller{
  my($peak_pkg, $params, $analysis, $signal_alignment, $control_alignment,
    $sam_ref_fai,$debug) = rearrange(
     [qw(peak_module peak_module_params analysis signal_alignment control_alignment
     sam_ref_fai debug)], @_);
     
  my (@params_array, $peak_module);

  if(defined $analysis){ #API mode

    if($peak_pkg){
      throw("Mututally exclusive paramters detected:\n\t-peak_module => $peak_pkg".
        "\n\t-analysis => $analysis(".$analysis->logic_name.')');
    }

    assert_ref($params, 'HASH', 'PeakCaller parameters');
    assert_ref($analysis, 'Bio::EnsEMBL::Analysis');

#     if(! defined $align_prefix){
#       throw('Must pass an -align_prefix in -analysis mode');
#     }
#     
    my $module = $analysis->module;

    $peak_module = validate_package_path($analysis->module);
    my $formats = $peak_module->input_formats;
    #my $filter_format = $self->param_silent('bam_filtered') ? undef : 'bam';
    #It is currently unsafe to filter here (control clash), so expect filtered file

    #The problem here is that we are returning a hash of files keys on the format
    #This conversion may cause clashes for fan job which share the same controls
    #(e.g. peak calling jobs if they require formats other than bam)
    #Collections jobs will be pre-processed/converted individually before submitting
    #the slice job.
    #So here we really only need the first available format

    #Restrict to bam (and bed explicitly for CCAT) for now
    my $format = 'bam';

    if($formats->[0] ne 'bam'){

      if($analysis->program eq 'CCAT'){

        warn 'Hardcoding CCAT format to bed until CollectionWriter can be made PeakCaller aware wrt to formats';
        $format = 'bed';
      }
      else{
        throw("It is currently unsafe to use any non-bam format at this point.\n".
          "This is due to the possibility of filtering/format conversion clashes between parallel\n".
          "jobs which share the same control files. Please implement PreprocessAlignments to\n".
          "group jobs by controls and set/handle FILTERING_CONTROL status");
      }
    }

    $formats = [$format];

    #May need to specify checksum param separately
    my $get_files_params = {debug              => $debug,
                            ref_fai            => $sam_ref_fai,  #Just in case we need to convert
                            skip_rmdups        => 1, #This will have been done in merge_bams
                            checksum           => undef,
                            #Specifying undef here turns on file based checksum generation/validation
                           };

#     my $align_file = get_files_by_formats($align_prefix,
#                                          $formats,
#                                          $get_files_params)->{$format};
# 
#     if(! defined $align_file){
#       throw("Failed to identify file using get_files_by_formats:\n\t".$align_prefix.".${format}");
#     }
# 
#     my $ctrl_file;
# 
#     if($control_prefix){
#       $ctrl_file = get_files_by_formats($control_prefix,
#                                        $formats,
#                                        $get_files_params)->{$format};
#     }

    $params->{-debug}          = $debug;
    $params->{-align_file}     = $signal_alignment;
    $params->{-control_file}   = $control_alignment;
    $params->{-parameters}     = $analysis->parameters;
    $params->{-program_file} ||= $analysis->program_file;
    #All the rest are passed in $params from RunPeaks
    #e.g.
    #-out_file_prefix
    #-out_dir
    #-convert_half_open
    #-is_half_open


    @params_array = %$params; #flatten hash
  }
  else{ #Script mode
    if(! defined $peak_pkg){
      throw('Must provide a -peak_module or -analysis parameter');
    }
    elsif($peak_pkg !~ /::/){
      #Might this be a path?
      warn "Full $peak_pkg namespace was not specified, defaulting to:\t".
        "Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::$peak_pkg\n";
      $peak_pkg = "Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller::$peak_pkg";
    }

    assert_ref($params, 'ARRAY', 'PeakCaller params');
    
    @params_array = @$params;
    
    $peak_module = validate_package_path($peak_pkg);
  }
  
#     use Data::Dumper;
#     print Dumper(\@params_array);
#     die;
  
  return $peak_module->new(@params_array);
}

sub _run_peak_caller{
  my($pcaller, $max_peaks, $file_types)  = rearrange([qw(peak_caller max_peaks file_types)], @_);
  assert_ref($pcaller, 'Bio::EnsEMBL::Funcgen::Sequencing::PeakCaller');

  #Do this first, so we fail early
  if( (defined $max_peaks) && (! $pcaller->can('filter_max_peaks')) ){
    throw(ref($pcaller).' cannot filter_max_peaks');
  }

  $pcaller->run; #eval in caller if required

  if(defined $max_peaks){

    foreach my $file_type(@$file_types){
      $pcaller->filter_max_peaks($max_peaks, $file_type)
    }
  }

  return;
}

sub run_peak_caller{
   my ($max_peaks, $file_types, $debug, $skip_qc) = rearrange([qw(max_peaks file_types debug skip_qc)], @_);
   my $peak_caller           = _init_peak_caller(@_);

   #Pass these explicitly, so we don't pass all the _init_peak_caller params too
   _run_peak_caller(-peak_caller => $peak_caller,
                    -max_peaks   => $max_peaks,
                    -file_types  => $file_types,
                    -skip_qc     => $skip_qc,
                    -debug       => $debug);



   #Now write/process/load as required
}

#Move peak caller post processing stuff in here too?
#This maybe overkill as it would require Set handling/defnition code


# Validate, count and define IDR threshold
#If you started with ~150 to 300K relaxed pre-IDR peaks for large genomes (human/mouse),
#then threshold of 0.01 or 0.02 generally works well.
#If you started with < 100K pre-IDR peaks for large genomes (human/mouse),
#then threshold of 0.05 is more appropriate.
#This is because the IDR sees a smaller noise component and the IDR scores get weaker.
#This is typically for use with peak callers that are unable to be adjusted to call large number of peaks (eg. PeakSeq or QuEST)
#What exactly are we counting here? Total number peaks across rep, average, or the max between reps?
#This also depends on the prevalence of the mark, it may be that a particular feature type genuinely does not have many genome wide hits

#idr_threshold is being skewed by the rep with the least peaks
#meaning more peaks will be identified from each replicate comparison
#This will however, still be low if we have had to truncate in line with
#the rep with the lowest number of pre-IDR peaks

sub pre_process_IDR{
  my $out_dir      = shift or throw('Must provide and out_dir argument');
  my $pre_idr_beds = shift;
  my $batch_name   = shift or throw('Must provide a batch_name argument');
  my $max_peaks_for_this_peak_caller    = shift or throw('Must provide a max_peaks_for_this_peak_caller argument');
  throw("out_dir does not exist:\t$out_dir") if ! -d $out_dir;
  assert_ref($pre_idr_beds, 'ARRAY');
  
  my $number_of_peaks_considered = $max_peaks_for_this_peak_caller;

  if(scalar(@$pre_idr_beds) < 2){
    throw("Cannot run IDR with less than 2 replicates:\n\t".$pre_idr_beds->[0]);
  }

  my ($mt_100k, $lt_100k, $log_txt);

  foreach my $bed_file(@$pre_idr_beds){

    if(! -f $bed_file){
      throw("Could not find pre-IDR bed file:\t$bed_file\n".
        'Need to dump bed from DB directly to required format!');
    }

    # Ignore comments and header
    my $cmd       = "grep -vE '(#|(^Region[[:space:]]+Start))' $bed_file | wc -l | awk '{print \$1}'";
    my $num_peaks = run_backtick_cmd($cmd);

    # $num_peaks can sometimes be 0. Throw here as it will just fail in the RunIDR step
    # We likely need to remove a replicate? 

    if($num_peaks == 0){
      throw("Found bed file with 0 peaks. This will cause IDR to fail, please remove or fix:\n\t$bed_file");
    }

    if($num_peaks > 100000){
      $mt_100k = 1;

      if($num_peaks < 150000){
        warn 'Pre-IDR peaks counts fall in threshold grey zone(100k-150k), defaulting to 0.01'.
          " but maybe consider 0.02 for:\n\t$bed_file\n";
      }
    }
    else{
      $lt_100k = 1;
    }

    if($num_peaks < $number_of_peaks_considered){
      # We take the lowest number of peaks, as we need comparable numbers of peaks across all inputs
      $number_of_peaks_considered = $num_peaks;
    }

    $log_txt .= $bed_file."\t".$num_peaks."\n";
  }

  # Note this does not yet support MACS yet, should prbably just ignore it as we filter to 100000
  my $idr_threshold   = ($number_of_peaks_considered < 100000) ? 0.05 : 0.01;
  # Could alternatively pass all thresholds back to the caller
  my $x_thresh_adjust = 0;

  if($lt_100k && $mt_100k){
    #$self->throw_no_retry("Identified different optimal thresholds due to pre-IDR peak counts spanning the 100k limit:\n".
    #  $log_txt);
    warn 'Identified different optimal thresholds due to pre-IDR peak counts spanning the 100k limit:'.
      "\n$log_txt\nDefaulting to lower threshold:\t$idr_threshold\n";
    $x_thresh_adjust = 1;
  }

  # TODO We need some mechanism to restart this job, to force the threshold, or by dropping 1/some of the replicates.

  my $cmd = "echo -e \"Pre-IDR File\tNum Peaks\n$log_txt\nIDR Threshold = $idr_threshold\"".
    "> ${out_dir}/${batch_name}-idr-stats.txt";
  run_system_cmd($cmd);

  # TODO parallelise the filtering and reformating to speed things up, so let's semphaore than to a simple CMD job.
  # can we even do this as we already have a semaphore on the RunIDR and
  # maybe with a job factory? I think this is not possible without another analysis
  # but we could use a dummy? which then submit the RunIDR and semaphores the PostProcessIDR
  # Just do here for now
  my @np_bed_files;

  foreach my $i(0..$#{$pre_idr_beds}){
    my $bed_file     = $pre_idr_beds->[$i];
    (my $np_bed_file = $bed_file) =~ s/\.txt$/.np_idr.txt/;
    # Never re-use np_idr output file in case it is truncated due to job failure/exit. 
    # Is there a danger of an np_idr file usage clash (pseudo/self consistency IDR?

    # SWEmbl output header::
    # Input  GSE30558_GSM758361_PrEC.sorted.bam
    # Reference      GSE30558_GSM758360_LNCaP.sorted.bam
    # Sequence length        0
    # Fragment length        0
    # Background     0.000000
    # Position Background    0.036383
    # Long Background        0.181917
    # Threshold      5.000000
    # Minimum count above bg 15
    # Penalty increase       70
    # Quality cutoff 0.000000
    # Result cutoff  0.000000
    # Penalty factor 0.552834

    # and fields:
    # Region        - Part of the genome build e.g. chromosome
    # Start pos.    - Base in region where peak starts
    # End pos.      - Base in region where peak ends
    # Count         - Number of reads in experimental sample in peak
    # Length        - Length of peak (distance between start and end pos.)
    # Unique pos.   - Number of unique bases within peak at which reads begin
    # Score         - The SWEMBL score, which is basically the count of filtered thresholded reads in the peak, minus the penalties (gap distances and reference sample reads).
    # Ref. count    - Number of reads in reference sample in peak
    # Max. Coverage - Depth of reads at summit
    # Summit        - Median position of highest read coverage

    # signalValue field was being set to (Count - Ref. Count)/min
    # Min is not really required here for ranking and simply add the header requirement
  
    # TODO Get header skipping regex from PeakCaller. Currently hardcoded for SWEmbl
    # TODO Handle pipes with perl pipe or IPC::open/run

    # Strip out the header and set signal.value to score, sort on score, filter based on $number_of_peaks_considered, resort based on position
    $cmd = 'awk \'BEGIN {OFS="\t"} { if($0 !~ /^(#|(Region[[:space:]]+Start))/) {print $1,$2,$3,".",$7,".",$7,-1,-1,int($9-$1)} }\' '.
      "$bed_file | sort -k 7nr,7nr | head -n $number_of_peaks_considered | sort -k 1,2n > ".$np_bed_file;
    run_system_cmd($cmd);

    # Sanity check we have the file with the correct number of lines
    $cmd = "wc -l $np_bed_file | awk '{print \$1}'";
    my $filtered_peaks = run_backtick_cmd($cmd);

    if($number_of_peaks_considered != $filtered_peaks){
      throw("Expected $number_of_peaks_considered in filtered pre-IDR bed file, but found $filtered_peaks:\n\t".$np_bed_file);
    }

    # Need to check this is != 0?


    # TODO check the feature_set_stat or states
    # such that we know the peak calling jobs has finished and passed QC!
    # Do this for each before we submit IDR jobs, as we may have to drop some reps
    push @np_bed_files, $np_bed_file;
  }

  return (\@np_bed_files, $idr_threshold, $x_thresh_adjust);
}

sub run_IDR {
  my ($out_dir, $output_prefix, $threshold, $bed_files, $batch_name) =
    rearrange([ qw(out_dir output_prefix threshold bed_files batch_name ) ], @_);

  defined $out_dir       or throw('Must provide an out_dir argument');
  defined $output_prefix or throw('Must provide an output_prefix argument');
  defined $threshold     or throw('Must provide an IDR threshold to count peaks');
  assert_ref($bed_files, 'ARRAY', 'bed_files');
  defined $batch_name    or throw('Must provide a batch name to log counts to idr-stats file');

  # Check we have 2 reps
  if(scalar (@$bed_files) != 2) {
    throw("run_IDR expect 2 replicate bed files:\t".join(' ', @$bed_files));
  }

  if($bed_files->[0] eq $bed_files->[1]) {
    throw("Pre-IDR peak files are identical:\t".join(' ', @$bed_files));
  }

  # Check output dir exists
  if(! -d $out_dir) {
    throw("Output directory does not exist:\t$out_dir");
  }

  my $idr_output_file = "${out_dir}/${output_prefix}-overlapped-peaks.new_idr.txt";
  
  my $cmd = "idr --idr-threshold $threshold --output-file $idr_output_file --plot --use-old-output-format --samples " . join(' ', @{$bed_files});

  print "Running:\n$cmd";
  run_system_cmd($cmd);
  
  my $png_file = $idr_output_file . '.png';
  if (! -e $png_file) {
    throw("Can't find expected png file $png_file for the idr!");
  }

  $cmd = "cat $idr_output_file | wc -l";
  my $num_peaks = run_backtick_cmd($cmd);

  # Warning: Parallelised appending to file!
  # This will also cause duplicate lines if the RunIDR jobs are rerun!
  $cmd = "echo -e \"IDR Comparison\tIDR Peaks\n$output_prefix\t$num_peaks\"".
    " >> ${out_dir}/${batch_name}-idr-stats.txt";
  run_system_cmd($cmd);

  return ($num_peaks, $png_file);
}

# TODO
# 1 Add more IDR based QC here
# 2 Rscript will currently submit and interactive job to farm if launched on a head node
#   else will run locally if already on a farm node

sub post_process_IDR{
  my $out_dir       = shift or throw('Must provide an out_dir argument');
  my $output_prefix = shift or throw('Must provide an output_prefix argument');
  my $idr_peaks     = shift;
  my $params        = shift || {};
  assert_ref($idr_peaks, 'ARRAY', 'idr_peaks');

  #add max_npairs, this may have to be a denominator of scalar(@$idr_comparison_files)?
  my ($idr_files, $npairs, $debug) =
    rearrange([ qw(idr_files npairs debug) ], %$params);

  my $max_peaks = int((sort {$a <=> $b} @{$idr_peaks})[-1]);#Take the highest!

  my $cmd = "echo -e \"IDR Max Peaks = $max_peaks\" >> ${out_dir}/${output_prefix}-idr-stats.txt";
  run_system_cmd($cmd);

  if($output_prefix && $idr_files){
    assert_ref($idr_files, 'ARRAY');#or just handle scalar arg here?
    $npairs ||= 1;

    if($npairs > scalar(@$idr_files)){
      throw("batch-consistency-plot.r will not handle npairs($npairs)".
        ' greater than the number of idr files specified('.scalar(@$idr_files).')');
    }

    my $script_path = which_path('batch-consistency-plot.r');
    $cmd = "Rscript $script_path 1 $out_dir/$output_prefix ".join(' ', map { $out_dir.'/'.$_ } @{$idr_files});
    run_system_cmd($cmd);

    warn "IDR plots are available here:\n\t$out_dir/${output_prefix}-plot.ps\n" if $debug;
  }
  elsif($output_prefix || $idr_files){
    throw('To run batch-consistency-plot.r both the output_prefix and idr_files parameters must be set');
  }

  warn "Final IDR Max peaks threshold:\t$max_peaks\n".
    "For more details, see ${out_dir}/${output_prefix}-idr-stats.txt\n" if $debug;

  return $max_peaks;
}

1;

