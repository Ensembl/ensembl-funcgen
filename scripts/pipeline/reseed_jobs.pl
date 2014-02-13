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



=head1 SYNOPSIS

  reseed_jobs.pl -url <String> (-job_ids <Int> ... | -logic_name <String>) \
    [ -ignore_retries -inlcude_ready -help -man ]

=head1 PARAMETERS

  Mandatory:
    -url              <String>    Hive DB URL
    -jobs_ids         <Int> ...   List of job IDs
    OR  
    -logic_name       <String>    Analysis logic_name
  Optional:
    
    WARNING - To be completely safe, all workers should be 
              stopped before using the following 'ignore' options.
    -ignore_retries   This will reseed jobs even when the max_retry_count has not been exceeded.
    -include_ready    This will also reseed jobs which are in the READY status
    
    -help             Prints a helpful message
    -man              Prints the man page

=head1 DESCRIPTION

This script appends to the input_id for a given FAILED job and reset the job
status to READY. This is primarily aimed at augmenting failing jobs which require
a slight update to the their input_ids to run successfully.

This is useful as it avoids having to seed new jobs, which can have various issues e.g.
  - Seeding new jobs may not be supported from the point of failure and reseeding upstream
    may cause more failures.
  - old FAILED jobs will remain in the system
  - new jobs may also fail if they are part of a previously successful
    branch down stream of the seed point, not necessarily part of the 
    FAILED branch.

Avoiding these FAILED job states makes managing the pipeline much easier, as there 
is only one status for a given process (i.e. we don't get a proliferation of input_ids
for the same process).

=cut

use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use Bio::EnsEMBL::Hive::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor;
use Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor;
use Bio::EnsEMBL::Hive::Utils      qw( destringify stringify );
use Bio::EnsEMBL::Utils::Exception qw( throw );

#Inject methods into the hive adaptor symbol tables/typeglobs

*Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor::fetch_dbID_by_data = sub {
  my ($self, $data) = @_;
 
  if(! defined $data){
    throw('Mandatory \'data\' argument not specified');  
  }
  
  stringify($data) if ref($data);
  my $sql = "SELECT analysis_data_id FROM analysis_data WHERE data = ?";
  my $sth = $self->prepare($sql);
  eval  { $sth->execute($data); };
  if($@){ throw("Failed to fetch_dbID_by data using:\t$sql ($data)\n$@"); }
  
  my ($data_id) = $sth->fetchrow_array();
  $sth->finish;
  
  return $data_id;
};

*Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor::update_if_needed = sub {
  my ($self, $data, $job_id) = @_;

  if(! (defined $data && defined $job_id)){
    throw('Mandatory \'data\' and/or \'job_id\' argument(s) not specified');  
  }

  #return 0 unless($data);
  #surely this never happens for sotre_if_needed as we were always 
  #testing for length > 255 in the caller?
  #if it did happen we would get an erroneous input_id of '_ext_input_analysis_data_id 0'
  #which doesn't exists and would probably die ungracefully
  #should probably throw or do the length test here instead?
  my $input_id;
  
  if(length($data) >= 255){   
    my $data_id = $self->fetch_dbID_by_data($data); 
  
    if($data_id) { #already stored      
      $input_id = '_ext_input_analysis_data_id '.$data_id;
    }
    else{  #else create a new one
      $input_id = '_ext_input_analysis_data_id '.$self->store($data);
    }
  }
  else{
    $input_id = $data; 
  }
   
  return $input_id;
};    
  
*Bio::EnsEMBL::Hive::DBSQL::AnalysisDataAdaptor::delete_by_dbID = sub {
  my ($self, $dbID) = @_;
  if(! defined $dbID){ throw('Mandatory dbID argument not met'); }
  
  eval  { $self->dbc->db_handle->do("DELETE from analysis_data where analysis_data_id=$dbID") };
  if($@){ throw("Failed to delete_by_dbID:\t$dbID\n$@"); }
  return;
};


#This is used in the update method below
*Bio::EnsEMBL::Hive::DBSQL::DBAdaptor::is_stored_and_valid = 
  \&Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor::is_stored_and_valid;

*Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor::update = sub {
  my ($self, $job) = @_;  
  $self->db->is_stored_and_valid('Bio::EnsEMBL::Hive::AnalysisJob', $job);

  #Get old job first and handle input_id updating
  my $job_id  = $job->dbID;
  my $old_job = $self->fetch_by_dbID($job_id);
  if(! defined $old_job){ throw("Could not fetch job with dbID:\t".$job_id); }  #Just in case
  
  my $input_id = $job->input_id;
  $input_id    = stringify($input_id) if ref($input_id);
  my $old_id   = $old_job->input_id;
  $old_id      = stringify($old_id)   if ref($old_id);
  my $adata_adaptor = $self->db->get_AnalysisDataAdaptor;
       
  if($old_id ne $input_id){  
    $input_id = $adata_adaptor->update_if_needed($input_id, $job_id);
    $job->input_id($input_id);
  }
  
  #Now actually do the update in the job table and return
  #We are updating the status here and below. Is this redundant
  #or is the designed to protect again doing this when running 
  #workers might pick up a half updated job?

  my $sql = "UPDATE job SET input_id='" . $job->input_id . "'".
    ",status='" . $job->status . "'".
    ",retry_count=" . $job->retry_count.
    ",semaphore_count=" . $job->semaphore_count;
    
  if (defined $job->semaphored_job_id){
    $sql .= ",semaphored_job_id=" . $job->semaphored_job_id;
  }
  
  $sql .= " WHERE job_id=" . $job->dbID;
  
  my $sth = $self->prepare($sql);
  eval  { $sth->execute }; 
  if($@){ 
    warn('Failed to update job(dbID='.$job->dbID.") using:\n\t".$sql."\n$@"); 
    return undef;
    #This is most likely due to an input_id_analysis key clash
    #do we want to throw here instead?   
  }
   
  
  $sth->finish;

  #This unless is questionable?
  #We might have manually updated the status to COMPLETED?
  #Also, why would we be updating the job if it is COMPLETED?
  #This should probably be removed or caught throw before the update above

  unless ($job->completed) {
    $self->update_status($job);
  }

  #Do any of the attributes in update_status need updating in the object?

  #Clean up analysis_data if required

  if($old_id ne $input_id){ #UPDATE job.input_id
    #Remove old analysis_data row if it is not linked to any other job
    my $old_ad_dbid = $adata_adaptor->fetch_dbID_by_data($old_id);
      
    if($old_ad_dbid){
      $old_id = '_ext_input_analysis_data_id '.$old_ad_dbid;     
      my @jobs = @{$self->fetch_all_by_input_id($old_id)};
      
      if(! scalar(@jobs)){ #No other jobs linked
        $adata_adaptor->delete_by_dbID($old_ad_dbid);        
      }
    }     
  }

  #Finally destringify input_id, so the returned job is useable
  $job->input_id(destringify($job->input_id));

  return $job;  
};



*Bio::EnsEMBL::Hive::DBSQL::AnalysisJobAdaptor::fetch_all_by_input_id = sub {
  my ($self, $input_id) = @_;
  
  if(! defined $input_id){
    throw('Mandatory \'input_id\' argument is not specified');  
  }
  
  $input_id = destringify($input_id) if ref($input_id);
  
  if(length($input_id) >= 255){
    my $data_id = $self->db->get_AnalysisDataAdaptor->fetch_dbID_by_data($input_id);
    
    if(defined $data_id){
      $input_id = '_ext_input_analysis_data_id '.$data_id;
    }
    else{
      undef $input_id;
      #as opposed to 0 which maybe a valid input_id??? 
    } 
  }
  
  return defined $input_id ? $self->generic_fetch("input_id=$input_id") : [];
};  



$| = 1;#for debug

my (%max_retries, @skipped_jobs);
my $updated_cnt     = 0;
my $skipped_cnt     = 0;
my $update_fail_cnt = 0;

sub main{
  my (@job_ids, $url, $hive_script_dir, $logic_name, $ignore_retries, $include_ready);
  my $append_string = '{recover => 1}';
  my @tmp_args = @ARGV;

  GetOptions (
              #Mandatory
              'url=s'             => \$url,
              'job_ids=i{,}'      => \@job_ids,
              'logic_name=s'      => \$logic_name,
     
              #Optional
              'append=s'           => \$append_string,
              'ignore_retries'     => \$ignore_retries, 
              'include_ready'      => \$include_ready,        
              'help'               => sub { pod2usage(-exitval => 0); }, 
              #removed ~? frm here as we don't want to exit with 0 for ?
              
              'man|m'              => sub { pod2usage(-exitval => 0, -verbose => 2); },
             ) or pod2usage(-exitval => 1, -message => "Specified parameters are:\t@tmp_args"); 
             
  #Do we want to catch the rest of ARGV here by using -- in the input and assign these to init_pipeline as extra args?
  
  if(defined $logic_name && (@job_ids) ){
    pod2usage(-exitval => 1,
              -message  => '-logic_name and -job_ids parameters are mutually exclusive, '.
                "please omit one:\n\t$0 @tmp_args");  
  }
  elsif(! ($logic_name || @job_ids)){
    pod2usage(-exitval => 1, -message  => 'You must specify either the -logic_name OR -job_ids parameters');
  }
  
  if(! defined $url){
    pod2usage(-exitval => 1, -message  => 'Mandatory -url parameter is missing.');
  }
  
  my $hive_dba   = Bio::EnsEMBL::Hive::DBSQL::DBAdaptor->new(-url => $url);
  my $job_a      = $hive_dba->get_AnalysisJobAdaptor;
  my $analysis_a = $hive_dba->get_AnalysisAdaptor;
  
  if(@job_ids){
    #Fetch all the jobs
    my @jobs = @{$job_a->fetch_all_by_dbID_list(\@job_ids)};
    
    foreach my $job(@jobs){
         
      if(! (($job->status eq 'FAILED') ||
          ($include_ready && ($job->status eq 'READY')) )){
        warn 'Skipping job_id '.$job->dbID.' as it is marked as '.$job->status."\n";
        $skipped_cnt ++;
        push @skipped_jobs, $job->dbID;
        next; #job       
      }
      
      &reseed_and_reset_job($job, $append_string, $job_a, $analysis_a, $ignore_retries);
    }
  }
  else{ #$logic_name
    #Fetch all failed jobs for given logic_name
    my $analysis = $analysis_a->fetch_by_logic_name($logic_name);
    my @jobs = @{$job_a->fetch_all_by_analysis_id_status($analysis->dbID, 'FAILED')};
    
    if($include_ready){
      push @jobs, @{$job_a->fetch_all_by_analysis_id_status($analysis->dbID, 'READY')};  
    }
     
    foreach my $job (@jobs){
      &reseed_and_reset_job($job, $append_string, $job_a, $analysis_a, $ignore_retries);
    }
  }

  print "\nUDPATE REPORT\n".
    "Update string:\t\t\t\t\t\t${append_string}\n".
    "Updated job count:\t\t\t\t\t$updated_cnt\n";
    
  print "Update failed count (likely due to duplicates):\t\t$update_fail_cnt\n" if $update_fail_cnt;
  print "Skipped job count (due to status or retry count):\t$skipped_cnt\n" if $skipped_cnt;
  print "Skipped/Failed update job dbIDs:\n@skipped_jobs\n" if @skipped_jobs;
  
  return;  
}# end of main

 


sub reseed_and_reset_job {
  my ($job, $append_string, $job_a, $analysis_a, $ignore_retries) = @_;
  
  if(! $ignore_retries){
    
    if(! exists $max_retries{$job->analysis_id}){
      $max_retries{$job->analysis_id} = $analysis_a->fetch_by_dbID($job->analysis_id)->max_retry_count;
    }
         
    if($job->retry_count < $max_retries{$job->analysis_id}){
      warn 'Skipping job_id '.$job->dbID.
        " as it has not yet exceeded the max retry count for it\'s analysis\n";  
      $skipped_cnt ++;
      push @skipped_jobs, $job->dbID;
      return; 
    }
  }
  
  
  $job->status('READY');
  $job->retry_count(1); 
   
  my $input_id = $job->input_id || {}; 
  #input_id is normally a string here
  $input_id    = destringify($input_id) if ! ref($input_id);  
  $job->input_id({ %{$input_id}, 
                   %{destringify( $append_string )} } );
 
  # DO THE UPDATE 
  if(! $job_a->update($job)){
    push @skipped_jobs, $job->dbID;
    $update_fail_cnt ++;    
  }
  else{
    $updated_cnt ++;
  } 
  
  return;
}# end of reseed_and_reset_job
  

main();

1;

