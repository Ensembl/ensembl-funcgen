
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

Bio::EnsEMBL::Funcgen::Hive::RunAligner

=head1 DESCRIPTION



=cut

package Bio::EnsEMBL::Funcgen::Hive::RunAligner;

use warnings;
use strict;
use Bio::EnsEMBL::Utils::Exception qw( throw );
use base qw( Bio::EnsEMBL::Funcgen::Hive::BaseDB );

#todo make this independant of ResultSet?


sub fetch_input {   # fetch parameters...
  my $self = shift;
  $self->SUPER::fetch_input();
  my $rset       = $self->fetch_Set_input('ResultSet');
  my $fastq_file = $self->get_param_method('fastq_file', 'required');



  #We need a list of Aligner specific param requirements which are not specified
  #in $analysis->parameters
  #bwa_index_root, gender, assembly
  #How are we going to genericise these?
  #We can't specify ResultSet->cell_type->gender
  #and we don't want to tie the aligner to use of a ResultSet
  #Can we make these co-optional in the constructor?
  #i.e. we can pass gender or a ResultSet?
  



  
  
  $self->get_param_method('output_dir', 'required'); 
  #This should have been set to a 'work' dir in MergeChunkResultSetFastqs
  
  #Do we even need this? The fastq chunks will already be in a work dir?
  
  
  
  my $analysis     = $rset->analysis;
  my $align_module = $self->validate_package_from_path($analysis->module);
  my $aligner      = $analysis->program || 
    throw('Aligner analysis cannot have an undef program attribute:'.$analysis->logic_name);
  $self->set_param_method('aligner', $aligner); 

  
  #my %aligner_reqs;  
  #foreach my $req(@{$align_module->required_parameters}){
   
  #  if(! $self->can($req)){
  #    throw("Could not call method for $aligner aligner creation:\t".$req);
  #  }
  #  else{
  #    $aligner_reqs{"-${req}"} = $self->$req;
  #  }    
  #}

  #validate program_file isn't a path?
  #would redefine bwa bin in the config for this analysis
  #this will change the bin_dir for everything else too, but we don't use bin_dir for anything else here 
  my $pfile_path = ( defined $self->bin_dir ) ? 
    $self->bin_dir.'/'.$analysis->program_file : $analysis->program_file;

 
  my $ref_fasta;
  
  #indexed_ref_fasta is batch flown (assuming all aligners will have indexed the ref fasta) 
  #indexes should be in the same locations, and this fasta will actually be a soft link 
  #to prevent  reundant/out of sync ref fasta files
  
  
   
  if(! defined $self->param_silent('indexed_ref_fasta')){ #This is batch flown
    my $gender         = $rset->cell_type->gender || 'male';
    my $species        = $self->species; 
  
    #TODO: Maybe check if the index file is really there? Eventually the bwa output will tell you though
  
    #index file suffix may change between aligners
    #best to pass just the target file, index root dir, species, gender
    #and let the Aligner construct the appropriate index file
  
  
    $ref_fasta = join('/', ($self->param_required('data_root_dir'),
                               $self->program.'_indexes',
                               $species,
                               $species.'_'.$gender.'_'.$self->assembly.'_unmasked.fasta'));
  }

  my $align_runnable = $align_module->new
   (
    -program_file      => $pfile_path,
    -parameters        => $analysis->parameters,
    -query_file        => $fastq_file,
    -reference_file    => $ref_fasta,
    #-output_dir        =>   
   # %aligner_reqs
   );
   
  $self->set_param_method('aligner', $align_runnable); 
  
  return;
}


sub run {
  shift->aligner->run;
  return;
}


sub write_output {  # Nothing to write
  #my $self = shift;
  return;
}



1;
