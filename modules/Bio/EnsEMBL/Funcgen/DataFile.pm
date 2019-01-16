=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2019] EMBL-European Bioinformatics Institute

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

=head1 DESCRIPTION

=cut

package Bio::EnsEMBL::Funcgen::DataFile;

use strict;

use Bio::EnsEMBL::Funcgen::GenericGetSetFunctionality qw(
  _generic_get_or_set
);

use Role::Tiny::With;
with 'Bio::EnsEMBL::Funcgen::GenericConstructor';

sub _constructor_parameters {
  return {
    dbID         => 'dbID',
    db           => 'db',
    table_id     => 'table_id',
    table_name   => 'table_name',
    path         => 'path',
    file_type    => 'file_type',
    md5sum       => 'md5sum',
  };
}

sub dbID         { return shift->_generic_get_or_set('dbID',         @_); }
sub db           { return shift->_generic_get_or_set('db',           @_); }
sub table_id     { return shift->_generic_get_or_set('table_id',     @_); }
sub table_name   { return shift->_generic_get_or_set('table_name',   @_); }

=head2 path

  Example    : my $file_name = $data_file->path;
               print "The file is here:$file_name\n";
  Description: Accessor for the path of the data file. This is the file name
               relative to the database file path.
  Returntype : String
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub path         { return shift->_generic_get_or_set('path',         @_); }

# Convenience method to get the path from the ftp site.
sub relative_ftp_site_path { 

    my $self = shift;
    my $path = $self->path;
    
    my $species = $self->db->db->species;
    
    use Bio::EnsEMBL::Funcgen::Utils::GoodUtils qw( create_species_assembly_path );
    my $species_assembly_path = create_species_assembly_path($species);
    
    #use File::Spec qw ( canonpath );
    return File::Spec->canonpath( $species_assembly_path . '/' . $path );
}

sub _delete_from_db {
    my $self = shift;
    $self->db->_delete($self);
}

=head2 file_type

  Example    : my $file_type = $data_file->file_type;
  Description: Accessor for the file type of the data file. It is one of 
               these:
                 - BAM,
                 - BAMCOV,
                 - BIGBED,
                 - BIGWIG,
                 - VCF,
                 - CRAM or
                 - DIR.
  Returntype : Enum(String)
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub file_type    { return shift->_generic_get_or_set('file_type',    @_); }

=head2 md5sum

  Example    : my $md5sum = $data_file->md5sum;
  Description: Accessor for the md5sum of the data file.
  Returntype : String
  Exceptions : None
  Caller     : general
  Status     : Stable

=cut
sub md5sum       { return shift->_generic_get_or_set('md5sum',       @_); }

=head2 summary_as_hash

  Example       : $summary = $peak_calling->summary_as_hash;
  Description   : Returns summary in a hash reference.
  Returns       : Hashref of descriptive strings
  Status        : Intended for internal use (REST)

=cut

sub summary_as_hash {
  my $self   = shift;
  
  return {
    file_type => $self->file_type,
    path      => $self->path,
    relative_ftp_site_path => $self->relative_ftp_site_path,
  };
}

1;
