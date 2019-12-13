package Bio::EnsEMBL::Funcgen::Utils::FtpUtils;
=head1 LICENSE

Copyright [1999-2015] Wellcome Trust Sanger Institute and the EMBL-European Bioinformatics Institute
Copyright [2016-2020] EMBL-European Bioinformatics Institute

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

     http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

=cut


=head1 CONTACT

  Please email comments or questions to the public Ensembl
  developers list at <http://lists.ensembl.org/mailman/listinfo/dev>.

  Questions may also be sent to the Ensembl help desk at
  <http://www.ensembl.org/Help/Contact>.

=cut
use strict;
use warnings;
use base qw( Exporter );
use vars qw( @EXPORT_OK );
use Carp;

@EXPORT_OK = qw(
  create_file_handle
  check_file_content_consistent
  get_all_regulatory_activity_file_infos_from_ftp
  fetch_epigenome_production_names
  get_expected_regulatory_activity_directory_names
);

sub create_epigenome_cache {

    my $funcgen_adaptor = shift;
    my $epigenome_adaptor = $funcgen_adaptor->get_EpigenomeAdaptor;
    my $all_epigenomes = $epigenome_adaptor->fetch_all;
    
    my $epigenome_cache = { 
        map  { $_->production_name => $_ } 
        grep { $_->production_name } 
        @$all_epigenomes 
    };
    return $epigenome_cache;
}

sub fetch_epigenome_production_names {

    my $dbc = shift;

    my $sql = 'select production_name from epigenome join regulatory_build_epigenome using (epigenome_id)';

    my $sth = $dbc->prepare($sql);
    $sth->execute;
    my $result = $sth->fetchall_arrayref;
    my $list_of_names = [ map { $_->[0] } @$result ];
    
    return $list_of_names;
}

sub build_regulatory_activity_directory_name {

    my $param = shift;
    
    my $dbc     = $param->{dbc};
    my $ftp_dir = $param->{ftp_dir};
    my $epigenome_production_name = $param->{epigenome_production_name};

    my $regulatory_activity_directory_name
        = $ftp_dir 
            . '/RegulatoryFeatureActivity/' 
            . $epigenome_production_name
    ;
    return $regulatory_activity_directory_name;
}

sub get_expected_regulatory_activity_directory_names {

    my $param = shift;
    
    my $dbc     = $param->{dbc};
    my $ftp_dir = $param->{ftp_dir};

    my $epigenome_production_names = fetch_epigenome_production_names($dbc);

    my $expected_regulatory_activity_directory_names = [
        map
            {
                {
                    dir_name 
                        => build_regulatory_activity_directory_name({
                            dbc                       => $dbc,
                            ftp_dir                   => $ftp_dir,
                            epigenome_production_name => $_,
                        }),
                    epigenome_production_name => $_
                }
            } @$epigenome_production_names
    ];
    return $expected_regulatory_activity_directory_names;
}

sub get_all_regulatory_activity_file_infos_from_ftp {

    my $param = shift;
    
    my $dbc     = $param->{dbc};
    my $ftp_dir = $param->{ftp_dir};
    my $species = $param->{species};
    
    if (! defined $species) {
      confess("Species has not been defined");
    }

    my $epigenome_production_names = fetch_epigenome_production_names($dbc);
    my $regulatory_activity_directory_names 
        = get_expected_regulatory_activity_directory_names($param);

    my $expected_regulatory_activity_directory_names = [
        map
            {
                $ftp_dir . '/RegulatoryFeatureActivity/' . $_->{dir_name}
            } @$regulatory_activity_directory_names
    ];
    
    my @regulatory_activity_file_names;
    foreach my $directory (@$regulatory_activity_directory_names) {
    
        my $dir_name = $directory->{dir_name};
        my $epigenome_production_name = $directory->{epigenome_production_name};
    
        # E.g.:
        # homo_sapiens.GRCh37.Fetal_Adrenal_Gland.Regulatory_Build.regulatory_activity.20180503.gff.gz
        my $pattern 
            = $dir_name
                . '/' . $species 
                . '.*.'
                . $epigenome_production_name
                . '.Regulatory_Build.regulatory_activity.*.gff.gz'
        ;
#         die(
#           "pattern = $pattern\n"
#           . "species = $species\n"
#           . "dir_name = $dir_name\n"
#         );
        my @expected_files = glob $pattern;
        
        if (@expected_files > 1) {
          confess("Found more than one file in directory! " . Dumper(\@expected_files));
        }
        
        my $expected_file  = $expected_files[0];
        push 
            @regulatory_activity_file_names, 
            {
                file_name                 => $expected_file,
                pattern                   => $pattern,
                epigenome_production_name => $epigenome_production_name,
            }
        ;
    }
    return \@regulatory_activity_file_names;
}

sub check_file_content_consistent {

    my $param = shift;
    
    my $funcgen_adaptor = $param->{funcgen_adaptor};
    my $ftp_dir = $param->{ftp_dir};
    my $species = $param->{species};
    my $check_these = $param->{check_these};
    
    my $dbc = $funcgen_adaptor->dbc;
    
    my $epigenome_cache = create_epigenome_cache($funcgen_adaptor);
    
    my $regulatory_feature_adaptor = $funcgen_adaptor->get_RegulatoryFeatureAdaptor;
    
    foreach my $regulatory_activity_file_name (@$check_these) {
    
        print("     Checking " . $regulatory_activity_file_name->{file_name} . "\n");
    
        my $file_name = $regulatory_activity_file_name->{file_name};
        my $epigenome_production_name
            = $regulatory_activity_file_name->{epigenome_production_name};
        
        my $epigenome = $epigenome_cache->{$epigenome_production_name};
        
        if (! defined $epigenome) {
            confess("Couldn't find epigenome $epigenome_production_name");
        }
        
        open my $fh, "zcat $file_name |" or die($file_name);
        
        my $number_of_lines_checked      =    0;
        my $max_number_of_lines_to_check = 1000;
        
        while (my $current_line = <$fh>) {
            chomp($current_line);
            
            $number_of_lines_checked++;
            
            if ($number_of_lines_checked >= $max_number_of_lines_to_check) {
              last;
            }
            
            my @f = split "\t", $current_line;
            
            my $attributes = $f[8];
            
            # Attributes look like this:
            #
            # activity=NA;bound_end=102119230;bound_start=102118695;description=TF binding site na in Lung;epigenome=Lung;feature_type=TF binding site;regulatory_feature_stable_id=ENSR00000368862
            #
            my %attribute_hash = map { split /=/ } split ';', $attributes;
            
            if (! exists $attribute_hash{regulatory_feature_stable_id}) {
                die("Couldn't find stable id in: $attributes");
            }
            my $stable_id = $attribute_hash{regulatory_feature_stable_id};
            
            my $regulatory_feature_from_file = {
                seq_region_name => $f[0],
                start           => $f[3],
                end             => $f[4],
                %attribute_hash
            };
            
            my $regulatory_feature_from_db = $regulatory_feature_adaptor->fetch_by_stable_id($regulatory_feature_from_file->{regulatory_feature_stable_id});
            
            my $start_ok        = $regulatory_feature_from_file->{start}           == $regulatory_feature_from_db->{start};
            my $end_ok          = $regulatory_feature_from_file->{end}             == $regulatory_feature_from_db->{end};
            my $feature_type_ok = $regulatory_feature_from_file->{feature_type}    eq $regulatory_feature_from_db->feature_type->name;
            my $seq_region_ok   = $regulatory_feature_from_file->{seq_region_name} eq $regulatory_feature_from_db->slice->seq_region_name;
            
            my $regulatory_activity_from_db   = $regulatory_feature_from_db->regulatory_activity_for_epigenome($epigenome)->activity;
            my $regulatory_activity_from_file = $regulatory_feature_from_file->{activity};
            
            my $activity_ok = $regulatory_feature_from_file->{activity} eq $regulatory_activity_from_db;
            
            if (! $seq_region_ok) {
                confess("Seq regions don't agree in $file_name! File: ".$regulatory_feature_from_file->{seq_region_name}." DB:" . $regulatory_feature_from_db->slice->name);
            }
            if (! $activity_ok) {
                confess("Regulatory activities don't agree: File: $regulatory_activity_from_file DB: $regulatory_activity_from_db!");
            }
            if (! $feature_type_ok) {
                confess("Feature types don't agree in $file_name File: ".$regulatory_feature_from_file->{feature_type}." DB: ".$regulatory_feature_from_db->feature_type->name."!");
            }
            if (! ($start_ok && $end_ok)) {
                confess("Coordinates don't agree in $file_name!");
            }
        }
        
        close($fh);
        print("Ok - File content for $file_name\n");
    }
    return;
}

sub create_file_handle {

  my $output_file = shift;

  my $output_fh;
  if ($output_file) {

    use File::Basename;
    my $ftp_dir = dirname($output_file);

    use File::Path qw(make_path);
    make_path($ftp_dir);

    use IO::File;
    $output_fh = IO::File->new(">$output_file");
  } else {
    $output_fh = *STDOUT;
  }
  return $output_fh;
}

1;