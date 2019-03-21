use strict;
package clsRegisterMetadata;

sub new {

    my ($class, $args) = @_;

	my $self = {
		epi_aceesion => $args->{epi_accession} || '',
		accession  => $args->{accession} || '',
		experiment_accession  => $args->{experiment_accession} || '',
		epigenome  => $args->{epigenome} || '',
		feature_type  => $args->{feature_type} || '',
		biological_replicate  => $args->{biological_replicate} || '',
		new_bio_replicate  => $args->{new_bio_replicate} || '',
		technical_replicate  => $args->{technical_replicate} || '',
		new_tech_replicate  => $args->{new_tech_replicate} || '',
		gender  => $args->{gender} || '',
		md5_checksum  => $args->{md5_checksum} || '',
		local_url  => $args->{local_url} || '',
		analysis  => $args->{analysis} || '',
		experimental_group  => $args->{experimental_group} || '',
		assay_xrefs  => $args->{assay_xrefs} || '',
		ontology_xrefs  => $args->{ontology_xrefs} || '',
		xrefs  => $args->{xrefs} || '',
		epigenome_description  => $args->{epigenome_description} || '',
		control_id  => $args->{control_id} || '',
		paired  => $args->{paired} || '',
		paired_end_tag  => $args->{paired_end_tag} || '',
		read_length  => $args->{read_length} || '',
		multiple  => $args->{multiple} || '',
		paired_with  => $args->{paired_with} || '',
		download_url  => $args->{download_url} || '',
		info  => $args->{info} || '',
		derived_from  => $args->{derived_from} || ''
	};
	return bless $self, $class;

}

sub set_epi_accession {
	my ($self, $epi_accession) = @_;
    $self->{epi_accession} = $epi_accession;	
}

sub set_accession {
	my ($self, $accession) = @_;
    $self->{accession} = $accession;	
}

sub set_experiment_accession {
	my ($self, $expAccession) = @_;
    $self->{experiment_accession} = $expAccession;	
}

sub set_epigenome {
	my ($self, $epigenome) = @_;
    $self->{epigenome} = $epigenome;	
}

sub set_feature_type {
	my ($self, $feature_type) = @_;
    $self->{feature_type} = $feature_type;	
}

sub set_biological_replicate {
	my ($self, $biological_replicate) = @_;
    $self->{biological_replicate} = $biological_replicate;	
}

sub set_new_bio_replicate {
	my ($self, $new_bio_replicate) = @_;
    $self->{new_bio_replicate} = $new_bio_replicate;	
}

sub set_technical_replicate {
	my ($self, $technical_replicate) = @_;
    $self->{technical_replicate} = $technical_replicate;	
}

sub set_new_tech_replicate {
	my ($self, $new_tech_replicate) = @_;
    $self->{new_tech_replicate} = $new_tech_replicate;	
}

sub set_gender {
	my ($self, $gender) = @_;
    $self->{gender} = $gender;	
}

sub set_md5_checksum {
	my ($self, $md5_checksum) = @_;
    $self->{md5_checksum} = $md5_checksum;	
}

sub set_local_url {
	my ($self, $local_url) = @_;
    $self->{local_url} = $local_url;	
}

sub set_analysis {
	my ($self, $analysis) = @_;
    $self->{analysis} = $analysis;	
}

sub set_experimental_group {
	my ($self, $experimental_group) = @_;
    $self->{experimental_group} = $experimental_group;	
}

sub set_assay_xrefs {
	my ($self, $assay_xrefs) = @_;
    $self->{assay_xrefs} = $assay_xrefs;	
}

sub set_ontology_xrefs {
	my ($self, $ontology_xrefs) = @_;
    $self->{ontology_xrefs} = $ontology_xrefs;	
}

sub set_xrefs {
	my ($self, $xrefs) = @_;
    $self->{xrefs} = $xrefs;	
}

sub set_epigenome_description {
	my ($self, $epigenome_description) = @_;
    $self->{epigenome_description} = $epigenome_description;	
}

sub set_control_id {
	my ($self, $control_id) = @_;
    $self->{control_id} = $control_id;	
}

sub set_paired {
	my ($self, $paired) = @_;
    $self->{paired} = $paired;	
}

sub set_paired_end_tag {
	my ($self, $paired_end_tag) = @_;
    $self->{paired_end_tag} = $paired_end_tag;	
}

sub set_read_length {
	my ($self, $read_length) = @_;
    $self->{read_length} = $read_length;	
}

sub set_multiple {
	my ($self, $multiple) = @_;
    $self->{multiple} = $multiple;	
}

sub set_paired_with {
	my ($self, $paired_with) = @_;
    $self->{paired_with} = $paired_with;	
}

sub set_download_url {
	my ($self, $download_url) = @_;
    $self->{download_url} = $download_url;	
}

sub set_info {
	my ($self, $info) = @_;
    $self->{info} = $info;	
}

sub set_derived_from {
	my ($self, $derived_from) = @_;
    $self->{derived_from} = $derived_from;	
}

sub get {
	my ($self, $field) = @_;
    return $self->{$field};
}

sub csv_row {
	my $self = shift;
	my $csv_row = join ("\t", $self->{epi_accession}, $self->{accession}, $self->{experiment_accession}, $self->{epigenome}, $self->{feature_type}, $self->{biological_replicate}, $self->{new_bio_replicate}, $self->{technical_replicate}, $self->{new_tech_replicate}, $self->{gender}, 
	$self->{md5_checksum}, $self->{local_url}, $self->{analysis}, $self->{experimental_group}, $self->{assay_xrefs}, $self->{ontology_xrefs}, $self->{xrefs}, $self->{epigenome_description}, 
	$self->{control_id}, $self->{paired}, $self->{paired_end_tag}, $self->{read_length}, $self->{multiple}, $self->{paired_with}, $self->{download_url}, $self->{info}, $self->{derived_from});
	
	return $csv_row;
}	


1;