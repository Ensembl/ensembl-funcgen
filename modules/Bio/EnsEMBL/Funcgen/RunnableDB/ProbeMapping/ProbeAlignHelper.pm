package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::ProbeAlignHelper;

use warnings;
use strict;
use Carp;
use Bio::EnsEMBL::Utils::Argument qw(rearrange);

sub new {
  my ($class, @args) = @_;

  my ( $funcgen_adaptor )
    = rearrange([
                 'funcgen_adaptor'
  ], @args );

  if (! defined $funcgen_adaptor) {
    confess('funcgen_adaptor is an obligatory parameter!');
  }  
  my $self = bless {},$class;  
  
  $self->funcgen_adaptor($funcgen_adaptor);
  
  return $self;
}

sub dbc {
  my $self = shift;
  return $self->{'funcgen_adaptor'}->dbc;
}

sub funcgen_adaptor {
  my ($self, $funcgen_adaptor) = @_;
  if ($funcgen_adaptor) {
    confess('Type error') unless ($funcgen_adaptor->isa('Bio::EnsEMBL::Funcgen::DBSQL::DBAdaptor'));    
    $self->{'funcgen_adaptor'} = $funcgen_adaptor;
  }
  return $self->{'funcgen_adaptor'};
}

sub analysis_adaptor {
  my $self = shift;
  if (! exists $self->{'analysis_adaptor'}) {
    $self->{'analysis_adaptor'} = $self->funcgen_adaptor->get_AnalysisAdaptor;
  }  
  return $self->{'analysis_adaptor'};
}

sub probe_adaptor {
  my $self = shift;
  if (! exists $self->{'probe_adaptor'}) {
    $self->{'probe_adaptor'} = $self->funcgen_adaptor->get_ProbeAdaptor;
  }  
  return $self->{'probe_adaptor'};
}

sub sql_helper {
  my $self = shift;
  
  if (! exists $self->{'sql_helper'}) {
    use Bio::EnsEMBL::Utils::SqlHelper;
    $self->{'sql_helper'} = Bio::EnsEMBL::Utils::SqlHelper->new(
      -DB_CONNECTION => $self->funcgen_adaptor->dbc
    );
  }  
  return $self->{'sql_helper'};
}

sub fetch_probe_and_analysis_for_probe_feature_linking {
    my $self                = shift;
    my $probe_seq_id        = shift;
    my $mapping_target_type = shift;
    
    my @probe = $self->fetch_probes_by_probe_seq_id($probe_seq_id);
    
    my @probe_analysis_pair;
    foreach my $current_probe (@probe) {
    
      my @analysis = $self->fetch_analyses_to_store_probes_features_as(
	$current_probe, 
	$mapping_target_type
      );
      
      push @probe_analysis_pair, {
	probe    => $current_probe,
	analysis => \@analysis,
      };
    }
    return @probe_analysis_pair;
}

sub fetch_probes_by_probe_seq_id {
  my $self = shift;
  my $probe_seq_id = shift;

  my $probe_id = $self->sql_helper->execute(
    -SQL      => 'select probe_id from probe where probe_seq_id=?',
    -PARAMS => [ $probe_seq_id ],
    # Callback flattens the result, so we get an arrayref of probe_ids from this.
    -CALLBACK => sub {
      my @row = @{ shift @_ };
      return $row[0];
    },
  );
  confess("No probe with probe_seq_id $probe_seq_id found!") unless (scalar @$probe_id);
  
  my @probe = map { $self->probe_adaptor->fetch_by_dbID($_) } @$probe_id;
  #confess('Type error') unless ($probe->isa('Bio::EnsEMBL::Funcgen::Probe'));
  
  return @probe;
}

sub fetch_analyses_to_store_probes_features_as {
  my $self                = shift;
  my $probe               = shift;
  my $mapping_target_type = shift;
  
  confess('Type error') unless ($probe->isa('Bio::EnsEMBL::Funcgen::Probe'));
  
  my $array = $probe->get_all_Arrays;
  my @analysis = map { 
    $self->fetch_analysis_to_store_array_probe_matches_as($_, $mapping_target_type) 
  } @$array;

  return @analysis;
}

sub fetch_analysis_to_store_array_probe_matches_as {
  my $self                = shift;
  my $array               = shift;
  my $mapping_target_type = shift;
  
  confess('Type error') unless ($array->isa('Bio::EnsEMBL::Funcgen::Array'));
  
  my $logic_name = array_to_logic_name($array, $mapping_target_type);      
  my $analysis = $self->analysis_adaptor->fetch_by_logic_name($logic_name);
  
  confess("No analysis with $logic_name found!") unless (defined $analysis);
  confess('Type error') unless ($analysis->isa('Bio::EnsEMBL::Analysis'));
  
  return $analysis;
}

sub array_to_logic_name {
  my $array               = shift;
  my $mapping_target_type = shift;      
  return array_class_to_logic_name($array->class, $mapping_target_type);
}

sub array_class_to_logic_name {
  my $array_class         = shift;
  my $mapping_target_type = shift;
  return "ProbeAlign_" . $array_class . "_" . $mapping_target_type;
}
  
1;

