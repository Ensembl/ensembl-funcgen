package Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::InsertAnalyses;

use strict;
use warnings;

use base ('Bio::EnsEMBL::Hive::Process');

sub run {   
    my $self = shift;
    
    my $tracking_dba_hash = $self->param('tracking_dba_hash');
    
    use Bio::EnsEMBL::Funcgen::RunnableDB::ProbeMapping::Utils qw( create_funcgen_adaptor );
    my $funcgen_adaptor = create_funcgen_adaptor($tracking_dba_hash);
    my $array_adaptor = $funcgen_adaptor->get_ArrayAdaptor;
    
    use List::MoreUtils qw( uniq );    
    my @array_class = uniq map { $_->class } @{$array_adaptor->generic_fetch};
  
    use Bio::EnsEMBL::DBSQL::DBAdaptor;
    use Bio::EnsEMBL::DBSQL::AnalysisAdaptor;
    
    my $dbc = Bio::EnsEMBL::DBSQL::DBConnection->new(%$tracking_dba_hash);
    my $analysis_adaptor = Bio::EnsEMBL::DBSQL::AnalysisAdaptor->new(
      Bio::EnsEMBL::DBSQL::DBAdaptor->new(-DBCONN => $dbc)
    );
    
    for my $current_array_class (@array_class) {
      for my $current_alignment_target ('transcript', 'genomic') {
	my $logic_name = join '_', 'ProbeAlign', $current_array_class, $current_alignment_target;
	print "$logic_name\n";
	
	my $obj = new Bio::EnsEMBL::Analysis(
	  -logic_name      => $logic_name,
	  -description     => 'some warm words about this analysis',
	  -display_label   => 'UniProt alignment',
	  -displayable     => '1',
	  -web_data        => 'web metadata info'
	);
	$analysis_adaptor->store($obj);
      }
    }    
}
1;
