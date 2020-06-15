Peak Calling
============

Peak calling is a computational method used to identify areas in a genome that have been enriched with aligned reads as a consequence of performing a ChIP-sequencing or MeDIP-seq experiment. 
These areas are those where a protein interacts with DNA.

Read alignment
^^^^^^^^^^^^^^
To maintain a standardised peak calling methodology, we start our analyses with raw reads from each experiment. We align reads (replicates are pooled) to the genome using BWA (with default parameters). All matches to the mitochondrial DNA are filtered out to avoid alignment anomalies due to similarities with autosomal regions.

Peak calling
^^^^^^^^^^^^
Peak calling is performed using two algorithms:

SWEMBL (S. Wilder et al., in preparation):
When replicates are available, permissive settings (-f 150 -R 0.0005 -d 150) are used, of which we choose the best peaks, so as to fit an IDR score of 0.01. This number of peaks is linearly adjusted to take into account replicate pairs with different numbers of peaks.
Otherwise, different stringent settings are used for very tight peaks (e.g. DNAse1) (-f 150 -R 0.0025) and larger ones (-f 150 -R 0.015).

CCAT (Xu et al., 2010)
The algorithm described by Xu et al, is specifically designed for peak calling of broad features, such as H3K36me3. The parameters were set empirically to (fragmentSize 200; slidingWinSize 150; movingStep 20; isStrandSensitiveMode 0; minCount 10; outputNum 100000; randomSeed 123456; minScore 4.0; bootstrapPass 50).
Peak filtering

The ENCODE DAC Blacklist regions are then used to filter the resultant peaks. These regions are stored in the core database as 'encode_excluded' MiscFeatures.
