 Samon & DESeq2 okay
https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#alignment-based-abundance-estimation-methods
From the very beginning


head Eve_samples.txt
B12 EveB12_1 /media/tertiary/transcriptome_raw/Eulimnogammarus_verrucosus/Sample_EveB12_1/Sample_EveB12_1._R1.fastq.gz /media/tertiary/transcriptome_raw/Eulimnogammarus_verrucosus/Sample_EveB12_1/Sample_EveB12_1._R2.fastq.gz
B12 EveB12_2 /media/tertiary/transcriptome_raw/Eulimnogammarus_verrucosus/Sample_EveB12_2/Sample_EveB12_2._R1.fastq.gz /media/tertiary/transcriptome_raw/Eulimnogammarus_verrucosus/Sample_EveB12_2/Sample_EveB12_2._R2.fastq.gz
B12 EveB12_3 /media/tertiary/transcriptome_raw/Eulimnogammarus_verrucosus/Sample_EveB12_3/Sample_EveB12_3._R1.fastq.gz /media/tertiary/transcriptome_raw/Eulimnogammarus_verrucosus/Sample_EveB12_3/Sample_EveB12_3._R2.fastq.gz

## prepare reference for salmon
$appdir/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_verrucosus/EveBCdTP1_ani.fasta --est_method salmon --trinity_mode --prep_reference
#stopped_here
$appdir/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_cyaneus/EcyBCdTP1_ani.fasta --est_method salmon --trinity_mode --prep_reference
$appdir/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/G_lacustris/GlaBCdTP1_ani.fasta --est_method salmon --trinity_mode --prep_reference

## map with salmon
$appdir/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_verrucosus/EveBCdTP1_ani.fasta --seqType fq --samples_file ../Eve_samples.txt --est_method salmon --trinity_mode --output_dir . --thread_count 10 --SS_lib_type FR

$appdir/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_cyaneus/EcyBCdTP1_ani.fasta --seqType fq --samples_file ../Ecy_samples.txt --est_method salmon --trinity_mode --output_dir . --thread_count 10 --SS_lib_type FR
#stopped_here

$appdir/trinityrnaseq-v2.10.0/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/G_lacustris/GlaBCdTP1_ani.fasta --seqType fq --samples_file ../Gla_samples.txt --est_method salmon --trinity_mode --output_dir . --thread_count 10 --SS_lib_type FR
## matrix
$appdir/trinityrnaseq-v2.10.0/util/abundance_estimates_to_matrix.pl --est_method salmon ../salmon/Eve*/quant.sf --gene_trans_map none --name_sample_by_basedir --out_prefix Eve
$appdir/trinityrnaseq-v2.10.0/util/abundance_estimates_to_matrix.pl --est_method salmon ../salmon/Ecy*/quant.sf --gene_trans_map none --name_sample_by_basedir --out_prefix Ecy
$appdir/trinityrnaseq-v2.10.0/util/abundance_estimates_to_matrix.pl --est_method salmon ../salmon/Gla*/quant.sf --gene_trans_map none --name_sample_by_basedir --out_prefix Gla

##differential expression
##salmon
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Eve.isoform.counts.matrix --method DESeq2 --samples_file ../Eve_samples.txt --contrasts contrasts_parallel.txt --o eve.salmon.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Ecy.isoform.counts.matrix --method DESeq2 --samples_file ../Ecy_samples.txt --contrasts contrasts_parallel.txt --o ecy.salmon.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Gla.isoform.counts.matrix --method DESeq2 --samples_file ../Gla_samples.txt --contrasts contrasts_parallel.txt --o gla.salmon.de.dir
## Differential expression
## edgeR
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Eve.isoform.counts.matrix --method edgeR --samples_file ../Eve_samples_upd.txt --contrasts ../contrasts.txt --o eve.rsem.edgeR.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Ecy.isoform.counts.matrix --method edgeR --samples_file ../Ecy_samples.txt --contrasts ../contrasts.txt --o ecy.rsem.edgeR.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Gla.isoform.counts.matrix --method edgeR --samples_file ../Gla_samples.txt --contrasts ../contrasts.txt --o gla.rsem.edgeR.de.dir
## salmon
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Eve.isoform.counts.matrix --method DESeq2 --samples_file ../Eve_samples_upd.txt --contrasts ../contrasts.txt --o eve.salmon.DESeq2.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Ecy.isoform.counts.matrix --method DESeq2 --samples_file ../Ecy_samples.txt --contrasts ../contrasts.txt --o ecy.salmon.DESeq2.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Gla.isoform.counts.matrix --method DESeq2 --samples_file ../Gla_samples.txt --contrasts ../contrasts.txt --o gla.salmon.DESeq2.de.dir

## be careful! Deleting data
drozdovapb@proteomics-asus:/media/main/sandbox/drozdovapb/DE/11-diff-DE/bowtie_rsem$ rm */bowtie.bam