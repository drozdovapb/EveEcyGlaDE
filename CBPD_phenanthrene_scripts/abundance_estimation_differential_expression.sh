## This code is based on the Trinity manual
## https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Transcript-Quantification#alignment-based-abundance-estimation-methods
## https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-Differential-Expression

=======================
## prepare reference
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_verrucosus/EveBCdTP1_ani.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_cyaneus/EcyBCdTP1_ani.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/G_lacustris/GlaBCdTP1_ani.fasta --est_method RSEM --aln_method bowtie --trinity_mode --prep_reference

## align (bowtie) & estimate abundance (rsem)
cd /media/main/sandbox/drozdovapb/DE/11-diff-DE/bowtie_rsem
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_verrucosus/EveBCdTP1_ani.fasta --seqType fq --samples_file ../Eve_samples_upd.txt --SS_lib_type FR  --est_method RSEM --aln_method bowtie --trinity_mode --output_dir . --thread_count 10
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_cyaneus/EcyBCdTP1_ani.fasta --seqType fq --samples_file ../Ecy_samples.txt --SS_lib_type FR  --est_method RSEM --aln_method bowtie --trinity_mode --output_dir . --thread_count 10
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/G_lacustris/GlaBCdTP1_ani.fasta --seqType fq --samples_file ../Gla_samples.txt --SS_lib_type FR  --est_method RSEM --aln_method bowtie --trinity_mode --output_dir . --thread_count 10
## matrix
$appdir/trinityrnaseq-2.8.5/util/abundance_estimates_to_matrix.pl --est_method RSEM  Eve*/*isoforms.results --gene_trans_map none --name_sample_by_basedir --out_prefix Eve 
$appdir/trinityrnaseq-2.8.5/util/abundance_estimates_to_matrix.pl --est_method RSEM  Ecy*/*isoforms.results --gene_trans_map none --name_sample_by_basedir --out_prefix Ecy
$appdir/trinityrnaseq-2.8.5/util/abundance_estimates_to_matrix.pl --est_method RSEM  Gla*/*isoforms.results --gene_trans_map none --name_sample_by_basedir --out_prefix Gla
##differential expression
## edgeR
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Eve.isoform.counts.matrix --method edgeR --samples_file ../Eve_samples_upd.txt --contrasts ../contrasts.txt --o eve.rsem.edgeR.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Ecy.isoform.counts.matrix --method edgeR --samples_file ../Ecy_samples.txt --contrasts ../contrasts.txt --o ecy.rsem.edgeR.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Gla.isoform.counts.matrix --method edgeR --samples_file ../Gla_samples.txt --contrasts ../contrasts.txt --o gla.rsem.edgeR.de.dir
## salmon
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Eve.isoform.counts.matrix --method DESeq2 --samples_file ../Eve_samples_upd.txt --contrasts ../contrasts.txt --o eve.salmon.DESeq2.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Ecy.isoform.counts.matrix --method DESeq2 --samples_file ../Ecy_samples.txt --contrasts ../contrasts.txt --o ecy.salmon.DESeq2.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Gla.isoform.counts.matrix --method DESeq2 --samples_file ../Gla_samples.txt --contrasts ../contrasts.txt --o gla.salmon.DESeq2.de.dir

## prepare reference for salmon
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_verrucosus/EveBCdTP1_ani.fasta --est_method salmon --trinity_mode --prep_reference
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_cyaneus/EcyBCdTP1_ani.fasta --est_method salmon --trinity_mode --prep_reference
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/G_lacustris/GlaBCdTP1_ani.fasta --est_method salmon --trinity_mode --prep_reference
## map with salmon
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_verrucosus/EveBCdTP1_ani.fasta --seqType fq --samples_file ../Eve_samples_upd.txt  --est_method salmon --trinity_mode --output_dir . --thread_count 10 --SS_lib_type FR
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/E_cyaneus/EcyBCdTP1_ani.fasta --seqType fq --samples_file ../Ecy_samples.txt  --est_method salmon --trinity_mode --output_dir . --thread_count 10 --SS_lib_type FR
$appdir/trinityrnaseq-2.8.5/util/align_and_estimate_abundance.pl --transcripts /media/main/transcriptome/assemblies/G_lacustris/GlaBCdTP1_ani.fasta --seqType fq --samples_file ../Gla_samples.txt  --est_method salmon --trinity_mode --output_dir . --thread_count 10 --SS_lib_type FR
## matrix
$appdir/trinityrnaseq-2.8.5/util/abundance_estimates_to_matrix.pl --est_method salmon  Eve*/quant.sf --gene_trans_map none --name_sample_by_basedir --out_prefix Eve
$appdir/trinityrnaseq-2.8.5/util/abundance_estimates_to_matrix.pl --est_method salmon  Ecy*/quant.sf --gene_trans_map none --name_sample_by_basedir --out_prefix Ecy
$appdir/trinityrnaseq-2.8.5/util/abundance_estimates_to_matrix.pl --est_method salmon  Gla*/quant.sf --gene_trans_map none --name_sample_by_basedir --out_prefix Gla



##differential expression

## contrasts.txt
##contrasts should be the other way round!
#PB03	B03
#Ph03	B03
#Ph03	PB03
#PB24	B24h
#Ph24	B24h
#Ph24	PB24



##salmon
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Eve.isoform.counts.matrix --method DESeq2 --samples_file ../Eve_samples_upd.txt --contrasts ../contrasts.txt --o eve.salmon..de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Ecy.isoform.counts.matrix --method DESeq2 --samples_file ../Ecy_samples.txt --contrasts ../contrasts.txt --o ecy.salmon.de.dir
$appdir/trinityrnaseq-2.8.5/Analysis/DifferentialExpression/run_DE_analysis.pl --matrix Gla.isoform.counts.matrix --method DESeq2 --samples_file ../Gla_samples.txt --contrasts ../contrasts.txt --o gla.salmon.de.dir
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


