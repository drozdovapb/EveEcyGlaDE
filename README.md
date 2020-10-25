# EveEcyGlaDE
Transcriptome sequencing of two Lake Baikal endemic (Eulimnogammarus verrucosus and E. cyaneus) and one Holarctic (Gammarus lacustris) transcriptome species from the Baikal region subjected to different stressors: scripts. 


Logical steps
* `read_cleanup_assembly` for read cleanup and transcriptome assembly.
* `salmon.sh` for read alignment. 
* `de_functions_new+decontamination.R` for DE from (functions based on the DESeq2 package). 
* `GO_enrichment_from_DESeq_results.R` for GO enrichment from DESeq results.

Other folders:
* `CBPD_phenanthrene_scripts`: the code used for transcriptomic data analysis for the Shatilina et al. manuscript published in CBPD (2020). Details inside.
* `gradual_temperature_increase`: the code used for transcriptomic data analysis for the Jakob et al. manuscript (submitted). Details inside.
* `plotting_scripts`: most of the useful code for figures.
* `misc_scripts`: other useful scripts for data analysis.
* `technical_report`: extensive story about the early stages of data analysis (mostly cleanup and assembly).
