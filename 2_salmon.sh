#example script for salmon run

cd /media/main/sandbox/drozdovapb/DE/4-salmon

/media/secondary/apps/Salmon-latest_linux_x86_64/bin/salmon index -t /media/main/transcriptome/assemblies/GlaBCdTP1_cor2.fasta -i GlaBCdTP1_cor2.index

mkdir gla2
mkdir gla2/gc

for sample in `ls /media/tertiary/transcriptome_raw/Sample_*/*R1.fastq.gz`; do base=$(basename  $sample _R1.fastq.gz); dir=$(dirname $sample); \
/media/secondary/apps/Salmon-latest_linux_x86_64/bin/salmon quant -i GlaBCdTP1_cor2.index/ -l ISF \
-1 $dir/${base}_R1.fastq.gz -2 $dir/${base}_R2.fastq.gz -p 11 -o gla2/gc/$base.quant.s --gcBias; done 

