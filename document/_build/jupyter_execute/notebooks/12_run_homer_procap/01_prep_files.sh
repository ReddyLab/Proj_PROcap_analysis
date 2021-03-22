# Check the file paths for the data

- PRO-cap
    - 0, 15, 30 min Dex treatment
- RNA-seq
    - short term Dex treatment 0, 5, 10, 15,..., 25 min
    - long term Dex treatment 0, 
    - /data/reddylab/projects/GGR/results/rna_seq/checkpoints/iter0/accepted_samples.txt
    - /data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0/*/STAR_2pass_featurecounts/Aligned.out.sorted.bam
    - 

FP_GEN=/data/reddylab/Kuei/annotation/gencode.v34.annotation.gtf
FD_CAP=/data/reddylab/Kuei/Dex_ProCap
FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads

FP_BAM_CAP_T00m=$FD_CAP/new_files/A549_control_alignments/A549_untreated_merged.bam
FP_BAM_CAP_T15m=$FD_CAP/new_files/A549_15min_Dexamethasone_alignments/A549_15min_Dexamethasone_merged.bam
FP_BAM_CAP_T60m=$FD_CAP/new_files/A549_1hr_Dexamethasone_alignments/A549_1h_Dexamethasone_merged.bam

FP_BAM_RNA_T00m_rep1=$FD_RNA_ALIGN/iter_short/t00m_rep1.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T00m_rep2=$FD_RNA_ALIGN/iter_short/t00m_rep2.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T00m_rep3=$FD_RNA_ALIGN/iter_short/t00m_rep3.star2.Aligned.out.sorted.bam

FP_BAM_RNA_T15m_rep1=$FD_RNA_ALIGN/iter_short/t15m_rep1.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T15m_rep2=$FD_RNA_ALIGN/iter_short/t15m_rep2.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T15m_rep3=$FD_RNA_ALIGN/iter_short/t15m_rep3.star2.Aligned.out.sorted.bam

FP_BAM_RNA_T60m_rep1=$FD_RNA_ALIGN/iter0/t1_rep1/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T60m_rep2=$FD_RNA_ALIGN/iter0/t1_rep2/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T60m_rep3=$FD_RNA_ALIGN/iter0/t1_rep3/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T60m_rep4=$FD_RNA_ALIGN/iter0/t1_rep4/STAR_2pass_featurecounts/Aligned.out.sorted.bam

ls $FP_GEN
echo ==================
ls $FD_CAP
echo ==================
ls $FD_GGR
echo ==================
ls $FD_RNA_ALIGN

ls $FP_BAM_CAP_T00m
ls $FP_BAM_CAP_T15m
ls $FP_BAM_CAP_T60m

ls $FP_BAM_RNA_T00m_rep1
ls $FP_BAM_RNA_T00m_rep2
ls $FP_BAM_RNA_T00m_rep3

ls $FP_BAM_RNA_T15m_rep1
ls $FP_BAM_RNA_T15m_rep2
ls $FP_BAM_RNA_T15m_rep3



/data/reddylab/projects/ggr/data/rna_seq/expression/A549.rnaseq.dex.featurecounts.genes.TPM.selected_samples.txt 

/data/reddylab/projects/ggr/data/rna_seq/expression/A549.rnaseq.dex.rsem.transcripts.TPM.txt

FD_ggr=/data/reddylab/projects/ggr
FD_GGR_RES=/data/reddylab/projects/GGR
FD_CHECKPOINTS=/data/reddylab/projects/GGR/results/rna_seq/checkpoints

/data/reddylab/projects/ggr/data/rna_seq/expression
/data/reddylab/projects/GGR/results/rna_seq/checkpoints/iter0/accepted_samples.txt
/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0
/data/reddylab/projects/GGR/results/rna_seq/checkpoints/iter_short/accepted_samples.txt
/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter_short

head -1 /data/reddylab/projects/ggr/data/rna_seq/expression/A549.rnaseq.dex.featurecounts.genes.TPM.selected_samples.txt 

ls $FD_CHECKPOINTS/iter0/accepted_samples.txt
ls $FD_CHECKPOINTS/iter_short/accepted_samples.txt

FP_SAMPLES_LONG=$FD_GGR/results/rna_seq/checkpoints/iter0/accepted_samples.txt
FP_SAMPLES_SHORT=$FD_GGR/results/rna_seq/checkpoints/iter_short/accepted_samples.txt

FD_CAP=/data/reddylab/Kuei/Dex_ProCap
ls $FD_CAP/new_files

ls $FD_CAP/new_files/A549_control_alignments

ls -1 $FD_CAP/new_files/A549_15min_Dexamethasone_alignments

ls -1 $FD_CAP/new_files/A549_1hr_Dexamethasone_alignments

/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0/t00_rep1/STAR_2pass_featurecounts/Aligned.out.sorted.bam

/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter_short/t00m_rep1.star2.Aligned.out.sorted.bam

FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads
ls $FD_RNA_ALIGN

/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter_short/t00m_rep1.star2.Aligned.out.sorted.bam
/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter_short/t00m_rep2.star2.Aligned.out.sorted.bam
/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter_short/t00m_rep3.star2.Aligned.out.sorted.bam

/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter_short/t15m_rep1.star2.Aligned.out.sorted.bam
/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter_short/t15m_rep2.star2.Aligned.out.sorted.bam
/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter_short/t15m_rep3.star2.Aligned.out.sorted.bam

FD_CAP=/data/reddylab/Kuei/Dex_ProCap
FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads

FD_CAP=/data/reddylab/Kuei/Dex_ProCap
FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads

FP_BAM_CAP_T00m=$FD_CAP/new_files/A549_control_alignments/A549_untreated_merged.bam
FP_BAM_CAP_T15m=$FD_CAP/new_files/A549_15min_Dexamethasone_alignments/A549_15min_Dexamethasone_merged.bam
FP_BAM_CAP_T60m=$FD_CAP/new_files/A549_1hr_Dexamethasone_aFD_RNA_LONG_ALIGN=data/rna_seq/mapped_reads/iter0lignments/A549_1h_Dexamethasone_merged.bam

FP_BAM_RNA_T00m_rep1=$FD_RNA_ALIGN/iter_short/t00m_rep1.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T00m_rep2=$FD_RNA_ALIGN/iter_short/t00m_rep2.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T00m_rep3=$FD_RNA_ALIGN/iter_short/t00m_rep3.star2.Aligned.out.sorted.bam

FP_BAM_RNA_T15m_rep1=$FD_RNA_ALIGN/iter_short/t15m_rep1.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T15m_rep2=$FD_RNA_ALIGN/iter_short/t15m_rep2.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T15m_rep3=$FD_RNA_ALIGN/iter_short/t15m_rep3.star2.Aligned.out.sorted.bam

FP_BAM_RNA_T60m_rep1=$FD_RNA_ALIGN/iter0/t1_rep1/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T60m_rep2=$FD_RNA_ALIGN/iter0/t1_rep2/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T60m_rep3=$FD_RNA_ALIGN/iter0/t1_rep3/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T60m_rep4=$FD_RNA_ALIGN/iter0/t1_rep4/STAR_2pass_featurecounts/Aligned.out.sorted.bam

FP_BAM_RNA_T60m_rep1=$FD_RNA_ALIGN/

/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0/t00_rep1/STAR_2pass_featurecounts/Aligned.out.sorted.bam

ls /data/reddylab/Kuei/Dex_ProCap

/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0/t1_rep1/STAR_2pass_featurecounts/Aligned.out.sorted.bam
/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0/t1_rep2/STAR_2pass_featurecounts/Aligned.out.sorted.bam
/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0/t1_rep3/STAR_2pass_featurecounts/Aligned.out.sorted.bam
/data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0/t1_rep4/STAR_2pass_featurecounts/Aligned.out.sorted.bam

## RNA-seq of long term Dex treatment

#cat /data/reddylab/projects/GGR/results/rna_seq/checkpoints/iter0/accepted_samples.txt
cat $FP_SAMPLES_LONG

ls -h /data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0/*/STAR_2pass_featurecounts/Aligned.out.sorted.bam

## RNA-seq of short term Dex treatment

#cat /data/reddylab/projects/GGR/results/rna_seq/checkpoints/iter_short/accepted_samples.txt
cat $FD_CHECKPOINTS/iter_short/accepted_samples.txt

ls /data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter_short/*.Aligned.out.sorted.bam

## Genome file

```
URL=ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna
FNAME=Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget $URL/$FNAME -P /gpfs/fs1/data/reddylab/Kuei/annotation
```

ls -l /gpfs/fs1/data/reddylab/Kuei/annotation

FD_GEN=/gpfs/fs1/data/reddylab/Kuei/annotation
gunzip -c $FD_GEN/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz > $FD_GEN/Homo_sapiens.GRCh38.dna.primary_assembly.fa

ls -l /gpfs/fs1/data/reddylab/Kuei/annotation

FD_GEN=/gpfs/fs1/data/reddylab/Kuei/annotation

head $FD_GEN/Homo_sapiens.GRCh38.dna.primary_assembly.fa

head -180 $FD_GEN/Homo_sapiens.GRCh38.dna.primary_assembly.fa | tail -30

tail $FD_GEN/Homo_sapiens.GRCh38.dna.primary_assembly.fa







%%bash
URL=ftp://ftp.ensembl.org/pub/release-102/fasta/homo_sapiens/dna
FNAME=Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget $URL/$FNAME -P /gpfs/fs1/data/reddylab/Kuei/annotation

