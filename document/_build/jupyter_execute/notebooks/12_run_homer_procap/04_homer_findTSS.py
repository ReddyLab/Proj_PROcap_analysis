# HOMER FindcsRNATSS

```

	findcsRNATSS.pl <csRNA tag directory> [options]

	Options:
		-o <prefix> 
		-i <csRNA input tag directory>
		-rna <RNAseq tag directory>
		-gtf <gtf file>
		-genome <genome>
		-cpu <#> (max CPUs)
		-minDistDiff <#> (default: 0.15)
		-defaultLog2Fold <#> (default: 1)
		-maxInputLog2Fold <#> (maximum log2 fold enrichment vs. input or RNA, default: 3)
		-maxRNALog2Fold <#> (maximum log2 fold enrichment vs. input or RNA, default: 3)
```

print(1)

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
findcsRNATSS.pl

test if I did not provide input

## FindcsRNATSS of time point 0 min

### Time point 00 min with RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 24G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findcsRNATSS_cap00_rna00.txt \
    <<'EOF'
#!/bin/bash

### Set Directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60
FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T00=$FD_OUT/tags/procap_t00
FD_TAG_RNA_T60=$FD_OUT/tags/rnaseq_t60
FD_TAG_RNA_T15=$FD_OUT/tags/rnaseq_t15
FD_TAG_RNA_T00=$FD_OUT/tags/rnaseq_t00


#findcsRNATSS.pl Exp1-csRNA-TagDir/ 
#    -o      outputPrefix 
#    -i      Exp1-input-TagDir/ 
#    -rna    Exp1-totalRNAseq-TagDir 
#    -gtf    genes.gtf 
#    -genome genome.fasta
PREFIX=$FD_OUT/out_findTSS/tss_cap00_rna00/out

findcsRNATSS.pl   $FD_TAG_CAP_T00/ \
        -rna      $FD_TAG_RNA_T00/ \
        -o        $PREFIX \
        -genome   hg38

EOF

### Time point 00 min without RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 24G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findcsRNATSS_cap00.txt \
    <<'EOF'
#!/bin/bash

### Set Directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60
FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T00=$FD_OUT/tags/procap_t00
FD_TAG_RNA_T60=$FD_OUT/tags/rnaseq_t60
FD_TAG_RNA_T15=$FD_OUT/tags/rnaseq_t15
FD_TAG_RNA_T00=$FD_OUT/tags/rnaseq_t00


#findcsRNATSS.pl Exp1-csRNA-TagDir/ 
#    -o      outputPrefix 
#    -i      Exp1-input-TagDir/ 
#    -rna    Exp1-totalRNAseq-TagDir 
#    -gtf    genes.gtf 
#    -genome genome.fasta
PREFIX=$FD_OUT/out_findTSS/tss_cap00/out

findcsRNATSS.pl   $FD_TAG_CAP_T00/ \
        -o        $PREFIX \
        -genome   hg38

EOF

## FindcsRNATSS of time point 15 min

### Time point 15 min with RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 24G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findcsRNATSS_cap15_rna15.txt \
    <<'EOF'
#!/bin/bash

### Set Directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60
FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T00=$FD_OUT/tags/procap_t00
FD_TAG_RNA_T60=$FD_OUT/tags/rnaseq_t60
FD_TAG_RNA_T15=$FD_OUT/tags/rnaseq_t15
FD_TAG_RNA_T00=$FD_OUT/tags/rnaseq_t00


#findcsRNATSS.pl Exp1-csRNA-TagDir/ 
#    -o      outputPrefix 
#    -i      Exp1-input-TagDir/ 
#    -rna    Exp1-totalRNAseq-TagDir 
#    -gtf    genes.gtf 
#    -genome genome.fasta
PREFIX=$FD_OUT/out_findTSS/tss_cap15_rna15/out

findcsRNATSS.pl   $FD_TAG_CAP_T15/ \
        -rna      $FD_TAG_RNA_T15/ \
        -o        $PREFIX \
        -genome   hg38

EOF

### Time point 15 min without RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 24G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findcsRNATSS_cap15.txt \
    <<'EOF'
#!/bin/bash

### Set Directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60
FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T00=$FD_OUT/tags/procap_t00
FD_TAG_RNA_T60=$FD_OUT/tags/rnaseq_t60
FD_TAG_RNA_T15=$FD_OUT/tags/rnaseq_t15
FD_TAG_RNA_T00=$FD_OUT/tags/rnaseq_t00


#findcsRNATSS.pl Exp1-csRNA-TagDir/ 
#    -o      outputPrefix 
#    -i      Exp1-input-TagDir/ 
#    -rna    Exp1-totalRNAseq-TagDir 
#    -gtf    genes.gtf 
#    -genome genome.fasta
PREFIX=$FD_OUT/out_findTSS/tss_cap15/out

findcsRNATSS.pl   $FD_TAG_CAP_T15/ \
        -o        $PREFIX \
        -genome   hg38

EOF

## FindcsRNATSS of time point 60 min

### Time point 60 min with RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 24G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findcsRNATSS_cap60_rna60.txt \
    <<'EOF'
#!/bin/bash

### Set Directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60
FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T00=$FD_OUT/tags/procap_t00
FD_TAG_RNA_T60=$FD_OUT/tags/rnaseq_t60
FD_TAG_RNA_T15=$FD_OUT/tags/rnaseq_t15
FD_TAG_RNA_T00=$FD_OUT/tags/rnaseq_t00


#findcsRNATSS.pl Exp1-csRNA-TagDir/ 
#    -o      outputPrefix 
#    -i      Exp1-input-TagDir/ 
#    -rna    Exp1-totalRNAseq-TagDir 
#    -gtf    genes.gtf 
#    -genome genome.fasta
PREFIX=$FD_OUT/out_findTSS/tss_cap60_rna60/out

findcsRNATSS.pl   $FD_TAG_CAP_T60/ \
        -rna      $FD_TAG_RNA_T60/ \
        -o        $PREFIX \
        -genome   hg38

EOF

### Time point 60 min without RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 24G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findcsRNATSS_cap60.txt \
    <<'EOF'
#!/bin/bash

### Set Directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60
FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T00=$FD_OUT/tags/procap_t00
FD_TAG_RNA_T60=$FD_OUT/tags/rnaseq_t60
FD_TAG_RNA_T15=$FD_OUT/tags/rnaseq_t15
FD_TAG_RNA_T00=$FD_OUT/tags/rnaseq_t00


#findcsRNATSS.pl Exp1-csRNA-TagDir/ 
#    -o      outputPrefix 
#    -i      Exp1-input-TagDir/ 
#    -rna    Exp1-totalRNAseq-TagDir 
#    -gtf    genes.gtf 
#    -genome genome.fasta
PREFIX=$FD_OUT/out_findTSS/tss_cap60/out

findcsRNATSS.pl   $FD_TAG_CAP_T60/ \
        -o        $PREFIX \
        -genome   hg38

EOF

-----

## Check the results

%%bash
for FP in $(ls /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findcsRNATSS_cap*.txt); do
    ls $FP
    #cat $FP | grep Difference
done

!head /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findcsRNATSS_cap60.txt

!tail /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findcsRNATSS_cap60.txt

-----

ls -d /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss_cap*

ls /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss_cap00/*tsv

ls /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss_cap00/*txt

-----

%%bash
ls /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss*/out.inputDistribution.txt

%%bash
head /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss*/out.inputDistribution.txt

<div class="alert-danger">
The problem I encounter here is that all out.inputDistribution.txt is empty. My guess: the reason might be: `Difference in distributions is too small`
</div>

%%bash
head /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss*/out.rnaDistribution.txt

-----

%%bash
head -3 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss*/out.stats.txt

cat /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss_t15v00_rna00_maxlen/out.stats.txt

%%bash
head -3 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss*/out.tss.txt

!head /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss_t60v00_rna00/out.tss.txt | cut -f11-15 | expand -t 20

-----

import numpy  as np
import pandas as pd
import os
import matplotlib.pyplot as plt

seq -> fastq -> trim/align -> BAM -> filtering (mapping quality)

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss_t00/out.tss.txt"
dat = pd.read_csv(fpath, sep="\t")
dat.head()

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss_t60v00_rna00/out.tss.txt"
dat = pd.read_csv(fpath, sep="\t")
dat.head()

Monday: 3pm

a lot of bioinformatics post analysis on Github 
BioSTAR
Google Group

Deeptools

dat['annotation'].unique()

dat.loc[1:3, ["csRNA", "csRNAinput", "rnaseq", "Log2Ratio vs. Input", "Log2Ratio vs. RNA"]]

np.log2(dat["csRNA"] / dat["csRNAinput"])

dat["Log2Ratio vs. Input"].hist(cumulative=False, density=1, bins=100)

dat["Log2Ratio vs. RNA"].hist(cumulative=True, density=1, bins=100)





