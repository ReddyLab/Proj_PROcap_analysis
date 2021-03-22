# Get the count using HOMER annotatePeaks.pl

print(1)

%%bash
FD_WRK=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
ls $FD_WRK/out_merge

%%bash
FD_WRK=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
ls $FD_WRK/out_merge/mergeTSS

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 24G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/annotateTSS.txt \
    <<'EOF'
#!/bin/bash

### set directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer

FD_TAG_CAP_T00=$FD_OUT/tags/procap_t00
FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60

FD_MERGE=$FD_OUT/out_merge/mergeTSS
FP_TSS_MERGE=${FD_MERGE}/tss_merge.txt

FD_ANNOT=$FD_OUT/out_annotate/annoTSS
FP_TSS_COUNT_RAW=$FD_ANNOT/tss_count_raw.txt
FP_TSS_COUNT_RLG=$FD_ANNOT/tss_count_rlg.txt

### 7. Quantifying TSS cluster strength and calculating differential expression
### When using csRNA-seq, be sure to specify "-strand + -fragLength 1" to ensure that 
### the program quantifies the reads based on their 5' ends (initiation sites).
#annotatePeaks.pl merged.tss.txt hg38 -strand + -fragLength 1 -raw -d Exp1-tagDir/ Exp2-tagDir/ > counts.txt
annotatePeaks.pl \
    $FP_TSS_MERGE \
    hg38 \
    -strand + -fragLength 1 \
    -raw -d $FD_TAG_CAP_T00/ $FD_TAG_CAP_T15/ $FD_TAG_CAP_T60/ > $FP_TSS_COUNT_RAW

annotatePeaks.pl \
    $FP_TSS_MERGE \
    hg38 \
    -strand + -fragLength 1 \
    -rlog -d $FD_TAG_CAP_T00/ $FD_TAG_CAP_T15/ $FD_TAG_CAP_T60/ > $FP_TSS_COUNT_RLG

EOF

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 24G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/annotateTSS_rna.txt \
    <<'EOF'
#!/bin/bash

### set directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer

FD_TAG_CAP_T00=$FD_OUT/tags/procap_t00
FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60

FD_MERGE=$FD_OUT/out_merge/mergeTSS
FP_TSS_MERGE=${FD_MERGE}/tss_merge_rna.txt

FD_ANNOT=$FD_OUT/out_annotate/annoTSS_rna
FP_TSS_COUNT_RAW=$FD_ANNOT/tss_count_raw.txt
FP_TSS_COUNT_RLG=$FD_ANNOT/tss_count_rlg.txt

### 7. Quantifying TSS cluster strength and calculating differential expression
### When using csRNA-seq, be sure to specify "-strand + -fragLength 1" to ensure that 
### the program quantifies the reads based on their 5' ends (initiation sites).
#annotatePeaks.pl merged.tss.txt hg38 -strand + -fragLength 1 -raw -d Exp1-tagDir/ Exp2-tagDir/ > counts.txt
annotatePeaks.pl \
    $FP_TSS_MERGE \
    hg38 \
    -strand + -fragLength 1 \
    -raw -d $FD_TAG_CAP_T00/ $FD_TAG_CAP_T15/ $FD_TAG_CAP_T60/ > $FP_TSS_COUNT_RAW

annotatePeaks.pl \
    $FP_TSS_MERGE \
    hg38 \
    -strand + -fragLength 1 \
    -rlog -d $FD_TAG_CAP_T00/ $FD_TAG_CAP_T15/ $FD_TAG_CAP_T60/ > $FP_TSS_COUNT_RLG

EOF



%%bash
ls -l /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_annotate/annoTSS/tss_count_raw.txt
head -2 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_annotate/annoTSS/tss_count_raw.txt

%%bash
ls -l /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_annotate/annoTSS/tss_count_rlg.txt
head -2 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_annotate/annoTSS/tss_count_rlg.txt

-----

import numpy  as np
import pandas as pd
import os
import matplotlib.pyplot as plt

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_annotate/annoTSS/tss_count_raw.txt"
dat_raw = pd.read_csv(fpath, sep="\t")
dat_raw.head(3)

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_annotate/annoTSS/tss_count_rlg.txt"
dat_rlg = pd.read_csv(fpath, sep="\t")
dat_rlg.head(3)

print(dat_raw.shape)
print(dat_rlg.shape)

dat = dat_raw
dat[dat["Gene Name"] == "PER1"].head(1)

dat.columns.values

dat = dat_raw
idx = ["Chr", "Start", "End", "Strand", "Peak Score", "Detailed Annotation", "Distance to TSS", "Nearest Refseq",
       '/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/tags/procap_t00/ Tag Count in given bp (9384553.0 Total, normalization factor = 1, effective total = 10000000)',
       '/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/tags/procap_t15/ Tag Count in given bp (9252351.5 Total, normalization factor = 1, effective total = 10000000)',
       '/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/tags/procap_t60/ Tag Count in given bp (13741782.5 Total, normalization factor = 1, effective total = 10000000)']
dat[dat["Gene Name"] == "PER1"][idx]



dat = dat_raw
dat = dat[dat["Gene Name"] == "PER1"]
dat["Chr"].astype(str) + ":" + dat["Start"].astype(str) + "-" + dat["End"].astype(str) + "(" + dat["Strand"].astype(str) + ")"

dat = dat_raw
idx = ['/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/tags/procap_t00/ Tag Count in given bp (9384553.0 Total, normalization factor = 1, effective total = 10000000)',
       '/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/tags/procap_t15/ Tag Count in given bp (9252351.5 Total, normalization factor = 1, effective total = 10000000)',
       '/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/tags/procap_t60/ Tag Count in given bp (13741782.5 Total, normalization factor = 1, effective total = 10000000)']

dat = dat[dat["Gene Name"] == "PER1"][idx].astype(str)
dat.iloc[:,0] + "-" + dat.iloc[:,1] + "-" + dat.iloc[:,2]

find region: chr17:8,154,824-8,159,095

dat = dat_raw
dat[(dat["Start"] > 8154824) & (dat["End"] < 8159095)]





idx = ["Chr", "Start", "End", "Strand"]
[idx]
dat["Chr"].astype(str) + ":" + dat["Start"].astype(str) + "-" + dat["End"].astype(str) + "(" + dat["Strand"].astype(str) + ")"





%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/annotateTSS_rna00.txt \
    <<'EOF'
#!/bin/bash

FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_IN=$FD_OUT/out_findTSS
FD_OT=$FD_OUT/out_merge/mergeTSS_rna00


cd /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/

FP_TSS_T15=tss_t15v00.tss.txt
FP_TSS_T60=tss_t60v00.tss.txt
FP_TSS_MERGE=tss_merge.tss.txt
FP_TSS_COUNT_RAW=tss_count_raw.tss.txt
FP_TSS_COUNT_RLG=tss_count_rlg.tss.txt

### 7. Quantifying TSS cluster strength and calculating differential expression
#annotatePeaks.pl merged.tss.txt hg38 -strand + -fragLength 1 -raw -d Exp1-tagDir/ Exp2-tagDir/ > counts.txt
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer
FD_GEN=/gpfs/fs1/data/reddylab/Kuei/annotation

FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60

annotatePeaks.pl \
    $FP_TSS_MERGE \
    $FD_GEN/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -gtf $FD_GEN/gencode.v34.annotation.gtf \
    -strand + -fragLength 1 \
    -raw -d $FD_TAG_CAP_T15/ $FD_TAG_CAP_T60/ > $FP_TSS_COUNT_RAW

annotatePeaks.pl \
    $FP_TSS_MERGE \
    $FD_GEN/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -gtf $FD_GEN/gencode.v34.annotation.gtf \
    -strand + -fragLength 1 \
    -rlog -d $FD_TAG_CAP_T15/ $FD_TAG_CAP_T60/ > $FP_TSS_COUNT_RLG

EOF

!ls /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss*.tss.txt | cut -f1-5 | expand -t 15

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_raw.tss.txt | tail -4

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_raw.tss.txt | tail -4 | cut -f1-5

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_raw.tss.txt | cut -f2-5 | expand -t 10

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_raw.tss.txt | cut -f6-8 | expand -t 30

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_raw.tss.txt | cut -f9-11 | expand -t 25

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_raw.tss.txt | cut -f12-19 | expand -t 15

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_raw.tss.txt | cut -f20 | expand -t 25

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_raw.tss.txt | cut -f21 | expand -t 25

-----

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_rlg.tss.txt | cut -f2-5 | expand -t 10

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_rlg.tss.txt | cut -f6-8 | expand -t 30

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_rlg.tss.txt | cut -f20 | expand -t 25

!head -5 /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_rlg.tss.txt | cut -f21 | expand -t 25



!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_raw.tss.txt

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count_rlg.tss.txt

