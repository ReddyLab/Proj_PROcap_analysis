# Merge Peak output from HOMER findPeaks 

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/mergePeak.txt \
    <<'EOF'
#!/bin/bash

### Set Directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_PEAK=$FD_OUT/out_findPeak
FD_MERGE=$FD_OUT/out_merge/mergePeak

FP_TSS_T00=${FD_PEAK}/tss_cap00.txt
FP_TSS_T15=${FD_PEAK}/tss_cap15.txt
FP_TSS_T60=${FD_PEAK}/tss_cap60.txt
FP_TSS_MERGE=${FD_MERGE}/tss_merge.txt
FP_TSS_SHIFT=${FD_MERGE}/tss_shift.txt
FP_TSS_SHIST=${FD_MERGE}/tss_shist.txt

### Merge TSS cluster positions from two separate experiments into 
### a single set of non-redundant TSS clusters
# mergePeaks exp1-rep1.tss.txt exp1-rep2.tss.txt exp2-rep1.tss.txt exp2-rep2.tss.txt -strand > merged.tss.txt
# mergePeaks exp1.tss.txt exp2.tss.txt -strand > merged.tss.txt
mergePeaks $FP_TSS_T00 $FP_TSS_T15 $FP_TSS_T60 -strand > $FP_TSS_MERGE

EOF

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/mergePeak_rna.txt \
    <<'EOF'
#!/bin/bash

### Set Directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_PEAK=$FD_OUT/out_findPeak
FD_MERGE=$FD_OUT/out_merge/mergePeak

FP_TSS_T00=${FD_PEAK}/tss_cap00_rna00.txt
FP_TSS_T15=${FD_PEAK}/tss_cap15_rna15.txt
FP_TSS_T60=${FD_PEAK}/tss_cap60_rna60.txt
FP_TSS_MERGE=${FD_MERGE}/tss_merge_rna.txt
FP_TSS_SHIFT=${FD_MERGE}/tss_shift_rna.txt
FP_TSS_SHIST=${FD_MERGE}/tss_shist_rna.txt

### Merge TSS cluster positions from two separate experiments into 
### a single set of non-redundant TSS clusters
# mergePeaks exp1-rep1.tss.txt exp1-rep2.tss.txt exp2-rep1.tss.txt exp2-rep2.tss.txt -strand > merged.tss.txt
# mergePeaks exp1.tss.txt exp2.tss.txt -strand > merged.tss.txt
mergePeaks $FP_TSS_T00 $FP_TSS_T15 $FP_TSS_T60 -strand > $FP_TSS_MERGE

EOF

-----

%%bash
wc -l /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_merge/mergePeak/tss_merge.txt
wc -l /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_merge/mergePeak/tss_merge_rna.txt

%%bash
head /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_merge/mergePeak/tss_merge.txt

%%bash
head /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_merge/mergePeak/tss_merge_rna.txt

-----

import numpy  as np
import pandas as pd
import os
import matplotlib.pyplot as plt

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_merge/mergePeak/tss_merge.txt"
dat = pd.read_csv(fpath, sep="\t")
dat.head(50)

dat.dropna()







