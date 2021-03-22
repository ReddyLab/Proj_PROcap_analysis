# Merge TSS output from HOMER findcsRNATSS 

%%bash
ls -1 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer

%%bash
ls -1 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_merge

%%bash
ls -1 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss_cap*/out.tss.txt

%%bash
head -3 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findTSS/tss_cap*/out.tss.txt

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/mergeTSS.txt \
    <<'EOF'
#!/bin/bash

### Set Directories
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_TSS=$FD_OUT/out_findTSS
FD_MERGE=$FD_OUT/out_merge/mergeTSS

FP_TSS_T00=${FD_TSS}/tss_cap00/out.tss.txt
FP_TSS_T15=${FD_TSS}/tss_cap15/out.tss.txt
FP_TSS_T60=${FD_TSS}/tss_cap60/out.tss.txt
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
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/mergeTSS_rna00.txt \
    <<'EOF'
#!/bin/bash

# mergePeaks exp1-rep1.tss.txt exp1-rep2.tss.txt exp2-rep1.tss.txt exp2-rep2.tss.txt -strand > merged.tss.txt
# mergePeaks exp1.tss.txt exp2.tss.txt -strand > merged.tss.txt

#FP_TSS1=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt
#FP_TSS2=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t60v00.tss.txt
#FP_TSSM=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/merge

FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_IN=$FD_OUT/out_findTSS
FD_OT=$FD_OUT/out_merge/mergeTSS_rna00

FP_TSS_T15=${FD_IN}/tss_t15v00_rna00/out.tss.txt
FP_TSS_T60=${FD_IN}/tss_t60v00_rna00/out.tss.txt
FP_TSS_MERGE=${FD_OT}/tss_merge.txt
FP_TSS_SHIFT=${FD_OT}/tss_shift.txt
FP_TSS_SHIST=${FD_OT}/tss_shist.txt

### Merge TSS cluster positions from two separate experiments into 
### a single set of non-redundant TSS clusters
mergePeaks $FP_TSS_T15 $FP_TSS_T60 -strand > $FP_TSS_MERGE

#annotatePeaks.pl <tss/peak/BED file> <genome> -p <tss/peak/BED file2> -pdist2 -strand + > outputFile.txt
#annotatePeaks.pl exp1.tss.txt hg38 -p exp2.tss.txt -pdist2 -strand + > output.txt
annotatePeaks.pl \
    $FP_TSS_T15 \
    hg38 \
    -p $FP_TSS_T60 \
    -pdist2 -strand + > $FP_TSS_SHIFT

#annotatePeaks.pl <tss/peak/BED file> <genome> -p <tss/peak/BED file2> -size <#> -hist <#> -strand + > outputFile.txt
#annotatePeaks.pl exp1.tss.txt hg38 -p exp2.tss.txt -size 1000 -hist 10 -strand + > outputFile.txt
annotatePeaks.pl \
    $FP_TSS_T15 \
    hg38 \
    -p $FP_TSS_T60 \
    -size 1000 -hist 10 -strand + > $FP_TSS_SHIST
EOF

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/mergeTSS_rna00_maxlen.txt \
    <<'EOF'
#!/bin/bash

# mergePeaks exp1-rep1.tss.txt exp1-rep2.tss.txt exp2-rep1.tss.txt exp2-rep2.tss.txt -strand > merged.tss.txt
# mergePeaks exp1.tss.txt exp2.tss.txt -strand > merged.tss.txt

#FP_TSS1=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt
#FP_TSS2=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t60v00.tss.txt
#FP_TSSM=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/merge

FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer
FD_IN=$FD_OUT/out_findTSS
FD_OT=$FD_OUT/out_merge/mergeTSS_rna00_maxlen

FP_TSS_T15=${FD_IN}/tss_t15v00_rna00_maxlen/out.tss.txt
FP_TSS_T60=${FD_IN}/tss_t60v00_rna00_maxlen/out.tss.txt
FP_TSS_MERGE=${FD_OT}/tss_merge.txt
FP_TSS_SHIFT=${FD_OT}/tss_shift.txt
FP_TSS_SHIST=${FD_OT}/tss_shist.txt

### Merge TSS cluster positions from two separate experiments into 
### a single set of non-redundant TSS clusters
mergePeaks $FP_TSS_T15 $FP_TSS_T60 -strand > $FP_TSS_MERGE

#annotatePeaks.pl <tss/peak/BED file> <genome> -p <tss/peak/BED file2> -pdist2 -strand + > outputFile.txt
#annotatePeaks.pl exp1.tss.txt hg38 -p exp2.tss.txt -pdist2 -strand + > output.txt
annotatePeaks.pl \
    $FP_TSS_T15 \
    hg38 \
    -p $FP_TSS_T60 \
    -pdist2 -strand + > $FP_TSS_SHIFT

#annotatePeaks.pl <tss/peak/BED file> <genome> -p <tss/peak/BED file2> -size <#> -hist <#> -strand + > outputFile.txt
#annotatePeaks.pl exp1.tss.txt hg38 -p exp2.tss.txt -size 1000 -hist 10 -strand + > outputFile.txt
annotatePeaks.pl \
    $FP_TSS_T15 \
    hg38 \
    -p $FP_TSS_T60 \
    -size 1000 -hist 10 -strand + > $FP_TSS_SHIST
EOF



!ls /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss*.tss.txt | cut -f1-5 | expand -t 15

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_merge.tss.txt

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_shift.tss.txt

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_shist.tss.txt

import numpy  as np
import pandas as pd
import os
import matplotlib.pyplot as plt

ls /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_merge/mergeTSS_rna00_maxlen

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_merge/mergeTSS_rna00_maxlen/tss_shift.txt"
dat = pd.read_csv(fpath, sep="\t")
dat.head()

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_merge/mergeTSS_rna00_maxlen/tss_shist.txt"
dat = pd.read_csv(fpath, sep="\t")
dat.head()

x = dat.iloc[:,0]
y = dat.iloc[:,1]
plt.plot(x, y)



%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
export PATH=/data/reddylab/software/homer/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/mergeTSS.txt \
    <<'EOF'
#!/bin/bash

# mergePeaks exp1-rep1.tss.txt exp1-rep2.tss.txt exp2-rep1.tss.txt exp2-rep2.tss.txt -strand > merged.tss.txt
# mergePeaks exp1.tss.txt exp2.tss.txt -strand > merged.tss.txt

#FP_TSS1=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt
#FP_TSS2=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t60v00.tss.txt
#FP_TSSM=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/merge

cd /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/

FP_TSS_T15=tss_t15v00.tss.txt
FP_TSS_T60=tss_t60v00.tss.txt
FP_TSS_MERGE=tss_merge.tss.txt
FP_TSS_COUNT=tss_count.tss.txt

mergePeaks $FP_TSS_T15 $FP_TSS_T60 -strand > $FP_TSS_MERGE

annotatePeaks.pl exp1.tss.txt hg38 -p exp2.tss.txt -size 1000 -hist 10 -strand + > outputFile.txt

#annotatePeaks.pl merged.tss.txt hg38 -strand + -fragLength 1 -raw -d Exp1-tagDir/ Exp2-tagDir/ > counts.txt
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer
FD_GEN=/gpfs/fs1/data/reddylab/Kuei/annotation

FD_TAG_CAP_T15=$FD_OUT/tags/procap_t15
FD_TAG_CAP_T60=$FD_OUT/tags/procap_t60

annotatePeaks.pl \
    $FP_TSS_MERGE \
    $FD_GEN/Homo_sapiens.GRCh38.dna.primary_assembly.fa \
    -gtf $FD_GEN/gencode.v34.annotation.gtf \
    -strand + -fragLength 1 -raw -d $FD_TAG_CAP_T15/ $FD_TAG_CAP_T60/ > $FP_TSS_COUNT


EOF

!ls /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss*.tss.txt | cut -f1-5 | expand -t 15

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_diffOutput.tss.txt

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_merge.tss.txt

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_count.tss.txt

ls -l /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/merge

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/merge

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt | cut -f1-5 | expand -t 15

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt | cut -f6-10 | expand -t 15

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt | cut -f11-15 | expand -t 15

!head /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt | cut -f16-20 | expand -t 15

