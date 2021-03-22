

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 32G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/findmotif.txt \
    <<'EOF'
#!/bin/bash

### 8. Analysis of Sequence Features and Motif Discovery at TSS
# annotatePeaks.pl <tss/peak/BED file> <genome> -size <#> -hist 1 -di > output.txt
# annotatePeaks.pl exp1.tss.txt hg38 -size 1000 -hist 1 -di > output.txt

#FP_TSS1=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt
#FP_TSS2=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t60v00.tss.txt
#FP_TSSM=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/merge


FD_IN=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS
FD_OT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_motif/tss_t15v00

FP_TSS_T15=${FD_IN}/tss_t15v00.tss.txt
FP_TSS_T60=${FD_IN}/tss_t60v00.tss.txt


#findMotifsGenome.pl exp1.tss.txt hg38 OutputDirectory/ -size -150,50
findMotifsGenome.pl $FP_TSS_T15 hg38 $FD_OT/ -size -150,50

EOF

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 32G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/findmotif_t60v00.txt \
    <<'EOF'
#!/bin/bash

### 8. Analysis of Sequence Features and Motif Discovery at TSS
# annotatePeaks.pl <tss/peak/BED file> <genome> -size <#> -hist 1 -di > output.txt
# annotatePeaks.pl exp1.tss.txt hg38 -size 1000 -hist 1 -di > output.txt

#FP_TSS1=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt
#FP_TSS2=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t60v00.tss.txt
#FP_TSSM=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/merge


FD_IN=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS
FD_OT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_motif/tss_t60v00

FP_TSS_T15=${FD_IN}/tss_t15v00.tss.txt
FP_TSS_T60=${FD_IN}/tss_t60v00.tss.txt

#findMotifsGenome.pl exp1.tss.txt hg38 OutputDirectory/ -size -150,50
findMotifsGenome.pl $FP_TSS_T60 hg38 $FD_OT/ -size -150,50

EOF

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 32G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/seqfreq.txt \
    <<'EOF'
#!/bin/bash

### 8. Analysis of Sequence Features and Motif Discovery at TSS
# annotatePeaks.pl <tss/peak/BED file> <genome> -size <#> -hist 1 -di > output.txt
# annotatePeaks.pl exp1.tss.txt hg38 -size 1000 -hist 1 -di > output.txt

#FP_TSS1=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt
#FP_TSS2=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t60v00.tss.txt
#FP_TSSM=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/merge


FD_IN=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS
FD_OT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_motif

FP_TSS_T15=${FD_IN}/tss_t15v00.tss.txt
FP_TSS_T60=${FD_IN}/tss_t60v00.tss.txt
FP_MOTIF_T15=${FD_OT}/motif_tss_t15v00.txt
FP_MOTIF_T15=${FD_OT}/motif_tss_t60v00.txt

FP_TSS_SEQFREQ_T15=tss_t15v00_seqfreq.txt

#annotatePeaks.pl exp1.tss.txt hg38 -size 1000 -hist 1 -di > output.txt
annotatePeaks.pl \
    $FP_TSS_T15 \
    hg38 \
    -size 1000 -hist 1 -di > $FP_MOTIF_T15

annotatePeaks.pl \
    $FP_TSS_T60 \
    hg38 \
    -size 1000 -hist 1 -di > $FP_MOTIF_T60

EOF

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 32G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/seqfreq.txt \
    <<'EOF'
#!/bin/bash

### 8. Analysis of Sequence Features and Motif Discovery at TSS
# annotatePeaks.pl <tss/peak/BED file> <genome> -size <#> -hist 1 -di > output.txt
# annotatePeaks.pl exp1.tss.txt hg38 -size 1000 -hist 1 -di > output.txt

#FP_TSS1=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t15v00.tss.txt
#FP_TSS2=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/tss_t60v00.tss.txt
#FP_TSSM=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/merge


FD_IN=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS
FD_OT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_motif

FP_TSS_T15=${FD_IN}/tss_t15v00.tss.txt
FP_TSS_T60=${FD_IN}/tss_t60v00.tss.txt
FP_MOTIF_T15=${FD_OT}/output.txt

cd /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/out_findTSS/

FP_TSS_T15=tss_t15v00.tss.txt
FP_TSS_T60=tss_t60v00.tss.txt
FP_TSS_MERGE=tss_merge.tss.txt
FP_TSS_COUNT_RAW=tss_count_raw.tss.txt
FP_TSS_COUNT_RLG=tss_count_rlg.tss.txt

FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer
FD_GEN=/gpfs/fs1/data/reddylab/Kuei/annotation

FP_TSS_SEQFREQ_T15=tss_t15v00_seqfreq.txt

#annotatePeaks.pl exp1.tss.txt hg38 -size 1000 -hist 1 -di > output.txt
annotatePeaks.pl \
    $FP_TSS_T15 \
    hg38 \
    -size 1000 -hist 1 -di > $FP_TSS_SEQFREQ_T15

EOF

cat /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/annotateTSS.txt

