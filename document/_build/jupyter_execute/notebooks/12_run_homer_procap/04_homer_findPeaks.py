# HOMER FindPeaks

```
General analysis options:
	-o <filename|auto> (file name for to output peaks, default: stdout)
		"-o auto" will send output to "<tag directory>/peaks.txt", ".../regions.txt",
		or ".../transcripts.txt" depending on the "-style" option
	-style <option> (Specialized options for specific analysis strategies)
		factor (transcription factor ChIP-Seq, uses -center, output: peaks.txt,  default)
		histone (histone modification ChIP-Seq, region based, uses -region -size 500 -L 0, regions.txt)
		groseq (de novo transcript identification from GroSeq data, transcripts.txt)
	tss (TSS identification from 5' RNA sequencing, tss.txt)
		dnase (Hypersensitivity [crawford style (nicking)], peaks.txt)
		super (Super Enhancers, superEnhancers.txt)
		superhistone (Super Enhancers from H3K27ac data, superEnhancers.txt)
		mC (Cytosine methylation (BS-seq/methylC-seq), regions.txt)
		damid (DamID enrichment from DpnI digestion, regions.txt)
		clip (CLIP-Seq enrichment, strand specific, peaks.txt)
```

print(1)

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
findPeaks 

## FindPeaks of time point 0 min

### Time point 00 min with RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findPeaks_cap00_rna00.txt \
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

#findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>
#findPeaks Exp1-csRNA-TagDir/ -i Exp1-input-TagDir/ -style tss > tssOutput.txt
#findPeaks IMR90-5GROseq/ -o auto -style tss -i IMR90-GROseq/
FP_OUT=$FD_OUT/out_findPeak/tss_cap00_rna00.txt
findPeaks $FD_TAG_CAP_T00/ -i $FD_TAG_RNA_T00/ -style tss > $FP_OUT

EOF

### Time point 00 min without RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findPeaks_cap00.txt \
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

#findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>
#findPeaks Exp1-csRNA-TagDir/ -i Exp1-input-TagDir/ -style tss > tssOutput.txt
#findPeaks IMR90-5GROseq/ -o auto -style tss -i IMR90-GROseq/
FP_OUT=$FD_OUT/out_findPeak/tss_cap00.txt
findPeaks $FD_TAG_CAP_T00/ -style tss > $FP_OUT

EOF

## FindPeaks of time point 15 min

### Time point 15 min with RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findPeaks_cap15_rna15.txt \
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

#findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>
#findPeaks Exp1-csRNA-TagDir/ -i Exp1-input-TagDir/ -style tss > tssOutput.txt
#findPeaks IMR90-5GROseq/ -o auto -style tss -i IMR90-GROseq/
FP_OUT=$FD_OUT/out_findPeak/tss_cap15_rna15.txt
findPeaks $FD_TAG_CAP_T15/ -i $FD_TAG_RNA_T15/ -style tss > $FP_OUT

EOF

### Time point 15 min without RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findPeaks_cap15.txt \
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

#findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>
#findPeaks Exp1-csRNA-TagDir/ -i Exp1-input-TagDir/ -style tss > tssOutput.txt
#findPeaks IMR90-5GROseq/ -o auto -style tss -i IMR90-GROseq/
FP_OUT=$FD_OUT/out_findPeak/tss_cap15.txt
findPeaks $FD_TAG_CAP_T15/ -style tss > $FP_OUT

EOF

## FindPeaks of time point 60 min

### Time point 60 min with RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findPeaks_cap60_rna60.txt \
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

#findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>
#findPeaks Exp1-csRNA-TagDir/ -i Exp1-input-TagDir/ -style tss > tssOutput.txt
#findPeaks IMR90-5GROseq/ -o auto -style tss -i IMR90-GROseq/
FP_OUT=$FD_OUT/out_findPeak/tss_cap60_rna60.txt
findPeaks $FD_TAG_CAP_T60/ -i $FD_TAG_RNA_T60/ -style tss > $FP_OUT

EOF

### Time point 60 min without RNA control

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/findPeaks_cap60.txt \
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

#findPeaks <tag directory> -style <factor | histone> -o auto -i <control tag directory>
#findPeaks Exp1-csRNA-TagDir/ -i Exp1-input-TagDir/ -style tss > tssOutput.txt
#findPeaks IMR90-5GROseq/ -o auto -style tss -i IMR90-GROseq/
FP_OUT=$FD_OUT/out_findPeak/tss_cap60.txt
findPeaks $FD_TAG_CAP_T60/ -style tss > $FP_OUT

EOF

-----

## Check output file

!head -40 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findPeak/tss_cap00_rna00.txt

!head -40 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findPeak/tss_cap00.txt

----

import numpy  as np
import pandas as pd
import os
import matplotlib.pyplot as plt

cname = "PeakID	chr	start	end	strand	Normalized Tag Count	focus ratio	findPeaks Score	Total Tags	Control Tags (normalized to IP Experiment)	Fold Change vs Control	p-value vs Control	Fold Change vs Local	p-value vs Local	Dispersion Ratio	Periodic Ratio"
cname = cname.split("\t")

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_findPeak/tss_t00.txt"
dat   = pd.read_csv(fpath, sep="\t", names=cname, comment="#")
dat.head()

dat[dat['chr'] == "chr1"]

dat["Control Tags (normalized to IP Experiment)"].hist(cumulative=False, density=1, bins=100)

np.log2(dat["Control Tags (normalized to IP Experiment)"] + 1).hist(cumulative=False, density=1, bins=100)

