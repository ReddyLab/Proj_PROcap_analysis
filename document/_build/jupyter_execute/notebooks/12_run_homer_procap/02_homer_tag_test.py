# Run HOMER makeTagDirectory of PRO-cap and RNA-seq data

The code is from the HOMER documentation: 

[HOMER: csRNA-seq Analysis Tutorial](http://homer.ucsd.edu/homer/ngs/csRNAseq/index.html)

%%bash
module load perl
module load gcc
export PATH=/data/reddylab/software/homer/bin/:$PATH
makeTagDirectory

## HOMER makeTagDirectory of single sample (Option: genome hg38)

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/maketag_t00_hg38.txt \
    <<'EOF'
#!/bin/bash

### set directory
FP_GEN=/data/reddylab/Kuei/annotation/gencode.v34.annotation.gtf
FD_CAP=/data/reddylab/Kuei/Dex_ProCap
FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads

### PRO-cap
FP_BAM_CAP_T00=$FD_CAP/new_files/A549_control_alignments/A549_untreated_merged.bam
FP_BAM_CAP_T15=$FD_CAP/new_files/A549_15min_Dexamethasone_alignments/A549_15min_Dexamethasone_merged.bam
FP_BAM_CAP_T60=$FD_CAP/new_files/A549_1hr_Dexamethasone_alignments/A549_1h_Dexamethasone_merged.bam

### make Tag dir for PRO-cap (control)
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer
FD_TAG=$FD_OUT/tags/procap_t00_hg38
FP_BAM=$FP_BAM_CAP_T00
makeTagDirectory $FD_TAG/ $FP_BAM -genome hg38 -checkGC

EOF

### Quality Control Checks

import numpy  as np
import pandas as pd
import os
import matplotlib.pyplot as plt

FD_TAG="/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/tags/procap_t00_hg38"

fdiry = FD_TAG
fname = "tagCountDistribution.txt"
fpath = os.path.join(fdiry, fname)

dat = pd.read_csv(fpath, sep="\t")
dat.head()

x=dat['Tags per tag position (Median = 0, tags per genomic bp = 0.003)']
y=dat['Fraction of Positions']

plt.plot(x, y, '-o', markersize=1)
plt.xlabel("Tag Position")
plt.ylabel("Fraction of Positions")
plt.xscale('log')
plt.yscale('log')
plt.show()

fdiry = FD_TAG
fname = "tagLengthDistribution.txt"
fpath = os.path.join(fdiry, fname)

dat = pd.read_csv(fpath, sep="\t")
dat.head()

x=dat['Tag Length (Average tag length = 61.291551)']
y=dat['Fraction of Tags']

plt.plot(x, y, '-o', markersize=1)
plt.xlabel("Tag Length")
plt.ylabel("Fraction of Tags")
plt.show()

fdiry = FD_TAG
fname = "tagFreqUniq.txt"
fpath = os.path.join(fdiry, fname)

dat = pd.read_csv(fpath, sep="\t")
dat.head()

dat.iloc[:,:5].set_index("Offset").plot.line()
plt.show()

fdiry = FD_TAG
fname = "tagAutocorrelation.txt"
fpath = os.path.join(fdiry, fname)

dat = pd.read_csv(fpath, sep="\t")
dat.head()

dat.set_index("Distance in bp(Fragment Length Estimate: 223)(Peak Width Estimate: 223)").plot.line()
plt.xlabel("Distance in bp")
plt.ylabel("Read counts relative to each 5' end")
plt.xlim(-600, 600)
plt.show()

dat.set_index("Distance in bp(Fragment Length Estimate: 223)(Peak Width Estimate: 223)").plot.line()
plt.xlabel("Distance in bp")
plt.ylabel("Read counts relative to each 5' end")
plt.xlim(-150, 150)
plt.show()

-----
## HOMER makeTagDirectory of single sample (Option: maxlen 65)

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/log/maketag_t00_maxlen.txt \
    <<'EOF'
#!/bin/bash

### set directory
FP_GEN=/data/reddylab/Kuei/annotation/gencode.v34.annotation.gtf
FD_CAP=/data/reddylab/Kuei/Dex_ProCap
FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads

### PRO-cap
FP_BAM_CAP_T00=$FD_CAP/new_files/A549_control_alignments/A549_untreated_merged.bam
FP_BAM_CAP_T15=$FD_CAP/new_files/A549_15min_Dexamethasone_alignments/A549_15min_Dexamethasone_merged.bam
FP_BAM_CAP_T60=$FD_CAP/new_files/A549_1hr_Dexamethasone_alignments/A549_1h_Dexamethasone_merged.bam

### make Tag dir for PRO-cap (control)
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer
FD_TAG=$FD_OUT/tags/procap_t00_maxlen
FP_BAM=$FP_BAM_CAP_T00
makeTagDirectory $FD_TAG/ $FP_BAM -genome hg38 -checkGC -maxlen 65

EOF

### Quality Control Checks:

import numpy  as np
import pandas as pd
import os
import matplotlib.pyplot as plt

FD_TAG="/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer/tags/procap_t00_maxlen/"

#### Read clonality/PCR duplicates

fdiry = FD_TAG
fname = "tagCountDistribution.txt"
fpath = os.path.join(fdiry, fname)

dat = pd.read_csv(fpath, sep="\t")
dat.head(3)

x=dat.loc[:, dat.columns.str.startswith('Tags per tag position')]
y=dat['Fraction of Positions']

plt.plot(x, y, '-o', markersize=1)
plt.xlabel("Tag Position")
plt.ylabel("Fraction of Positions")
plt.xscale('log')
plt.yscale('log')
plt.show()

#### Read length distribution

fdiry = FD_TAG
fname = "tagLengthDistribution.txt"
fpath = os.path.join(fdiry, fname)

dat = pd.read_csv(fpath, sep="\t")
dat.head(3)

x=dat.loc[:, dat.columns.str.startswith('Tag Length')]
y=dat['Fraction of Tags']

plt.plot(x, y, '-o', markersize=1)
plt.xlabel("Tag Length")
plt.ylabel("Fraction of Tags")
plt.ylim(0, 0.05)
plt.show()

#### Nucleotide preferences near the 5' end of the read

fdiry = FD_TAG
fname = "tagFreqUniq.txt"
fpath = os.path.join(fdiry, fname)

dat = pd.read_csv(fpath, sep="\t")
dat.head(3)

dat.iloc[:,:5].set_index("Offset").plot.line()
plt.show()

#### Distribution of reads relative to one another

fdiry = FD_TAG
fname = "tagAutocorrelation.txt"
fpath = os.path.join(fdiry, fname)

dat = pd.read_csv(fpath, sep="\t")
dat.head(3)

val = dat.columns.values.astype(str)
idx = np.char.startswith(val, "Distance")

dat.set_index(val[idx][0]).plot.line()
plt.xlabel("Distance in bp")
plt.ylabel("Read counts relative to each 5' end")
plt.xlim(-150, 150)
plt.show()

-----

## Run makeTagDirectory for all sample

```
### set directory
FP_GEN=/data/reddylab/Kuei/annotation/gencode.v34.annotation.gtf
FD_CAP=/data/reddylab/Kuei/Dex_ProCap
FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads

### PRO-cap
FP_BAM_CAP_T00=$FD_CAP/new_files/A549_control_alignments/A549_untreated_merged.bam
FP_BAM_CAP_T15=$FD_CAP/new_files/A549_15min_Dexamethasone_alignments/A549_15min_Dexamethasone_merged.bam
FP_BAM_CAP_T60=$FD_CAP/new_files/A549_1hr_Dexamethasone_alignments/A549_1h_Dexamethasone_merged.bam

### make Tag dir for PRO-cap (control)
FD_OUT=/gpfs/fs1/data/reddylab/Kuei/Dex_ProCap/run_homer
FD_TAG=$FD_OUT/tags/procap_t00_hg38
FP_BAM=$FP_BAM_CAP_T00
makeTagDirectory $FD_TAG/ $FP_BAM -genome hg38 -checkGC
```

Creating Tag Directories

`makeTagDirectory <Output Directory Name> [options] <alignment file1> [alignment file 2] ...`

### PRO-cap (all samples)

ls /data/reddylab/Kuei/Dex_PROcap

ls -1 /data/reddylab/Kuei/Dex_PROcap/run_homer/tags

ls -1 /data/reddylab/Kuei/Dex_PROcap/new_files/*/*bam

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/maketag_cap_all.txt \
    <<'EOF'
#!/bin/bash

### set directory
FP_GEN=/data/reddylab/Kuei/annotation/gencode.v34.annotation.gtf
FD_CAP=/data/reddylab/Kuei/Dex_PROcap
FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads

### PRO-cap
FP_BAM_CAP_T00=$FD_CAP/new_files/A549_control_alignments/A549_untreated_merged.bam
FP_BAM_CAP_T15=$FD_CAP/new_files/A549_15min_Dexamethasone_alignments/A549_15min_Dexamethasone_merged.bam
FP_BAM_CAP_T60=$FD_CAP/new_files/A549_1hr_Dexamethasone_alignments/A549_1h_Dexamethasone_merged.bam

### make Tag dir for PRO-cap (control)
FD_OUT=$FD_CAP/run_homer
FD_TAG=$FD_OUT/tags/procap_t00
FP_BAM=$FP_BAM_CAP_T00
makeTagDirectory $FD_TAG/ $FP_BAM -genome hg38 -checkGC

### make Tag dir for PRO-cap (15m)
FD_OUT=$FD_CAP/run_homer
FD_TAG=$FD_OUT/tags/procap_t15
FP_BAM=$FP_BAM_CAP_T15
makeTagDirectory $FD_TAG/ $FP_BAM -genome hg38 -checkGC

### make Tag dir for PRO-cap (1hr)
FD_OUT=$FD_CAP/run_homer
FD_TAG=$FD_OUT/tags/procap_t60
FP_BAM=$FP_BAM_CAP_T60
makeTagDirectory $FD_TAG/ $FP_BAM -genome hg38 -checkGC

EOF

### PRO-cap (all samples; options: maxlen)

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/maketag_cap_all_maxlen.txt \
    <<'EOF'
#!/bin/bash

### set directory
FP_GEN=/data/reddylab/Kuei/annotation/gencode.v34.annotation.gtf
FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads
FD_CAP=/data/reddylab/Kuei/Dex_PROcap
FD_OUT=/data/reddylab/Kuei/Dex_PROcap/run_homer

### PRO-cap
FP_BAM_CAP_T00=$FD_CAP/new_files/A549_control_alignments/A549_untreated_merged.bam
FP_BAM_CAP_T15=$FD_CAP/new_files/A549_15min_Dexamethasone_alignments/A549_15min_Dexamethasone_merged.bam
FP_BAM_CAP_T60=$FD_CAP/new_files/A549_1hr_Dexamethasone_alignments/A549_1h_Dexamethasone_merged.bam

### make Tag dir for PRO-cap (control)
FD_TAG=$FD_OUT/tags/procap_t00_maxlen
FP_BAM=$FP_BAM_CAP_T00
makeTagDirectory $FD_TAG/ $FP_BAM -genome hg38 -checkGC -maxlen 65

### make Tag dir for PRO-cap (15m)
FD_TAG=$FD_OUT/tags/procap_t15_maxlen
FP_BAM=$FP_BAM_CAP_T15
makeTagDirectory $FD_TAG/ $FP_BAM -genome hg38 -checkGC -maxlen 65

### make Tag dir for PRO-cap (1hr)
FD_TAG=$FD_OUT/tags/procap_t60_maxlen
FP_BAM=$FP_BAM_CAP_T60
makeTagDirectory $FD_TAG/ $FP_BAM -genome hg38 -checkGC -maxlen 65

EOF

### RNA-seq (3 time points)

ls /data/reddylab/projects/GGR/data/rna_seq/mapped_reads

!ls /data/reddylab/projects/GGR/data/rna_seq/mapped_reads/iter0/t00_rep?/

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/maketag_rna_all.txt \
    <<'EOF'
#!/bin/bash

### set directory
FP_GEN=/data/reddylab/Kuei/annotation/gencode.v34.annotation.gtf
FD_GGR=/data/reddylab/projects/GGR/
FD_RNA_ALIGN=$FD_GGR/data/rna_seq/mapped_reads
FD_CAP=/data/reddylab/Kuei/Dex_PROcap
FD_OUT=/data/reddylab/Kuei/Dex_PROcap/run_homer

### RNA-seq
FP_BAM_RNA_T00_short_rep1=$FD_RNA_ALIGN/iter_short/t00m_rep1.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T00_short_rep2=$FD_RNA_ALIGN/iter_short/t00m_rep2.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T00_short_rep3=$FD_RNA_ALIGN/iter_short/t00m_rep3.star2.Aligned.out.sorted.bam

FP_BAM_RNA_T00_long_rep1=$FD_RNA_ALIGN/iter0/t00_rep1/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T00_long_rep2=$FD_RNA_ALIGN/iter0/t00_rep2/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T00_long_rep3=$FD_RNA_ALIGN/iter0/t00_rep3/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T00_long_rep4=$FD_RNA_ALIGN/iter0/t00_rep4/STAR_2pass_featurecounts/Aligned.out.sorted.bam

FP_BAM_RNA_T15_short_rep1=$FD_RNA_ALIGN/iter_short/t15m_rep1.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T15_short_rep2=$FD_RNA_ALIGN/iter_short/t15m_rep2.star2.Aligned.out.sorted.bam
FP_BAM_RNA_T15_short_rep3=$FD_RNA_ALIGN/iter_short/t15m_rep3.star2.Aligned.out.sorted.bam

FP_BAM_RNA_T60_long_rep1=$FD_RNA_ALIGN/iter0/t1_rep1/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T60_long_rep2=$FD_RNA_ALIGN/iter0/t1_rep2/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T60_long_rep3=$FD_RNA_ALIGN/iter0/t1_rep3/STAR_2pass_featurecounts/Aligned.out.sorted.bam
FP_BAM_RNA_T60_long_rep4=$FD_RNA_ALIGN/iter0/t1_rep4/STAR_2pass_featurecounts/Aligned.out.sorted.bam

### make Tag dir for RNA-seq (control)
FD_TAG=$FD_OUT/tags/rnaseq_t00
FP_BAM1=$FP_BAM_RNA_T00_short_rep1
FP_BAM2=$FP_BAM_RNA_T00_short_rep2
FP_BAM3=$FP_BAM_RNA_T00_short_rep3
makeTagDirectory $FD_TAG/ $FP_BAM1 $FP_BAM2 $FP_BAM3 -genome hg38 -checkGC

### make Tag dir for RNA-seq (15m)
FD_TAG=$FD_OUT/tags/rnaseq_t15
FP_BAM1=$FP_BAM_RNA_T15_short_rep1
FP_BAM2=$FP_BAM_RNA_T15_short_rep2
FP_BAM3=$FP_BAM_RNA_T15_short_rep3
makeTagDirectory $FD_TAG/ $FP_BAM1 $FP_BAM2 $FP_BAM3 -genome hg38 -checkGC

### make Tag dir for RNA-seq (1hr)
FD_TAG=$FD_OUT/tags/rnaseq_t60
FP_BAM1=$FP_BAM_RNA_T60_long_rep1
FP_BAM2=$FP_BAM_RNA_T60_long_rep2
FP_BAM3=$FP_BAM_RNA_T60_long_rep3
FP_BAM4=$FP_BAM_RNA_T60_long_rep4
makeTagDirectory $FD_TAG/ $FP_BAM1 $FP_BAM2 $FP_BAM3 $FP_BAM4 -genome hg38 -checkGC

EOF

-----

makeTagDirectory basically parses through the alignment file and splits the tags into separate files based on their chromosome.  As a result, several *.tags.tsv files are created in the output directory.  These are made to very efficiently return to the data during downstream analysis.  This also helps speed up the analysis of very large data sets without running out of memory.

In the end, your output directory will contain several *.tags.tsv files, as well as a file named "tagInfo.txt".  

This file contains information about your sequencing run, including the total number of tags considered.  This file is used by later programs to quickly reference information about the experiment, and can be manually modified to set certain parameters for analysis.

!head /data/reddylab/Kuei/Dex_PROcap/run_homer/tags/procap_t00/tagInfo.txt

!head /data/reddylab/Kuei/Dex_PROcap/run_homer/tags/procap_t00/chr1.tags.tsv

-----

## Check output message

ls -1 /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/tags

cat /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/maketag_cap_all.txt | grep Avg

cat /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/maketag_cap_all_maxlen.txt | grep Avg

cat /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/maketag_rna_all.txt | grep Avg

## Quality Control Checks

import numpy  as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import seaborn as sns

FD_OUT="/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/tags"
os.listdir(FD_OUT)

FDIRYS_PROCAP = [
    'procap_t00',
    'procap_t15',
    'procap_t60']

FDIRYS_PROCAP_MAXLEN = [
    'procap_t00_maxlen',
    'procap_t15_maxlen',
    'procap_t60_maxlen']

FDIRYS_RNASEQ = [
    'rnaseq_t00',
    'rnaseq_t15',
    'rnaseq_t60']

FDIRYS = np.r_[FDIRYS_PROCAP, FDIRYS_PROCAP_MAXLEN, FDIRYS_RNASEQ]

### Read clonality/PCR duplicates

### set file and column names
cname = ["Tags per tag position", "Fraction of Positions"]
fname = "tagCountDistribution.txt"

### read in the results (procap)
lst = []
for fdiry in FDIRYS_PROCAP:
    fpath = os.path.join(FD_OUT, fdiry, fname)
    dat   = pd.read_csv(fpath, sep="\t", header=0, names=cname)
    dat   = dat.assign(sample=lambda x: fdiry)
    lst.append(dat)
dat_merge_procap = pd.concat(lst)

### read in the results (procap maxlen)
lst = []
for fdiry in FDIRYS_PROCAP_MAXLEN:
    fpath = os.path.join(FD_OUT, fdiry, fname)
    dat   = pd.read_csv(fpath, sep="\t", header=0, names=cname)
    dat   = dat.assign(sample=lambda x: fdiry)
    lst.append(dat)
dat_merge_procap_maxlen = pd.concat(lst)

### read in the results (rnaseq)
lst = []
for fdiry in FDIRYS_RNASEQ:
    fpath = os.path.join(FD_OUT, fdiry, fname)
    dat   = pd.read_csv(fpath, sep="\t", header=0, names=cname)
    dat   = dat.assign(sample=lambda x: fdiry)
    lst.append(dat)
dat_merge_rnaseq = pd.concat(lst)
    
###
display(dat_merge_procap.head(3))
display(dat_merge_procap_maxlen.head(3))
display(dat_merge_rnaseq.head(3))

fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(12, 3))
plt.subplots_adjust(hspace=0.5, wspace=0.3)

ax = axes[0]
sns.lineplot(data=dat_merge_procap, ax=ax, x=cname[0], y=cname[1], hue="sample")
ax.set_xscale('log')
ax.set_yscale('log')

ax = axes[1]
sns.lineplot(data=dat_merge_procap_maxlen, ax=ax, x=cname[0], y=cname[1], hue="sample")
ax.set_xscale('log')
ax.set_yscale('log')

ax = axes[2]
sns.lineplot(data=dat_merge_rnaseq, ax=ax, x=cname[0], y=cname[1], hue="sample")
ax.set_xscale('log')
ax.set_yscale('log')

plt.show()

### 
cname = ["Tags per tag position", "Fraction of Positions"]
fname = "tagCountDistribution.txt"

###
lst = []
for fdiry in FDIRYS:
    fpath = os.path.join(FD_OUT, fdiry, fname)
    dat   = pd.read_csv(fpath, sep="\t", header=0, names=cname)
    dat   = dat.assign(sample=lambda x: fdiry)
    lst.append(dat)
dat_merge = pd.concat(lst)

fig, ax = plt.subplots(figsize=(5, 3))
sns.lineplot(
    data = dat_merge,
    ax   = ax,
    x    = "Tags per tag position", 
    y    = "Fraction of Positions", 
    hue  = "sample")
plt.xscale('log')
plt.yscale('log')
plt.show()

### Read length distribution

### 
cname = ["Tag Length", "Fraction of Tags"]
fname = "tagLengthDistribution.txt"

### read in the results (procap)
lst = []
for fdiry in FDIRYS_PROCAP:
    fpath = os.path.join(FD_OUT, fdiry, fname)
    dat   = pd.read_csv(fpath, sep="\t", header=0, names=cname)
    dat   = dat.assign(sample=lambda x: fdiry)
    lst.append(dat)
dat_merge_procap = pd.concat(lst)

### read in the results (procap maxlen)
lst = []
for fdiry in FDIRYS_PROCAP_MAXLEN:
    fpath = os.path.join(FD_OUT, fdiry, fname)
    dat   = pd.read_csv(fpath, sep="\t", header=0, names=cname)
    dat   = dat.assign(sample=lambda x: fdiry)
    lst.append(dat)
dat_merge_procap_maxlen = pd.concat(lst)

### read in the results (rnaseq)
lst = []
for fdiry in FDIRYS_RNASEQ:
    fpath = os.path.join(FD_OUT, fdiry, fname)
    dat   = pd.read_csv(fpath, sep="\t", header=0, names=cname)
    dat   = dat.assign(sample=lambda x: fdiry)
    lst.append(dat)
dat_merge_rnaseq = pd.concat(lst)
    
###
display(dat_merge_procap.head(3))
display(dat_merge_procap_maxlen.head(3))
display(dat_merge_rnaseq.head(3))

fig, axes = plt.subplots(ncols=3, nrows=1, figsize=(12, 3))
plt.subplots_adjust(hspace=0.5, wspace=0.3)

sns.lineplot(data=dat_merge_procap,        ax=axes[0], x=cname[0], y=cname[1], hue="sample")
sns.lineplot(data=dat_merge_procap_maxlen, ax=axes[1], x=cname[0], y=cname[1], hue="sample")
sns.lineplot(data=dat_merge_rnaseq,        ax=axes[2], x=cname[0], y=cname[1], hue="sample")
plt.show()

According to the document, most csRNA-seq reads are between 20-55 nt in length. In the PRO-cap data I have, most reads are around 140bp

### Nucleotide preferences

fname = "tagFreqUniq.txt"
lst = []
for fdiry in FDIRYS:
    fpath = os.path.join(FD_OUT, fdiry, fname)
    dat   = pd.read_csv(fpath, sep="\t")
    dat   = dat.assign(sample=lambda x: fdiry)
    lst.append(dat)

fig, axes = plt.subplots(ncols=3, nrows=3, figsize=(12, 9))
plt.subplots_adjust(hspace=0.5, wspace=0.5)

axes = axes.ravel()
for ax, dat in zip(axes, lst):
    dat.iloc[:,:5].set_index("Offset").plot.line(ax=ax)
    ax.set_title(np.unique(dat['sample']))
plt.show()

According to the document:  `Most csRNA-seq data should have a strong preference for Initiator-like sequences near the 5' end (usually C-1A0). Most species also show evidence for a TATA box at approx. -30 bp (e.g. spike in A/T content)`

### Distribution of reads relative to one another

fname = "tagAutocorrelation.txt"
lst = []
for fdiry in FDIRYS:
    fpath = os.path.join(FD_OUT, fdiry, fname)
    dat   = pd.read_csv(fpath, sep="\t")
    dat   = dat.assign(sample=lambda x: fdiry)
    lst.append(dat)

fig, axes = plt.subplots(ncols=3, nrows=3, figsize=(12, 9))
plt.subplots_adjust(hspace=0.5, wspace=0.5)

axes = axes.ravel()
for ax, dat in zip(axes, lst):
    dat.set_index(dat.columns[0]).plot.line(ax=ax, legend=False)
    ax.set_xlabel("Distance in bp")
    ax.set_ylabel("Read counts relative to each 5' end")
    ax.set_xlim(-600, 600)
    ax.set_title(np.unique(dat['sample']))
    
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))    
plt.show()

According to the document: `In the case of transcription initiation data, reads from TSS will typically cluster on the same strand within 50-100nt from each other. Due to the [normally] bidirectional nature of transcription from regulatory elements, reads will also typically appear upstream of TSS on the opposite strand (i.e. bidirectional transcription).`

Below I have zoomed in further to check if the peak of the opposite strand of PRO-cap is within 50-100nt.

fig, axes = plt.subplots(ncols=3, nrows=3, figsize=(12, 9))
plt.subplots_adjust(hspace=0.5, wspace=0.5)

axes = axes.ravel()
for ax, dat in zip(axes, lst):
    dat.set_index(dat.columns[0]).plot.line(ax=ax, legend=False)
    ax.set_xlabel("Distance in bp")
    ax.set_ylabel("Read counts relative to each 5' end")
    ax.set_xlim(-150, 150)
    ax.set_title(np.unique(dat['sample']))
    
plt.legend(loc='center left', bbox_to_anchor=(1.0, 0.5))    
plt.show()











