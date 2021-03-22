# Create Genome Browser Visualization Files

%%bash
module load perl
module load gcc
source /data/reddylab/software/miniconda2/bin/activate alex_dev
export PATH=/data/reddylab/software/homer/bin/:$PATH
sbatch -pnew,all \
    --mem 16G \
    -o /gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/log/makeUCSCfile_bedgraph.txt \
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

###
#makeUCSCfile Exp1-csRNA-tagDir/ -style tss -strand +  > exp1.posStrand.bedGraph
#makeUCSCfile Exp1-csRNA-tagDir/ -style tss -strand - -neg  > exp1.negStrand.bedGraph
PREFIX_BED_CAP_T00=$FD_OUT/out_bedgraph/procap_t00
PREFIX_BED_CAP_T15=$FD_OUT/out_bedgraph/procap_t15
PREFIX_BED_CAP_T60=$FD_OUT/out_bedgraph/procap_t60

echo "Generate BedGraph for T00 sample"
makeUCSCfile $FD_TAG_CAP_T00/ -style tss -strand +       > ${PREFIX_BED_CAP_T00}.plus.bedGraph
makeUCSCfile $FD_TAG_CAP_T00/ -style tss -strand - -neg  > ${PREFIX_BED_CAP_T00}.minus.bedGraph

echo "Generate BedGraph for T15 sample"
makeUCSCfile $FD_TAG_CAP_T15/ -style tss -strand +       > ${PREFIX_BED_CAP_T15}.plus.bedGraph
makeUCSCfile $FD_TAG_CAP_T15/ -style tss -strand - -neg  > ${PREFIX_BED_CAP_T15}.minus.bedGraph

echo "Generate BedGraph for T60 sample"
makeUCSCfile $FD_TAG_CAP_T60/ -style tss -strand +       > ${PREFIX_BED_CAP_T60}.plus.bedGraph
makeUCSCfile $FD_TAG_CAP_T60/ -style tss -strand - -neg  > ${PREFIX_BED_CAP_T60}.minus.bedGraph

EOF

makeUCSCfile Exp1-csRNA-tagDir/ -style tss -strand +  > exp1.posStrand.bedGraph
makeUCSCfile Exp1-csRNA-tagDir/ -style tss -strand - -neg  > exp1.negStrand.bedGraph