###
library("tidyverse")
library("vroom")
library("RColorBrewer")

###
library("scales")
library("gridExtra")
library("reshape2")

###
library("pheatmap")

###
#library("corrplot")
#library(GGally)

###
FD_DAT="/data/reddylab/gjohnson/whole_genome_STARRseq/wgss3/expression_analyses"
FN_EXP="A549.GR.compounds.rnaseq.all.counts.txt"
FP_EXP=file.path(FD_DAT, FN_EXP)

### Bio annot.
library("biomaRt")
library("DOSE")
library("clusterProfiler")
library("org.Hs.eg.db")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
ensembl = useMart("ensembl", dataset="hsapiens_gene_ensembl")
grch38  = useMart("ensembl", dataset="hsapiens_gene_ensembl")

#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
fun_map = function(
    val, 
    fil = "ensembl_gene_id", 
    att=c("chromosome_name","ensembl_gene_id","entrezgene_id","external_gene_name")){
    att = unique(c(fil, att))
    res = getBM(
        attributes = att,
        filters    = fil,
        values     = val,
        mart       = grch38)
    return(res)
}