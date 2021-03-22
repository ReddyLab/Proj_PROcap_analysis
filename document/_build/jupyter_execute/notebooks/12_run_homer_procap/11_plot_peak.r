suppressMessages(suppressWarnings(source("config.R")))

gene_sel = c(
    "CIDEC", "PER1", "ANGPTL4", "TSC22D3", "STOM", "ZFP36", "NFKBIA",
    "JUN", "FOSL1",  "CCL2", "IL11")
length(gene_sel)

gene_pos = c(
    "CIDEC", 
    "PER1", 
    "EKI2", "ETNK2",
    "ANGPTL4",
    "TSC22D3",
    "RGS2",
    "CDC42EP3",
    "STOM",
    "THBD",
    "DUSP1",
    "GADD45A",
    "KLF4",
    "ZFP36",
    "ERRFI1",
    "SDPR",
    "CTGF",
    "NFKBIA",
    "BIPARP",
    "BIRC3",
    "BCL6",
    "ELL2",
    "CEBPB",
    "CEBPD",
    "CITED2",
    "ARL4D",
    "GALNT4",
    "KLF6",
    "ZC3H12A")

gene_neg = c(
    "IL11",
    "RND1",
    "CCL2",
    "MAP3K14",
    "SNKPLK2",
    "GEM",
    "IER3",
    "IER2",
    "ELF3",
    "IER5L",
    "FOSL1",
    "JUN",
    "HES1",
    "ENC1",
    "SERTAD3",
    "BHLHE40",
    "BCL3",
    "KLF10"
)

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_annotate/annoTSS/tss_count_raw.txt"
dat_raw = read_tsv(fpath)
head(dat_raw, 3)

fpath = "/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer/out_annotate/annoTSS/tss_count_rlg.txt"
dat_rlg = read_tsv(fpath)
head(dat_rlg, 3)



dat = dat_raw %>% dplyr::filter(`Gene Name`=="PER1", grepl(x = `Detailed Annotation`, pattern="promoter"))
head(dat)



dat = dat_raw
dat = dat %>% dplyr::filter(`Gene Name` %in% gene_pos, grepl(x = `Detailed Annotation`, pattern="promoter"))
dat = dat %>% dplyr::select(
    Chr, Start, End, Strand, `Gene Name`, 
    starts_with("/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer")) %>%
    unite("loc", Chr:Strand)
colnames(dat) = c("loc", "gene", "t00", "t15", "t60")

dat = dat %>% gather(time_point, count, -gene, -loc)
dat$gene = factor(dat$gene, levels=gene_pos)
head(dat)

options(repr.plot.height=5, repr.plot.width=7)
gpt = ggplot(dat, aes(x=time_point, y=count, color=gene, group=loc)) + 
    geom_line() +
    theme_bw() + 
    labs(x="Time Point", y="Peak Count (raw)")
    
print(gpt)

dat = dat_rlg
dat = dat %>% dplyr::filter(`Gene Name` %in% gene_pos, grepl(x = `Detailed Annotation`, pattern="promoter"))
dat = dat %>% dplyr::select(
    Chr, Start, End, Strand, `Gene Name`, 
    starts_with("/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer")) %>%
    unite("loc", Chr:Strand)
colnames(dat) = c("loc", "gene", "t00", "t15", "t60")

dat = dat %>% gather(time_point, count, -gene, -loc)
dat$gene = factor(dat$gene, levels=gene_pos)
head(dat)

options(repr.plot.height=5, repr.plot.width=7)
gpt = ggplot(dat, aes(x=time_point, y=count, color=gene, group=loc)) + 
    geom_line() +
    theme_bw() + 
    labs(x="Time Point", y="Peak Count (rlg)")
print(gpt)

gene = gene_neg
dat = dat_raw

dat = dat %>% dplyr::filter(`Gene Name` %in% gene, grepl(x = `Detailed Annotation`, pattern="promoter"))
dat = dat %>% dplyr::select(
    Chr, Start, End, Strand, `Gene Name`, 
    starts_with("/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer")) %>%
    unite("loc", Chr:Strand)
colnames(dat) = c("loc", "gene", "t00", "t15", "t60")

dat = dat %>% gather(time_point, count, -gene, -loc)
dat$gene = factor(dat$gene, levels=gene)
head(dat)

options(repr.plot.height=5, repr.plot.width=6)
gpt = ggplot(dat, aes(x=time_point, y=count, color=gene, group=loc)) + 
    geom_line() +
    theme_bw() + 
    labs(x="Time Point", y="Peak Count (raw)")
    
print(gpt)

gene = gene_neg
dat = dat_rlg

dat = dat %>% dplyr::filter(`Gene Name` %in% gene, grepl(x = `Detailed Annotation`, pattern="promoter"))
dat = dat %>% dplyr::select(
    Chr, Start, End, Strand, `Gene Name`, 
    starts_with("/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer")) %>%
    unite("loc", Chr:Strand)
colnames(dat) = c("loc", "gene", "t00", "t15", "t60")

dat = dat %>% gather(time_point, count, -gene, -loc)
dat$gene = factor(dat$gene, levels=gene)
head(dat)

options(repr.plot.height=5, repr.plot.width=6)
gpt = ggplot(dat, aes(x=time_point, y=count, color=gene, group=loc)) + 
    geom_line() +
    theme_bw() + 
    labs(x="Time Point", y="Peak Count (rlg)")
print(gpt)

