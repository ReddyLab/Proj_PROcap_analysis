# Summarize and visualize peaks from HOMER pipeline

suppressMessages(suppressWarnings(source("config.R")))

FD_OUT

###
fdiry = file.path(FD_OUT, "Dex_PROcap", "run_homer", "out_annotate", "annoTSS")
fname = "tss_count_raw.txt"
fpath = file.path(fdiry, fname)

###
dat_raw = read_tsv(fpath)
head(dat_raw, 3)

###
fdiry = file.path(FD_OUT, "Dex_PROcap", "run_homer", "out_annotate", "annoTSS")
fname = "tss_count_rlg.txt"
fpath = file.path(fdiry, fname)

###
dat_rlg = read_tsv(fpath)
head(dat_rlg, 3)

## Summarize Annotation across the peaks

options(repr.plot.width=5, repr.plot.height=3)
dat = dat_rlg
txt = dat$Annotation
txt = str_remove(string = txt, pattern = "\\(.*")
qplot(txt) + theme_bw() + theme(axis.text.x = element_text(hjust=1, vjust=0.5, angle=90, size=10))

## Visualize

dat = dat_rlg
dat = dat %>% mutate(annotation = str_remove(string = Annotation, pattern = " \\(.*"))
dat = dat %>% dplyr::select(
    Chr, Start, End, Strand, `Gene Name`, annotation,
    starts_with("/gpfs/fs1/data/reddylab/Kuei/Dex_PROcap/run_homer")) %>%
    unite("grange", Start:End,  sep = "-") %>%
    unite("loc",    Chr:grange, sep = ":")

colnames(dat) = c("loc", "strand", "gene", "annotation", "t00", "t15", "t60")
dat = dat %>% mutate(d15 = t15 - t00, d60 = t60 - t00, avg=mean(t00 + t15 + t60))
dat = dat %>% mutate(label = paste0(gene, " (", annotation, ")"))
head(dat)

options(repr.plot.width=5, repr.plot.height=3)
ggplot(dat, aes(x=d15, y=d60, color=annotation)) + 
    geom_point(size=0.1, alpha=0.5) + 
    labs(x="t15 - t00", y="t60 - t00") +
    theme_bw()

options(repr.plot.width=5, repr.plot.height=3)
tmp = dat %>% dplyr::filter(annotation=="promoter-TSS")
ggplot(tmp, aes(x=d15, y=d60, color=annotation)) + 
    geom_point(size=0.1, alpha=0.5) + 
    labs(x="t15 - t00", y="t60 - t00") +
    theme_bw()

### Interaction

fig = plot_ly(
    data = tmp, 
    x = ~d15, y = ~d60, 
    type   = "scatter", 
    mode   = 'markers',
    marker = list(size=3, opacity=0.5),
    hoverinfo = 'text',
    text = ~label,
    width  = 500, 
    height = 500)
fig

chrom      = "chr5"
chromStart = 143275931
chromEnd   = 143437512

dat = dat_raw
dat %>% dplyr::filter(Chr=="chr5", Start > 143275931, End < 143437512)

tmp %>% dplyr::filter(d15 > 0, d60 > 0) %>% na.omit %>% head

tmp %>% dplyr::filter(gene == "NR3C1") %>% na.omit %>% head