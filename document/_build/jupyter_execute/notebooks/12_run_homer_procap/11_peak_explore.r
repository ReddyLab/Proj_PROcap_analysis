# Exploring PRO-cap Peaks

suppressMessages(suppressWarnings(source("config_genome.R")))



###
fdiry = file.path(FD_OUT, "Dex_PROcap", "run_homer", "out_annotate", "annoTSS")
fname = "tss_count_raw.txt"
fpath = file.path(fdiry, fname)

###
dat_raw = read_tsv(fpath)
head(dat_raw, 3)

dat = dat_raw
dat = dat %>% mutate(annotation = str_remove(string = Annotation, pattern = " \\(.*"))
dat = dat %>% dplyr::select(Chr, Start, End, Strand, `Peak Score`, annotation)
colnames(dat) = c("Chr", "Start", "End", "Strand", "Score", "Annotation")
head(dat)

grg = GRanges(
    seqnames = dat$Chr,
    ranges = IRanges(
        start = dat$Start,
        end   = dat$End
    ),
    strand = dat$Strand,
    score  = dat$Score,
    annotation = dat$Annotation
)

grg

mcols(grg)

autoplot(grg)

autoplot(grg, aes(fill=as.factor(annotation), color=as.factor(annotation)))

tmp = disjoin(grg)
tmp

