# Data preparation script: import Borzoi model predictions into misha tracks.
# This script documents the steps used to prepare data and is NOT meant to be re-run directly.
# It references local paths to bigWig prediction outputs from the Borzoi training pipeline.
# The resulting tracks are provided in the data/mm10/tracks/ directory.
#
# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .R
#       format_name: light
#       format_version: '1.5'
#       jupytext_version: 1.16.4
#   kernelspec:
#     display_name: R 4.4.1
#     language: R
#     name: r_4.4.1
# ---

# + vscode={"languageId": "r"}
library(here)
library("misha")
library("misha.ext")
library("zoo")
library("tglkmeans")
library("misha.ext")
library(tidyverse)
library(tgstat)
library(prego)
library(tgutil)
#gdb.reload()
options(gmax.data.size = 1e10)
options(gmultitasking = FALSE)


link_dir = link_dir =  paste0(here(),'/')
setwd(here())

gsetroot(paste0(link_dir,'data/mm10/'))

gdb.reload()

source(paste0(link_dir,'code/seq2epi_utils.r'))
source(paste0(link_dir,'code/fig_fun.r'))



mod = init_pipe()
mod$epi_tss = mod$tss

# + vscode={"languageId": "r"}
compute_crop_bp <- function(seq_len = 65536, bin_size_bp = 32, prediction_fraction = 0.5) {
  # total number of bins Borzoi outputs
  total_bins <- seq_len / bin_size_bp
  
  # number of bins we keep in the center
  pred_bins  <- total_bins * prediction_fraction
  
  # bins to trim each side
  crop_left_bins  <- floor((total_bins - pred_bins) / 2)
  crop_right_bins <- (total_bins - pred_bins) - crop_left_bins
  
  # convert to base pairs
  crop_left_bp  <- crop_left_bins  * bin_size_bp
  crop_right_bp <- crop_right_bins * bin_size_bp
  
  list(
    total_bp      = seq_len,
    total_bins    = as.integer(total_bins),
    kept_bins     = as.integer(pred_bins),
    kept_bp       = as.integer(pred_bins * bin_size_bp),
    crop_left_bp  = as.integer(crop_left_bp),
    crop_right_bp = as.integer(crop_right_bp)
  )
}
compute_crop_bp()

# + vscode={"languageId": "r"}

align_to_bin_size <- function(regions, bin_size = 32) {
    regions %>%
        mutate(
            start = floor(start / bin_size) * bin_size,
            end = floor(end / bin_size) * bin_size
        )
}

perpare_regions <- function(regions, seq_len = 65536, bin_size = 32) {
    crops <- compute_crop_bp(seq_len, bin_size)
    cli::cli_alert("Kept: {.val {crops$kept_bp}}, total: {.val {crops$total_bp}}, crop_left: {.val {crops$crop_left_bp}}, crop_right: {.val {crops$crop_right_bp}}")
    cli::cli_alert_info("initial number of regions: {.val {nrow(regions)}}")
    
    regions <- regions %>% gintervals.normalize(seq_len) 
    new_regions <- NULL

    while(is.null(new_regions) || nrow(regions) != nrow(new_regions)) {
        if (!is.null(new_regions)) {
            regions <- new_regions
        }
        new_regions <- regions %>% 
            gintervals.normalize(seq_len) %>%
            gintervals.canonic()        
        cli::cli_alert_info("New number of regions: {.val {nrow(new_regions)}}")
    }
    regions <- regions %>% 
        gintervals.normalize(seq_len) 
    regions <- align_to_bin_size(regions, bin_size)      
    stopifnot(all((regions$start %% 32) == 0))
    stopifnot(all((regions$end %% 32) == 0))
    
    return(regions)
}

set.seed(60427)
random_regs <- gintervals.random(size = 2e3, n = 5e3, filter = mod$cgdom_ann)
initial_regs <- bind_rows(mod$cgdom_ann %>% select(chrom, start, end), random_regs)
regs <- perpare_regions(initial_regs) 
cli::cli_alert_info("Covered: {.val {scales::comma(gintervals.covered_bp(regs))}bp}")
dir.create(here("output/data-for-borzoi"), showWarnings = FALSE)
fwrite(regs, here("output/data-for-borzoi/borzoi_regions.bed"), sep = "\t", col.names = FALSE, quote = FALSE, scipen = 9999)

# + vscode={"languageId": "r"}
sps(7, 7)
regs %>% mutate(l = end - start) %>% arrange(l) %>% head
quantile(regs$end - regs$start, 0:20/20)
regs %>% 
    ggplot(aes(x = end - start, y = 1-after_stat(y))) + stat_ecdf() 

# + vscode={"languageId": "r"}
gen_k27_vt(mod)
gen_k4_vt(mod)
gvtrack.create("e65_epi_atac", "mut_multi.WT_E65.epi", func = "sum")
gvtrack.iterator("e65_epi_atac", sshift = -140, eshift = 140)
gvtrack.create("e70_epi_atac", "mut_multi.WT_E70.epi", func = "sum")
gvtrack.iterator("e70_epi_atac", sshift = -140, eshift = 140)
gvtrack.create("e75_epi_atac", "mut_multi.WT_E75.epi", func = "sum")
gvtrack.iterator("e75_epi_atac", sshift = -140, eshift = 140)
borzoi_vt <- c('ES_cnt','ES2i_cnt','EB4_cnt','e75_ecto_cnt','e75_emeso_cnt', "EB4_cnt_k4", "e65_epi_atac", "e70_epi_atac", "e75_epi_atac")
borzoi_data <- gextract.left_join(borzoi_vt, intervals = regs, iterator = 32) %>% 
    mutate_at(borzoi_vt, ~ ifelse(is.na(.), 0, .)) %>%
    group_by(chrom) %>%     
    mutate_at(borzoi_vt, ~ iceqream::norm01(log2(1 + .))) %>%     
    ungroup() %>% 
    mutate(region = paste0(chrom1, ":", start1, "-", end1)) %>% 
    select(-chrom1, -start1, -end1) %>% 
    select(chrom, start, end, region, everything()) %>% 
    arrange(chrom, start) %fcache_df% here("output/data-for-borzoi/borzoi_data.tsv")

# + vscode={"languageId": "r"}
arrow::write_parquet(borzoi_data, here("output/data-for-borzoi/borzoi_data.parquet"), compression = "snappy")

# + vscode={"languageId": "r"}
borzoi_data_no_norm <- gextract.left_join(borzoi_vt, intervals = regs, iterator = 32) %>% 
    mutate_at(borzoi_vt, ~ ifelse(is.na(.), 0, .)) %>%
    group_by(chrom) %>%     
    # mutate_at(borzoi_vt, ~ log2(1 + .)) %>%     
    ungroup() %>% 
    mutate(region = paste0(chrom1, ":", start1, "-", end1)) %>% 
    select(-chrom1, -start1, -end1) %>% 
    select(chrom, start, end, region, everything()) %>% 
    arrange(chrom, start) %fcache_df% here("output/data-for-borzoi/borzoi_data_no_norm.tsv")
arrow::write_parquet(borzoi_data_no_norm, here("output/data-for-borzoi/borzoi_data_no_norm.parquet"), compression = "snappy")

# + vscode={"languageId": "r"}
borzoi_data %>% count(region) %>% count(n, name = "num_of_regs")

# + vscode={"languageId": "r"}
browser()

# + vscode={"languageId": "r"}
regs_500k <- perpare_regions(initial_regs, 524288)
cli::cli_alert_info("Covered: {.val {scales::comma(gintervals.covered_bp(regs_500k))}bp}")
dir.create(here("output/data-for-borzoi"), showWarnings = FALSE)
fwrite(regs_500k, here("output/data-for-borzoi/borzoi_regions_500k.bed"), sep = "\t", col.names = FALSE, quote = FALSE, scipen = 9999)

# + vscode={"languageId": "r"}
gen_k27_vt(mod)
gen_k4_vt(mod)
gvtrack.create("e65_epi_atac", "mut_multi.WT_E65.epi", func = "sum")
gvtrack.iterator("e65_epi_atac", sshift = -140, eshift = 140)
gvtrack.create("e70_epi_atac", "mut_multi.WT_E70.epi", func = "sum")
gvtrack.iterator("e70_epi_atac", sshift = -140, eshift = 140)
gvtrack.create("e75_epi_atac", "mut_multi.WT_E75.epi", func = "sum")
gvtrack.iterator("e75_epi_atac", sshift = -140, eshift = 140)
borzoi_vt <- c('ES_cnt','ES2i_cnt','EB4_cnt','e75_ecto_cnt','e75_emeso_cnt', "EB4_cnt_k4", "e65_epi_atac", "e70_epi_atac", "e75_epi_atac")
borzoi_data_no_norm <- gextract.left_join(borzoi_vt, intervals = regs_500k, iterator = 32) %>% 
    mutate_at(borzoi_vt, ~ ifelse(is.na(.), 0, .)) %>%
    group_by(chrom) %>%     
    # mutate_at(borzoi_vt, ~ log2(1 + .)) %>%     
    ungroup() %>% 
    mutate(region = paste0(chrom1, ":", start1, "-", end1)) %>% 
    select(-chrom1, -start1, -end1) %>% 
    select(chrom, start, end, region, everything()) %>% 
    arrange(chrom, start) %fcache_df% here("output/data-for-borzoi/borzoi_data_no_norm_500k.tsv")
arrow::write_parquet(borzoi_data_no_norm, here("output/data-for-borzoi/borzoi_data_no_norm_500k.parquet"), compression = "snappy")

# + vscode={"languageId": "r"}


# + vscode={"languageId": "r"}
gdir.create("borzoi", showWarnings = FALSE)
gtrack.import(file = "/home/aviezerl/proj/ebpcg/seq2epi_paper_raw/output/data-for-borzoi/predictions/64kb_context_pytroch/EB4_cnt.bw", track = "borzoi.finetune_64k_EB4_cnt", description = "", binsize = 20)

# + vscode={"languageId": "r"}
gdir.create("borzoi", showWarnings = FALSE)
gtrack.import(file = "/home/aviezerl/proj/ebpcg/seq2epi_paper_raw/output/data-for-borzoi/predictions/64kb_context_no_norm_pytroch/EB4_cnt.bw", track = "borzoi.finetune_64k_no_norm_EB4_cnt", description = "", binsize = 20)

# + vscode={"languageId": "r"}


# + vscode={"languageId": "r"}
pred_eb4_cnt <- rtracklayer::import("/home/aviezerl/proj/ebpcg/seq2epi_paper_raw/output/data-for-borzoi/predictions/64kb_context_pytroch/EB4_cnt.bw") %>% 
    as.data.frame() %>% 
    rename(chrom = seqnames) %>% 
    mutate(start = start - 1, end = end) %>% 
    select(-strand)

# + vscode={"languageId": "r"}
val_chroms <- c("chr8", "chr10")
test_chroms <- c("chr9", "chr18")
a <- pred_eb4_cnt %>% left_join(borzoi_data) %>% filter(!is.na(EB4_cnt)) %>% mutate(type = case_when(chrom %in% val_chroms ~ "val", chrom %in% test_chroms ~ "test", TRUE ~ "train"))

# + vscode={"languageId": "r"}
a %>% group_by(type) %>% summarise(r2 = cor(EB4_cnt, score)^2)

# + vscode={"languageId": "r"}


# + vscode={"languageId": "r"}
create_bw <- function(df, filename, column) {
    df <- as_tibble(df)
    score <- df[[column]]
    # score <- iceqream::norm01(df[[column]])

    gr <- GenomicRanges::GRanges(
        seqnames = df$chrom,
        # We convert to 1-based here so GRanges is happy.
        # rtracklayer will convert it back to 0-based when writing the file.
        ranges = IRanges::IRanges(start = df$start + 1, end = df$end),
        score = score
    )
    sl <- gintervals.all() %>%
        select(chrom, end) %>%
        tibble::deframe()
    sl <- sl[as.character(levels(GenomeInfoDb::seqnames(gr)))]
    GenomeInfoDb::seqlengths(gr) <- sl

    rtracklayer::export.bw(gr, filename)
}
dir.create(here("output/data-for-borzoi/borzoi_bw"), showWarnings = FALSE)
purrr::walk(borzoi_vt, ~ create_bw(borzoi_data, here("output/data-for-borzoi/borzoi_bw", paste0(.x, ".bw")), .x))

# + vscode={"languageId": "r"}


# + vscode={"languageId": "r"}
pred_eb4_cnt <- rtracklayer::import("/net/mraid20/ifs/wisdom/tanay_lab/tgdata/users/evghenic/Proj/polycomb/seq2epi_paper_raw/output/borzoi-output/it7_64kb_context/track_2.bw") %>% 
    as.data.frame() %>% 
    rename(chrom = seqnames) %>% 
    mutate(start = start - 1, end = end - 1) %>% 
    select(-strand)

# + vscode={"languageId": "r"}
# mod$gw$cg_trace = readRDS('./data/cg_trace_mm10.rds')

# + vscode={"languageId": "r"}
gvtrack.create("d_ltr", "intervs.global.rmsk_ltr", "distance")
gvtrack.create("d_line", "intervs.global.rmsk_line", "distance")
gvtrack.create("d_sine", "intervs.global.rmsk_sine", "distance")
regions <- gscreen("seq.CG_500_mean_new > 0.02 & d_ltr != 0 & d_line != 0")
set.seed(60427)
regions_e <- bind_rows(mod$cgdom_ann  %>% select(chrom, start, end), regions) %>%
    gintervals.canonic()
rand_genome <- gintervals.random(size = 2e3, n = 2e4, filter = regions_e)
reg32 <- giterator.intervals(intervals = gintervals.canonic(bind_rows(regions_e, rand_genome)), iterator = 32) %>% 
    mutate(start = floor(start / 32) * 32, end = start + 32)
reg32 <- gextract(c("d_ltr", "d_line", "d_sine"), intervals = reg32, iterator = reg32) %>% 
    filter(d_ltr != 0 & d_line != 0)  %>% 
    arrange(intervalID) %>% 
    select(-intervalID) %>% 
    select(chrom, start, end)
large_regions <- giterator.intervals(intervals = gintervals.canonic(reg32), iterator = 524288)
regs <- reg32 %>% gintervals.neighbors(large_regions) %>%
    distinct(chrom1, start1, end1) %>%
    select(chrom = chrom1, start = start1, end = end1)
reg32 <- giterator.intervals(intervals = regs, iterator = 32) %>% 
    mutate(start = floor(start / 32) * 32, end = start + 32)

dim(reg32)
scales::comma(gintervals.covered_bp(reg32))

# + vscode={"languageId": "r"}
gen_k27_vt(mod)
gen_k4_vt(mod)

# + vscode={"languageId": "r"}
borzoi_vt <- c('ES_cnt','ES2i_cnt','EB4_cnt','e75_ecto_cnt','e75_emeso_cnt', "EB4_cnt_k4")
data_borzoi <- {
    d <- gextract(borzoi_vt, intervals = reg32, iterator = reg32)
    d <- d %>% 
        arrange(intervalID) %>% 
        select(-intervalID) %>% 
        group_by(chrom) %>% 
        mutate_at(borzoi_vt, ~rollmean(., 10, f='e')) %>% 
        ungroup()
} %cache_df% here("output/data_borzoi.tsv")
dim(data_borzoi)

# + vscode={"languageId": "r"}
bin_size <- 32L
seq_len  <- 65536L
n_bins   <- seq_len / bin_size  # 2048


# + vscode={"languageId": "r"}


# + vscode={"languageId": "r"}
seq_len <- 65536
d %>% gintervals.neighbors(giterator.intervals(intervals = gintervals.canonic(d), iterator=seq_len))

# + vscode={"languageId": "r"}
	k27_1k = rollmean(cg_trace$k27, 5, f='e')##on 200bp bins
	cg_trace$k27_1k = k27_1k
	cg_trace$lk27_1k = log2(2+k27_1k)

# + vscode={"languageId": "r"}
100/32

# + vscode={"languageId": "r"}
gen_k27_vt

# + vscode={"languageId": "r"}
gen_k4_vt

# + vscode={"languageId": "r"}
sum(regions$end - regions$start)

# + vscode={"languageId": "r"}
gintervals

# + vscode={"languageId": "r"}
cgd_e = mod$cgdom_ann %>% mutate(start = start - 2e6, end = end + 2e6)%>% 
    gintervals.force_range() 

# + vscode={"languageId": "r"}
regions_e <- bind_rows(cgd_e  %>% select(chrom, start, end), regions) %>% 
    gintervals.canonic() 

# + vscode={"languageId": "r"}
dim(regions_e)
sum(regions_e$end - regions_e$start)

# + vscode={"languageId": "r"}


# + vscode={"languageId": "r"}
gext = gextract(ndx$short_name,iterator = cg_trace_f,intervals = cg_trace_f)

gext[is.na(gext)] = 0




dim(cg_trace_f)
dim(cg_trace)


mat = gext[ ,ndx$short_name]

#mat_n = t(t(mat)/apply(mat, 2, median))
mat_n = (t(t(mat)/ndx$cov))*2e7
#mat_n2 = mat_n/rowMeans(mat_n)
mat_n2 = mat_n 
mat_n2l = log2(mat_n2+1)



df = cbind(gext[,c('chrom','start','end')],mat_n2l )

# + vscode={"languageId": "r"}
create_bw <- function(df, filename, column) {
    df <- as_tibble(df)
    score <- iceqream::norm01(df[[column]])
    #score[score==0]=1e-5
    #score <- (df[[column]])
    gr <- GRanges(
        seqnames = df$chrom,
        ranges = IRanges(start = df$start, end = df$end),
        score = score
    )
    sl <- gintervals.all() %>%
        select(chrom, end) %>%
        tibble::deframe()
    sl <- sl[as.character(levels(seqnames(gr)))]
    seqlengths(gr) <- sl

    export.bw(gr, filename)
}


# + vscode={"languageId": "r"}
dir.create(("/net/mraid20/export/tgdata/users/evghenic/Proj/ml_course/proj/data/it7/cgd-bw"), showWarnings = FALSE)

# + vscode={"languageId": "r"}



df$end = df$end-1 ##otherwise can't create bigwigs

df_f = df #%>% filter(EB4_cnt>3 & e75_emeso_cnt>3)

# + vscode={"languageId": "r"}
for (i in colnames(df_f[,-c(1,2,3)])){
    path = paste0("/net/mraid20/export/tgdata/users/evghenic/Proj/ml_course/proj/data/it7/cgd-bw/",i,'.bw')
    print(path)
    create_bw(df_f, path, i)

}

# + vscode={"languageId": "r"}


# + vscode={"languageId": "r"}


# + vscode={"languageId": "r"}

