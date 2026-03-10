
library(iceqream)
library(tidyverse)
library(misha)
library(misha.ext)
library(tgutil)
library(prego)
library(here)
traj_data <- readr::read_rds("/home/aviezerl/proj/motif_reg/enhflow/data/traj_data_with_dinucs.rds")
all_intervals <- traj_data %>% select(chrom:end, peak_name)

gsetroot("/home/aviezerl/mm10")
cell_types <- gsub("or_rk_sf.", "", gtrack.ls("^or_rk_sf"))
cell_types <- cell_types[cell_types != "all"]
cell_types <- grep("marginal", cell_types, value = TRUE, invert = TRUE)
cell_types <- grep("emb.", cell_types, value = TRUE, invert = TRUE)
atac_data <- preprocess_data(project_name = "or_rk_sf", cell_types = cell_types, anchor_cell_type = "epi", peaks = all_intervals, figures_dir = here("figures/ap-atac-cell-types"), const_quantile = 0.8)

atac <- atac_data$peaks %>% 
        mutate(epi = atac_data$atac_norm_prob[, "epi"], exe_meso = atac_data$atac_norm_prob[, "exe_meso"])
fwrite(atac, "data/atac_ap.csv")

mod<- readr::read_rds("data/cache_mod.RDS")
doms <- cgdd$cgdom_ann


source(here("code/seq2epi_utils.r"))
source(here("code/fig_fun.r"))
gsetroot("data/mm10")
gen_k27_vt(mod)
diff_mat <- gen_diff_mat(mod)
doms_diff <- cbind(doms, diff_mat)

atac %>% gintervals.neighbors(doms_diff) %>%
        mutate(dist = cut(abs(dist), breaks = c(0, 1e3, 2e3, 4e3, 1e4, 2e4, 5e4, 1e5, 2e5, 5e6), include.lowest = TRUE)) %>%
        ggplot(aes(x=dist, y=exe_meso - epi)) + geom_boxplot()


# energies <- readr::read_rds("/home/aviezerl/proj/motif_reg/output/after_oct7/epi_to_nas_revision/new_inter_energies_all_n.rds")

mdb <- create_motif_db(iceqream::motif_db, prior = 0.01)
pwm <- gextract_pwm(atac, dataset = mdb, bidirect = TRUE) %cache_rds% here("output/pwm_ap.rds")

mat <- pwm %>% select(-(chrom:end), -const, -tss_dist, -epi, -exe_meso) %>%
    remove_rownames() %>%
    column_to_rownames("peak_name") %>%
    as.matrix()

pwm_norm <- norm_energy_matrix(mat, mat, q = 0.995) %cache_rds% here("output/pwm_norm_ap.rds")


atac <- atac %>% mutate(type = ifelse(chrom %in% c("chr2", "chr8", "chr12", "chr18"), "test", "train")) %>% 
    gintervals.normalize(500)

norm_intervals <- atac %>% select(chrom:end, peak_name)

traj_model <- iq_regression(
        peak_intervals = atac %>% select(chrom:end, peak_name),
        atac_scores = atac %>% select(epi, exe_meso),
        motif_energies = pwm_norm,    
        normalize_energies = FALSE,
        norm_intervals = norm_intervals,
        seed = 60427,
        test_idxs = which(atac$type == "test"),
        train_idxs = which(atac$type == "train"),
        n_prego_motifs = 5,      
        peaks_size = 500,
        max_motif_num = 50,                
        energy_norm_quantile = 0.995,        
        distill_on_diff = TRUE,
        include_interactions = TRUE,
        pssm_db = iceqream::motif_db,
        max_n_interactions = 250,
        output_dir = here("output/iq-ap-epi-exe-meso-model")
    )

traj_model_f <- readr::read_rds(here("output/iq-ap-epi-exe-meso-model/iq_regression_filtered_model.rds"))
iq_model <- create_iq_model(traj_model_f)
atac_pred <- predict(iq_model, intervals = atac)
atac <- atac %>% mutate(pred = atac_pred, obs = exe_meso - epi)

atac %>% 
   group_by(type) %>%
   summarise(r2 = cor(obs, pred)^2)
fwrite(atac, here("data/atac_ap_pred.csv"))


motif_db <- purrr::imap_dfr(traj_model_f@motif_models, ~ .x$pssm %>% mutate(motif = .y)) %>%
        select(motif, pos, A:T)

fwrite(motif_db, here("data/motif_db_ap.tsv"), sep = "\t")

