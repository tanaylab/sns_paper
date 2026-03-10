

fig6_diff_cres = function(mod,nm1,nm2,nm1_a,nm2_a){
#####

#lm seeds
seeds = readRDS('./data/lm_10test_noX_seeds_logist_dinucs_logist_fig1.rds')


f_pcg = seeds$modality=='pcg'
f_txg = seeds$modality=='txg'
f_mix = seeds$modality=='mix'
seeds$pred_seed_k27 = seeds$k27_pred
seeds$pred_seed_k4 = seeds$k4_pred
seeds$k27_resp = seeds$eb4_k27_mean##!!!

seeds$modality_temp = ifelse(seeds$pred_seed_k27>6.4,'pcg','mix')
seeds$modality_temp  = ifelse(seeds$pred_seed_k4 >= 6 & seeds$pred_seed_k27 < 5.3,
                               'txg',seeds$modality_temp )
table(seeds$modality_temp)


f_x = mod$cgdom_ann$chrom == 'chrX'
mod$cgdom_ann = mod$cgdom_ann[!f_x,]
mod$cgdom_ann$pred_seed_k27 = seeds$pred_seed_k27
mod$cgdom_ann$pred_seed_k4 = seeds$pred_seed_k4
mod$cgdom_ann$pred = seeds$pred_seed_k4 - seeds$pred_seed_k27
mod$seqmod_loc_pred = mod$cgdom_ann$pred


mod$epi_tracks = mod$epi_tracks_all


gen_k4_vt(mod)
gen_track_vext()

gvtrack.create('EB5_epi_cnt_n', 'jk.epipcg.pcg.CRJK_0293_k27me3_j1_eb_d5_epi_norm', 
        "avg")
gvtrack.iterator('EB5_epi_cnt_n', sshift = -500, eshift = 500)


gvtrack.create('EB5_meso_cnt_n', 'jk.epipcg.pcg.CRJK_0294_k27me3_j1_eb_d5_meso_norm', 
        "avg")
gvtrack.iterator('EB5_meso_cnt_n', sshift = -500, eshift = 500)


gvtrack.create('e75_ecto_cnt_n', 'jk.epipcg.pcg.CRJK_0321_k27me3_e75_ecto_norm', 
        "avg")
gvtrack.iterator('e75_ecto_cnt_n', sshift = -500, eshift = 500)


gvtrack.create('e75_emeso_cnt_n', 'jk.epipcg.pcg.CRJK_0324_k27me3_e75_e_meso_norm', 
        "avg")
gvtrack.iterator('e75_emeso_cnt_n', sshift = -500, eshift = 500)

names_k27 = c('EB5_epi_cnt_n','EB5_meso_cnt_n','e75_ecto_cnt_n','e75_emeso_cnt_n')

peaks_eb_cen = peaks_eb()

dim(peaks_eb_cen)

peaks_invivo = as.data.frame(fread('./data/atac_ap.csv'))

tracks_atac = c('jk.epipcg.multieb24.mcEpiblast','jk.epipcg.multieb24.mcn_meso',
                'jk.epipcg.multieb24.mcGut',
                'jk.epipcg.atac.CRJK_0389_atac_wt_to_wt_eb_d3',
                'jk.epipcg.atac.CRJK_0373_atac_wt_to_wt_epi_eb_d5',
'jk.epipcg.atac.CRJK_0377_atac_wt_to_wt_meso_eb_d5',
'jk.epipcg.atac.CRJK_0376_atac_wt_to_wt_endo_eb_d5',
               'mut_multi.WT_E75.epi',
              'mut_multi.WT_E75.decto',
              'mut_multi.WT_E75.exe_meso' )
names(tracks_atac) = c('multeb_epi','multeb_meso','multeb_endo','atac_eb3','atac_epi_d5','atac_meso_d5','atac_endo_d5','atac_e75_epi','atac_e75_decto','atac_e75_exe_meso')


####
    
    
    
    
    
    
    
if (nm1 == 'e75_ecto_cnt_n' | nm2 =='e75_ecto_cnt_n'){
if (file.exists('./data/atac_mat_invivo_fig6.rds')) {
        atac_mat = readRDS('./data/atac_mat_invivo_fig6.rds')
       
    }else {
atac_mat = gen_atac_mat(tracks_atac =tracks_atac ,atac_peaks = peaks_invivo)

saveRDS(atac_mat,'./data/atac_mat_invivo_fig6.rds')
        }

message('./data/atac_mat_invivo_fig6.rds')
}else{

if (file.exists('./data/atac_mat_meeb_fig6.rds')) {
        atac_mat = readRDS('./data/atac_mat_meeb_fig6.rds')
       
    }else {
atac_mat = gen_atac_mat(tracks_atac =tracks_atac ,atac_peaks = peaks_eb_cen)

saveRDS(atac_mat,'./data/atac_mat_meeb_fig6.rds')
        }

message('./data/atac_mat_meeb_fig6.rds')
}




if (file.exists('./data/ribo_cov_fig6.rds')) {
        ribo_cov = readRDS('./data/ribo_cov_fig6.rds')
       
    }else {

ribo_cov = fig5_ribo_cov_atac(mod=mod,tracks_atac=tracks_atac)
saveRDS(ribo_cov,'./data/ribo_cov_fig6.rds')
        }


cov_vect = fig5_track_cov()

if (file.exists('./data/all_seed_k27_fig6.rds')) {
        all_seed_k27 = readRDS('./data/all_seed_k27_fig6.rds')
       
    }else {
all_seed_k27 = fig5_gext_all_tracks(mod = mod)
saveRDS(all_seed_k27,'./data/all_seed_k27_fig6.rds')
        }




if (file.exists('./data/quant_fig6.rds')) {
        quant = readRDS('./data/quant_fig6.rds')
       
    }else {
cg_trace <- readRDS(here("data/cg_trace_mm10.rds"))

cg_trace <- gextract.left_join("mapab.umap_k100", intervals = cg_trace, iterator = cg_trace)



cg_trace$dist_black =  cg_trace %>%
        select(chrom, start, end) %>%
        gintervals.neighbors("mapab.Encode_blacklist_v2") %>% select(dist) %>% pull()

f_norp <- cg_trace$d_ltr != 0 &
              cg_trace$d_line != 0 &            
              cg_trace$start > 3e6 &
              !is.na(cg_trace$mapab.umap_k100) &
              cg_trace$mapab.umap_k100 > 0.9 &
               cg_trace$dist_black !=0 
cg_trace_f = cg_trace[ f_norp,]

quant = data.frame()
for (tr in names_k27){
    q = gquantiles(tr,c(0.2,.5,.98,.99,.995),iterator = 20,intervals = cg_trace_f)
    quant = rbind(quant,q)
}

rownames(quant) = names_k27
colnames(quant) = c(0.2,.5,.98,.99,.995)
saveRDS(quant,'./data/quant_fig6.rds')
}



plt = fig5_create_plt_mat(all_seed_k27=all_seed_k27,nm1 = nm1,nm2 = nm2,quant =quant )

#emeso
gg = plt %>% filter(col2 %in% c('k4',paste0(nm2,"_loss"),'cons')) %>% ggplot(aes(x = col2,y=pred_seed_k27)) +
geom_boxplot(outlier.shape = NA)+ theme_bw()+
coord_cartesian(ylim = c(4.7,7))
print(gg)
#save_baseR_to_ppt(plot(1,1),'./figs/fig5_iq_pred_vs_pcg_loss_emeso.pptx')
#save_gg_to_ppt(gg=gg,link = './figs/fig5_iq_pred_vs_pcg_loss_emeso.pptx')

#meso no X final
options( repr.plot.width=8,repr.plot.height=8)
gg = plt %>%
  ggplot(aes(x = track1, y = track2,col=col)) +
  geom_point(size = 0.8,alpha=1) +
  labs(
    title = paste(
      nm2, "loss =", as.numeric(table(plt$col)[1]), "|",
      nm1, "loss =", as.numeric(table(plt$col)[2]),
        ' cor = ',cor(plt[,'track1'],plt[,'track2'],m='p')
    ),
    x = nm1,
    y = nm2) +
#geom_abline(slope = 1)+
 coord_cartesian(xlim = c(3.5,9),ylim = c(3.5,9))+
theme_bw()+scale_color_manual(values = c('gray','turquoise','gray'))+
  theme(legend.position = "none")
print(gg)

colnames(plt) = c(nm1,nm2,'pred','pred_seed_k27','pred_seed_k4','col','col2')

seeds=cbind(all_seed_k27[, c('chrom','start','end')],plt)

#write.csv(seeds,'./data/fig5_meeb_meso_epi_seeds.csv')



atac_normed = fig5_norm_atac(atac_mat = atac_mat,nm1_a = nm1_a,nm2_a = nm2_a,ribo_cov = ribo_cov)


#options( repr.plot.width=8,repr.plot.height=8)
#atac_normed %>%
#  ggplot(aes(x = track1_a, y = track2_a)) +
# coord_cartesian(xlim = c(2,12),ylim = c(2,12))+
#  geom_point(size = 0.7) +
#  labs(
   
 #   x = nm1_a,
 #   y = nm2_a) +
#geom_abline(slope = 1)+

#theme_bw()+scale_color_manual(values = c('gold','darkgreen','black'))+
#  theme(legend.position = "none")


#write.csv(seeds %>% filter(col=='e75_ecto_cnt_n_loss'),'./data/tableS3_e75_ecto_loss.csv')

#write.csv(seeds %>% filter(col=='EB5_epi_cnt_n_loss'),'./data/tableS3_EBecto_loss.csv')

#emeso
gg_list = cres_near_pcg_seed(seeds = seeds,nm1 =nm1 ,nm2 = nm2,atac_normed = atac_normed,
                        more = more,delta=delta,nm1_a = nm1_a ,nm2_a = nm2_a)

gg = gg_list$gg
save_baseR_to_ppt(plot_func = plot(1,1),link_ppt = paste0('./figs/cres_near_',nm1,'_to_',nm2,
                                                          '_atac_',nm1_a,'_',nm2_a,'.pptx'))
save_gg_to_ppt(gg = gg, link = paste0('./figs/cres_near_',nm1,'_to_',nm2,
                                                          '_atac_',nm1_a,'_',nm2_a,'.pptx'))
 print(gg)   
}

plot_pwm_cluster_heatmap_final <- function(
  hc,
  cl,
  nhit,
  d2,
  ehit,
  other_motifs_clean,
  figure_motifs_clean
) {

  other_motifs_clean = final_figure_motifs_clean
  set.seed(42)
  sps(24,24)
  ord <- hc$order
  dist_mat <- (nhit - d2)
  colnames(d2) <- colnames(ehit)
  rownames(d2) <- colnames(ehit)
  colnames(dist_mat) <- colnames(ehit)
  rownames(dist_mat) <- colnames(ehit)

  dist_mat_ord <- dist_mat[ord, ord, drop = FALSE]

  ann <- data.frame(cluster = factor(cl))
  rownames(ann) <- names(cl)
  ann_ord <- ann[ord, , drop = FALSE]
  nr <- nrow(dist_mat_ord)
  nc <- ncol(dist_mat_ord)

  motifs20 = readRDS('./data/motifs20_fig1.rds')

  rep_idx <- tapply(seq_len(nr), ann_ord$cluster, function(ix) {
    hits <- ix[ rownames(dist_mat_ord)[ix] %in% motifs20 ]
    if (length(hits)) hits[1] else ix[1]
  })
  rep_idx <- as.integer(rep_idx)
  labels_row <- rep("", nr)
  labels_col <- rep("", nc)

  labels_row[rep_idx] <- rownames(dist_mat_ord)[rep_idx]

  hit_c <- which(colnames(dist_mat_ord) %in% figure_motifs_clean)
  labels_col[hit_c] <- colnames(dist_mat_ord)[hit_c]

  cl_order = unique(cl[hc$order])

  if (any(hit_c)) labels_col[hit_c] <- paste0("bold(\"", labels_col[hit_c], "\")")
  labels_row <- parse(text = ifelse(labels_row == "", "''", labels_row))
  labels_col <- parse(text = ifelse(labels_col == "", "''", labels_col))

  clust_cols <- grDevices::hcl.colors(
    n = 20,
    palette = "Dark 3"
  )
  set.seed(42)
  clust_cols = sample(clust_cols)
  names(clust_cols) <- cl_order
  names(clust_cols) <- levels(ann_ord$cluster)
  ann_colors <- list(
    cluster = clust_cols
  )

  maxv <- max(nhit - as.matrix(dist_mat_ord), na.rm = TRUE)

  white_frac <- 0.5
  white_max  <- maxv * white_frac

  bk <- c(
    seq(0, white_max, length.out = 50),
    seq(white_max, maxv, length.out = 151)
  )

  bk <- unique(bk)

  cols <- colorRampPalette(c("white", "red", "darkred", "yellow"))(length(bk) - 1)

  ann_ord$cluster <- factor(as.character(ann_ord$cluster), levels = cl_order)

  mat = nhit - as.matrix(dist_mat_ord)
  gap_w <- 5

  gaps <- which(ann_ord$cluster[-1] != ann_ord$cluster[-nrow(mat)])

  for (g in rev(gaps)) {
    ins <- matrix(NA_real_, gap_w, ncol(mat)); rownames(ins) <- rep("", gap_w)
    mat <- rbind(mat[1:g, , drop=FALSE], ins, mat[(g+1):nrow(mat), , drop=FALSE])

    ins2 <- matrix(NA_real_, nrow(mat), gap_w); colnames(ins2) <- rep("", gap_w)
    mat <- cbind(mat[, 1:g, drop=FALSE], ins2, mat[, (g+1):ncol(mat), drop=FALSE])
  }

  p <- pheatmap(mat, na_col="white", cluster_rows=FALSE, cluster_cols=FALSE, border_color=NA,
    breaks = bk, color = cols,
    labels_row = labels_col,
    labels_col = FALSE,
    fontsize_row = 24,fontsize_col = 24,
    annotation_row = ann_ord,
    annotation_col = ann_ord,
    annotation_colors = ann_colors,
    main = sprintf("Distance used for clustering: nhit - overlap (h=%.2f)", nhit*0.6))

  return(p)
}

eval_linear_models = function(test_lm,resp){
#BEST big
modelbig_best = readRDS('./data/lm_10test_noX_seeds_logist_dinucs_logist_fig1.rds')
pred_big_best = modelbig_best$k4_pred - modelbig_best$k27_pred
#message(cor(pred_big_best[test_lm], resp[test_lm],m="p")**2) 
#message(cor(pred_big_best[!test_lm], resp[!test_lm],m="p")**2)
resp

#BEST 20 motifs
model20_best = readRDS('./data/lm_10test_noX_seeds_20mtfs_logist_sqrdi_fig1.rds')
pred_20_best = model20_best$k4_pred - model20_best$k27_pred
#message(cor(pred_20_best[test_lm], resp[test_lm],m="p")**2) 
#message(cor(pred_20_best[!test_lm], resp[!test_lm],m="p")**2)

#20 motifs center
model20_cent = readRDS('./data/lm_10test_noX_seeds_20mtfs_logist_sqrdi_center_fig1.rds')
pred_20_cent = model20_cent$k4_pred - model20_cent$k27_pred
#message(cor(pred_20_cent[test_lm], resp[test_lm],m="p")**2) 
#message(cor(pred_20_cent[!test_lm], resp[!test_lm],m="p")**2)

rsqr_te = c(cor(pred_big_best[test_lm], resp[test_lm],m="p")**2,
           cor(pred_20_best[test_lm], resp[test_lm],m="p")**2,
           cor(pred_20_cent[test_lm], resp[test_lm],m="p")**2)
rsqr_tr = c(cor(pred_big_best[!test_lm], resp[!test_lm],m="p")**2,
           cor(pred_20_best[!test_lm], resp[!test_lm],m="p")**2,
           cor(pred_20_cent[!test_lm], resp[!test_lm],m="p")**2)

df = data.frame(rsqr_te=rsqr_te,rsqr_tr=rsqr_tr, model = c('big','mtfs20','mtfs20_cent'))

df_long <- df %>%
  pivot_longer(cols = c(rsqr_tr, rsqr_te),
               names_to = "split", values_to = "rsq") %>%
  mutate(split = recode(split, rsqr_tr = "train", rsqr_te = "test"))

library(ggplot2)
sps(8,8)
pd <- position_dodge(width = 0.75)

gg = ggplot(df_long, aes(x = model, y = rsq, fill = split)) +
  geom_col(position = pd, width = 0.7) +
  geom_text(
    aes(label = sprintf("%.2f", rsq)),
    position = pd,
    vjust = -0.4,
    size = 5
  ) +
  labs(
    x = NULL,
    y = expression(R^2),
    fill = "split",
    title = "Train vs test R² by model"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    plot.title   = element_text(size = 20),
    axis.title.y = element_text(size = 18),
    axis.text.x  = element_text(size = 16, angle = 30, hjust = 1),
    axis.text.y  = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text  = element_text(size = 14)
  ) +
  # ensures labels above bars aren't clipped
  coord_cartesian(clip = "off") +
  theme(plot.margin = margin(10, 10, 10, 10))
print(gg)
#save_gg_to_pptx(gg = gg,path = './figs/rmodels_performances_sup1.pptx')
}

gen_cg_corplot_gg = function(mod) {
  stopifnot(requireNamespace("ggplot2", quietly = TRUE))
  # reshape helper; you can switch to tidyr if you prefer
  if (!requireNamespace("reshape2", quietly = TRUE)) {
    stop("Please install either 'reshape2' (or edit to use 'tidyr::pivot_longer').")
  }

  set.seed(42)
  ndx = mod$epi_tracks
  ndx = ndx[ !ndx$track_k27 %in% c('jk.epipcg.pcg.CRJK_0361_k27me3_e7_meso_a70ls',
                                   'jk.epipcg.lit.encode.es.ENCFF595SIA_H3K27me3_es'), ]

  gvtrack.create("CG2k", "seq.CG_500_mean_new", "avg")
  gvtrack.iterator("CG2k", sshift=-1000, eshift=1000)
  gvtrack.create("CG1k", "seq.CG_500_mean_new", "avg")
  gvtrack.iterator("CG1k", sshift=-500, eshift=500)

  doms = mod$k27_doms
  doms = doms[ !doms$chrom %in% c('chrX','chrM','chrY'), ]
  samp_d = doms[sample(1:nrow(doms), 1000), ]
  samp_d$start = samp_d$start - 3000
  samp_d$end   = samp_d$end   + 3000

  prf = gextract(c(ndx$short_name, "seq.CG_500_mean_new", "CG1k", "CG2k"),
                 iterator = 20, intervals = samp_d)

  prf_v  = prf[, ndx$short_name, drop = FALSE]
  lprf_v = log2(.1 + t(t(prf_v) / colMeans(prf_v)))#0.1

  # CG class masks from the raw CG profile
  f_hcg = prf$seq.CG_500_mean_new > 0.05
  f_lcg = prf$seq.CG_500_mean_new < 0.02
  f_icg = !f_hcg & !f_lcg

  # Build a single factor for class
  cg_class = ifelse(f_hcg, "High CG",
                    ifelse(f_lcg, "Low CG", "Intermediate CG"))
  cg_class = factor(cg_class, levels = c("Low CG", "Intermediate CG", "High CG"))

  # Long format for ggplot
  ldf = as.data.frame(lprf_v)
  ldf$CG_class = cg_class
  ldf_long = reshape2::melt(ldf, id.vars = "CG_class",
                            variable.name = "Track", value.name = "value")

  # Plot with ggplot2
p = ggplot2::ggplot(ldf_long, ggplot2::aes(x = value, color = CG_class)) +
  ggplot2::geom_density() +
  ggplot2::scale_color_manual(values = c("Low CG" = "black",
                                         "Intermediate CG" = "orange",
                                         "High CG" = "darkred")) +
  ggplot2::coord_cartesian(xlim = c(-3, 6.5), ylim = c(0, 0.3)) +#.3 #6.5
  ggplot2::facet_wrap(~ Track, nrow = 1, scales = "fixed") +
  ggplot2::labs(x = NULL, y = NULL,
                title = "Density of log2-normalized signal by CG class") +
  ggplot2::theme_bw(base_size = 10) +
  ggplot2::theme(
    panel.grid.major = ggplot2::element_line(size = 0.2),
    panel.grid.minor = ggplot2::element_blank(),
    legend.title = ggplot2::element_blank(),
    legend.position = "top",
    # --- Facet title styling ---
    strip.background = ggplot2::element_blank(),   # removes the box
    strip.text = ggplot2::element_text(size = 9, face = "bold"),
    panel.spacing = grid::unit(0, "lines")         # no gaps between plots
  )

 p = p + ggplot2::theme(panel.spacing = grid::unit(0, "lines"))
  print(p)
  #invisible(p)
}


generate_cg_cor = function(mod){
gvtrack.create("CG15k", "seq.CG_500_mean_new", "avg")
gvtrack.iterator("CG15k", sshift=-750,eshift=750)
gvtrack.create("CG2k", "seq.CG_500_mean_new", "avg")
gvtrack.iterator("CG2k", sshift=-1000,eshift=1000)
gvtrack.create("CG1k", "seq.CG_500_mean_new", "avg")
gvtrack.iterator("CG1k", sshift=-500,eshift=500)
gvtrack.create("CG500", "seq.CG_500_mean_new", "avg")

cg = mod$cgdom_ann

cg$dist_cgd = misha.ext::gintervals.neighbors1(cg,cg,maxneighbors = 1, mindist = 1)$dist



genome = gintervals.all()
#foc_chroms = genome[ genome$chrom %in% c('chr1','chr2','chr18','chr10','chr11'),]

foc_chroms = genome[ genome$chrom %in% c('chr1'),]

stats = plt_cg_histone_trends(foc_chroms=foc_chroms,mod=mod)

colnames(stats) = c('c_k27', 'c_k272', 'r2_k27', 'c_k4', 'c_k42', 'r2_k4', 'c_g_k27', 'c_g_k272', 'c_g_k4', 'c_g_k42')

rownames(stats) = seq(0,5000,50)

gg = plot_cor_to_cg(vars = c('c_k27','c_k4'),stats = stats)
print(gg)
#save_gg_to_ppt(gg,"./figs/cg_cor.pptx")
}

gen_stats_fig1a = function(){
gvtrack.create("CG15k", "seq.CG_500_mean_new", "avg")
gvtrack.iterator("CG15k", sshift=-750,eshift=750)
gvtrack.create("CG2k", "seq.CG_500_mean_new", "avg")
gvtrack.iterator("CG2k", sshift=-1000,eshift=1000)
gvtrack.create("CG1k", "seq.CG_500_mean_new", "avg")
gvtrack.iterator("CG1k", sshift=-500,eshift=500)
gvtrack.create("CG500", "seq.CG_500_mean_new", "avg")

cg = mod$cgdom_ann

cg$dist_cgd = misha.ext::gintervals.neighbors1(cg,cg,maxneighbors = 1, mindist = 1)$dist



genome = gintervals.all()
#foc_chroms = genome[ genome$chrom %in% c('chr1','chr2','chr18','chr10','chr11'),]

foc_chroms = genome[ genome$chrom %in% c('chr1'),]

stats = plt_cg_histone_trends(foc_chroms=foc_chroms,mod=mod)

colnames(stats) = c('c_k27', 'c_k272', 'r2_k27', 'c_k4', 'c_k42', 'r2_k4', 'c_g_k27', 'c_g_k272', 'c_g_k4', 'c_g_k42')

rownames(stats) = seq(0,5000,50)
    return(stats)
    }
plot_cor_to_cg = function(vars = c('c_g_k27','c_g_k4'),stats){

stats_gg = stats %>% as.data.frame()%>%rownames_to_column(var = 'scale')

stats_gg$scale = factor(levels = stats_gg$scale,stats_gg$scale)



options( repr.plot.width=15,repr.plot.height=10)
gg = stats_gg[,c('scale',vars)] %>% pivot_longer(-1) %>% ggplot(aes(x=scale, y = (value), group = name,col=name))+
coord_cartesian(ylim = c(0.1,.62))+#0.25,.62
scale_color_manual(values = c('darkblue','darkred'))+
geom_line(size=.7)+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                                         panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())
    print(gg)
    return(gg)
    }

plt_cg_histone_trends = function(foc_chroms,mod)
{
    
   # tracks = c(
    #'jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm',
    #'jk.epipcg.pcg.CRJK_0403_k27me3_wt_to_wt_eb_d3',
    #'jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3')
	#tracks = c('EB4_cnt','EB4_cnt_k4')
    short_nms = c('k27','k4')
	#for(i in tracks) {
	#	gvtrack.create(gsub('jk.epipcg.','',i), i, "sum") 
	#	gvtrack.iterator(gsub('jk.epipcg.','',i), sshift = -140, eshift= 140)
	#}

    short_nms = c('k27','k4')

		gvtrack.create('k27', 'jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm', "avg") 
		gvtrack.iterator('k27', sshift = -500, eshift= 500)

		gvtrack.create('k4', 'jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3', "sum") 
		gvtrack.iterator('k4', sshift = -140, eshift= 140)
	gvtrack.create("CGh", "seq.CG", "sum")
	gvtrack.create("GCh", "seq.G_or_C", "sum")

	dst = gextract(short_nms, iterator=200, intervals=foc_chroms,colnames = c('k27','k4'))
	dst = dst[dst$start > 3e+6,] #remove prefix of NNNN
	lk4 = log2(6+dst$k4)#2
	lk27 = dst$k27
	stats = matrix(nrow=0,ncol=10)


	for(h in seq(0,5000,50)) { 
		gvtrack.iterator("CGh", sshift=-h, eshift=h); 
		gvtrack.iterator("GCh", sshift=-h, eshift=h); 
		cg_prf = gextract(c("CGh","GCh"), iterator=dst, intervals=dst); 
		cg_prf$CGh[is.na(cg_prf$CGh)]=0; 
		c_k27 = cor(lk27, cg_prf$CGh, u="p")
		c_k272 = cor(lk27, cg_prf$CGh**2, u="p")
		c_k4 = cor(lk4, cg_prf$CGh, u="p")
		c_k42 = cor(lk4, cg_prf$CGh**2, u="p")

		c_g_k27 = cor(lk27, cg_prf$GCh, u="p")
		c_g_k272 = cor(lk27, cg_prf$GCh**2, u="p")
		c_g_k4 = cor(lk4, cg_prf$GCh, u="p")
		c_g_k42 = cor(lk4, cg_prf$GCh**2, u="p")

		cg_prf[is.na(cg_prf)] = 0
	
		cg = cg_prf$CGh
		cg2 = cg**2	
		cg3 = cg**3	
		r2_k4 = summary(lm(lk4 ~ cg + cg2 +cg3))$r.squared
		r2_k27 = summary(lm(lk27 ~ cg + cg2 +cg3))$r.squared

		stats = rbind(stats, c(c_k27, c_k272, r2_k27, c_k4, c_k42, r2_k4, c_g_k27, c_g_k272, c_g_k4, c_g_k42))
		message("at ", h)
	}
    return(stats)
}

domain_cov_stat = function(cgd,len_min,len_max){

cgd = mod$cgdom_ann
cgd = cgd[ cgd$l > len_min  & cgd$l < len_max,]
lk27 = cgd$eb4_k27_mean
lk4 = cgd$l6_eb4_k4_mean
fgood = lk4 > 5 | lk27 > 6
cgd = cgd[ fgood,]
n_all = nrow(cgd)
n = nrow(cgd[ cgd$eb4_k27_cov_1>0 & cgd$eb4_k4_cov_1 ==0 & cgd$eb4_biv_cov_1==0,])
message('k27 n = ',n,' % = ', (n/n_all)*100)
n = nrow(cgd[ cgd$eb4_k27_cov_1>0 & cgd$eb4_k4_cov_1 ==0 & cgd$eb4_biv_cov_1>0,])
message('biv k27 n = ',n,' % = ', (n/n_all)*100)
n = nrow(cgd[ cgd$eb4_k27_cov_1==0 & cgd$eb4_k4_cov_1 ==0 & cgd$eb4_biv_cov_1>0,])
message('biv  n = ',n,' % = ', (n/n_all)*100)
n = nrow(cgd[ cgd$eb4_k27_cov_1==0 & cgd$eb4_k4_cov_1 > 0 & cgd$eb4_biv_cov_1>0,])
message('k4 biv  n = ',n,' % = ', (n/n_all)*100)
n = nrow(cgd[ cgd$eb4_k27_cov_1>0 & cgd$eb4_k4_cov_1 > 0 & cgd$eb4_biv_cov_1>0,])
message('k4 k27 biv  n = ',n,' % = ', (n/n_all)*100)
n = nrow(cgd[ cgd$eb4_k27_cov_1==0 & cgd$eb4_k4_cov_1 > 0 & cgd$eb4_biv_cov_1==0,])
message('k4   n = ',n,' % = ', (n/n_all)*100)
n = nrow(cgd[ cgd$eb4_k27_cov_1==0 & cgd$eb4_k4_cov_1 == 0 & cgd$eb4_biv_cov_1==0,])
message('else   n = ',n,' % = ', (n/n_all)*100)
    }

pcg_plot_energ_heats_legends = 
function (cgdd, len_min, len_max, base_dir = "figs/") 
{
    cgdd = cgdd[cgdd$l > len_min & cgdd$l < len_max, ]
	lk27 = cgdd$eb4_k27_mean
	lk4 = cgdd$l6_eb4_k4_mean
    fgood = lk4 > 5 | lk27 > 6
    grad_score = (lk27 - lk4)[fgood]
    bad_score = (lk27 - lk4)[!fgood]
    shades = colorRampPalette(c("white", "gray", "red", "black"))(1000)
    dst = cgdd[fgood, c("eb4_k27_cov_1", "eb4_biv_cov_1", "eb4_k4_cov_1")]
	
    dst_n = dst/rowSums(dst)
	dst_n0 = dst_n
	dst_n0[is.na(dst_n0)]=0
	grad_score = (1e-3+dst_n0[,1]+dst_n0[,2])/(1e-3+dst_n0[,3]+dst_n0[,2])
	
    dst_n = dst_n[order(grad_score, decreasing = F), ]
    dst_n_sm = apply(dst_n, 2, zoo::rollmean, 20, na.rm = T, f = "e")
    n = nrow(dst_n)
    plot(-0, -0, xlim = c(0, 1), ylim = c(1, nrow(dst_n)), xaxt = "n", 
        yaxt = "n")
    polygon(x = c(dst_n_sm[, 3], rep(0, t = n)), y = c(1:n, n:1), 
        col = "darkred", border = NA)
    polygon(x = c(dst_n_sm[, 3], rev(dst_n_sm[, 2] + dst_n_sm[, 
        3])), y = c(1:n, n:1), col = "black", border = NA)
    polygon(x = c(rep(1, t = n), rev(dst_n_sm[, 3] + dst_n_sm[, 
        2])), y = c(1:n, n:1), col = "darkblue", NA)
    par(mar = c(0, 0, 0, 0))
    barplot(zoo::rollmean(lk27[fgood][order(grad_score, decreasing = F)], 
        150, f = "e") - 2, border = "darkblue", col = "darkblue", 
        horiz = T, xaxt = "n", yaxt = "n")
    abline(v = 2, lty = 2, lwd = 2)
    abline(v = 4, lty = 2, lwd = 2)
}



generate_len_plots = function(mod){
cgd = mod$cgdom_ann
cgdr = mod$cgd_rpt
cgd$l = cgd$l+500
cgdr$l = cgdr$l + 500

plt_cgd_repeats = function(cgd,cgdr){
plot(density(log2(cgd$l)), lwd=1, ylim=c(0,3.2), xlim=c(9.5,12), col="darkgreen")
	lines(density(log2(cgdr$l[cgdr$line>0.5])), lwd=.7, col="darkgray")
	lines(density(log2(cgdr$l[cgdr$ltr>0.5])), lwd=.7, col="darkgray",lty=2)
    }

plt_cgd_repeats(cgd = cgd,cgdr = cgdr)

save_baseR_to_ppt(plt_cgd_repeats(cgd = cgd,cgdr = cgdr),'./figs/cgdds_len_repeats.pptx')

lbin = as.numeric(cut(cgd$l,c(0,seq(800,3000,200))))

save_baseR_to_ppt(boxplot(split(cgd$cg_max*4/(cgd$gc_max**2), lbin), outline=FALSE,col="gray",ylim = c(0.4,1.1))
                  ,'./figs/cgdd_len_cg_gc.pptx')
boxplot(split(cgd$cg_max*4/(cgd$gc_max**2), lbin), outline=FALSE,col="gray",ylim = c(0.4,1.1))
        
save_baseR_to_ppt(boxplot(split(log2(100+abs(cgd$tss_dist)), lbin),outline=FALSE, col="gold",ylim = c(6.6,19)),
                  './figs/cgdd_len_tss_dist.pptx')
boxplot(split(log2(100+abs(cgd$tss_dist)), lbin),outline=FALSE, col="gold",ylim = c(6.6,19))
save_baseR_to_ppt(boxplot(split(cgd$eb4_k27_mean, lbin),outline=FALSE,ylim = c(4,8), col="darkblue"),
                  './figs/cgdd_len_k27.pptx')
boxplot(split(cgd$eb4_k27_mean, lbin),outline=FALSE, col="darkblue",ylim = c(4,8))
#boxplot(split(cgd$k27_f1, lbin),outline=FALSE, col="pink")
#boxplot(split(log2(1+cgd$k27_chip_f1), lbin),outline=FALSE, col="violet")
save_baseR_to_ppt(boxplot(split(log2(6+cgd$eb4_k4_mean), lbin),outline=FALSE, col="darkred"),
                  './figs/cgdd_len_k4.pptx')
boxplot(split(log2(6+cgd$eb4_k4_mean), lbin),outline=FALSE, col="darkred")
    }

generate_len_plots_hg = function(mod){
cgd = mod$cgdom_ann
cgdr = mod$cgd_rpt
cgd$l = cgd$l+500
cgdr$l = cgdr$l + 500
plt_cgd_repeats = function(cgd,cgdr){
plot(density(log2(cgd$l)), lwd=1, ylim=c(0,3.2), xlim=c(9.5,12), col="darkgreen")
	lines(density(log2(cgdr$l[cgdr$line>0.5])), lwd=.7, col="darkgray")
	lines(density(log2(cgdr$l[cgdr$ltr>0.5])), lwd=.7, col="darkgray",lty=2)
    }

plt_cgd_repeats(cgd = cgd,cgdr = cgdr)

save_baseR_to_ppt(plt_cgd_repeats(cgd = cgd,cgdr = cgdr),'./figs/cgdds_len_repeats_hg.pptx')

lbin = as.numeric(cut(cgd$l,c(0,seq(800,3000,200))))

save_baseR_to_ppt(boxplot(split(cgd$cg_max*4/(cgd$gc_max**2), lbin), outline=FALSE,col="gray",ylim = c(0.4,1.1))
                  ,'./figs/cgdd_len_cg_gc_hg.pptx')
boxplot(split(cgd$cg_max*4/(cgd$gc_max**2), lbin), outline=FALSE,col="gray",ylim = c(0.4,1.1))
        
save_baseR_to_ppt(boxplot(split(log2(100+abs(cgd$tss_dist)), lbin),outline=FALSE, col="gold",ylim = c(6.6,19)),
                  './figs/cgdd_len_tss_dist_hg.pptx')
boxplot(split(log2(100+abs(cgd$tss_dist)), lbin),outline=FALSE, col="gold",ylim = c(6.6,19))
save_baseR_to_ppt(boxplot(split(cgd$l10_h1_k27_sum, lbin),outline=FALSE, col="darkblue"),
                  './figs/cgdd_len_k27_hg.pptx')
boxplot(split(cgd$l10_h1_k27_sum, lbin),outline=FALSE, col="darkblue")
save_baseR_to_ppt(boxplot(split(log2(6+cgd$h1_k4_sum), lbin),outline=FALSE, col="darkred"),
                  './figs/cgdd_len_k4_hg.pptx')
boxplot(split(log2(6+cgd$h1_k4_sum), lbin),outline=FALSE, col="darkred")
    }

plt_k27_k4_cov = function(){
	cgd = mod$cgdom_ann
	cgdr = mod$cgd_rpt


lk27 = log2(1+cgd$eb4_k27_max)
	lk4 = log2(1+cgd$eb4_k4_max)
	k27_rat = cgd$eb4_k27_rat
	k4_rat = cgd$eb4_k4_rat
	#l_class = as.numeric(cut(cgd$l, c(-1,500,750,1000,1e+6)))
    l_class = as.numeric(cut(cgd$l, c(-1,700,1000,1e+6)))

k27cov = cgd$eb4_k27_cov_1
			k4cov = cgd$eb4_k4_cov_1
			bivcov = cgd$eb4_biv_cov_1
			dst = cgd[,c("eb4_k27_cov_1", "eb4_biv_cov_1", "eb4_k4_cov_1")]




layout(matrix(1:1, ncol=1))
		f_out = lk27 < 6 & lk4 < 6
		dst_n = dst/rowSums(dst)
		for(li in 2:2) {
			par(mar = c(2,2,1,1))
			f = l_class==li & !f_out

			dst_i = dst_n[f,]

			n = nrow(dst_i)
			ord = order((1e-3+dst_i[,1]+dst_i[,2])/(1e-3+dst_i[,3]+dst_i[,2]))
			dst_i = dst_i[ord,]
			plot(-0,-0, xlim=c(0,1), ylim=c(1,nrow(dst_i)), xaxt='n',yaxt='n')#
			polygon(x = c(dst_i[,3], rep(0,t=n)), y=c(1:n, n:1), col="darkred", border=NA)
			polygon(x = c(dst_i[,3], rev(dst_i[,2]+ dst_i[,3])), y=c(1:n, n:1), col="black", border=NA)
			polygon(x = c(rep(1,t=n), rev(dst_i[,3]+dst_i[,2])), y=c(1:n, n:1), col="darkblue", NA)
		}
   
    }


plt_k27_k4_cov_995q = function(){
	cgd = mod$cgdom_ann
	cgdr = mod$cgd_rpt


lk27 = log2(1+cgd$eb4_k27_max)
	lk4 = log2(1+cgd$eb4_k4_max)
	k27_rat = cgd$eb4_k27_rat
	k4_rat = cgd$eb4_k4_rat
	#l_class = as.numeric(cut(cgd$l, c(-1,500,750,1000,1e+6)))
    l_class = as.numeric(cut(cgd$l, c(-1,700,1000,1e+6)))

k27cov = cgd$eb4_k27_cov
			k4cov = cgd$eb4_k4_cov
			bivcov = cgd$eb4_biv_cov
			dst = cgd[,c("eb4_k27_cov", "eb4_biv_cov", "eb4_k4_cov")]




layout(matrix(1:1, ncol=1))
		f_out = lk27 < 6 & lk4 < 6
		dst_n = dst/rowSums(dst)
		for(li in 2:2) {
			par(mar = c(2,2,1,1))
			f = l_class==li & !f_out

			dst_i = dst_n[f,]

			n = nrow(dst_i)
			ord = order((1e-3+dst_i[,1]+dst_i[,2])/(1e-3+dst_i[,3]+dst_i[,2]))
			dst_i = dst_i[ord,]
			plot(-0,-0, xlim=c(0,1), ylim=c(1,nrow(dst_i)), xaxt='n',yaxt='n')#
			polygon(x = c(dst_i[,3], rep(0,t=n)), y=c(1:n, n:1), col="darkred", border=NA)
			polygon(x = c(dst_i[,3], rev(dst_i[,2]+ dst_i[,3])), y=c(1:n, n:1), col="black", border=NA)
			polygon(x = c(rep(1,t=n), rev(dst_i[,3]+dst_i[,2])), y=c(1:n, n:1), col="darkblue", NA)
		}
   
    }


plt_regionAtac_k27_k4 = function(mod, chr, st, end, marg5=2e+3, marg3=2e+3, mark=c(),  k4_max_fact=1){
#'jk.epipcg.multieb24.mcEpiblast'
 #   'jk.epipcg.atac.CRJK_0389_atac_wt_to_wt_eb_d3'
tracks = c(
    'jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a',
#'jk.epipcg.pcg.CRJK_0403_k27me3_wt_to_wt_eb_d3_a70ls',
'jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3')
short_nms = c('k27','k4')
	for(i in tracks) {
		gvtrack.create(gsub('jk.epipcg.','',i), i, "sum") 
		gvtrack.iterator(gsub('jk.epipcg.','',i), sshift = -140, eshift= 140)
	}
prof = gextract(gsub('jk.epipcg.','',tracks), intervals=gintervals(chr, st-marg5, end+marg3), iterator=20,colnames=short_nms)
prof[is.na(prof)] = 0 

#ds_rat = 
#k27 cov6487273.4
#k4 coverage5481174

k27_vs = as.matrix(prof[, 'k27'])# * tmax/tmax_k27
        message("into dsamp k27",  " targ ", floor(sum(k27_vs)*0.72))
        post_ds = downsample_matrix(k27_vs,  target_n = floor(sum(k27_vs)*0.72), seed=19)
        k27_vs = as.numeric(post_ds)
    
k4_vs = as.matrix(prof[, 'k4'])# * tmax/tmax_k4
       message("into dsamp k4",  " targ ", floor(sum(k4_vs)*1))
        post_ds = downsample_matrix(k4_vs,  target_n = floor(sum(k4_vs)*1), seed=19)
        k4_vs = as.numeric(post_ds)    
tmax = max(as.numeric(c(k27_vs,k4_vs)))
#tmax =  max(prof[, 'k27'])
# Set up the plot area with equal-sized panels and no spacing between them
# Important: call this BEFORE any plots are created
par(mfrow=c(2,1), oma=c(4,4,2,4), mar=c(0,0,0,0), mgp=c(2,1,0))



# Plot 2: K27 signal (scaled)
tmax_k27 = max(prof[, 'k27'])# * k4_max_fact


#message("scaled k27 max: ", tmax_k27, " (original max: ", max(prof[, 'k27']), ")")
plot(prof$start, pmin(k27_vs, tmax),
     type="l", lwd=1, col="darkblue",
     ylim=c(0, tmax),
     xaxt='n',  # No x-axis
     yaxt='s',  # Show y-axis
     main="", xlab="", ylab="")
polygon(c(prof$start, rev(prof$start)), c(pmin(k27_vs, tmax), rep(1, length(pmin(k27_vs, tmax)))), 
        col="darkblue", border=NA)
mtext("K27 Signal", side=2, line=2, cex=0.9)
		abline(v=st)
		abline(v=end)
		if(!is.null(mark) & length(mark)>0) {
			for(x in mark) {
				abline(v=x, col="blue", lwd=2)
			}
		}
# Plot 1: ATAC signal
#tmax_atac = max(prof[, 'atac'])
#atac_vs = prof[, 'atac'] * tmax/tmax_atac
#message("scaled atac max: ", max(atac_vs), " (original max: ", max(prof[, 'atac']), ")")
#plot(prof$start, pmin(atac_vs, tmax), 
#     type="l", col='darkred', lwd=3,
#     ylim=c(0, tmax), 
#     #xaxt='n',  # Hide x-axis
#     yaxt='s',  # Show y-axis
#     main="",   # No individual plot title
#     xlab="", ylab="")
#mtext("Genomic Position", side=1, line=2.5, outer=TRUE)
#mtext("ATAC Signal", side=2, line=2, cex=0.9)
#		abline(v=st)
#		abline(v=end)
#		if(!is.null(mark) & length(mark)>0) {
#			for(x in mark) {
#				abline(v=x, col="blue", lwd=2)
#			}
#		}
# Plot 3: K4 signal (scaled)
tmax_k4 = max(as.numeric(prof[, 'k4']))

#message("scaled k4 max: ", max(k4_vs), " (original max: ", max(prof[, 'k4']), ")")
plot(prof$start, pmin(k4_vs, tmax),
     type="l", lwd=1, col="darkred",
     ylim=c(0, tmax),
     yaxt='s',  # Show y-axis
     main="", xlab="", ylab="")
polygon(c(prof$start, rev(prof$start)), c(pmin(k4_vs, tmax), rep(1, length(pmin(k4_vs, tmax)))), 
        col="darkred", border=NA)
mtext("K4 Signal", side=2, line=2, cex=0.9)


# Add the x-axis label to the bottom of the entire figure
mtext("Genomic Position", side=1, line=2.5, outer=TRUE)
		abline(v=st)
		abline(v=end)
		if(!is.null(mark) & length(mark)>0) {
			for(x in mark) {
				abline(v=x, col="blue", lwd=2)
			}
		}
   }
plot_dense_scatter = function (data, x, y, lim = c(8.5)) 
{
    options(repr.plot.width = 10, repr.plot.height = 10)
    temp = data[, c(x, y)]
    colnames(temp) = c("x", "y")
    p_coldens = densCols(x = (temp[, "x"]), y = (temp[, "y"]), 
        colramp = colorRampPalette(c("darkgray", "blue3", "red", 
            "yellow")))
    temp$col = p_coldens
    gg = temp %>% ggplot(aes(x = (x), y = y, col = col)) + geom_point(alpha = 1, 
        size = 0.8) + labs(x = x, y = y) + theme(axis.title.y = element_text(angle = 0, 
        vjust = 0.5), axis.text.x = element_text(angle = -45, 
        hjust = 0)) + theme_bw() + scale_color_identity() + xlim(c(0, 
        lim))
    print(gg)
}
compute_k27_doms = function(mod, th_percent="X0.99")
{
	ndx = mod$epi_tracks
	for(i in 1:nrow(ndx)) {
#		message("running ", ndx$short_name[i])
		gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "sum") 
		gvtrack.iterator(ndx$short_name[i], sshift = -140, eshift= 140)
	}
	
	gvtrack.create("cg_max", "seq.CG_500_mean_new", "max")
	gvtrack.create("gc_max", "seq.GC500_bin20", "max")
    th_mat = mod$k27_track_thresh
	cond_or = paste(paste(ndx$short_name, th_mat[,th_percent], sep=">"), collapse=" | ")
	doms = gscreen(cond_or)
	doms = doms[doms$chrom != "chrM" & doms$chrom != "chrY",]
	doms$l = doms$end - doms$start

	doms1k = doms[doms$l>1000,]

   doms1k_cg = gextract(c("seq.GC500_bin20", "gc_max", "cg_max"), iterator=doms1k, intervals=doms1k, colnames=c("GC","GC_max", "CG_max"))
	doms1k_tss = gintervals.neighbors(doms1k, mod$tss[,c("chrom","start","end","strand","geneSymbol")])
	
	doms1k = cbind(doms1k, doms1k_cg[,c("GC", "GC_max", "CG_max")])
	doms1k = cbind(doms1k, doms1k_tss[,c("dist","geneSymbol")])

	mod$k27_doms = doms1k
	mod$k27_all_doms = doms

	#mat = gextract(ndx$short_name, intervals=doms1k, iterator=doms1k)
	#mat_n = t(t(mat[,4:(ncol(mat)-1)])/apply(mat[,4:(ncol(mat)-1)], 2, median))
	#mat_n2 = mat_n / rowMeans(mat_n)
	#mod$k27_mat = mat
	#mod$k27_mat_n2 = mat_n2

	return(mod)
}
compute_k27_doms_ES = function(mod, th_percent="X0.99")
{
	
    gvtrack.create('ES_chip','jk.epipcg.lit.encode.es.ENCFF595SIA_H3K27me3_es', "sum") 
		gvtrack.iterator('ES_chip', sshift = -140, eshift= 140)

gvtrack.create('ES_cnt_n','jk.epipcg.pcg.CRJK_0211_k27me3_es_50k_norm', "avg") 
		gvtrack.iterator('ES_cnt_n', sshift = -500, eshift= 500)
short_nms = c('ES_cnt_n','ES_chip')
    

	gvtrack.create("cg_max", "seq.CG_500_mean_new", "max")
	gvtrack.create("gc_max", "seq.GC500_bin20", "max")
 
    cond_or = 'ES_cnt_n > 6.05 | ES_ch > 117.2'
	doms = gscreen(cond_or)
	doms = doms[doms$chrom != "chrM" & doms$chrom != "chrY",]
	doms$l = doms$end - doms$start

	doms1k = doms[doms$l>300,]

   doms1k_cg = gextract(c("seq.GC500_bin20", "gc_max", "cg_max"), iterator=doms1k, intervals=doms1k, colnames=c("GC","GC_max", "CG_max"))
	doms1k_tss = gintervals.neighbors(doms1k, mod$tss[,c("chrom","start","end","strand","geneSymbol")])
	
	doms1k = cbind(doms1k, doms1k_cg[,c("GC", "GC_max", "CG_max")])
	doms1k = cbind(doms1k, doms1k_tss[,c("dist","geneSymbol")])

	mod$k27_doms = doms1k
	mod$k27_all_doms = doms

	#mat = gextract(ndx$short_name, intervals=doms1k, iterator=doms1k)
	#mat_n = t(t(mat[,4:(ncol(mat)-1)])/apply(mat[,4:(ncol(mat)-1)], 2, median))
	#mat_n2 = mat_n / rowMeans(mat_n)
	#mod$k27_mat = mat
	#mod$k27_mat_n2 = mat_n2

	return(mod)
}
pcg_report_locmod_iq_meth_pptx = function(mod, fit_type)
{
#prediction - f_out, f_ambig, f_train, f_test
#vs. lk27, k4_cov_1
#scatter and box
    mod$test_chroms = c("chr2", "chr8", "chr12", "chr18")  
	cgd = mod$cgdom_ann
	f_500 = cgd$l>300
	k27_k4_l10  =  log2(mod$cgdom_ann$es_k27_max + mod$cgdom_ann$es_k4_max+10)
    k4l10 = log2( mod$cgdom_ann$es_k4_max+10)

	#k27_cov_1 = cgd$eb4_k27_cov_1
	#f_out = lk27<4.5 & lk4<4.5
	#f_ambig = mod$cgdom_ann$f_ambig
    
	f_train = !cgd$chrom %in% mod$test_chroms
	f_test = cgd$chrom %in% mod$test_chroms
	f_tss = TRUE#cgd$tss_dist == 0
        resp4 = k27_k4_l10
        resp10 = ifelse(resp4 > 7, 1, ifelse(resp4 < 4,
            0, NA))
    print(table(resp10))
    f_ambig = resp4 > 4 & resp4 < 7
    to_fit = resp10
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)
	pred = mod$cgdom_ann$pred5mc
#	to_fit = mod$locmod[[fit_type]]$to_fit
	
	#to_fit = ifelse(cgd$eb4_biv_cov_1==0 & cgd$eb4_k27_cov==0,1,ifelse(cgd$eb4_k4_cov_1==0, 0, NA))
    
	plt1 = function(){
    layout(matrix(1:4,nrow=2))
	for(tss in c(1)) {
		if(tss == 1) { 
			filt = f_tss
		} else {
			filt = !f_tss
		}
		f1 = filt & f_500 & !f_ambig  & f_train
		f2 = filt & f_500 & !f_ambig  & f_test
		#f3 = filt & f_500 & f_ambig 
		#f4 = filt & f_500 & f_out
	for(i in 1:2) {
		if(i == 1) {
			stat = k27_k4_l10
		} else {
			stat = k4l10 
		}
		tr_cor = round(cor(stat[f1], pred[f1],m="s"),3)
		te_cor = round(cor(stat[f2], pred[f2], m="s"),3)
		
		message("comp ranges")
		#ranges = c(min(pred),seq(ceiling(quantile(pred,0.05)), floor(quantile(pred,0.95)), 0.5), max(pred))
        ranges = seq(0, 1, 0.125)
        ranges =c(-1,.2,.4,.6,.8,1.1)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred[f1],ranges)), boxwex=0.15, col="blue",outline=FALSE, 
								main=sprintf("%s, tss = %s, cor = %s,%s", fit_type, tss, tr_cor, te_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col="gold", add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}}
    save_baseR_to_ppt(plt1(),'./figs/bxplt_iq_meth.pptx')

	plt_auc = function(){
	f1 = f_tss & f_500 & !f_ambig  & f_train
	f2 = f_tss & f_500 & !f_ambig  & f_test
	#f3 = !f_tss & f_500 & !f_ambig & !f_out & f_train
	#f4 = !f_tss& f_500 & !f_ambig & !f_out & f_test
	roc1 = comp_auc(pred, to_fit, f1)
	roc2 = comp_auc(pred, to_fit, f2)
	#roc3 = comp_auc(pred, to_fit, f3)
	#roc4 = comp_auc(pred, to_fit, f4)

	plot(roc1$fn, roc1$tp, t="l", lwd=2, col="blue", xlab="fn", ylab="tp",
					 main=sprintf("auc %s,%s, train np = %s, na = %s, nn = %s", 
										round(roc1$auc,3), round(roc2$auc,3),
										sum(f1 & !is.na(to_fit) & to_fit),
										sum(f1 & is.na(to_fit)),
										sum(f1 & !is.na(to_fit) & !to_fit)))
	lines(roc2$fn, roc2$tp, t="l", lwd=2, col="gold")
	#lines(roc3$fn, roc3$tp, t="l", lwd=2, "blue", lty=2)
	#lines(roc4$fn, roc4$tp, t="l", lwd=2, col="gold", lty=2)
	abline(a=0,b=1)
    }
    save_baseR_to_ppt(plt_auc(),'./figs/auc_iq_meth.pptx')

}



pcg_report_locmod_iq_pptx_hg = function(mod,fit_type)
{
#prediction - f_out, f_ambig, f_train, f_test
#vs. lk27, k4_cov_1
#scatter and box
    mod$test_chroms = c("chr2", "chr8", "chr12", "chr18")  
	cgd = mod$cgdom_ann
	f_500 = cgd$l>500
	lk27 = cgd$l10_h1_k27_sum
	lk4 = cgd$l10_h1_k4_sum
	k4_cov_1 = cgd$es_k4_cov_1
	#k27_cov_1 = cgd$eb4_k27_cov_1
	f_out = cgd$l10_h1_k27_sum< 6.5& cgd$l10_h1_k4_sum < 6
	#f_ambig = mod$cgdom_ann$f_ambig
    f_ambig = cgd$l10_h1_k27_sum> 4.5& cgd$l10_h1_k27_sum < 6.5
	f_train = !cgd$chrom %in% mod$test_chroms
	f_test = cgd$chrom %in% mod$test_chroms
	f_tss = cgd$tss_dist == 0
    resp = lk27#log2(10+k27)
    resp4 = lk4#log2(10+k4)
    resp10 = ifelse(resp<=4 & resp4 > 7, 1, ifelse(resp>=7, 0, NA))
    print(table(resp10))
    to_fit = resp10
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)
	pred = mod$cgdom_ann$pred
#	to_fit = mod$locmod[[fit_type]]$to_fit
	
	#to_fit = ifelse(cgd$eb4_biv_cov_1==0 & cgd$eb4_k27_cov==0,1,ifelse(cgd$eb4_k4_cov_1==0, 0, NA))
    
	plt1 = function(){
    layout(matrix(1:4,nrow=2))
	for(tss in c(1,0)) {
		if(tss == 1) { 
			filt = f_tss
		} else {
			filt = !f_tss
		}
		f1 = filt & f_500 & !f_ambig & !f_out & f_train
		f2 = filt & f_500 & !f_ambig & !f_out & f_test
		f3 = filt & f_500 & f_ambig 
		f4 = filt & f_500 & f_out
	for(i in 1:2) {
		if(i == 1) {
			stat = lk27
		} else {
			stat = lk4 - lk27
		}
		tr_cor = round(cor(stat[f1], pred[f1],m="s"),3)
		te_cor = round(cor(stat[f2], pred[f2], m="s"),3)
		
		message("comp ranges")
		#ranges = c(min(pred),seq(ceiling(quantile(pred,0.05)), floor(quantile(pred,0.95)), 0.5), max(pred))
        #ranges = seq(0, 1, 0.125)
        #ranges =c(-1,.2,.4,.6,.8,1.1)
		ranges = c(1.5,4.5,5,5.5,6,6.5,9)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred[f1],ranges)), boxwex=0.15, col="blue",outline=FALSE, 
								main=sprintf("%s, tss = %s, cor train = %s,cor test %s", fit_type, tss, tr_cor, te_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col="gold", add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}}
    save_baseR_to_ppt(plt1(),'./figs/bxplt_iq_hg.pptx')

	plt_auc = function(){
	f1 = f_tss & f_500 & !f_ambig & !f_out & f_train
	f2 = f_tss & f_500 & !f_ambig & !f_out & f_test
	f3 = !f_tss & f_500 & !f_ambig & !f_out & f_train
	f4 = !f_tss& f_500 & !f_ambig & !f_out & f_test
	roc1 = comp_auc(-pred, to_fit, f1)
	roc2 = comp_auc(-pred, to_fit, f2)
	roc3 = comp_auc(-pred, to_fit, f3)
	roc4 = comp_auc(-pred, to_fit, f4)

	plot(roc1$fn, roc1$tp, t="l", lwd=2, col="blue", xlab="fn", ylab="tp",
					 main=sprintf("auc train %s, test %s, train np = %s, na = %s, nn = %s", 
										round(roc1$auc,3), round(roc2$auc,3),
										sum(f1 & !is.na(to_fit) & to_fit),
										sum(f1 & is.na(to_fit)),
										sum(f1 & !is.na(to_fit) & !to_fit)))
	lines(roc2$fn, roc2$tp, t="l", lwd=2, col="gold")
	#lines(roc3$fn, roc3$tp, t="l", lwd=2, "blue", lty=2)
	#lines(roc4$fn, roc4$tp, t="l", lwd=2, col="gold", lty=2)
	abline(a=0,b=1)
    }
    save_baseR_to_ppt(plt_auc(),'./figs/auc_iq_hg.pptx')

}

pcg_report_locmod_iq_pptx = function(mod, fit_type='IQ')
{
#prediction - f_out, f_ambig, f_train, f_test
#vs. lk27, k4_cov_1
#scatter and box
	iq = as.data.frame(fread('./output/iq-mm10-model_preds-diff-preds.tsv'))
    iq = iq %>% arrange(chrom,start)

    
    cgd = mod$cgdom_ann
	f_500 = cgd$l>500
	lk27 = log2(10+cgd$eb4_k27_max)
	lk4 = log2(10+cgd$eb4_k4_max)
	k4_cov_1 = cgd$eb4_k4_cov_1
	k27_cov_1 = cgd$eb4_k27_cov_1
	f_out = lk27<6 & lk4<6
	f_ambig = k4_cov_1 > 0 & k27_cov_1 > 0
	f_train = !cgd$chrom %in% mod$test_chroms
	f_test = cgd$chrom %in% mod$test_chroms
	f_tss = cgd$tss_dist == 0
	
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)
	pred = iq$pred
#	to_fit = mod$locmod[[fit_type]]$to_fit
	
	to_fit = ifelse(cgd$eb4_biv_cov_1==0 & cgd$eb4_k27_cov==0,1,ifelse(cgd$eb4_k4_cov_1==0, 0, NA))

	#png(sprintf("figs/model_%s_preds.png", fit_type), w=1000,h=800)
	plt1 = function(){
    layout(matrix(1:4,nrow=2))
	for(tss in c(1,0)) {
		if(tss == 1) { 
			filt = f_tss
		} else {
			filt = !f_tss
		}
		f1 = filt & f_500 & !f_ambig & !f_out & f_train
		f2 = filt & f_500 & !f_ambig & !f_out & f_test
		f3 = filt & f_500 & f_ambig 
		f4 = filt & f_500 & f_out
	for(i in 1:2) {
		if(i == 1) {
			stat = lk27
		} else {
			stat = k4_cov_1
		}
		tr_cor = round(cor(stat[f1], pred[f1],m="s"),3)
		te_cor = round(cor(stat[f2], pred[f2], m="s"),3)
		
		message("comp ranges")
		#ranges = c(min(pred),seq(ceiling(quantile(pred,0.05)), floor(quantile(pred,0.95)), 0.5), max(pred))
        ranges = seq(0, 1, 0.125)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred[f1],ranges)), boxwex=0.15, col="blue",outline=FALSE, 
								main=sprintf("%s, tss = %s, cor = %s,%s", fit_type, tss, tr_cor, te_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col="gold", add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}}
    save_baseR_to_ppt(plt1(),'./figs/bxplt_iq.pptx')
	print(plt1())
	#dev.off()

	#pdf(sprintf("figs/model_roc_%s.pdf", fit_type), w=6, h=6)
    plt_auc = function(){
	f1 = f_tss & f_500 & !f_ambig & !f_out & f_train
	f2 = f_tss & f_500 & !f_ambig & !f_out & f_test
	f3 = !f_tss & f_500 & !f_ambig & !f_out & f_train
	f4 = !f_tss& f_500 & !f_ambig & !f_out & f_test
	roc1 = comp_auc(pred, to_fit, f1)
	roc2 = comp_auc(pred, to_fit, f2)
	roc3 = comp_auc(pred, to_fit, f3)
	roc4 = comp_auc(pred, to_fit, f4)

	plot(roc1$fn, roc1$tp, t="l", lwd=2, col="blue", xlab="fn", ylab="tp",
					 main=sprintf("auc %s,%s, train np = %s, na = %s, nn = %s", 
										round(roc1$auc,3), round(roc2$auc,3),
										sum(f1 & !is.na(to_fit) & to_fit),
										sum(f1 & is.na(to_fit)),
										sum(f1 & !is.na(to_fit) & !to_fit)))
	lines(roc2$fn, roc2$tp, t="l", lwd=2, col="gold")
	#lines(roc3$fn, roc3$tp, t="l", lwd=2, "blue", lty=2)
	#lines(roc4$fn, roc4$tp, t="l", lwd=2, col="gold", lty=2)
	abline(a=0,b=1)
    }
    save_baseR_to_ppt(plt_auc(),'./figs/auc_iq.pptx')
	print(plt_auc())
	#dev.off()

}

comp_auc = function (pred, out, f = NULL) 
{
    if (is.null(f)) {
        f = rep(T, length(pred))
    }
    pred = pred[!is.na(out)]
    f = f[!is.na(out)]
    out = out[!is.na(out)]
    res = data.frame(pred[f], ifelse(out[f] == 0, 1, ifelse(out[f] == 
        1, 0, NA)))
    res = res[!is.na(res[, 2]), ]
    fn = cumsum(res[order(res[, 1]), 2] == 0)
    tp = cumsum(res[order(res[, 1]), 2] == 1)
    fn = fn/max(fn)
    tp = tp/max(tp)
    dlt_fn = c(fn, 1) - c(0, fn)
    auc = sum(dlt_fn * c(tp, 0))
    return(list(auc = auc, fn = fn, tp = tp))
}

pcg_spat_energ_dist = function (mod, cgdd, mot_db) 
{
    cgdd = cgdd[cgdd$l > 500, ]
    gen_track_vext()
    interv_cgd = cgdd[, c("chrom", "start", "end", "tss_strand")]
    colnames(interv_cgd) = c("chrom", "start", "end", "strand")
    interv_cgd$ID = 1:nrow(interv_cgd)
    gvtrack.create("d_cgdd", interv_cgd, "distance", 800)
    interv_cgd_1k = interv_cgd
    interv_cgd_1k$start = interv_cgd_1k$start - 1000
    interv_cgd_1k$end = interv_cgd_1k$end + 1000
    interv_cgd_3k = interv_cgd
    interv_cgd_3k$start = interv_cgd_3k$start - 3000
    interv_cgd_3k$end = interv_cgd_3k$end + 3000
    bins = giterator.intervals(c("1"), intervals = interv_cgd_1k, 
        iterator = 20)
    cgd_d = gextract(c("d_cgdd"), intervals = bins, iterator = bins)
    bins_assoc = gintervals.neighbors(cgd_d, interv_cgd)
    cgd_d$ID = bins_assoc$ID
    cgd_d$k4_cov = cgdd[cgd_d$ID, "v_k4_cov"]
    cgd_d$k27_cov = cgdd[cgd_d$ID, "v_k27_cov"]
    bins3k = giterator.intervals(c("1"), intervals = interv_cgd_3k, 
        iterator = 20)
    cgd_d3k = gextract(c("d_cgdd"), intervals = bins3k, iterator = bins3k)
    bins_assoc3k = gintervals.neighbors(cgd_d3k, interv_cgd)
    cgd_d3k$ID = bins_assoc3k$ID
    cgd_d3k$k4_cov = cgdd[cgd_d3k$ID, "v_k4_cov"]
    cgd_d3k$k27_cov = cgdd[cgd_d3k$ID, "v_k27_cov"]
    mots = unique(mot_db$motif)
    for (mot in mots) {
        gvtrack.create(sprintf("PWM%s", mot), func = "pwm", params = list(pssm = mot_db[mot_db$motif == 
            mot, ]))
    }
    mots_vt = paste("PWM", mots, sep = "")
    gvtrack.create("GCcount", "seq.G_or_C", "sum")
    gvtrack.iterator("GCcount", sshift = 0, eshift = 15)
    gvtrack.create("CGcount", "seq.CG", "sum")
    gvtrack.iterator("CGcount", sshift = 0, eshift = 15)
    energ = gextract(mots_vt, intervals = bins, iterator = bins)
    gc = gextract(c("atac_ext", "GCcount", "seq.CG_500_mean_new"), 
        intervals = bins, iterator = bins)
    gc3k = gextract(c("atac_ext", "GCcount", "seq.CG_500_mean_new", 
        "CGcount"), intervals = bins3k, iterator = bins3k)
    work = list(energ = energ, gc = gc, gc3k = gc3k, cgd_d = cgd_d, 
        cgd_d3k = cgd_d3k, cgdd = cgdd)
    return(work)
}

pcg_spat_hg_energ_dist = function (mod, cgdd, mot_db) 
{
    hpcg_gen_atac_vt()
    interv_cgd = cgdd[, c("chrom", "start", "end", "tss_strand")]
    colnames(interv_cgd) = c("chrom", "start", "end", "strand")
    interv_cgd$ID = 1:nrow(interv_cgd)
    gvtrack.create("d_cgdd", interv_cgd, "distance", 800)
    interv_cgd_1k = interv_cgd
    interv_cgd_1k$start = interv_cgd_1k$start - 1000
    interv_cgd_1k$end = interv_cgd_1k$end + 1000
    interv_cgd_3k = interv_cgd
    interv_cgd_3k$start = interv_cgd_3k$start - 3000
    interv_cgd_3k$end = interv_cgd_3k$end + 3000
    bins = giterator.intervals(c("1"), intervals = interv_cgd_1k, 
        iterator = 20)
    cgd_d = gextract(c("d_cgdd"), intervals = bins, iterator = bins)
    bins_assoc = gintervals.neighbors(cgd_d, interv_cgd)
    cgd_d$ID = bins_assoc$ID
    cgd_d$k4_cov = cgdd[cgd_d$ID, "es_k4_cov"]
    cgd_d$k27_cov = cgdd[cgd_d$ID, "es_k27_cov"]
    bins3k = giterator.intervals(c("1"), intervals = interv_cgd_3k, 
        iterator = 20)
    cgd_d3k = gextract(c("d_cgdd"), intervals = bins3k, iterator = bins3k)
    bins_assoc3k = gintervals.neighbors(cgd_d3k, interv_cgd)
    cgd_d3k$ID = bins_assoc3k$ID
    cgd_d3k$k4_cov = cgdd[cgd_d3k$ID, "es_k4_cov"]
    cgd_d3k$k27_cov = cgdd[cgd_d3k$ID, "es_k27_cov"]
    mots = unique(mot_db$motif)
    for (mot in mots) {
        gvtrack.create(sprintf("PWM%s", mot), func = "pwm", params = list(pssm = mot_db[mot_db$motif == 
            mot, ]))
    }
    mots_vt = paste("PWM", mots, sep = "")
    gvtrack.create("GCcount", "seq.G_or_C", "sum")
    gvtrack.iterator("GCcount", sshift = 0, eshift = 15)
    gvtrack.create("CGcount", "seq.CG", "sum")
    gvtrack.iterator("CGcount", sshift = 0, eshift = 15)
    energ = gextract(mots_vt, intervals = bins, iterator = bins)
    gc = gextract(c("atac_ext", "GCcount", "seq.CG_500_mean"), 
        intervals = bins, iterator = bins)
    gc3k = gextract(c("atac_ext", "GCcount", "seq.CG_500_mean", 
        "CGcount"), intervals = bins3k, iterator = bins3k)
    work = list(energ = energ, gc = gc, gc3k = gc3k, cgd_d = cgd_d, 
        cgd_d3k = cgd_d3k, cgdd = cgdd)
    return(work)
}


pcg_plot_energ_heats = function (work, base_dir = "figs/cgd_mod",mm10=TRUE,legend=FALSE) 
{


	if ( !mm10 ){base_dir = paste0(base_dir,'_hg19')}
    lk27 = work$cgdd$v_k27_max
    lk4 = work$cgdd$v_k4_max
    fgood = lk4 > 5 | lk27 > 6
    grad_score = (lk27 - lk4)[fgood]
    bad_score = (lk27 - lk4)[!fgood]
    shades = colorRampPalette(c("white", "gray", "red", "black"))(1000)
    for (i in 4:(ncol(work$energ) - 1)) {
        tf_name = colnames(work$energ)[i]
        tf_name = sub("PWM", "", tf_name)
        message("plot ", tf_name)
        tst = cbind(work$cgd_d, tf = work$energ[, i])
        tst$d_cgdd_b = round(tst$d_cgdd/10)
        emat = tst %>% reshape2::dcast(d_cgdd_b ~ ID, value.var = "tf", 
            fun.aggregate = max)
        edf = as.data.frame(t(emat[, -1]))
        f_na = edf == -Inf
        e_T_low = quantile(edf[!f_na], 0.1)
        e_T_high = quantile(edf[!f_na], 0.995)
        edf[f_na] = e_T_low
        edf_good = edf[fgood, ]
        emat_ord = edf_good[order(grad_score, decreasing = T), 
            ]
        emat_ord = pmax(emat_ord, e_T_low)
        emat_150 = t(apply(emat_ord > e_T_high, 1, zoo::rollmean, 
            15, na.rm = T, f = "e"))
        emat_150_150 = apply(emat_150, 2, zoo::rollmean, 150, na.rm = T, 
            f = "e")
        breaks = c(-1e-04, seq(0, max(emat_150_150) * 1.0001, 
            l = 1000))
        ymax = round(1000 * max(emat_150_150))
        pheatmap::pheatmap(emat_150_150, cluster_cols = F, cluster_rows = F, 
            col = shades, filename = sprintf("%s/heat_%s_%s.png", 
                base_dir, tf_name, ymax), legend = legend, w = 5, 
            h = 12, breaks = breaks, show_rownames = F)
        emat_bad = edf[!fgood, ]
        emat_bad = emat_bad[order(bad_score, decreasing = T), 
            ]
        emat_bad = pmax(emat_bad, e_T_low)
        emat_bad_150 = t(apply(emat_bad > e_T_high, 1, zoo::rollmean, 
            15, na.rm = T, f = "e"))
        emat_bad_150_150 = apply(emat_bad_150, 2, zoo::rollmean, 150, 
            na.rm = T, f = "e")
        pheatmap::pheatmap(emat_bad_150_150, cluster_cols = F, 
            cluster_rows = F, col = shades, filename = sprintf("%s/bad_heat_%s.png", 
                base_dir, tf_name), legend = legend, w = 5, h = 3, 
            breaks = breaks, show_rownames = F)
    }
    dst = work$cgdd[fgood, c("v_k27_cov", "v_biv_cov", "v_k4_cov")]
    dst_n = dst/rowSums(dst)
    dst_n = dst_n[order(grad_score, decreasing = T), ]
    dst_n_sm = apply(dst_n, 2, zoo::rollmean, 20, na.rm = T, f = "e")
    n = nrow(dst_n)
    png(sprintf("%s/dom_key.png", base_dir), h = 1000, w = 300)
    plot(-0, -0, xlim = c(0, 1), ylim = c(1, nrow(dst_n)), xaxt = "n", 
        yaxt = "n")
    polygon(x = c(dst_n_sm[, 3], rep(0, t = n)), y = c(1:n, n:1), 
        col = "darkred", border = NA)
    polygon(x = c(dst_n_sm[, 3], rev(dst_n_sm[, 2] + dst_n_sm[, 
        3])), y = c(1:n, n:1), col = "black", border = NA)
    polygon(x = c(rep(1, t = n), rev(dst_n_sm[, 3] + dst_n_sm[, 
        2])), y = c(1:n, n:1), col = "darkblue", NA)
    dev.off()
    png(sprintf("%s/dom_k4k27.png", base_dir), h = 1000, w = 250)
    par(mar = c(0, 0, 0, 0))
    barplot(zoo::rollmean(lk27[fgood][order(grad_score, decreasing = T)], 
        150, f = "e") - 2, border = "darkblue", col = "darkblue", 
        horiz = T, xaxt = "n", yaxt = "n")
    abline(v = 2, lty = 2, lwd = 2)
    abline(v = 4, lty = 2, lwd = 2)
    dev.off()
}
pcg_plot_energ_heats_hg = function (work, base_dir = "figs/cgd_mod_hg19",legend = F) 
    
{
    
    
    lk27 = log2(4 + work$cgdd$h1_k27_sum)
    lk4 = log2(4 + work$cgdd$h1_k4_sum)
    fgood = lk4 > 5 | lk27 > 5
    grad_score = (lk27 - lk4)[fgood]
    bad_score = (lk27 - lk4)[!fgood]
    shades = colorRampPalette(c("white", "gray", "red", "black"))(1000)
    for (i in 4:(ncol(work$energ) - 1)) {
        tf_name = colnames(work$energ)[i]
        tf_name = sub("PWM", "", tf_name)
        message("plot ", tf_name)
        tst = cbind(work$cgd_d, tf = work$energ[, i])
        tst$d_cgdd_b = round(tst$d_cgdd/10)
        emat = tst %>% reshape2::dcast(d_cgdd_b ~ ID, value.var = "tf", 
            fun.aggregate = max)
        edf = as.data.frame(t(emat[, -1]))
        f_na = edf == -Inf
        e_T_low = quantile(edf[!f_na], 0.1)
        e_T_high = quantile(edf[!f_na], 0.995)
        edf[f_na] = e_T_low
        edf_good = edf[fgood, ]
        emat_ord = edf_good[order(grad_score, decreasing = T), 
            ]
        emat_ord = pmax(emat_ord, e_T_low)
        emat_150 = t(apply(emat_ord > e_T_high, 1, zoo::rollmean, 
            15, na.rm = T, f = "e"))
        emat_150_150 = apply(emat_150, 2, zoo::rollmean, 150, na.rm = T, 
            f = "e")
        breaks = c(-1e-04, seq(0, max(emat_150_150) * 1.0001, 
            l = 1000))
        ymax = round(1000 * max(emat_150_150))
        pheatmap::pheatmap(emat_150_150, cluster_cols = F, cluster_rows = F, 
            col = shades, filename = sprintf("%s/heat_%s_%s.png", 
                base_dir, tf_name, ymax), legend = legend, w = 5, 
            h = 12, breaks = breaks, show_rownames = F)
        emat_bad = edf[!fgood, ]
        emat_bad = emat_bad[order(bad_score, decreasing = T), 
            ]
        emat_bad = pmax(emat_bad, e_T_low)
        emat_bad_150 = t(apply(emat_bad > e_T_high, 1, zoo::rollmean, 
            15, na.rm = T, f = "e"))
        emat_bad_150_150 = apply(emat_bad_150, 2, zoo::rollmean, 150, 
            na.rm = T, f = "e")
        pheatmap::pheatmap(emat_bad_150_150, cluster_cols = F, 
            cluster_rows = F, col = shades, filename = sprintf("%s/bad_heat_%s.png", 
                base_dir, tf_name), legend = legend, w = 5, h = 3, 
            breaks = breaks, show_rownames = F)
    }
    dst = work$cgdd[fgood, c("es_k27_cov", "es_biv_cov", "es_k4_cov")]
    dst_n = dst/rowSums(dst)
    dst_n = dst_n[order(grad_score, decreasing = T), ]
    dst_n_sm = apply(dst_n, 2, zoo::rollmean, 20, na.rm = T, f = "e")
    n = nrow(dst_n)
    png(sprintf("%s/dom_key.png", base_dir), h = 1000, w = 300)
    plot(-0, -0, xlim = c(0, 1), ylim = c(1, nrow(dst_n)), xaxt = "n", 
        yaxt = "n")
    polygon(x = c(dst_n_sm[, 3], rep(0, t = n)), y = c(1:n, n:1), 
        col = "darkred", border = NA)
    polygon(x = c(dst_n_sm[, 3], rev(dst_n_sm[, 2] + dst_n_sm[, 
        3])), y = c(1:n, n:1), col = "black", border = NA)
    polygon(x = c(rep(1, t = n), rev(dst_n_sm[, 3] + dst_n_sm[, 
        2])), y = c(1:n, n:1), col = "darkblue", NA)
    dev.off()
    png(sprintf("%s/dom_k4k27.png", base_dir), h = 1000, w = 250)
    par(mar = c(0, 0, 0, 0))
    barplot(zoo::rollmean(lk27[fgood][order(grad_score, decreasing = T)], 
        150, f = "e") - 2, border = "darkblue", col = "darkblue", 
        horiz = T, xaxt = "n", yaxt = "n")
    abline(v = 2, lty = 2, lwd = 2)
    abline(v = 4, lty = 2, lwd = 2)
    dev.off()
}
####EDF1
rna_emb_vs_eb = function(){
mat = readRDS('data/mat_mc_rna_emb_vs_eb_fig1.rds')

bad_gene_names = c(grep('Eif',rownames(mat),v=T),
grep("^Gm[0-9]",rownames(mat),v=T),
grep("Rps",rownames(mat),v=T),
grep("Rpl",rownames(mat),v=T),
grep("Rik",rownames(mat),v=T),
grep("AK[0-9]",rownames(mat),v=T),                   
grep("mKIA",rownames(mat),v=T),
grep("Hsp",rownames(mat),v=T)                  
                   
                  )

length(bad_gene_names)

f_bad = rownames(mat) %in% bad_gene_names

mat_l = log2(mat[ !f_bad,]+3e-5)

dim(mat_l)
dim(mat)

c = cor(mat_l)

c_f = c[ grepl('eb_',colnames(c)),grepl('emb_',colnames(c))]

top_corr <- apply(c_f, 2, function(x) {
  sort_x <- sort(x, decreasing = TRUE)      # sort correlations high → low
  #sort_x <- sort_x[names(sort_x) != names(x)[1]] # remove self-correlation
  head(sort_x, 10)                          # take top 10
})

# 3. Check result
 

best_emb = names(sort(colMeans(top_corr),decreasing=T)[1:10])

top_corr_best <- apply(c_f[ ,best_emb], 2, function(x) {
  sort_x <- sort(x, decreasing = TRUE)      # sort correlations high → low
  #sort_x <- sort_x[names(sort_x) != names(x)[1]] # remove self-correlation
  head(sort_x, 10)   
    names(sort_x[1:10])
})

# 3. Check result


best_emb_mc = colnames(top_corr_best)

best_eb_mc <- unique(as.vector(as.matrix(top_corr_best)))


#options(repr.plot.width = 15, repr.plot.height = 10)
#pheatmap::pheatmap(top_corr)

#options(repr.plot.width = 10, repr.plot.height = 10)
#pheatmap::pheatmap(c_f)


rna_eb = rowMeans(mat_l[, best_eb_mc])
rna_emb = rowMeans(mat_l[, best_emb_mc])

cor(rna_eb,rna_emb)
cor(rna_eb,rna_emb,m='s')
cor(rna_eb,rna_emb)^2

col = ((rna_emb - rna_eb)>2)

col2 = ((rna_emb - rna_eb)< -2 )

rna_emb[col]

rna_emb[col2]
    

df = data.frame(rna_eb,rna_emb)
options(repr.plot.width = 8, repr.plot.height = 8)
gg = ggplot(df,aes( x = rna_eb, y= rna_emb)) + geom_point(size=.8,alpha =1)+ theme_bw()+
    geom_abline(slope = 1,intercept = 0)+
    geom_abline(slope = 1,intercept = -2)+
    geom_abline(slope = 1,intercept = 2)
print(gg)
}

compare_atac_emb_vs_eb = function(){
tracks = c('jk.epipcg.multieb.epi','mut_multi.WT_E65.epi')

ds_a = TRUE

p1 = readRDS('data/peaks_multieb_epi.rds')
p1 = p1$peaks
p2 = readRDS('data/peaks_mut_multi_WT_E65_epi.rds')
p2 = p2$peaks
dim(p2)
peaks = gintervals.rbind(p1,p2)

dim(peaks)

peaks = (gintervals.canonic(peaks))

dim(peaks)

atac_meta = mod$track_atac_ribo
atac_meta = atac_meta[ atac_meta$tr_atac %in% tracks,]

atac_meta = atac_meta[match(tracks,atac_meta$tr_atac),]
short_names = atac_meta$short_name

ribo = atac_meta$ribo_c

tr_dsamp_r = min(ribo)/ribo
names(tr_dsamp_r) = short_names

for (tr in tracks){
gvtrack.create(atac_meta[ atac_meta$tr_atac==tr,]$short_name, tr, "sum")
gvtrack.iterator(atac_meta[ atac_meta$tr_atac==tr,]$short_name, sshift=-140, eshift=140)
}





peaks$start = floor(peaks$start + (peaks$end - peaks$start)/2) - 150 
peaks$end = peaks$start + 300

peaks = peaks %>% arrange(chrom,start)

gext = gextract(short_names,intervals = peaks,iterator = peaks,colnames = short_names)

gext[is.na(gext)] = 0

cs = colSums(gext[,short_names])

cs

prof_rel = gext

if (ds_a == TRUE){
for (nm in short_names) {
  if(!is.null(tr_dsamp_r[nm])) {
    to_ds = matrix(prof_rel[,nm],ncol=1)
    r_ds = tr_dsamp_r[nm]
    message("into dsamp, ", nm, " targ ", floor(sum(to_ds)*r_ds))
    post_ds = downsample_matrix(to_ds,  target_n = floor(sum(to_ds)*r_ds), seed=19)
    prof_rel[,nm] = post_ds
    }
}
}

gext = prof_rel

#mat  = as.data.frame((t(t(gext[,short_names])/as.numeric(cs)))*40e6)


mat  = as.data.frame(gext[,short_names])
mat = log2(7+mat)

temp = cbind(gext[,c('chrom','start','end')],mat)

cg = mod$cgdom_ann

mat2 = cbind(gext[,c('chrom','start','end','intervalID')],mat)

mat3 = misha.ext::gintervals.neighbors1(mat2,cg)

mat3$type_seed = paste0(mat3$type,mat3$type_50k)

mat3$dist[is.na(mat3$dist)] = 100e6

mat3$type2 = ifelse(abs(mat3$dist)==0, mat3$type,'cg_low')

library('RColorBrewer')

cols = brewer.pal(5, 'Set1')

cols = rev(c('black','#440154FF','#31688EFF','darkgreen','brown'))



data = mat3 %>% filter(type2!='cg_low')

x= short_names[1]
y = short_names[2]
data$geneSymbol2 = paste0(data$tss_gene,'_',data$tss_dist)

options(repr.plot.width = 8, repr.plot.height = 8)
    library(ggrepel)
    temp = data[, c(x, y,'geneSymbol2','type')]
    colnames(temp) = c("x", "y",'geneSymbol','type')
  
    temp$col2 = ifelse(abs(temp$x - temp$y) > 2 & temp$type %in%c('cl1','cl2','cl3','cl4')  ,"darkblue","gray")
  

gg =  temp %>% ggplot(aes(x = (x), y = y)) + geom_point(alpha = 1, 
        size = 2)+
labs(x = x, y = y,subtitle = paste0('cor = ',sprintf("%.2f", cor(temp$x,temp$y)))) +
    #coord_cartesian(ylim = c(-4,4))+
theme(plot.subtitle = element_text(hjust = 0.5) )+
    geom_text_repel(data = filter(temp, temp$col2!="gray"),size=5,max.overlaps = 20,
                    aes(label=geneSymbol))+geom_abline(slope = 1,intercept = 0)#+
    #scale_color_manual(values = rev(c('gray','#440154FF','#31688EFF','darkgreen','brown')))
print(gg)

data = mat3 %>% filter(type2!='cg_low')

x= short_names[1]
y = short_names[2]
data$geneSymbol2 = paste0(data$tss_gene,'_',data$tss_dist)

options(repr.plot.width =8, repr.plot.height = 8)
    library(ggrepel)
    temp = data[, c(x, y,'geneSymbol2','type')]
    colnames(temp) = c("x", "y",'geneSymbol','type')
  
    temp$col2 = ifelse(abs(temp$x - temp$y) > 2 & temp$type %in%c('cl1','cl2','cl3','cl4')  ,"darkblue","gray")
  

gg =  temp %>% ggplot(aes(x = (x), y = y)) + geom_point(alpha = 1, 
        size = .8)+
#labs(x = x, y = y,subtitle = paste0('cor = ',sprintf("%.2f", cor(temp$x,temp$y)))) +
    #coord_cartesian(ylim = c(-4,4))+
theme(plot.subtitle = element_text(hjust = 0.5) )+coord_cartesian(xlim = c(5,12.5),ylim = c(5,12.5))+
   # geom_text_repel(data = filter(temp, temp$col2!="gray"),size=5,max.overlaps = 20,
                   # aes(label=geneSymbol))+geom_abline(slope = 1,intercept = 0)+
geom_abline(slope = 1,intercept = 2)+
geom_abline(slope = 1,intercept = 0)+theme_bw()+
geom_abline(slope = 1,intercept = -2)
    #scale_color_manual(values = rev(c('gray','#440154FF','#31688EFF','darkgreen','brown')))
print(gg)
print('r squared')
print(cor(temp$x,temp$y)^2)
}



cnt_vs_chip = function(){
gvtrack.create("ES_chip", "jk.epipcg.lit.encode.es.ENCFF595SIA_H3K27me3_es", 
        "sum")
    gvtrack.iterator("ES_chip", sshift = -140, eshift = 140)
    gvtrack.create("ES_cnt_n", "jk.epipcg.pcg.CRJK_0211_k27me3_es_50k_norm", 
        "avg")
    gvtrack.iterator("ES_cnt_n", sshift = -500, eshift = 500)
    cnt_doms_a = gscreen("ES_cnt_n > 6.059")#.99 percentile
    cnt_doms_a$l = cnt_doms_a$end - cnt_doms_a$start
    cnt_doms = cnt_doms_a[cnt_doms_a$l > 1000, ]
    chip_doms_a = gscreen("ES_chip > 117.201529")#.99 percentile
    chip_doms_a$l = chip_doms_a$end - chip_doms_a$start
    chip_doms = chip_doms_a[chip_doms_a$l > 1000, ]
    cnt_doms = cnt_doms %>% arrange(chrom, start)
    chip_doms = chip_doms %>% arrange(chrom, start)
doms = rbind(cnt_doms,chip_doms) %>% gintervals.canonic()
    prof = gextract(c(c("ES_chip", "ES_cnt_n"), "seq.CG_500_mean_new", 
        "seq.GC500_bin20"), intervals = doms, iterator = 20)

    prof[is.na(prof)] = 0
    
    mat = doms
    mat$ES_cnt_n = tapply(prof$ES_cnt_n,prof$intervalID,mean)
mat$ES_chip = tapply(prof$ES_chip,prof$intervalID,mean)



    mat2 = mat
    mat2$ES_cnt = mat2$ES_cnt_n
    mat2$ES_ch = log2(7 + mat2$ES_ch)
    numeric_vector = mat2$ES_ch
    max_value <- max(numeric_vector)
    ranges <- list(paste(8, "to", max_value), "7 to 8", "6 to 7", 
        "5 to 6", paste("< 5"))
 
binned_data <- cut(numeric_vector, breaks = c(-Inf,7,8,8.5, max_value), right = TRUE, include.lowest = TRUE)
    mat2$bin = binned_data
    mat2$fill = "cnt"
    ggcnt = mat2 %>% ggplot(aes(x = bin, y = ES_cnt, fill = fill)) + 
        geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(6, 
        7)) + theme_bw() + scale_fill_manual(values = "darkviolet")
    print(ggcnt)
table(mat2$bin)

    numeric_vector = mat2$ES_cnt
    max_value <- max(numeric_vector)
    ranges <- list(paste(8, "to", max_value), "7 to 8", "6 to 7", 
        "5 to 6", paste("< 5"))
    binned_data <- cut(numeric_vector, breaks = c(-Inf,6.5,6.8,7, max_value), right = TRUE, include.lowest = TRUE)
    mat2$bin = binned_data
    mat2$fill = "chip"
    ggchip = mat2 %>% ggplot(aes(x = bin, y = ES_ch, fill = fill)) + 
        geom_boxplot(outlier.shape = NA) + coord_cartesian(ylim = c(6, 
        9)) + theme_bw() + scale_fill_manual(values = "darkgreen")
    print(ggchip)
table(mat2$bin)
return(list(ggchip=ggchip,ggcnt=ggcnt))
}


plt_region_cnt_vs_chip = function(mod, chr, st, end, marg5=2e+3, marg3=2e+3, mark=c(),  k4_max_fact=1){
gvtrack.create('ES_chip','jk.epipcg.lit.encode.es.ENCFF595SIA_H3K27me3_es', "sum") 
		gvtrack.iterator('ES_chip', sshift = -140, eshift= 140)

gvtrack.create('ES_cnt_n','jk.epipcg.pcg.CRJK_0211_k27me3_es_50k_norm', "avg") 
		gvtrack.iterator('ES_cnt_n', sshift = -500, eshift= 500)
short_nms = c('k27','k4')
prof = gextract(c('ES_cnt_n','ES_chip'), intervals=gintervals(chr, st-marg5, end+marg3), iterator=20,colnames=short_nms)
prof[is.na(prof)] = 0 

#ds_rat = 
#k27 cov6487273.4
#k4 coverage5481174

k27_vs = as.matrix(prof[, 'k27'])# * tmax/tmax_k27
        #message("into dsamp k27",  " targ ", floor(sum(k27_vs)*1))
       # post_ds = downsample_matrix(k27_vs,  target_n = floor(sum(k27_vs)*1), seed=19)
       # k27_vs = as.numeric(post_ds)
    
k4_vs = as.matrix(prof[, 'k4'])# * tmax/tmax_k4
       #message("into dsamp k4",  " targ ", floor(sum(k4_vs)*0.1))
       # post_ds = downsample_matrix(k4_vs,  target_n = floor(sum(k4_vs)*0.1), seed=19)
        #k4_vs = as.numeric(post_ds)

k27_vs = zoo::rollmean(k27_vs,k=50,f='e')
k4_vs = zoo::rollmean(k4_vs,k=50,f='e')	
k27_vs[is.na(k27_vs)] = 0
k4_vs[is.na(k4_vs)] = 0	
#k27_vs = (2^(k27_vs))-6
k4_vs = log2(6+k4_vs)
tmax = max(as.numeric(c(k27_vs,k4_vs)))
#tmax =  max(prof[, 'k27'])
# Set up the plot area with equal-sized panels and no spacing between them
# Important: call this BEFORE any plots are created
par(mfrow=c(2,1), oma=c(4,4,2,4), mar=c(0,0,0,0), mgp=c(2,1,0))



# Plot 2: K27 signal (scaled)
tmax_k27 = max(prof[, 'k27'])# * k4_max_fact


#message("scaled k27 max: ", tmax_k27, " (original max: ", max(prof[, 'k27']), ")")
plot(prof$start, pmin(k27_vs, tmax),
     type="l", lwd=.3, col="darkviolet",
     ylim=c(5, 8.1),
     xaxt='n',  # No x-axis
     yaxt='s',  # Show y-axis
     main="", xlab="", ylab="")
polygon(c(prof$start, rev(prof$start)), c(pmin(k27_vs, tmax), rep(1, length(pmin(k27_vs, tmax)))), 
        col="darkviolet", border=NA)
mtext("CnT Signal", side=2, line=2, cex=0.9)
		abline(v=st)
		abline(v=end)
		if(!is.null(mark) & length(mark)>0) {
			for(x in mark) {
				abline(v=x, col="blue", lwd=2)
			}
		}
# Plot 1: ATAC signal
#tmax_atac = max(prof[, 'atac'])
#atac_vs = prof[, 'atac'] * tmax/tmax_atac
#message("scaled atac max: ", max(atac_vs), " (original max: ", max(prof[, 'atac']), ")")
#plot(prof$start, pmin(atac_vs, tmax), 
#     type="l", col='darkred', lwd=3,
#     ylim=c(0, tmax), 
#     #xaxt='n',  # Hide x-axis
#     yaxt='s',  # Show y-axis
#     main="",   # No individual plot title
#     xlab="", ylab="")
#mtext("Genomic Position", side=1, line=2.5, outer=TRUE)
#mtext("ATAC Signal", side=2, line=2, cex=0.9)
#		abline(v=st)
#		abline(v=end)
#		if(!is.null(mark) & length(mark)>0) {
#			for(x in mark) {
#				abline(v=x, col="blue", lwd=2)
#			}
#		}
# Plot 3: K4 signal (scaled)
tmax_k4 = max(as.numeric(prof[, 'k4']))

#message("scaled k4 max: ", max(k4_vs), " (original max: ", max(prof[, 'k4']), ")")
plot(prof$start, pmin(k4_vs, tmax),
     type="l", lwd=.3, col="darkgreen",
     ylim=c(2.5, 10),
     yaxt='s',  # Show y-axis
     main="", xlab="", ylab="")
mtext("ChIP Signal", side=2, line=2, cex=0.9)
polygon(c(prof$start, rev(prof$start)), c(pmin(k4_vs, tmax), rep(1, length(pmin(k4_vs, tmax)))), 
        col="darkgreen", border=NA)

# Add the x-axis label to the bottom of the entire figure
mtext("Genomic Position", side=1, line=2.5, outer=TRUE)
		abline(v=st)
		abline(v=end)
		if(!is.null(mark) & length(mark)>0) {
			for(x in mark) {
				abline(v=x, col="blue", lwd=2)
			}
		}
   }
pcg_gen_cnt_comparisons = function(mod)
{
	set.seed(42)
	#tr_es = mod$epi_tracks[grepl("ES_", mod$epi_tracks$short_name),]
	samp_d = mod$k27_doms[sample(1:nrow(mod$k27_doms), 2000),]
	samp_d$start = samp_d$start - 3000
	samp_d$end = samp_d$end + 3000
    gvtrack.create('ES_ch','jk.epipcg.lit.encode.es.ENCFF595SIA_H3K27me3_es', "sum") 
		gvtrack.iterator('ES_ch', sshift = -140, eshift= 140)

gvtrack.create('ES_cnt_n','jk.epipcg.pcg.CRJK_0211_k27me3_es_50k_norm', "avg") 
		#gvtrack.iterator('ES_cnt_n', sshift = -500, eshift= 500)
short_nms = c('ES_cnt_n','ES_ch')
	#horiz = 140
	#for(i in 1:nrow(tr_es)) {
	#	gvtrack.create(tr_es$short_name[i], tr_es$track_k27[i], "sum") 
	#	gvtrack.iterator(tr_es$short_name[i], sshift = -horiz, eshift= horiz)
	#}
	#cnt_doms_a = gscreen("EB4_cnt > 18.615989")
	cnt_doms_a = gscreen("ES_cnt_n > 6.17")#.99
	cnt_doms_a$l = cnt_doms_a$end - cnt_doms_a$start
	cnt_doms = cnt_doms_a[cnt_doms_a$l>300,]
	cnt_bord5 = cnt_doms
	cnt_bord5$end = cnt_bord5$start+1
	cnt_bord5$strand = 1
	cnt_bord3 = cnt_doms
	cnt_bord3$start = cnt_bord3$end-1
	cnt_bord3$strand = -1
	cnt_bord = rbind(cnt_bord3, cnt_bord5)
	gvtrack.create("cnt_dom_dist", cnt_bord, "distance")

	chip_doms_a = gscreen("ES_ch > 117.201529")#.99
	chip_doms_a$l = chip_doms_a$end - chip_doms_a$start
	chip_doms = chip_doms_a[chip_doms_a$l>300,]
	chip_bord5 = chip_doms
	chip_bord5$end = chip_bord5$start+1
	chip_bord5$strand = 1
	chip_bord3 = chip_doms
	chip_bord3$start = chip_bord3$end-1
	chip_bord3$strand = -1
	chip_bord = rbind(chip_bord3, chip_bord5)
	gvtrack.create("chip_dom_dist", chip_bord, "distance")


	prof = gextract(c('2^(ES_cnt_n)-14','ES_chip',"cnt_dom_dist", "chip_dom_dist", "seq.CG_500_mean_new", "seq.GC500_bin20"),
                    intervals=samp_d, iterator=20,colnames = c(short_nms,"cnt_dom_dist", "chip_dom_dist", "seq.CG_500_mean_new", "seq.GC500_bin20"))
	prof = prof[rowSums(is.na(prof))==0,]
	prf_v  = prof[,short_nms]
	lprf_v = log2(0.05+t(t(prf_v)/colMeans(prf_v)))

	cnt_dst_bin = 20*pmin(pmax(floor(prof$cnt_dom_dist/20),-250),250)
	chip_dst_bin = 20*pmin(pmax(floor(prof$chip_dom_dist/20),-250),250)
	
	cnt_cnt_trend = tapply(lprf_v[,"ES_cnt_n"], cnt_dst_bin, mean)
	chip_cnt_trend = tapply(lprf_v[,"ES_ch"], cnt_dst_bin, mean)

	cnt_chip_trend = tapply(lprf_v[,"ES_cnt_n"], chip_dst_bin, mean)
	chip_chip_trend = tapply(lprf_v[,"ES_ch"], chip_dst_bin, mean)

	#pdf("figs/cnt_chip_borders.pdf", w=12,h=6)
	#layout(matrix(1:2,nrow=1))
    plt1 = function(){
	plot(names(cnt_cnt_trend), cnt_cnt_trend, xlim=c(-4000,2000), type="l", col="darkviolet",
         lwd=3, main="at CnT domain borders", ylim=c(-2.5,2), xlab="dist from border")
	lines(names(chip_cnt_trend), chip_cnt_trend, lwd=3, col="darkgreen")
	#grid()
	abline(v=0)}
    save_baseR_to_ppt(plot_func = plt1(),link_ppt = 'figs/cnt_chip_borders_atCnTbord.pptx' )
	plt2 = function(){
    plot(names(cnt_chip_trend), cnt_chip_trend, xlim=c(-4000,2000), type="l",  col="darkviolet",
         lwd=3, main="at ChIP domain borders", ylim=c(-2.5,2), xlab="dist from border")
	lines(names(chip_chip_trend), chip_chip_trend, lwd=3, col="darkgreen")
	#grid()
	abline(v=0)}
    save_baseR_to_ppt(plot_func = plt2(),link_ppt = 'figs/cnt_chip_borders_atChIPbord.pptx' )
	#dev.off()
	print(plt1())
	print(plt2())
}


####Fig2
generate_cooperativity = function(iq,t_k4=7,t_k27=8,t_k4_on_k27 = 6.4){
    
k4 = iq[ iq$pred_seed_k4 > t_k4 & iq$pred_seed_k27 < t_k4_on_k27,]
k27 = iq[ iq$pred_seed_k27 > t_k27 ,]

iq$d_k4 =  misha.ext::gintervals.neighbors1(iq,k4,mindist = 1)$dist

iq$d_k27 =  misha.ext::gintervals.neighbors1(iq,k27,mindist = 1)$dist



#max pred split to near k4 neark27
#min pred to near k4 neark27

iq$bin = ifelse(iq$pred_seed_k27>t_k27  & iq$d_k4<10000 & iq$d_k27>10000,'pred_k27_near_k4','other')
iq$bin = ifelse(iq$pred_seed_k27>t_k27  & iq$d_k4>10000 & iq$d_k27<10000,'pred_k27_near_k27',iq$bin)

iq$bin = ifelse(iq$pred_seed_k4 > t_k4 & iq$pred_seed_k27 < t_k4_on_k27  & iq$d_k4<10000 & iq$d_k27>10000,'pred_k4_near_k4',iq$bin)
iq$bin = ifelse(iq$pred_seed_k4 > t_k4 & iq$pred_seed_k27 < t_k4_on_k27  & iq$d_k4>10000 & iq$d_k27<10000,'pred_k4_near_k27',iq$bin)

print(table(iq$bin))

gg1 = iq %>%filter(bin%in%c('pred_k27_near_k27',  'pred_k27_near_k4'))%>% ggplot(aes(bin,k27_resp))+ 
geom_boxplot(outlier.shape = NA)+theme_bw()+
coord_cartesian(ylim = c(4,9))
print(gg1)

gg2 = iq %>%filter(bin%in%c('pred_k4_near_k27',  'pred_k4_near_k4'))%>% ggplot(aes(bin,k27_resp)) + 
geom_boxplot(outlier.shape = NA)+theme_bw()+
coord_cartesian(ylim = c(4,9))
print(gg2)  
save_gg_to_ppt(gg1,"./figs/cooperativity.pptx" )      
save_gg_to_ppt(gg2,"./figs/cooperativity2.pptx" )    
}

generate_decay = function(){
cg_trace = mod$gw$cg_trace

f_txg =  cg_trace$d_txg < cg_trace$d_pcg & cg_trace$d_txg < cg_trace$d_mix

f_mix =  cg_trace$d_mix < cg_trace$d_pcg & cg_trace$d_mix < cg_trace$d_txg
f_pcg = cg_trace$d_pcg < cg_trace$d_txg & cg_trace$d_pcg < cg_trace$d_mix

gg_pcg = signal_decay(f=f_pcg,mod=mod,color='darkblue',type = 'd_pcg')
print(gg_pcg)

gg_txg = signal_decay(f=f_txg,mod=mod,color='red',type='d_txg')
print(gg_txg)

gg_mix = signal_decay(f=f_mix,mod=mod,color='gold',type='d_mix')
print(gg_mix)

save_gg_to_ppt(gg_mix,'./figs/decay1.pptx')

save_gg_to_ppt(gg_txg,'./figs/decay2.pptx')

save_gg_to_ppt(gg_pcg,'./figs/decay3.pptx')
    }

signal_decay = function(f,mod,color,type){

cg_trace = mod$gw$cg_trace
trace_pcg = cg_trace[f,]

trace_pcg$bin = cut(trace_pcg$CG,c(-1,.005,0.01,.015,0.035))

trace_pcg = trace_pcg[! is.na(trace_pcg$bin),]

trace_pcg$d_bin = floor(log2(100+trace_pcg[,type]))

trace_gg = trace_pcg %>% group_by(bin,d_bin) %>% summarise(mean(k27_1k))

colnames(trace_gg) = c('cg_bin','d_bin','k271k')

shades=colorRampPalette(c("gray", color))

gg = trace_gg%>% ggplot(aes(x= d_bin,y = log2(k271k),col = cg_bin))+ geom_line()+
theme_bw()+xlim(c(8,19))+ylim(c(-2,6.5))+scale_color_manual(values = shades(4))+
geom_hline(yintercept = -1)+geom_hline(yintercept = 1)+geom_hline(yintercept = 3)+
geom_hline(yintercept = 5)+
theme(panel.grid = element_blank())+
geom_vline(xintercept = log2(500+100))+
geom_vline(xintercept = log2(1000+100))+
geom_vline(xintercept = log2(2000+100))+
geom_vline(xintercept = log2(4000+100))+
geom_vline(xintercept = log2(8000+100))+
geom_vline(xintercept = log2(16000+100))+
geom_vline(xintercept = log2(32000+100))+
geom_vline(xintercept = log2(64000+100))+
geom_vline(xintercept = log2(128000+100))+
geom_vline(xintercept = log2(256000+100))+
geom_vline(xintercept = log2(512000+100))+
geom_abline(slope = -0.85,intercept = 14)
return(gg)
}
generate_decay_new = function(){
cg_trace = mod$gw$cg_trace

cg_trace <- gextract.left_join("mapab.umap_k100", intervals = cg_trace, iterator = cg_trace)


cg_trace$dist_black =  cg_trace %>%
        select(chrom, start, end) %>%
        gintervals.neighbors("mapab.Encode_blacklist_v2") %>% select(dist) %>% pull()

f_norp <- cg_trace$d_ltr != 0 &
              cg_trace$d_line != 0 &            
              cg_trace$start > 3e6 &
              !is.na(cg_trace$mapab.umap_k100) &
              cg_trace$mapab.umap_k100 > 0.9 &
               cg_trace$dist_black !=0 

cg_trace = cg_trace[f_norp,]
cg_trace = cg_trace %>% arrange(chrom,start)
f_txg =  cg_trace$d_txg < cg_trace$d_pcg & cg_trace$d_txg < cg_trace$d_mix

f_mix =  cg_trace$d_mix < cg_trace$d_pcg & cg_trace$d_mix < cg_trace$d_txg
f_pcg = cg_trace$d_pcg < cg_trace$d_txg & cg_trace$d_pcg < cg_trace$d_mix

gg_pcg = signal_decay_new(f=f_pcg,mod=mod,color='darkblue',type = 'd_pcg')
print(gg_pcg)

gg_txg = signal_decay_new(f=f_txg,mod=mod,color='red',type='d_txg')
print(gg_txg)

gg_mix = signal_decay_new(f=f_mix,mod=mod,color='gold',type='d_mix')
print(gg_mix)

save_gg_to_pptx(gg_mix,'./figs/decay1.pptx')

save_gg_to_pptx(gg_txg,'./figs/decay2.pptx')

save_gg_to_pptx(gg_pcg,'./figs/decay3.pptx')
    }

signal_decay_new = function(f,mod,color,type){

cg_trace = mod$gw$cg_trace
#gvtrack.create("k27_eb4_sum", "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a", "sum")  
gvtrack.create("k27_eb4_sum", "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm", "avg")  
f_norp = cg_trace$d_ltr !=0 & cg_trace$d_line!=0 &  cg_trace$d_sine != 0

cg_trace = cg_trace %>% arrange(chrom,start)

if ( "k27_eb4_sum" %in% colnames(cg_trace)) {
  message('runing from precomputed cg_trace')
} else {
  gext = gextract("k27_eb4_sum", intervals = cg_trace, iterator = cg_trace)
    gext[is.na(gext)] = 0
    cg_trace$k27_eb4_sum = (2^gext$k27_eb4_sum) - 6
    message('saving cg_trace')
    saveRDS(cg_trace,'data/cg_trace_mm10.rds')
}
#gext = gextract("k27_eb4_sum",intervals =cg_trace,iterator = cg_trace)
#gext[is.na(gext)] = 0
#cg_trace$k27_eb4_sum = (2^gext$k27_eb4_sum) - 6
cg_trace = cg_trace[f_norp,]
trace_pcg = cg_trace[f,]

trace_pcg$bin = cut(trace_pcg$CG,c(-1,.005,0.01,.015,0.035))

trace_pcg = trace_pcg[! is.na(trace_pcg$bin),]

trace_pcg$d_bin = floor(log2(100+trace_pcg[,type]))

trace_gg = trace_pcg %>% group_by(bin,d_bin) %>% summarise(mean(k27_eb4_sum))

colnames(trace_gg) = c('cg_bin','d_bin','k27_eb4_sum')

shades=colorRampPalette(c("gray", color))

gg = trace_gg%>% ggplot(aes(x= d_bin,y = log2(k27_eb4_sum),col = cg_bin))+ geom_line()+
theme_bw()+
    xlim(c(8,19))+ylim(c(4,7.2))+
    scale_color_manual(values = shades(4))+
geom_hline(yintercept = 4)+geom_hline(yintercept = 5)+geom_hline(yintercept = 6)+
geom_hline(yintercept = 7)+
theme(panel.grid = element_blank())+
#geom_vline(xintercept = log2(500+100))+
geom_vline(xintercept = log2(1000+100))+
#geom_vline(xintercept = log2(2000+100))+
geom_vline(xintercept = log2(4000+100))+
#geom_vline(xintercept = log2(8000+100))+
geom_vline(xintercept = log2(16000+100))+
#geom_vline(xintercept = log2(32000+100))+
geom_vline(xintercept = log2(64000+100))+
#geom_vline(xintercept = log2(128000+100))+
geom_vline(xintercept = log2(256000+100))+
#geom_vline(xintercept = log2(512000+100))+
geom_abline(slope = -.85,intercept = 15)
#return(trace_gg)
return(gg)
}

pcg_report_locmod_gw_pptx = function(mod, fit_type='IQ',perc_genome=0)
{
#prediction - f_out, f_ambig, f_train, f_test
#vs. lk27, k4_cov_1
#scatter and box
	set.seed(42)
	gw = mod$gw
	cg_trace = gw$cg_trace
cg_trace <- gextract.left_join("mapab.umap_k100", intervals = cg_trace, iterator = cg_trace)
cg_trace$dist_black =  cg_trace %>%
        select(chrom, start, end) %>%
        gintervals.neighbors("mapab.Encode_blacklist_v2") %>% select(dist) %>% pull()
f_norp <- cg_trace$d_ltr != 0 &
              cg_trace$d_line != 0 &            
              cg_trace$start > 3e6 &
              !is.na(cg_trace$mapab.umap_k100) &
              cg_trace$mapab.umap_k100 > 0.9 &
               cg_trace$dist_black !=0 
cg_trace_f = cg_trace[f_norp,]
	p = runif(n=nrow(gw$cg_trace))
	f_sub = f_norp & (p >= perc_genome | cg_trace$CG>0.02) 
#& cg_trace$cg_all!=1

	#ptrain = runif(n=nrow(gw$cg_trace))
	f_test_chr = cg_trace$chrom %in% c('chr14','chr10','chr15','chr4')#c("chr1", "chr4", "chr6","chr10","chr13", "chr16", "chr19")
	f_train_chr =!f_test_chr

	f_test = f_test_chr & f_sub 
	f_test_all = f_test_chr & f_norp
	f_train = f_train_chr & f_sub 
	f_train_all = f_train_chr & f_norp
    lk27 = cg_trace$lk27_1k
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)
	pred = cg_trace$pred
#	to_fit = mod$locmod[[fit_type]]$to_fit
	
	to_fit = ifelse(cg_trace$lk27_1k> 6 ,1,0)

	#png(sprintf("figs/model_%s_preds.png", fit_type), w=1000,h=800)
	plt1 = function(){
    layout(matrix(1:4,nrow=2))
	for(tss in c(1)) {
		#if(tss == 1) { 
		#	filt = f_tss
		#} else {
		#	filt = !f_tss
		#}
		f1 = f_train
		f2 = f_test
		#f3 = filt & f_500 & f_ambig 
		#f4 = filt & f_500 & f_out
	for(i in 1) {
		if(i == 1) {
			stat = lk27
		} else {
			stat = k4_cov_1
		}
		tr_cor = round((cor(stat[f1], pred[f1],m="p")^2),3)
		te_cor = round((cor(stat[f2], pred[f2], m="p")^2),3)
		
		message("comp ranges")
		#ranges = c(min(pred),seq(ceiling(quantile(pred,0.05)), floor(quantile(pred,0.95)), 0.5), max(pred))
        ranges = c(4,5.5,6,6.5,7,7.5,8,8.5,9)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred[f1],ranges)), boxwex=0.15, col="blue",outline=FALSE, 
								main=sprintf("rsqr train = %s , rsqr test = %s", tr_cor, te_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col="gold", add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}}
	plt1()
    save_baseR_to_ppt(plt1(),'./figs/bxplt_gw.pptx')
	#dev.off()
	#don't plot auc
	#pdf(sprintf("figs/model_roc_%s.pdf", fit_type), w=6, h=6)
    plt_auc = function(){
	f1 = f_train
	f2 = f_test
	#f3 = !f_tss & f_500 & !f_ambig & !f_out & f_train
	#f4 = !f_tss& f_500 & !f_ambig & !f_out & f_test
	roc1 = comp_auc(pred, to_fit, f1)
	roc2 = comp_auc(pred, to_fit, f2)
	#roc3 = comp_auc(pred, to_fit, f3)
	#roc4 = comp_auc(pred, to_fit, f4)

	plot(roc1$fn, roc1$tp, t="l", lwd=2, col="blue", xlab="fn", ylab="tp",
					 main=sprintf("auc %s,%s, train np = %s, na = %s, nn = %s", 
										round(roc1$auc,3), round(roc2$auc,3),
										sum(f1 & !is.na(to_fit) & to_fit),
										sum(f1 & is.na(to_fit)),
										sum(f1 & !is.na(to_fit) & !to_fit)))
	lines(roc2$fn, roc2$tp, t="l", lwd=2, col="gold")
	#lines(roc3$fn, roc3$tp, t="l", lwd=2, "blue", lty=2)
	#lines(roc4$fn, roc4$tp, t="l", lwd=2, col="gold", lty=2)
	abline(a=0,b=1)
    }
    #save_baseR_to_ppt(plt_auc(),'./figs/auc_gw.pptx')
	#dev.off()

}

plt_genome_pred_ppt = function(mod, g = NA, off5 = NA, off3 = NA,mark_reg=NA,
			 chrom = NULL, locus = NULL, win_smooth = 1, show_center=T, 
			plot_pred_doms = F,
			label_tss = F,
			pred_lwd=1, plot_base_pred = F,
			 fn=NULL, fn_w = 15, fn_h=10, more_tracks=c())
{
	#if(!is.null(fn)) {
	#	if(grepl("pdf", fn)) {
	#		pdf(fn, w=fn_w, h=fn_h)
	#	} else {
	#		png(fn, w=fn_w, h=fn_h)
	#	}
	#}
	cg_trace = mod$gw$cg_trace
	if(is.null(chrom)) {
		f = mod$tss$geneSymbol==g
		if(sum(f) == 0) {
			return
		}
		hits = mod$tss[f,]
		locus = mean(mod$tss$start[f])
		chrom = hits[1,"chrom"]
	}
	f = as.character(cg_trace$chrom) == chrom & (cg_trace$start > locus + off5 & cg_trace$start < locus + off3)

	tss = mod$epi_tss[mod$epi_tss$chrom == chrom & (mod$epi_tss$start > locus + off5 & mod$epi_tss$start < locus + off3),]

	add_n = length(more_tracks)
	layout(matrix(1:(2+add_n), ncol=1),h=c(rep(1, 3+add_n),1.3))
	par(mar = c(0, 3,2, 3))
	maxy = max(max(cg_trace$lk27_1k[f]),7)
	pred_sm = cg_trace$pred[f]
	pred_segs = which(pred_sm[-1]>4 & pred_sm[-length(pred_sm)] < 4)
	pred_segs = c(pred_segs, which(pred_sm[-1]<4 & pred_sm[-length(pred_sm)] > 4))
	pred_segs_xs = cg_trace$start[f][pred_segs]
	obs = cg_trace$lk27_1k[f]
	r2 = cor(pmax(obs,2), pmax(pred_sm,2))**2
	obs_given_pred = round(sum(obs > 6.5 & pred_sm > 5.88)/(1+sum(pred_sm > 5.88)),3)
	pred_given_obs = round(sum(obs > 6.5 & pred_sm > 5.88)/(1+sum(obs > 6.5)),3)
#quant0.2 background
	plot(cg_trace$start[f], cg_trace$lk27_1k[f], 
				t="l", ylim=c(4.74,maxy), lwd=1,col="blue", 
				xaxt='n', 
				main=sprintf("EB4 locus r2 = %s, obs|pred %s, pred|obs %s", 
								round(r2,3), obs_given_pred, pred_given_obs), 
				xlab=sprintf("chr %s, %s", cg_trace$chrom[f][1], g))
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(cg_trace$lk27_1k[f], rep(1, length(cg_trace$lk27_1k[f]))), 
        col="blue", border=NA)
	#if(nrow(tss) > 0) {
	#	points(tss$start, pmin(1.5+tss$epi_rna-log2(1e-5),7.5), pch=19, col="red",cex=1.5)
	#}
	
	abline(h=6.85, lty=2)#.98 quant
	
	if(plot_pred_doms) {
		segments(x0 = pred_segs_xs, x1=pred_segs_xs, 
			y0 = rep(0,length(pred_segs_xs)),
			y1 = rep(maxy,length(pred_segs_xs)), lty=1)
	}
	
	if(show_center) {
		abline(v=locus, col="black", lwd=1)
	}
    if(!is.na(mark_reg[1])) {
        for (mark in 1:length(mark_reg)){
		abline(v=mark_reg[mark], col="black", lwd=1)}
	}
	#par(mar = c(0, 3,1, 3))
    par(mar = c(2,3,1, 3))
	plot(cg_trace$start[f], pred_sm, t="l", col="darkblue", ylim=c(4.83,maxy), xlab=cg_trace$chrom[f][1],
         lwd=pred_lwd, main="prediction")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(pred_sm, rep(1, length(pred_sm))), 
        col="darkblue", border=NA)
		
		abline(h=6.24, lty=2)##.98 quant
		
    if(!is.na(mark_reg[1])) {
        for (mark in 1:length(mark_reg)){
		abline(v=mark_reg[mark], col="black", lwd=1)}
        }
	if(plot_base_pred) {
		pred_sm_b = zoo::rollmean(cg_trace$pred_base[f], win_smooth, f='e')
		lines(cg_trace$start[f], pred_sm_b, t="l", col="cyan", ylim=c(1,maxy), lwd=pred_lwd, xlab=cg_trace$chrom[f][1])
	}
	if(label_tss) {
		if(nrow(tss) > 0) {
			par(srt=-45)
			text(x = tss$start, y = rep(maxy-2,t=nrow(tss)), labels=tss$geneSymbol, cex=1.5)
		}
	}

	if(plot_pred_doms) {
		segments(x0 = pred_segs_xs, x1=pred_segs_xs, 
			y0 = rep(0,length(pred_segs_xs)),
			y1 = rep(maxy,length(pred_segs_xs)), lty=1)
	}
	#abline(h=log2(19), lty=2)
	if(show_center) {
		abline(v=locus, col="black",lwd=2)
	}
	if(length(more_tracks) != 0) {
		profs = gextract(more_tracks, intervals=cg_trace[f,1:3], iterator=cg_trace[f,1:3])
		for(i in 1:length(more_tracks)) {
			sm_prf = zoo::rollmean(profs[,more_tracks[i]], 5, f='e')
			plot(cg_trace$start[f], log2(2+sm_prf), t="l", col="gray", ylim=c(1,maxy), lwd=pred_lwd, main=more_tracks[[i]])
		}
	}

	#par(mar = c(1,3,1, 3))
	#plot(cg_trace$start[f], pmin(log2(10+cg_trace$atac[f]),9), t="l", col="red", ylim=c(3,9), xaxt='n', main="atac")
	#par(mar = c(5,3,1,3))
	#plot(cg_trace$start[f], pmin(cg_trace$CG[f],0.1), t="l", col="darkgreen", ylim=c(0,0.1), 
    #xlab=cg_trace$chrom[f][1], main="CG")
	#if(plot_pred_doms) {
	#	segments(x0 = pred_segs_xs, x1=pred_segs_xs, 
	#		y0 = rep(0,length(pred_segs_xs)),
	#		y1 = rep(maxy,length(pred_segs_xs)), lty=1)
	#}
	if(!is.null(fn)) {
		#dev.off()
	}
}



pcg_plot_atac_heats = function(work, base_dir = "figs/cgd_mod",legend = FALSE)
{
	shades = colorRampPalette(c("white","gray","blue","red","black"))(1000)

	lk27 =work$cgdd$v_k27_max
	lk4 = work$cgdd$v_k4_max

	fgood = lk4 > 5 | lk27 > 6
	grad_score = (lk27 - lk4)[fgood]
	bad_score = (lk27 - lk4)[!fgood]

	tst = cbind(work$cgd_d, tf = work$gc$atac_ext)
	tst$d_cgdd_b = round(tst$d_cgdd/10)

	emat = tst %>% reshape2::dcast(d_cgdd_b ~ ID, value.var = "tf", 
													fun.aggregate = max)
	f = is.na(emat) | emat == -Inf
	emat[f] = 0
   edf = as.data.frame(t(emat[,-1])) 
   e_T_high = 1024
	e_T_low = 0

   edf_good = edf[fgood,]
	emat_ord = edf_good[order(grad_score, decreasing = T),]

	emat_ord = pmax(emat_ord, e_T_low)
	emat_150 = t(apply(emat_ord,1,zoo::rollmean, 3, na.rm=T, f='e'))
   emat_150_150 = apply(emat_150, 2, zoo::rollmean, 150, na.rm=T, f='e')

   breaks = seq(20,1800,l=1001)
   pheatmap::pheatmap(emat_150_150, cluster_cols= F, cluster_rows =F,
               col=shades,
               filename = sprintf("%s/heat_atac.png", base_dir),
               legend=legend, w=5,h=12, breaks = breaks, show_rownames=F)

   emat_bad = edf[!fgood,]
   emat_bad = emat_bad[order(bad_score, decreasing = T),]
   emat_bad = pmax(emat_bad, e_T_low)
   emat_bad_150 = t(apply(emat_bad, 1, zoo::rollmean, 3, na.rm=T, f='e'))
   emat_bad_150_150 = apply(emat_bad_150, 2, zoo::rollmean, 150, na.rm=T, f='e')
   pheatmap::pheatmap(emat_bad_150_150, cluster_cols= F, cluster_rows =F, col=shades, filename = sprintf("%s/bad_heat_atac.png", base_dir),
            legend=legend, w=5,h=3, breaks = breaks, show_rownames=F)
}
plt_region_per_tr = function (mod,take,track_names , chr, st, end,k=1, ds_k27 = TRUE, marg5 = 2000, marg3 = 2000,
                              mark_seeds=F,
    mark = c(), wk4 = F, k4_max_fact = 1) 
{
take = track_names
#ndx = mod$epi_tracks %>% filter(short_name%in% take)
 #   for (i in 1:nrow(ndx)) {
 #       message("running ", ndx$short_name[i])
 #       gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "sum")
   #     gvtrack.iterator(ndx$short_name[i], sshift = -140, eshift = 140)
  #  }
#coverage = readRDS('./data/tracks_cov.rds')

#coverage_v = coverage$cov_vect
#names(coverage_v) = rownames(coverage)
 #   ndx = ndx[match(take, ndx$short_name), ]
 #   tracks_k27 = ndx$track_k27
 #   short_namesk27 = ndx$short_name
    tracks_k27 = take
   short_namesk27 = take
    
    #####
    if (mark_seeds==T){
    seeds = mod$diff_doms$a_doms_e$EB4_cnt

    seeds = seeds[seeds$l> 300,]

    seeds$strand = 1

    seed_nei = gintervals.neighbors1(gintervals(chrom = chr,start = st,end = end,strand=1),seeds,maxneighbors = 4,maxdist = 4e5,mindist = -4e5)


    if (nrow(seed_nei[!is.na(seed_nei$dist),])>0){
        seed_nei = seed_nei[!is.na(seed_nei$dist),]
        mark_seed_st = seed_nei$start1
        mark_seed_en = seed_nei$end1
    
    
    } else {
        mark_seed_st = NULL
        mark_seed_en = NULL
    }
    }
    #####
    
    
    cov = data.frame(row.names = short_namesk27,cov  = rep(NA,length(short_namesk27)))
    
 
 #   hoxs = grep("Hox", mod$tss$geneSymbol, v = T)
 #   hox_int = mod$tss[mod$tss$geneSymbol %in% hoxs, ]
 #   hox_int$start = hox_int$start - 500
 #   hox_int$end = hox_int$end + 500
 #   hox_int = (gintervals.canonic(hox_int))
 #   h = gextract(c(short_namesk27), intervals = hox_int, iterator = hox_int)
 #   hs = colSums(h[, short_namesk27])
 #   hs_df = as.data.frame(hs)
 #   cov$hox = hs_df$hs
    prof = gextract(short_namesk27, intervals = gintervals(chr, 
        st - marg5, end + marg3), iterator = 20)
    prof[is.na(prof)] = 0

	

ylim <-( max(apply(prof[, short_namesk27], 2, function(x) {
  max(zoo::rollmean(x, k = k, fill = 0))
})))
    print(ylim)
    #gvtrack.create("WT_Epi_rep1", "jk.epipcg.meth.meissN20.WT_Epi_rep1", 
    #    "avg")
    #gvtrack.iterator("WT_Epi_rep1", sshift = -100, eshift = 100)

    prof$intervalID = NULL
    layout(matrix(1:(length(short_namesk27)), ncol = 1), h = c(1.4, rep(1, length(short_namesk27) - 
        2), 1.4))
    for (i in 1:length(short_namesk27)) {
        if (i == length(short_namesk27)) {
            par(mar = c(2, 4, 0.8, 4))
        }
        else if (i == 1) {
            par(mar = c(0, 4, 2, 4))
        }
        else {
            par(mar = c(0, 4, 0.8, 4))
        }
        
            c((prof[, 3 + i]), rep(0, length(prof$start)))
            plot(prof$start, (zoo::rollmean(prof[, 3 + i],k = k,fill = 0)), pch = 19, type = "l", 
                lwd = 3, ylim = c(4.9, ylim), xaxt = ifelse(i == 
                  length(short_namesk27), "s", "n"), main = short_namesk27[i], col = "darkblue", 
                ylab = NA)
            polygon(c(prof$start, rev(prof$start)), c((zoo::rollmean(prof[, 
                3 + i],k = k,fill = 0)), rep(0, length(prof$start))), col = "darkblue", 
                border = NA)
       

        abline(v = st)
        abline(v = end)
        if (!is.null(mark) & length(mark) > 0) {
            for (x in mark) {
                abline(v = x, col = "blue", lwd = 2)
            }
        }
        if (mark_seeds==T){
        if (!is.null(mark_seed_st) ) {
            for (x in mark_seed_st) {
                abline(v = x, col = "green", lwd = 2)
            }
        } 
        if (!is.null(mark_seed_en) ) {
            for (x in mark_seed_en) {
                abline(v = x, col = "red", lwd = 2)
            }
        }
            }
    }
}
plot_atac_cgd_bin_distrib = function(mod)
{
	gen_track_vext()
	cgdd_all = mod$cgdom_ann
	cgdd = cgdd_all[cgdd_all$l>500,]
	interv_cgd = cgdd[,c("chrom","start","end","tss_strand")]
	interv_cgd_all = cgdd_all[,c("chrom","start","end","tss_strand")]
	colnames(interv_cgd) = c("chrom", "start", "end", "strand")
	interv_cgd$ID = 1:nrow(interv_cgd)
	gvtrack.create("d_cgdd", interv_cgd, "distance", 800)
	gvtrack.create("d_cgdd_all", interv_cgd, "distance", 800)
	interv_cgd_1k = interv_cgd
	interv_cgd_1k$start = interv_cgd_1k$start - 1000
	interv_cgd_1k$end = interv_cgd_1k$end + 1000

	bins = giterator.intervals(c("1"), intervals = interv_cgd_1k, iterator=20)
	print('gext')
	atac_bins = gextract(c("atac_ext", "seq.CG_500_mean_new", "d_cgdd"), intervals = bins ,iterator=bins, colnames=c("atac", "CG", "d_cgdd"))

	bins_assoc = gintervals.neighbors(atac_bins, interv_cgd)
	atac_bins$ID = bins_assoc$ID

	latac = log2(4+atac_bins$atac)
	#lk27 = log2(4+cgdd$eb4_k27_sum)[atac_bins$ID]
	#lk4 = log2(4+cgdd$eb4_k4_sum)[atac_bins$ID]
	lk27 = cgdd$eb4_k27_mean[atac_bins$ID]
	lk4 = cgdd$l6_eb4_k4_mean[atac_bins$ID]
	fnz = (lk4 != log2(6) | lk27 != log2(6)) & !is.na(latac)
	f = fnz & abs(atac_bins$d_cgdd) < 400

	#pdf("figs/cgdd_atac.pdf", w=8,h=5)

	plt1 = function(){
    plot(density(latac[f & lk27>7 & lk4  > 7], na.rm=T), col="black", 
										xlim=c(3,13), ylim=c(0,0.5), lwd=3, main=NA, xlab=NA, ylab=NA)
#skip plotting the "bad" interacls - it look like background, don't add much
#	lines(density(latac[f & lk27 < 4 & lk4  < 4], na.rm=T), col="darkgray", lwd=3)
	lines(density(latac[f & lk27 > 7 & lk4  < 6.5], na.rm=T), col="blue", lwd=3)
	lines(density(latac[f & lk27 < 6 & lk4  > 7], na.rm=T), col="red", lwd=3)
	lines(density(latac[!f], na.rm=T), col="gray", lty=2,lwd=3)
        }
    print(plt1())
    save_baseR_to_ppt(plot_func = plt1(),link_ppt = "figs/cgdd_atac.pptx" )
	#dev.off()


	#k27_hits4 = gscreen("EB4_cnt > 4.83")	
	##k27_bins = gextract(c("EB4_cnt", "atac_ext", "seq.CG_500_mean_new", "seq.GC200_bin20", "d_cgdd_all"), iterator = 20, intervals = k27_hits4, colnames=c("k27","atac","cg", "gc", "d_cgdd_all"))
	#k27_bins = k27_bins[k27_bins$chrom != "chrM" & k27_bins$chrom != "chrY",]
	#latac = log2(4+k27_bins$atac)
	#
	#latac = pmin(latac, 13)
	#
	#f_far = abs(k27_bins$d_cgdd_all) > 1400
	#f_k27_4 = k27_bins$k27 < 4.83 #q.5 all genome
	#f_k27_5 = k27_bins$k27 < 6.49 & ! f_k27_4#q98 all genome
	#f_k27_6 = k27_bins$k27 > 7.06 #q .99 all genome
	#
	#cgbin = cut(pmin(k27_bins$cg,0.06),seq(0,0.06,0.005))
	#gcbin = cut(pmax(pmin(k27_bins$gc,0.7),0.36),seq(0.35,0.7,0.05))
	#
	##png("figs/gc_k27_atac_far.png", w=800,h=400)
    #plt2 = function(){
	#f = f_k27_4 & f_far; 
	#boxplot(split(latac[f], gcbin[f]), boxwex=0.2, las=2, ylim=c(3,13), pch=19, cex=0.3,outline = FALSE)
	#f = f_k27_5 & f_far; 
	#boxplot(split(latac[f], gcbin[f]), boxwex=0.2, las=2, col="lightblue", at=0.3+1:7, add=T, pch=19, cex=0.3,outline = FALSE)
	#f = f_k27_6 & f_far; 
	#boxplot(split(latac[f], gcbin[f]), boxwex=0.2, las=2, col="darkblue", at=0.6+1:7, add=T, pch=19, cex=0.3,outline = FALSE)
	#abline(h=6)
	#abline(h=8)
	#abline(h=10)
    #    }
	##dev.off()
    #print(plt2())
    #save_baseR_to_ppt(plot_func = plt2(),link_ppt = "figs/gc_k27_atac_far.pptx" )
	#
	##png("figs/gc_k27_atac_near.png", w=800,h=400)
    #plt3 = function(){
	#f = f_k27_4 & !f_far; 
	#boxplot(split(latac[f], gcbin[f]), boxwex=0.2, las=2, ylim=c(3,13), pch=19, cex=0.3,outline = FALSE)
	#f = f_k27_5 & !f_far; 
	#boxplot(split(latac[f], gcbin[f]), boxwex=0.2, las=2, col="lightblue", at=0.3+1:7, add=T, pch=19, cex=0.3,outline = FALSE)
	#f = f_k27_6 & !f_far; 
	#boxplot(split(latac[f], gcbin[f]), boxwex=0.2, las=2, col="darkblue", at=0.6+1:7, add=T, pch=19, cex=0.3,outline = FALSE)
	#abline(h=6)
	#abline(h=8)
	#abline(h=10)
    #    }
    #print(plt3())
    #save_baseR_to_ppt(plot_func = plt3(),link_ppt = "figs/gc_k27_atac_near.pptx" )
	##dev.off()
	#
	##png("figs/cg_k27_atac_far.png", w=800,h=400)
    #plt4 = function(){
	#f = f_k27_4 & f_far; 
	#boxplot(split(latac[f], cgbin[f]), boxwex=0.2, las=2, ylim=c(3,13), pch=19, cex=0.3,outline = FALSE)
	#f = f_k27_5 & f_far; 
	#boxplot(split(latac[f], cgbin[f]), boxwex=0.2, las=2, col="lightblue", at=0.3+1:12, add=T, pch=19, cex=0.3,outline = FALSE)
	#f = f_k27_6 & f_far; 
	#boxplot(split(latac[f], cgbin[f]), boxwex=0.2, las=2, col="darkblue", at=0.6+1:12, add=T, pch=19, cex=0.3,outline = FALSE)
	#abline(h=6)
	#abline(h=8)
	#abline(h=10)
    #    }
    #print(plt4())
    #save_baseR_to_ppt(plot_func = plt4(),link_ppt = "figs/cg_k27_atac_far.pptx" )
	##dev.off()
	#
	##png("figs/cg_k27_atac_near.png", w=800, h=400)
    #plt5 = function(){
	#f = f_k27_4 & !f_far; 
	#boxplot(split(latac[f], cgbin[f]), boxwex=0.2, las=2, ylim=c(3,13), pch=19, cex=0.3,outline = FALSE)
	#f = f_k27_5 & !f_far; 
	#boxplot(split(latac[f], cgbin[f]), boxwex=0.2, las=2, col="lightblue", at=0.3+1:12, add=T, pch=19, cex=0.3,outline = FALSE)
	#f = f_k27_6 & !f_far; 
	#boxplot(split(latac[f], cgbin[f]), boxwex=0.2, las=2, col="darkblue", at=0.6+1:12, add=T, pch=19, cex=0.3,outline = FALSE)
	#abline(h=6)
	#abline(h=8)
	#abline(h=10)
    #    }
    #print(plt5())
    #save_baseR_to_ppt(plot_func = plt5(),link_ppt = "figs/cg_k27_atac_near.pptx" )
	##dev.off()

}

#####


save_gg_to_ppt = function(gg,link){
library(officer)
library(rvg)
if(file.exists( link)) {
doc <- officer::read_pptx(link) 

editable_graph <- dml(ggobj = gg,fonts = 'Arial')



doc <- ph_with(x = doc, editable_graph,
   location = ph_location_type(type = "body") )
print(doc ,target = link)

doc <- add_slide(doc, layout = "Title and Content", master = "Office Theme")
    }
    else{message('create an emty presentation file')}
    }

save_baseR_to_ppt = function(plot_func,link_ppt ){
library(officer)
library(rvg)
ppt <- read_pptx()
ppt <- add_slide(ppt, layout = "Title and Content", master = "Office Theme")
ppt <- ph_with(ppt, value = "My R Plot", 
               location = ph_location_type(type = "title"))
ppt <- ph_with(ppt, value = dml(code =plot_func), 
               location = ph_location_type(type = "body"))
print(ppt, target = link_ppt)
}



pcg_build_gw_feats_hg = function (mod, quick_mode = F) 
{
    set.seed(42)
	all_c = gintervals.all()
    all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"), ]
    options(gmax.data.size = 1e+09)
    gen_track_vext()
    gen_k27_vt(mod)
    f_pcg = mod$seqmod_loc_pred5mc>.5 & mod$seqmod_loc_pred > 5
    f_txg = mod$seqmod_loc_pred5mc>.5 & mod$seqmod_loc_pred <3.8
    f_mix = mod$seqmod_loc_pred5mc>.5 & !f_pcg & !f_txg
    f_5mc = mod$seqmod_loc_pred5mc <= .5
    cgd_pcg = mod$cgdom_ann[f_pcg, 1:3]
    cgd_txg = mod$cgdom_ann[f_txg, 1:3]
    cgd_mix = mod$cgdom_ann[f_mix, 1:3]
    cgd_5mc = mod$cgdom_ann[f_5mc, 1:3]
    gvtrack.create("d_cgd_pcg", cgd_pcg, "distance")
    gvtrack.create("d_cgd_txg", cgd_txg, "distance")
    gvtrack.create("d_cgd_mix", cgd_mix, "distance")
    gvtrack.create("d_cgd_5mc", cgd_5mc, "distance")
    cg_trace = gextract(c("d_cgd_pcg", "d_cgd_txg", "d_cgd_mix", 'd_cgd_5mc',
        "atac_ext", "h1_1", "h1_1", "h1_1", "h1_1", "h1_1", "seq.CG_500_mean"), 
        intervals = all_c, iterator = 200, colnames = c("d_pcg", 
            "d_txg", "d_mix",'d_5mc', "atac", "k27", "k27_eb3", "k27_es", 
            "k27_ecto", "k27_emeso", "cg_all"))
    cg_trace$cg_all[is.na(cg_trace$cg_all)] = 0
    cg_trace$cg_all = ifelse(cg_trace$cg_all > 0.04, 1, 0)
    cgd_pcg$strand = 1
    cgd_txg$strand = 1
    cgd_mix$strand = 1
    cgd_5mc$strand = 1
    print('cgd_5mc')
    print(head(cgd_5mc))
    pcg_5 = gintervals.neighbors(cg_trace, cgd_pcg, mindist = 0, 
        na = T)
    pcg_3 = gintervals.neighbors(cg_trace, cgd_pcg, maxdist = 0, 
        na = T)
    mix_5 = gintervals.neighbors(cg_trace, cgd_mix, mindist = 0, 
        na = T)
    mix_3 = gintervals.neighbors(cg_trace, cgd_mix, maxdist = 0, 
        na = T)
    txg_5 = gintervals.neighbors(cg_trace, cgd_txg, mindist = 0, 
        na = T)
    txg_3 = gintervals.neighbors(cg_trace, cgd_txg, maxdist = 0, 
        na = T)
    mc_5 = gintervals.neighbors(cg_trace, cgd_5mc, mindist = 0, 
        na = T)
    mc_3 = gintervals.neighbors(cg_trace, cgd_5mc, maxdist = 0, 
        na = T)    
    cg_trace$d_pcg5 = ifelse(is.na(pcg_5$dist), 4e+06, pmin(pcg_5$dist, 
        4e+06))
    cg_trace$d_pcg3 = ifelse(is.na(pcg_3$dist), -4e+06, pmax(pcg_3$dist, 
        -4e+06))
    cg_trace$d_mix5 = ifelse(is.na(mix_5$dist), 4e+06, pmin(mix_5$dist, 
        4e+06))
    cg_trace$d_mix3 = ifelse(is.na(mix_3$dist), -4e+06, pmax(mix_3$dist, 
        -4e+06))
    cg_trace$d_txg5 = ifelse(is.na(txg_5$dist), 4e+06, pmin(txg_5$dist, 
        4e+06))
    cg_trace$d_txg3 = ifelse(is.na(txg_3$dist), -4e+06, pmax(txg_3$dist, 
        -4e+06))
    cg_trace$d_5mc5 = ifelse(is.na(mc_5$dist), 4e+06, pmin(mc_5$dist, 
        4e+06))
    cg_trace$d_5mc3 = ifelse(is.na(mc_3$dist), -4e+06, pmax(mc_3$dist, 
        -4e+06))    
    gc_trace = gextract(c("seq.GC_500_mean", "seq.CG_500_mean", 
        "seq.GC_500_mean", "tn5bias_mean200bp", "seq.GC_500_mean"), 
        intervals = all_c, iterator = 200, colnames = c("GC", 
            "CG", "IQbase", "tn5bias", "deep_base"))
    f = is.infinite(gc_trace$GC) | is.na(gc_trace$GC)
    gc_traceGC = 0.3
    gc_trace$GC[is.na(gc_trace$GC)] = min(cg_trace$GC, na.rm = T)
    gc_trace$CG[is.na(gc_trace$CG)] = min(cg_trace$CG, na.rm = T)
    cg_trace$GC = gc_trace$GC
    cg_trace$CG = gc_trace$CG
    cg_trace$tn5bias = gc_trace$tn5bias
    gvtrack.create("blacklist", "ENCODE.blacklist", "distance")
    gvtrack.create("bgc", "intervs.global.bgc", "distance")
    gvtrack.create("d_ltr", "intervs.global.rmsk_ltr", "distance")
    gvtrack.create("d_line", "intervs.global.rmsk_line", "distance")
    gvtrack.create("d_sine", "intervs.global.rmsk_sine", "distance")
	gvtrack.create("simp_d", "intervs.global.rmsk_simple_repeat", "distance")
	gvtrack.create("lowcomplex_d", "intervs.global.rmsk_low_complexity", "distance")	
	
	
	
    rpt_trace = gextract(c("d_line", "d_ltr", "d_sine", "bgc", 
        "blacklist", "mapab.length_50", "mapab.wgEncodeDukeMapabilityUniqueness20bp","simp_d","lowcomplex_d"), 
        intervals = all_c, iterator = 200, colnames = c("d_line", 
            "d_ltr", "d_sine", "d_bgc", "d_blacklist", "mapab_50", 
            "Duke20bp","simp_d","lowcomplex_d"))
    cg_trace$d_ltr = rpt_trace$d_ltr
    cg_trace$d_line = rpt_trace$d_line
    cg_trace$d_sine = rpt_trace$d_sine
    cg_trace$d_bgc = rpt_trace$d_bgc
    cg_trace$d_blacklist = rpt_trace$d_blacklist
    cg_trace$mapab_50 = rpt_trace$mapab_50
    cg_trace$Duke20bp = rpt_trace$Duke20bp
	cg_trace$simp_d = rpt_trace$simp_d
	cg_trace$lowcomplex_d = rpt_trace$lowcomplex_d	
    cg_trace$k27[is.na(cg_trace$k27)] = 0
    cg_trace$atac[is.na(cg_trace$atac)] = 0
    is_pcg = ifelse(cg_trace$d_pcg == 0, 1, 0)
    is_txg = ifelse(cg_trace$d_txg == 0, 1, 0)
    is_mix = ifelse(cg_trace$d_mix == 0, 1, 0)
    is_5mc = ifelse(cg_trace$d_5mc == 0, 1, 0)
    k27_1k = zoo::rollmean(cg_trace$k27, 5, f = "e")
    cg_trace$k27_1k = k27_1k
    cg_trace$lk27_1k = log2(2 + k27_1k)
    cg_sc = matrix(cg_trace$cg_all, ncol = 1)
    cgpcg_sc = matrix(is_pcg, ncol = 1)
    cgtxg_sc = matrix(is_txg, ncol = 1)
    cgmix_sc = matrix(is_mix, ncol = 1)
    cg5mc_sc = matrix(is_5mc, ncol = 1)
    print(('cg5mc_sc'))
    print(head(cg5mc_sc))
    for (scale in 2^(1:12)) {
        message("scale ", scale)
        cg_sc = cbind(cg_sc, zoo::rollmean(cg_trace$cg_all, 1 + scale, 
            f = "e"))
        cgpcg_sc = cbind(cgpcg_sc, zoo::rollmean(is_pcg, 1 + scale, 
            f = "e"))
        cgtxg_sc = cbind(cgtxg_sc, zoo::rollmean(is_txg, 1 + scale, 
            f = "e"))
        cgmix_sc = cbind(cgmix_sc, zoo::rollmean(is_mix, 1 + scale, 
            f = "e"))
        cg5mc_sc = cbind(cg5mc_sc, zoo::rollmean(is_5mc, 1 + scale, 
            f = "e"))        
        
    }
    colnames(cgmix_sc) = paste(rep("cgmix", 13), 1:13, sep = "_")
    colnames(cgpcg_sc) = paste(rep("cgpcg", 13), 1:13, sep = "_")
    colnames(cgtxg_sc) = paste(rep("cgtxg", 13), 1:13, sep = "_")
    colnames(cg5mc_sc) = paste(rep("cg5mc", 13), 1:13, sep = "_")    
    n = nrow(cg_trace)
    for (i in 6:12) {
        scale = 2^(i - 1)
        half_win = scale/2
        pcg_pad3 = c(cgpcg_sc[-(1:half_win), i], rep(cgpcg_sc[n, 
            i], half_win))
        pcg_pad5 = c(rep(cgpcg_sc[1, i], half_win), cgpcg_sc[-((n - 
            half_win + 1):n), i])
        txg_pad3 = c(cgtxg_sc[-(1:half_win), i], rep(cgtxg_sc[n, 
            i], half_win))
        txg_pad5 = c(rep(cgtxg_sc[1, i], half_win), cgtxg_sc[-((n - 
            half_win + 1):n), i])
        mix_pad3 = c(cgmix_sc[-(1:half_win), i], rep(cgmix_sc[n, 
            i], half_win))
        mix_pad5 = c(rep(cgmix_sc[1, i], half_win), cgmix_sc[-((n - 
            half_win + 1):n), i])
        mc_pad3 = c(cg5mc_sc[-(1:half_win), i], rep(cg5mc_sc[n, 
            i], half_win))
        mc_pad5 = c(rep(cg5mc_sc[1, i], half_win), cg5mc_sc[-((n - 
            half_win + 1):n), i])        
        if (i == 6) {
            cgpcg_3 = matrix(pcg_pad3, ncol = 1)
            cgpcg_5 = matrix(pcg_pad5, ncol = 1)
            cgtxg_3 = matrix(txg_pad3, ncol = 1)
            cgtxg_5 = matrix(txg_pad5, ncol = 1)
            cgmix_3 = matrix(mix_pad3, ncol = 1)
            cgmix_5 = matrix(mix_pad5, ncol = 1)
            cg5mc_3 = matrix(mc_pad3, ncol = 1)
            cg5mc_5 = matrix(mc_pad5, ncol = 1)
        }
        else {
            cgpcg_3 = cbind(cgpcg_3, pcg_pad3)
            cgpcg_5 = cbind(cgpcg_5, pcg_pad5)
            cgtxg_3 = cbind(cgtxg_3, txg_pad3)
            cgtxg_5 = cbind(cgtxg_5, txg_pad5)
            cgmix_3 = cbind(cgmix_3, mix_pad3)
            cgmix_5 = cbind(cgmix_5, mix_pad5)
            cg5mc_3 = cbind(cg5mc_3, mc_pad3)
            cg5mc_5 = cbind(cg5mc_5, mc_pad5)            
        }
    }
    colnames(cgmix_3) = paste(rep("cgmix3", 7), 1:7, sep = "_")
    colnames(cgpcg_3) = paste(rep("cgpcg3", 7), 1:7, sep = "_")
    colnames(cgtxg_3) = paste(rep("cgtxg3", 7), 1:7, sep = "_")
    colnames(cg5mc_3) = paste(rep("cg5mc3", 7), 1:7, sep = "_")    
    colnames(cgmix_5) = paste(rep("cgmix5", 7), 1:7, sep = "_")
    colnames(cgpcg_5) = paste(rep("cgpcg5", 7), 1:7, sep = "_")
    colnames(cgtxg_5) = paste(rep("cgtxg5", 7), 1:7, sep = "_")
    colnames(cg5mc_5) = paste(rep("cg5mc5", 7), 1:7, sep = "_")    
    inv_d_txg = 1/(cg_trace$d_txg + 100)
    inv_d_mix = 1/(cg_trace$d_mix + 100)
    inv_d_pcg = 1/(cg_trace$d_pcg + 100)
    inv_d_5mc = 1/(cg_trace$d_5mc + 100)    
    min_d_cg = pmax(pmax(inv_d_txg, inv_d_mix), inv_d_pcg)
    min_d_cg5mc = pmax(min_d_cg,inv_d_5mc)    
    iq_score = misha.ext::gintervals.neighbors1(cg_trace[,c('chrom','start','end')],mod$cgdom_ann[,c('chrom','start','end','pred')],maxneighbors = 1)
    iq_score5mc = misha.ext::gintervals.neighbors1(cg_trace[,c('chrom','start','end')],mod$cgdom_ann_orig[,c('chrom','start','end','pred5mc')],maxneighbors = 1)

    
    feats = cbind(ld_mix = log2(inv_d_mix), ld_pcg = log2(inv_d_pcg), 
        ld_txg = log2(inv_d_txg), min_dcg = min_d_cg, cgmix_sc[, 
            1:11], cgtxg_sc[, 1:11], cgpcg_sc[, 1:11], GC = gc_trace$GC, 
        CG = gc_trace$CG, tn5b = gc_trace$tn5bias, dltr = abs(rpt_trace$d_ltr), 
        dline = abs(rpt_trace$d_line), ldltr = log2(1 + abs(rpt_trace$d_ltr)), 
        ldline = log2(1 + abs(rpt_trace$d_line)),
		#min_d5mc = min_d_cg5mc, 
		cg5mc_sc[, 1:11],
          ld_5mc = log2(inv_d_5mc), dblck =  abs(rpt_trace$d_blacklist), ldblck = log2(1 + abs(rpt_trace$d_blacklist)),  duke20 = rpt_trace$Duke20bp ,
		  dsine =  abs(rpt_trace$d_sine), ldsine = log2(1 + abs(rpt_trace$d_sine)),
		  dlowc =  abs(rpt_trace$lowcomplex_d), ldlowc = log2(1 + abs(rpt_trace$lowcomplex_d)),		  
		  dsimpr =  abs(rpt_trace$simp_d), ldsimr = log2(1 + abs(rpt_trace$simp_d)),	
		  dbgc =  abs(rpt_trace$d_bgc), ldbgc = log2(1 + abs(rpt_trace$d_bgc)),		  
                  iq_score = as.numeric(iq_score$pred),iq_score5mc = as.numeric(iq_score5mc$pred5mc))
                 
                 
                 
    feats35 = cbind(cgmix_3, cgmix_5, cgtxg_3, cgtxg_5, cgpcg_3, 
        cgpcg_5,cg5mc_3,cg5mc_5)
    mod$gw = list()
    mod$gw$cg_trace = cg_trace
    mod$gw$cgpcg_sc = cgpcg_sc
    mod$gw$cgtxg_sc = cgtxg_sc
    mod$gw$cgmix_sc = cgmix_sc
    mod$gw$feats = feats
    mod$gw$feats35 = feats35
    f = is.infinite(mod$gw$feats[, "GC"]) | is.na(mod$gw$feats[, 
        "GC"])
    mod$gw$feats[f, "GC"] = 0.3
    f = is.infinite(mod$gw$feats[, "CG"]) | is.na(mod$gw$feats[, 
        "CG"])
    mod$gw$feats[f, "CG"] = 0.001
    return(mod)
}


pcg_build_gw_pred_hg = function(mod,  rebuild_iqdn=F, rebuild_base=F, rebuild_atac=F,perc_genome=0)
{
	gw = mod$gw
	cg_trace = gw$cg_trace
	cg_trace$Duke20bp[is.na(cg_trace$Duke20bp)] = 0
	cg_trace$d_bgc[is.na(cg_trace$d_bgc)] = 1e6
	f_norp = cg_trace$d_ltr !=0 & cg_trace$d_line!=0 & cg_trace$simp_d!=0 & cg_trace$d_sine!=0 & cg_trace$lowcomplex_d!=0 & cg_trace$d_blacklist!=0 & 
	rowSums(is.na(cg_trace))==0 & rowSums(is.na(gw$feats))==0  & cg_trace$Duke20bp>.5 #& cg_trace$start > 3e+6
	p = runif(n=nrow(gw$cg_trace))
	f_sub = f_norp & (p>= perc_genome | cg_trace$CG>0.02) 
#& cg_trace$cg_all!=1

	ptrain = runif(n=nrow(gw$cg_trace))
	#f_test_chr = cg_trace$chrom %in% c("chr1", "chr4", "chr6","chr10","chr13", "chr16", "chr19")
	f_test_chr = cg_trace$chrom %in% c("chr2", "chr3", "chr6","chr9","chr13", "chr16", "chr18")
	f_train_chr =!f_test_chr

	f_test = f_test_chr & f_sub 
	f_test_all = f_test_chr & f_norp
	f_train = f_train_chr & f_sub 
	f_train_all = f_train_chr & f_norp

	feats = gw$feats

	mix = rowSums(feats[,paste("cgmix_", 4:11, sep="")])
	pcg = rowSums(feats[,paste("cgpcg_", 4:11, sep="")])
	mix2 = mix**2
	pcg2 = pcg**2
	feats = cbind(feats, mix2=mix2, pcg2=pcg2)
	
	if(rebuild_iqdn) {
		feats_iqdn = cbind(feats, gw$feats_iqdn)
		formu2iq = paste("resp ~ ", paste(colnames(feats2iq),collapse="+"))
		formu2iq = paste(formu2iq," + GC*CG",sep="")
		mod = pcg_build_mod_from_feats_hg(mod, resp = cg_trace$lk27_1k, 
								formu = formu2iq, feats= feats2iq, 
								f_train = f_train, f_test = f_test, 
								f_train_all = f_train_all, f_test_all = f_test_all, tag = "base2iqdn")
	}

	if(rebuild_base) {
		formu = paste("resp ~ ", paste(colnames(feats),collapse="+"))
		formu = paste(formu," + GC*CG",sep="")
		mod = pcg_build_mod_from_feats_hg(mod, resp = cg_trace$lk27_1k, 
								formu = formu, feats= feats, 
								f_train = f_train, f_test = f_test, 
								f_train_all = f_train_all, f_test_all = f_test_all, tag = "base")
	}

	if(rebuild_atac) {
		feats_atac = cbind(feats, gw$feats_iqdn, atac = log2(10+cg_trace$atac))
		formu_atac = paste("resp ~ ", paste(colnames(feats_atac),collapse="+"))
		formu_atac = paste(formu_atac, " + GC*CG",sep="")
		mod = pcg_build_mod_from_feats_hg(mod, resp = cg_trace$lk27_1k, 
								formu = formu_atac, feats= feats_atac,
								f_train = f_train, f_test = f_test, 
								f_train_all = f_train_all, f_test_all = f_test_all, tag = "base_atac")
	}
	return(mod)
}

pcg_build_mod_from_feats_hg = function(mod, resp, formu, feats, f_train, f_test, f_train_all, f_test_all, tag)
{
	#extract features from formu
	fnms = colnames(feats)
	fnms_f = c()
	for(fnm in fnms) {
		if(grepl(fnm, formu)) {
			fnms_f = c(fnms_f, fnm)
		}
	}
	if(length(fnms_f)==0) {
		message("not overlap of features and formula ", formu)
		return
	}

	feats = feats[,fnms_f]
	stats = c(formu = formu)

#train linear and extract QC

	feats_r = as.data.frame(feats)
	feats_r$resp = resp
	lmod = lm(formu, feats_r[f_train,])
	
	lm_pred = predict(lmod, feats_r)
	lm_pred = pmax(lm_pred, 1)
	lm_stat = c(lm_r2_test = cor(resp[f_test], lm_pred[f_test])**2,
					lm_r2_train = cor(resp[f_train], lm_pred[f_train])**2,
					lm_r2_test_all = cor(resp[f_test_all], lm_pred[f_test_all])**2,
					lm_r2_train_all = cor(resp[f_train_all], lm_pred[f_train_all])**2,
					lm_auc_999_train = comp_auc(lm_pred, ifelse(resp > 7.16,1,0), f_train_all)$auc,
					lm_auc_999_test = comp_auc(lm_pred, ifelse(resp > 7.16,1,0), f_test_all)$auc,
					lm_auc_998_train = comp_auc(lm_pred, ifelse(resp > 6.4,1,0), f_train_all)$auc,
					lm_auc_998_test = comp_auc(lm_pred, ifelse(resp > 6.4,1,0), f_test_all)$auc,
					lm_auc_997_train = comp_auc(lm_pred, ifelse(resp > 5.8,1,0), f_train_all)$auc,
					lm_auc_997_test = comp_auc(lm_pred, ifelse(resp > 5.8,1,0), f_test_all)$auc,
					lm_auc_994_train = comp_auc(lm_pred, ifelse(resp > 4.3,1,0), f_train_all)$auc,
					lm_auc_994_test = comp_auc(lm_pred, ifelse(resp > 4.3,1,0), f_test_all)$auc)

	dtrain = xgb.DMatrix(as.matrix(feats)[f_train,], label = resp[f_train])
	dtest = xgb.DMatrix(as.matrix(feats)[f_test,], label = resp[f_test])
set.seed(42)
#cv_model <- xgb.cv(
 #eta = 0.05, max_depth=6,
#  data = dtrain,
 # nrounds = 1000,
 # nfold = 5,
 # early_stopping_rounds = 10,
 #   nthread = 56,
  
 # metrics = "rmse",
 # verbose = 1
#)

#best_nrounds <- cv_model$best_iteration
#print(cat("Optimal number of rounds:", best_nrounds, "\n"))
	message("start training ")
	xgb_mod <- xgb.train(
 	 	  	data = dtrain,
			verbose=0,
    		nrounds = 500, objective = "reg:squarederror", 
			early_stopping_rounds = 3,
	 		eta = 0.05, max_depth=4,
			#watchlist = list(train = dtrain, test = dtest)
			  evals = list(train = dtrain, test = dtest),  # ← use evals, not watchlist
  eval_metric = "rmse",                         # ← recommended for early stopping
			)

	xg_pred = predict(xgb_mod, feats)
	xg_stat = c(xg_r2_test = cor(resp[f_test], xg_pred[f_test])**2,
					xg_r2_train = cor(resp[f_train], xg_pred[f_train])**2,
					xg_r2_test_all = cor(resp[f_test_all], xg_pred[f_test_all])**2,
					xg_r2_train_all = cor(resp[f_train_all], xg_pred[f_train_all])**2,
					xg_auc_999_train = comp_auc(xg_pred, ifelse(resp > 7.16,1,0), f_train_all)$auc,
					xg_auc_999_test = comp_auc(xg_pred, ifelse(resp > 7.16,1,0), f_test_all)$auc,					
					xg_auc_998_train = comp_auc(xg_pred, ifelse(resp > 6.4,1,0), f_train_all)$auc,
					xg_auc_998_test = comp_auc(xg_pred, ifelse(resp > 6.4,1,0), f_test_all)$auc,
					xg_auc_997_train = comp_auc(xg_pred, ifelse(resp > 5.8,1,0), f_train_all)$auc,
					xg_auc_997_test = comp_auc(xg_pred, ifelse(resp > 5.8,1,0), f_test_all)$auc,
					xg_auc_994_train = comp_auc(xg_pred, ifelse(resp > 4.3,1,0), f_train_all)$auc,
					xg_auc_994_test = comp_auc(xg_pred, ifelse(resp > 4.3,1,0), f_test_all)$auc,
					xg_n_under4 = sum((xg_pred < 2 & resp > 6)[f_test_all | f_train_all]),
					xg_n_overr4 = sum((xg_pred > 6 & resp < 2)[f_test_all | f_train_all]))

	mod$gw$xg_mods[[tag]] = xgb_mod
	mod$gw$lm_mods[[tag]] = lmod
	mod$gw$xg_pred[[tag]] = xg_pred
	mod$gw$lm_pred[[tag]] = lm_pred
	mod$gw$mods_stat[[tag]] = list(formu=formu, lm_stat = lm_stat, xg_stat = xg_stat)
	return(mod)

}


generate_cooperativity_hg = function(iq,t_k4=3.8,t_k27=5){
    
  k4 = iq[ iq$pred5mc>.5& iq$pred<t_k4 ,]
k27 = iq[ iq$pred5mc>.5&iq$pred>t_k27 ,]

iq$d_k4 =  misha.ext::gintervals.neighbors1(iq,k4,mindist = 1)$dist

iq$d_k27 =  misha.ext::gintervals.neighbors1(iq,k27,mindist = 1)$dist



#max pred split to near k4 neark27
#min pred to near k4 neark27

iq$bin = ifelse(iq$pred>t_k27  & iq$d_k4<10000 & iq$d_k27>10000,'pred_max_near_k4','other')
iq$bin = ifelse(iq$pred>t_k27  & iq$d_k4>10000 & iq$d_k27<10000,'pred_max_near_k27',iq$bin)

iq$bin = ifelse(iq$pred<t_k4  & iq$d_k4<10000 & iq$d_k27>10000,'pred_min_near_k4',iq$bin)
iq$bin = ifelse(iq$pred<t_k4  & iq$d_k4>10000 & iq$d_k27<10000,'pred_min_near_k27',iq$bin)


print(table(iq$bin))

gg1 = iq %>%filter(bin%in%c('pred_max_near_k27',  'pred_max_near_k4'))%>% ggplot(aes(bin,log2(4+es_k27_max)))+ 
geom_boxplot(outlier.shape = NA)+theme_bw()+
coord_cartesian(ylim = c(2,8))
print(gg1)

gg2 = iq %>%filter(bin%in%c('pred_min_near_k27',  'pred_min_near_k4'))%>% ggplot(aes(bin,log2(4+es_k27_max))) + 
geom_boxplot(outlier.shape = NA)+theme_bw()+
coord_cartesian(ylim = c(2,8))
print(gg2)  
save_gg_to_ppt(gg1,"./figs/cooperativity_hg.pptx" )      
save_gg_to_ppt(gg2,"./figs/cooperativity2_hg.pptx" )    
}


generate_decay_hg = function(){
cg_trace = mod$gw$cg_trace

f_txg =  cg_trace$d_txg < cg_trace$d_pcg & cg_trace$d_txg < cg_trace$d_mix

f_mix =  cg_trace$d_mix < cg_trace$d_pcg & cg_trace$d_mix < cg_trace$d_txg
f_pcg = cg_trace$d_pcg < cg_trace$d_txg & cg_trace$d_pcg < cg_trace$d_mix

gg_pcg = signal_decay_hg(f=f_pcg,mod=mod,color='darkblue',type = 'd_pcg')
print(gg_pcg)

gg_txg = signal_decay_hg(f=f_txg,mod=mod,color='red',type='d_txg')
print(gg_txg)

gg_mix = signal_decay_hg(f=f_mix,mod=mod,color='gold',type='d_mix')
print(gg_mix)

save_gg_to_ppt(gg_mix,'./figs/decay1_hg.pptx')

save_gg_to_ppt(gg_txg,'./figs/decay2_hg.pptx')

save_gg_to_ppt(gg_pcg,'./figs/decay3_hg.pptx')
    }

signal_decay_hg = function(f,mod,color,type){

cg_trace = mod$gw$cg_trace
trace_pcg = cg_trace[f,]

trace_pcg$bin = cut(trace_pcg$CG,c(-1,.005,0.01,.015,0.035))

trace_pcg = trace_pcg[! is.na(trace_pcg$bin),]

trace_pcg$d_bin = floor(log2(100+trace_pcg[,type]))

trace_gg = trace_pcg %>% group_by(bin,d_bin) %>% summarise(mean(k27_1k))

colnames(trace_gg) = c('cg_bin','d_bin','k271k')

shades=colorRampPalette(c("gray", color))

gg = trace_gg%>% ggplot(aes(x= d_bin,y = log2(k271k),col = cg_bin))+ geom_line()+
theme_bw()+xlim(c(8,19))+ylim(c(-2,6.5))+scale_color_manual(values = shades(4))+
geom_abline(slope = -0.85,intercept = 12)
return(gg)
}

generate_decay_hg_new = function(){
cg_trace = mod$gw$cg_trace
f_norp = cg_trace$d_ltr !=0 & cg_trace$d_line!=0 &  cg_trace$d_sine != 0
cg_trace = cg_trace[f_norp,]
cg_trace = cg_trace %>% arrange(chrom,start)
f_txg =  cg_trace$d_txg < cg_trace$d_pcg & cg_trace$d_txg < cg_trace$d_mix

f_mix =  cg_trace$d_mix < cg_trace$d_pcg & cg_trace$d_mix < cg_trace$d_txg
f_pcg = cg_trace$d_pcg < cg_trace$d_txg & cg_trace$d_pcg < cg_trace$d_mix

gg_pcg = signal_decay_hg_new(f=f_pcg,mod=mod,color='darkblue',type = 'd_pcg')
print(gg_pcg)

gg_txg = signal_decay_hg_new(f=f_txg,mod=mod,color='red',type='d_txg')
print(gg_txg)

gg_mix = signal_decay_hg_new(f=f_mix,mod=mod,color='gold',type='d_mix')
print(gg_mix)

save_gg_to_ppt(gg_mix,'./figs/decay1_hg.pptx')

save_gg_to_ppt(gg_txg,'./figs/decay2_hg.pptx')

save_gg_to_ppt(gg_pcg,'./figs/decay3_hg.pptx')
    }

signal_decay_hg_new = function(f,mod,color,type){

cg_trace = mod$gw$cg_trace
gvtrack.create("k27_h1_sum", "jk.henikoff_n19.h1_h3k27me3_2reps", "sum")    
f_norp = cg_trace$d_ltr !=0 & cg_trace$d_line!=0 &  cg_trace$d_sine != 0
cg_trace = cg_trace[f_norp,]
cg_trace = cg_trace %>% arrange(chrom,start)
gext = gextract("k27_h1_sum",intervals =cg_trace,iterator = cg_trace)
gext[is.na(gext)] = 0
cg_trace$k27_h1_sum = gext$k27_h1_sum
trace_pcg = cg_trace[f,]

trace_pcg$bin = cut(trace_pcg$CG,c(-1,.005,0.01,.015,0.035))

trace_pcg = trace_pcg[! is.na(trace_pcg$bin),]

trace_pcg$d_bin = floor(log2(100+trace_pcg[,type]))

trace_gg = trace_pcg %>% group_by(bin,d_bin) %>% summarise(mean(k27_h1_sum))

colnames(trace_gg) = c('cg_bin','d_bin','k27_h1_sum')

shades=colorRampPalette(c("gray", color))

gg = trace_gg%>% ggplot(aes(x= d_bin,y = log2(k27_h1_sum),col = cg_bin))+ geom_line()+
theme_bw()+xlim(c(8,19))+ylim(c(-2,5.5))+scale_color_manual(values = shades(4))+
geom_hline(yintercept = -1)+geom_hline(yintercept = 1)+geom_hline(yintercept = 3)+
geom_hline(yintercept = 5)+
theme(panel.grid = element_blank())+
#geom_vline(xintercept = log2(500+100))+
geom_vline(xintercept = log2(1000+100))+
#geom_vline(xintercept = log2(2000+100))+
geom_vline(xintercept = log2(4000+100))+
#geom_vline(xintercept = log2(8000+100))+
geom_vline(xintercept = log2(16000+100))+
#geom_vline(xintercept = log2(32000+100))+
geom_vline(xintercept = log2(64000+100))+
#geom_vline(xintercept = log2(128000+100))+
geom_vline(xintercept = log2(256000+100))+
#geom_vline(xintercept = log2(512000+100))+
geom_abline(slope = -.85,intercept = 12)
geom_abline(slope = -0.85,intercept = 12)
return(gg)
}
pcg_report_locmod_iq_pptx_gw_hg = function(mod, fit_type='IQ')
{
#prediction - f_out, f_ambig, f_train, f_test
#vs. lk27, k4_cov_1
#scatter and box
	set.seed(42)
	gw = mod$gw
	cg_trace = gw$cg_trace
	
    f_norp = cg_trace$d_ltr !=0 & cg_trace$d_line!=0 & cg_trace$simp_d!=0 & cg_trace$d_sine!=0 & cg_trace$lowcomplex_d!=0 & cg_trace$d_blacklist!=0 & 
	rowSums(is.na(cg_trace))==0 & rowSums(is.na(gw$feats))==0  & cg_trace$Duke20bp>.5 #& cg_trace$start > 3e+6
	p = runif(n=nrow(gw$cg_trace))
	f_sub = f_norp & (p>= 0 | cg_trace$CG>0.02) 
#& cg_trace$cg_all!=1

	ptrain = runif(n=nrow(gw$cg_trace))
	f_test_chr = cg_trace$chrom %in% c("chr2", "chr3", "chr6","chr9","chr13", "chr16", "chr18")
	f_train_chr =!f_test_chr

	f_test = f_test_chr & f_sub 
	f_test_all = f_test_chr & f_norp
	f_train = f_train_chr & f_sub 
	f_train_all = f_train_chr & f_norp
    lk27 = cg_trace$lk27_1k
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)
	pred = cg_trace$pred
#	to_fit = mod$locmod[[fit_type]]$to_fit
	
	to_fit = ifelse(cg_trace$lk27_1k> 6 ,1,0)

	#png(sprintf("figs/model_%s_preds.png", fit_type), w=1000,h=800)
	plt1 = function(){
    layout(matrix(1,nrow=1))
	for(tss in c(1)) {
		#if(tss == 1) { 
		#	filt = f_tss
		#} else {
		#	filt = !f_tss
		#}
		f1 = f_train
		f2 = f_test
		#f3 = filt & f_500 & f_ambig 
		#f4 = filt & f_500 & f_out
	for(i in 1) {
		if(i == 1) {
			stat = lk27
		} else {
			stat = k4_cov_1
		}
		tr_cor = round(cor(stat[f1], pred[f1],m="p")^2,3)
		te_cor = round(cor(stat[f2], pred[f2], m="p")^2,3)
		
		message("comp ranges")
		#ranges = c(min(pred),seq(ceiling(quantile(pred,0.05)), floor(quantile(pred,0.95)), 0.5), max(pred))
        ranges = c(0,1,2,3,4,5,6,7,9)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred[f1],ranges)), boxwex=0.15, col="blue",outline=FALSE, 
								main=sprintf("%s, tss = %s, rsqrd tr = %s, rsqrd te = %s", fit_type, tss, tr_cor, te_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col="gold", add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}}
    save_baseR_to_ppt(plt1(),'./figs/bxplt_gw_hg.pptx')
	#dev.off()

	#pdf(sprintf("figs/model_roc_%s.pdf", fit_type), w=6, h=6)
    plt_auc = function(){
	f1 = f_train
	f2 = f_test
	#f3 = !f_tss & f_500 & !f_ambig & !f_out & f_train
	#f4 = !f_tss& f_500 & !f_ambig & !f_out & f_test
	roc1 = comp_auc(pred, to_fit, f1)
	roc2 = comp_auc(pred, to_fit, f2)
	#roc3 = comp_auc(pred, to_fit, f3)
	#roc4 = comp_auc(pred, to_fit, f4)

	plot(roc1$fn, roc1$tp, t="l", lwd=2, col="blue", xlab="fn", ylab="tp",
					 main=sprintf("auc %s,%s, train np = %s, na = %s, nn = %s", 
										round(roc1$auc,3), round(roc2$auc,3),
										sum(f1 & !is.na(to_fit) & to_fit),
										sum(f1 & is.na(to_fit)),
										sum(f1 & !is.na(to_fit) & !to_fit)))
	lines(roc2$fn, roc2$tp, t="l", lwd=2, col="gold")
	#lines(roc3$fn, roc3$tp, t="l", lwd=2, "blue", lty=2)
	#lines(roc4$fn, roc4$tp, t="l", lwd=2, col="gold", lty=2)
	abline(a=0,b=1)
    }
    #save_baseR_to_ppt(plt_auc(),'./figs/auc_gw_hg.pptx')
	#dev.off()

}

plt_genome_pred_ppt_hg = function(mod, g = NA, off5 = NA, off3 = NA,
			 chrom = NULL, locus = NULL, win_smooth = 1, show_center=T, 
			plot_pred_doms = F,
			label_tss = F,
			pred_lwd=1, plot_base_pred = F,
			 fn=NULL, fn_w = 15, fn_h=10, more_tracks=c())
{
#	if(!is.null(fn)) {
#		if(grepl("pdf", fn)) {
#			pdf(fn, w=fn_w, h=fn_h)
#		} else {
#			png(fn, w=fn_w, h=fn_h)
#		}
#	}
	cg_trace = mod$gw$cg_trace
	if(is.null(chrom)) {
		f = mod$tss$geneSymbol==g
		if(sum(f) == 0) {
			return
		}
		hits = mod$tss[f,]
		locus = mean(mod$tss$start[f])
		chrom = hits[1,"chrom"]
	}
	f = as.character(cg_trace$chrom) == chrom & (cg_trace$start > locus + off5 & cg_trace$start < locus + off3)

	tss = mod$epi_tss[mod$epi_tss$chrom == chrom & (mod$epi_tss$start > locus + off5 & mod$epi_tss$start < locus + off3),]
    #####
    
    f_pcg = mod$seqmod_loc_pred5mc>.25 & mod$seqmod_loc_pred < 0.4
    f_txg = mod$seqmod_loc_pred5mc>.25 & mod$seqmod_loc_pred > 0.95
    f_mix = mod$seqmod_loc_pred5mc>.25 & !f_pcg & !f_txg
    f_5mc = mod$seqmod_loc_pred5mc <= .25
	cgd_pcg = mod$cgdom_ann[f_pcg,1:3]
	cgd_txg = mod$cgdom_ann[f_txg,1:3]
	cgd_mix = mod$cgdom_ann[f_mix,1:3]
    cgd_5mc = mod$cgdom_ann[f_5mc,1:3]
    rpt = gintervals.diff(mod$cgd_rpt,mod$cgdom_ann)
    cgd_pcg = cgd_pcg[cgd_pcg$chrom == chrom & (cgd_pcg$start > locus + off5 & cgd_pcg$start < locus + off3),]
    cgd_txg = cgd_txg[cgd_txg$chrom == chrom & (cgd_txg$start > locus + off5 & cgd_txg$start < locus + off3),]
    cgd_mix = cgd_mix[cgd_mix$chrom == chrom & (cgd_mix$start > locus + off5 & cgd_mix$start < locus + off3),]
    cgd_5mc = cgd_5mc[cgd_5mc$chrom == chrom & (cgd_5mc$start > locus + off5 & cgd_5mc$start < locus + off3),] 
    rpt = rpt[rpt$chrom == chrom & (rpt$start > locus + off5 & rpt$start < locus + off3),]
    ######
	add_n = length(more_tracks)
    if (length(more_tracks)==0){
      layout(matrix(1:(2+add_n), ncol=1), heights=c(.95,1.3))  
    }else{
    layout(matrix(1:(2+add_n), ncol=1), heights=c(1.05,.98 ,1.4))}
	#layout(matrix(1:(2+add_n), ncol=1),h=c(rep(1, 2+add_n),1.3))
	par(mar = c(0, 3,2, 3))
	maxy = max(max(cg_trace$lk27_1k[f]),7)
	pred_sm = zoo::rollmean(cg_trace$pred[f], win_smooth, f='e')
	pred_segs = which(pred_sm[-1]>4 & pred_sm[-length(pred_sm)] < 4)
	pred_segs = c(pred_segs, which(pred_sm[-1]<4 & pred_sm[-length(pred_sm)] > 4))
	pred_segs_xs = cg_trace$start[f][pred_segs]
	obs = cg_trace$lk27_1k[f]
	#r2 = cor(pmax(obs,2), pmax(pred_sm,2))**2
    r2 = cor((obs), (pred_sm))**2
	obs_given_pred = round(sum(obs > 4 & pred_sm > 4)/(1+sum(pred_sm > 4)),3)
	pred_given_obs = round(sum(obs > 4 & pred_sm > 4)/(1+sum(obs > 4)),3)


    
	plot(cg_trace$start[f], cg_trace$lk27_1k[f], 
				t="l", ylim=c(1,maxy), lwd=1,col="blue", 
				xaxt='n', 
				main=sprintf("EB4 locus r2 = %s, obs|pred %s, pred|obs %s", 
								round(r2,3), obs_given_pred, pred_given_obs), 
				xlab=sprintf("chr %s, %s", cg_trace$chrom[f][1], g))
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(cg_trace$lk27_1k[f], rep(1, length(cg_trace$lk27_1k[f]))), 
        col="blue", border=NA)
	#if(nrow(tss) > 0) {
	#	points(tss$start, rep(5,nrow(tss)), pch=19, col="darkgreen",cex=1.5)
	#}
    
   # if(nrow(cgd_pcg) > 0) {
	#	points(cgd_pcg$start,rep(2,nrow(cgd_pcg)) , pch=17, col="darkblue",cex=3)
	#}
   # if(nrow(cgd_txg) > 0) {
	#	points(cgd_txg$start,rep(4.2,nrow(cgd_txg)) , pch=17, col="darkred",cex=3)
	#}
    #if(nrow(cgd_mix) > 0) {
	#	points(cgd_mix$start,rep(3,nrow(cgd_mix)) , pch=17, col="violet",cex=3)
	#}    
    #if(nrow(cgd_5mc) > 0) {
	#	points(cgd_5mc$start,rep(1,nrow(cgd_5mc)) , pch=17, col="black",cex=3)
	#}   
    #if(nrow(rpt) > 0) {
	#	points(rpt$start,rep(1,nrow(rpt)) , pch=13, col="black",cex=3)
	#}    
	abline(h=log2(19), lty=2)
	if(plot_pred_doms) {
		segments(x0 = pred_segs_xs, x1=pred_segs_xs, 
			y0 = rep(0,length(pred_segs_xs)),
			y1 = rep(maxy,length(pred_segs_xs)), lty=1)
	}
	
	if(show_center) {
		abline(v=locus, col="black", lwd=2)
	}
    if ( length(more_tracks)==0){
    par(mar = c(5,3,1,3))}
	else{ par(mar = c(0, 3,1, 3))}
    if ( length(more_tracks)==0){
    plot(cg_trace$start[f], pred_sm, t="l", col="darkblue", ylim=c(1,maxy), xlab=(cg_trace$chrom[f][1]), lwd=pred_lwd, main="prediction")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(pred_sm, rep(1, length(pred_sm))), 
        col="darkblue", border=NA)
    }else{
    
	plot(cg_trace$start[f], pred_sm, t="l", col="darkblue", ylim=c(1,maxy), xaxt='n', lwd=pred_lwd, main="prediction")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(pred_sm, rep(1, length(pred_sm))), 
        col="darkblue", border=NA)
        }
	if(plot_base_pred) {
		pred_sm_b = zoo::rollmean(cg_trace$pred_base[f], win_smooth, f='e')
		lines(cg_trace$start[f], pred_sm_b, t="h", col="cyan", ylim=c(1,maxy), xaxt='n', lwd=pred_lwd)
	}
	if(label_tss) {
		if(nrow(tss) > 0) {
			par(srt=-45)
			text(x = tss$start, y = rep(maxy-2,t=nrow(tss)), labels=tss$geneSymbol, cex=1.5)
		}
	}

	if(plot_pred_doms) {
		segments(x0 = pred_segs_xs, x1=pred_segs_xs, 
			y0 = rep(0,length(pred_segs_xs)),
			y1 = rep(maxy,length(pred_segs_xs)), lty=1)
	}
	abline(h=log2(19), lty=2)
	if(show_center) {
		abline(v=locus, col="black",lwd=2)
	}
	if(length(more_tracks) != 0) {
		profs = gextract(more_tracks, intervals=cg_trace[f,1:3], iterator=cg_trace[f,1:3])
		for(i in 1:length(more_tracks)) {
			#sm_prf = zoo::rollmean(profs[,more_tracks[i]], 5, f='e')
            sm_prf = (profs[,more_tracks[i]])
            sm_prf[is.na(sm_prf)] = 0
            par(mar = c(5,3,1,3))
            
            #r2_d = cor(pmax(obs,2), pmax(log2(2+sm_prf),4))**2
            r2_d = cor((obs), (log2(2+sm_prf)))**2
            obs_given_pred_d = round(sum(obs > 6 & log2(2+sm_prf) > 6)/(1+sum(log2(2+sm_prf) > 6)),3)
            pred_given_obs_d = round(sum(obs > 6 & log2(2+sm_prf) > 6)/(1+sum(obs > 6)),3)
    
            
			plot(cg_trace$start[f], log2(2+sm_prf), t="l", col="gray", ylim=c(1,maxy),
                 xlab=(cg_trace$chrom[f][1]), lwd=pred_lwd, 
                 #main=more_tracks[[i]],
              main=sprintf("Enformer_celltyping locus r2 = %s, obs|pred %s, pred|obs %s", 
								round(r2_d,3), obs_given_pred_d, pred_given_obs_d))
            polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(log2(2+sm_prf), rep(1, length(log2(2+sm_prf)))), 
                col="gray", border=NA)
		}
	}

	#par(mar = c(1,3,1, 3))
	#plot(cg_trace$start[f], pmin(log2(10+cg_trace$atac[f]),9), t="l", col="red", ylim=c(3,9), xaxt='n', main="atac")
	#par(mar = c(5,3,1,3))
	#plot(cg_trace$start[f], pmin(cg_trace$CG[f],0.1), t="l", col="darkgreen", ylim=c(0,0.1), xlab=cg_trace$chrom[f][1], main="CG")
	#if(plot_pred_doms) {
	#	segments(x0 = pred_segs_xs, x1=pred_segs_xs, 
	#		y0 = rep(0,length(pred_segs_xs)),
	#		y1 = rep(maxy,length(pred_segs_xs)), lty=1)
	#}
	if(!is.null(fn)) {
		#dev.off()
	}
}



####Fig3EDF
pcg_report_locmod_borzoi_pptx = function(mod, fit_type='IQ',perc_genome=0)
{
#prediction - f_out, f_ambig, f_train, f_test
#vs. lk27, k4_cov_1
#scatter and box
	set.seed(42)
	gw = mod$gw
	 
	cg_trace = gw$cg_trace
	  f_norp = cg_trace$d_ltr != 0 & cg_trace$d_line != 0 & rowSums(is.na(cg_trace)) == 
        0 & rowSums(is.na(gw$feats)) == 0 & rowSums(is.na(gw$feats_iqdn)) == 
        0 & cg_trace$start > 3e+06
    p = runif(n = nrow(gw$cg_trace))
	f_sub = f_norp & (p >= perc_genome | cg_trace$CG > 0.02)
	#f_norp = cg_trace$d_ltr !=0 & cg_trace$d_line!=0 #& rowSums(is.na(cg_trace))==0 & 
    #rowSums(is.na(gw$feats))==0 & rowSums(is.na(gw$feats_iqdn)) == 0 & cg_trace$start > 3e+6
	
	#f_sub = f_norp 
#& cg_trace$cg_all!=1

	#ptrain = runif(n=nrow(gw$cg_trace))

	f_test_chr = cg_trace$chrom %in% c('chr9', 'chr18') & cg_trace$type =='test'
    #f_val_chr = cg_trace$chrom %in% c('chr8', 'chr10')  & cg_trace$type =='val'
	f_val_chr = cg_trace$chrom %in% c('chr8', 'chr10','chr9', 'chr18')  & cg_trace$type %in%c('val','test')
	f_train_chr =!f_test_chr & !f_val_chr

	f_test = f_test_chr & f_sub 

	f_train = f_train_chr & f_sub 
    f_val = f_val_chr & f_sub 
    lk27 = cg_trace$lk27_1k
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)
	pred = cg_trace$pred_borzoi


	#png(sprintf("figs/model_%s_preds.png", fit_type), w=1000,h=800)
	plt1 = function(){
    layout(matrix(1:4,nrow=2))
	for(tss in c(1)) {
		#if(tss == 1) { 
		#	filt = f_tss
		#} else {
		#	filt = !f_tss
		#}
		f1 = f_train
		f2 = f_test
        f3 = f_val
		#f3 = filt & f_500 & f_ambig 
		#f4 = filt & f_500 & f_out	
		stat = lk27
		
		tr_cor = round((cor(stat[f1], pred[f1],m="p")^2),3)
		te_cor = round((cor(stat[f2], pred[f2], m="p")^2),3)
		val_cor = round((cor(stat[f3], pred[f3], m="p")^2),3)
		message("comp ranges")
		#ranges = c(min(pred),seq(ceiling(quantile(pred,0.05)), floor(quantile(pred,0.95)), 0.5), max(pred))
        ranges = c(4,5.5,6,6.5,7,7.5,8,8.5,9)
        #ranges = c(-1,0.05,.1,.15,.2,.4,.6,.7,1.1)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred[f1],ranges)), boxwex=0.15, col="blue",outline=FALSE, 
								main=sprintf("%s, tss = %s, rsqr tr, te, val = %s,%s,%s", fit_type, tss, tr_cor, te_cor,val_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col="gold", add=T,xaxt='n')
        boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.4+1:nr,outline=FALSE,
												boxwex=0.15, col="#008080", add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}
    save_baseR_to_ppt(plt1(),'./figs/bxplt_gw_borzoi.pptx')
	plt1()
	#dev.off()



}

pcg_report_locmod_borzoi_pptx2 = function(mod, fit_type='IQ')
{
#prediction - f_out, f_ambig, f_train, f_test
#vs. lk27, k4_cov_1
#scatter and box
	set.seed(42)
	gw = mod$gw
	cg_trace = gw$cg_trace
	
	f_norp = cg_trace$d_ltr !=0 & cg_trace$d_line!=0 & rowSums(is.na(cg_trace))==0 & 
    rowSums(is.na(gw$feats))==0 & rowSums(is.na(gw$feats_iqdn)) == 0 & cg_trace$start > 3e+6
	p = runif(n=nrow(gw$cg_trace))
	f_sub = f_norp & (p>0.1 | cg_trace$CG>0.02) 
#& cg_trace$cg_all!=1

	#ptrain = runif(n=nrow(gw$cg_trace))

	f_test_chr = cg_trace$chrom %in% c('chr9', 'chr18')
    f_val_chr = cg_trace$chrom %in% c('chr8', 'chr10') 
	f_train_chr =!f_test_chr & !f_val_chr

	f_test = f_test_chr & f_sub 

	f_train = f_train_chr & f_sub 
    f_val = f_val_chr & f_sub 
    lk27 = cg_trace$truth_xgb
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)
	pred = cg_trace$pred_borzoi


	#png(sprintf("figs/model_%s_preds.png", fit_type), w=1000,h=800)
	plt1 = function(){
    layout(matrix(1:4,nrow=2))
	for(tss in c(1)) {
		#if(tss == 1) { 
		#	filt = f_tss
		#} else {
		#	filt = !f_tss
		#}
		f1 = f_train
		f2 = f_test
        f3 = f_val
		#f3 = filt & f_500 & f_ambig 
		#f4 = filt & f_500 & f_out
	for(i in 1) {
		if(i == 1) {
			stat = lk27
		} else {
			stat = k4_cov_1
		}
		tr_cor = round((cor(stat[f1], pred[f1],m="p")^2),3)
		te_cor = round((cor(stat[f2], pred[f2], m="p")^2),3)
		val_cor = round((cor(stat[f3], pred[f3], m="p")^2),3)
		message("comp ranges")
		#ranges = c(min(pred),seq(ceiling(quantile(pred,0.05)), floor(quantile(pred,0.95)), 0.5), max(pred))
        #ranges = c(0,1,2,3,4,5,6,7,9)
        ranges = c(-1,0.05,.1,.15,.2,.4,.6,.7,1.1)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred[f1],ranges)), boxwex=0.15, col="blue",outline=FALSE, 
								main=sprintf("%s, tss = %s, cor = %s,%s,%s", fit_type, tss, tr_cor, te_cor,val_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col="gold", add=T,xaxt='n')
        boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.4+1:nr,outline=FALSE,
												boxwex=0.15, col="#008080", add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}}
    save_baseR_to_ppt(plt1(),'./figs/bxplt_gw_borzoi2.pptx')
	#dev.off()



}



plt_genome_pred_ppt_borzoi = function(mod, g = NA, off5 = NA, off3 = NA,mark_reg=NA,
			 chrom = NULL, locus = NULL, win_smooth = 1, show_center=T, 
			plot_pred_doms = F,
			label_tss = F,
			pred_lwd=1, plot_base_pred = F,
			 fn=NULL, fn_w = 15, fn_h=10, more_tracks=c())
{
	#if(!is.null(fn)) {
	#	if(grepl("pdf", fn)) {
	#		pdf(fn, w=fn_w, h=fn_h)
	#	} else {
	#		png(fn, w=fn_w, h=fn_h)
	#	}
	#}
	cg_trace = cg_trace
	if(is.null(chrom)) {
		f = mod$tss$geneSymbol==g
		if(sum(f) == 0) {
			return
		}
		hits = mod$tss[f,]
		locus = mean(mod$tss$start[f])
		chrom = hits[1,"chrom"]
	}
	f = as.character(cg_trace$chrom) == chrom & (cg_trace$start > locus + off5 & cg_trace$start < locus + off3)

	tss = mod$epi_tss[mod$epi_tss$chrom == chrom & (mod$epi_tss$start > locus + off5 & mod$epi_tss$start < locus + off3),]

	add_n = length(more_tracks)
	#layout(matrix(1:(2+add_n), ncol=1),h=c(rep(1, 3+add_n),1.3))
    layout(matrix(1:(2+add_n), ncol=1), heights=c(1.1, .95, 1.35))

	par(mar = c(0, 3, 2, 3))
	maxy = max(max(cg_trace$lk27_1k[f]),7)
	pred_sm = cg_trace$pred_sns_brz2k[f]
	pred_segs = which(pred_sm[-1]>4 & pred_sm[-length(pred_sm)] < 4)
	pred_segs = c(pred_segs, which(pred_sm[-1]<4 & pred_sm[-length(pred_sm)] > 4))
	pred_segs_xs = cg_trace$start[f][pred_segs]
	obs = cg_trace$lk27_1k[f]
	r2 = cor(pmax(obs,2), pmax(pred_sm,2))**2
	obs_given_pred = round(sum(obs > 4 & pred_sm > 4)/(1+sum(pred_sm > 4)),3)
	pred_given_obs = round(sum(obs > 4 & pred_sm > 4)/(1+sum(obs > 4)),3)

	plot(cg_trace$start[f], cg_trace$lk27_1k[f], 
				t="l", ylim=c(1,maxy), lwd=1,col="blue", 
				xaxt='n', 
				main=sprintf("EB4 locus r2 = %s, obs|pred %s, pred|obs %s", 
								round(r2,3), obs_given_pred, pred_given_obs), 
				xlab=sprintf("chr %s, %s", cg_trace$chrom[f][1], g))
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(cg_trace$lk27_1k[f], rep(1, length(cg_trace$lk27_1k[f]))), 
        col="blue", border=NA)
	#if(nrow(tss) > 0) {
	#	points(tss$start, pmin(1.5+tss$epi_rna-log2(1e-5),7.5), pch=19, col="red",cex=1.5)
	#}
	abline(h=log2(19), lty=2)
	if(plot_pred_doms) {
		segments(x0 = pred_segs_xs, x1=pred_segs_xs, 
			y0 = rep(0,length(pred_segs_xs)),
			y1 = rep(maxy,length(pred_segs_xs)), lty=1)
	}
	
	if(show_center) {
		abline(v=locus, col="black", lwd=1)
	}
    if(!is.na(mark_reg[1])) {
        for (mark in 1:length(mark_reg)){
		abline(v=mark_reg[mark], col="black", lwd=1)}
	}
	par(mar = c(0, 3, 1, 3))
    #par(mar = c(2,3,1, 3))
	plot(cg_trace$start[f], pred_sm, t="l", col="darkblue",xaxt='n', ylim=c(1,maxy), xlab=cg_trace$chrom[f][1],
         lwd=pred_lwd, main="prediction")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(pred_sm, rep(1, length(pred_sm))), 
        col="darkblue", border=NA)
    if(!is.na(mark_reg[1])) {
        for (mark in 1:length(mark_reg)){
		abline(v=mark_reg[mark], col="black", lwd=1)}
        }
	if(plot_base_pred) {
		pred_sm_b = zoo::rollmean(cg_trace$pred_base[f], win_smooth, f='e')
		lines(cg_trace$start[f], pred_sm_b, t="l", col="cyan", ylim=c(1,maxy), lwd=pred_lwd, xlab=cg_trace$chrom[f][1])
	}
	if(label_tss) {
		if(nrow(tss) > 0) {
			par(srt=-45)
			text(x = tss$start, y = rep(maxy-2,t=nrow(tss)), labels=tss$geneSymbol, cex=1.5)
		}
	}

	if(plot_pred_doms) {
		segments(x0 = pred_segs_xs, x1=pred_segs_xs, 
			y0 = rep(0,length(pred_segs_xs)),
			y1 = rep(maxy,length(pred_segs_xs)), lty=1)
	}
	abline(h=log2(19), lty=2)
	if(show_center) {
		abline(v=locus, col="black",lwd=2)
	}
	if(length(more_tracks) != 0) {
		profs = gextract(more_tracks, intervals=cg_trace[f,1:3], iterator=cg_trace[f,1:3])
        profs[is.na(profs)] = 0 
		for(i in 1:length(more_tracks)) {
            par(mar = c(4, 3, 1, 3))
			sm_prf = zoo::rollmean(profs[,more_tracks[i]], 5, f='e')
            #sm_prf = profs[,more_tracks[i]]
            obs = cg_trace$lk27_1k[f]
            r2 = cor(obs, sm_prf)**2
            obs_given_pred = round(sum(obs > 4 & sm_prf > 0.1)/(1+sum(sm_prf > 0.1)),3)
            pred_given_obs = round(sum(obs > 4 & sm_prf > 0.1)/(1+sum(obs > 4)),3)
			plot(cg_trace$start[f], (sm_prf), t="l", col="darkgreen", ylim=c(0,.6), lwd=pred_lwd,
                				main=sprintf("%s borzoi locus r2 = %s, obs|pred %s, pred|obs %s", 
								more_tracks[[i]],round(r2,3), obs_given_pred, pred_given_obs),)
            polygon(
                  c(cg_trace$start[f], rev(cg_trace$start[f])), 
                  c(sm_prf, rep(0, length(sm_prf))), 
                  col="darkgreen", border=NA
)
		}
	}

	#par(mar = c(1,3,1, 3))
	#plot(cg_trace$start[f], pmin(log2(10+cg_trace$atac[f]),9), t="l", col="red", ylim=c(3,9), xaxt='n', main="atac")
	#par(mar = c(5,3,1,3))
	#plot(cg_trace$start[f], pmin(cg_trace$CG[f],0.1), t="l", col="darkgreen", ylim=c(0,0.1), 
    #xlab=cg_trace$chrom[f][1], main="CG")
	#if(plot_pred_doms) {
	#	segments(x0 = pred_segs_xs, x1=pred_segs_xs, 
	#		y0 = rep(0,length(pred_segs_xs)),
	#		y1 = rep(maxy,length(pred_segs_xs)), lty=1)
	#}
	if(!is.null(fn)) {
		#dev.off()
	}
}
pcg_pred_misses = function(mod, gap_T = 4)
{

	options(gmax.data.size=1e+9)
	cgt = mod$gw$cg_trace

cg_trace <- gextract.left_join("mapab.umap_k100", intervals = cgt, iterator = cgt)



cg_trace$dist_black =  cg_trace %>%
        select(chrom, start, end) %>%
        gintervals.neighbors("mapab.Encode_blacklist_v2") %>% select(dist) %>% pull()

f_norp <- cg_trace$d_ltr != 0 &
              cg_trace$d_line != 0 &            
              cg_trace$start > 3e6 &
              !is.na(cg_trace$mapab.umap_k100) &
              cg_trace$mapab.umap_k100 > 0.9 &
               cg_trace$dist_black !=0 
	
	
	
	
f_rp =	!f_norp



	
	dlt = (cgt$lk27_1k-cgt$pred_sns_lm)
	dlt[f_rp] = 0

	gap_under200k = zoo::rollmean(!f_rp & dlt < -4, 1000, f='e')
	gap_over200k = zoo::rollmean(!f_rp & dlt > 4, 1000, f='e')

	cgt$gap_over200k = gap_over200k
	cgt$gap_under200k = gap_under200k
	
	cgt_gene = gintervals.neighbors(cgt, mod$tss)

	cgt$gene = cgt_gene$gene
	cgt$tss_dist = cgt_gene$dist

# TODO plot outlier regions

	mod$gw$cg_trace = cgt
	dlt4k = cgt$lk27_1k-cgt$pred_sns_lm
	dlt4k[f_rp] = 0
	dlt4k = zoo::rollmean(dlt4k, 10,f='e')
	
	n = length(dlt4k)
	up = dlt4k[-n] < gap_T & dlt4k[-1] >= gap_T 
	down = dlt4k[-n] >= gap_T & dlt4k[-1] < gap_T 

	over_doms = data.frame(chrom = cgt$chrom[-n][up],
							start = cgt$start[-n][up],
							end = cgt$start[-n][down])
	over_doms$l = over_doms$end-over_doms$start
	a = gintervals.neighbors(over_doms, mod$tss)
	over_doms$g = a$gene
	over_doms$g_dist = a$dist
	
	up = dlt4k[-n] > -gap_T & dlt4k[-1] <= -gap_T 
	down = dlt4k[-n] <= -gap_T & dlt4k[-1] > -gap_T 

	under_doms = data.frame(chrom = cgt$chrom[-n][up],
							start = cgt$start[-n][up],
							end = cgt$start[-n][down])
	under_doms$l = under_doms$end-under_doms$start

	a = gintervals.neighbors(under_doms, mod$tss)
	under_doms$g = a$gene
	under_doms$g_dist = a$dist

	mod$gw$over_doms = over_doms
	mod$gw$under_doms = under_doms 

	return(mod)
}
pcg_pred_misses_borzoi = function (mod, gap_T = .6) 
{

    options(gmax.data.size = 1e+09)
    cgt = mod$gw$cg_trace
	
	cg_trace <- gextract.left_join("mapab.umap_k100", intervals = cgt, iterator = cgt)



cg_trace$dist_black =  cg_trace %>%
        select(chrom, start, end) %>%
        gintervals.neighbors("mapab.Encode_blacklist_v2") %>% select(dist) %>% pull()

f_norp <- cg_trace$d_ltr != 0 &
              cg_trace$d_line != 0 &            
              cg_trace$start > 3e6 &
              !is.na(cg_trace$mapab.umap_k100) &
              cg_trace$mapab.umap_k100 > 0.9 &
               cg_trace$dist_black !=0 
	
	
	
	
f_rp =	!f_norp
    #f_rp = cgt$d_ltr == 0 | cgt$d_line == 0
    dlt = (cgt$lk27_1k - cgt$pred_borzoi)
    dlt[f_rp] = 0
    gap_under200k = zoo::rollmean(!f_rp & dlt < -4, 1000, f = "e")
    gap_over200k = zoo::rollmean(!f_rp & dlt > 4, 1000, f = "e")
    cgt$gap_over200k = gap_over200k
    cgt$gap_under200k = gap_under200k
    cgt_gene = gintervals.neighbors(cgt, mod$tss)
    cgt$gene = cgt_gene$gene
    cgt$tss_dist = cgt_gene$dist
    mod$gw$cg_trace = cgt
	dlt4k = cgt$lk27_1k - cgt$pred_borzoi
	dlt4k[f_rp] = 0
	dlt4k = zoo::rollmean(dlt4k, 10,f='e')
    n = length(dlt4k)
    up = dlt4k[-n] < gap_T & dlt4k[-1] >= gap_T
    down = dlt4k[-n] >= gap_T & dlt4k[-1] < gap_T
    over_doms = data.frame(chrom = cgt$chrom[-n][up], start = cgt$start[-n][up], 
        end = cgt$start[-n][down])
    over_doms$l = over_doms$end - over_doms$start
    a = gintervals.neighbors(over_doms, mod$tss)
    over_doms$g = a$gene
    over_doms$g_dist = a$dist
    up = dlt4k[-n] > -gap_T & dlt4k[-1] <= -gap_T
    down = dlt4k[-n] <= -gap_T & dlt4k[-1] > -gap_T
    under_doms = data.frame(chrom = cgt$chrom[-n][up], start = cgt$start[-n][up], 
        end = cgt$start[-n][down])
    under_doms$l = under_doms$end - under_doms$start
    a = gintervals.neighbors(under_doms, mod$tss)
    under_doms$g = a$gene
    under_doms$g_dist = a$dist
    mod$gw$over_doms_bor = over_doms
    mod$gw$under_doms_bor = under_doms
    return(mod)
}

plt_overdoms_in_xgb = function(){
#over_doms = (mod$gw$over_doms)
test_chroms = c('chr4','chr10','chr14','chr15')
f_test3 = mod$gw$over_doms$chrom %in% test_chroms
print(table(f_test3))
over_doms = (mod$gw$over_doms[f_test3,])
gext = gextract(c('seq.IQ.pcg.flashzoi.mm10.rf524k_EB4_cnt','EB4_cnt',
                  'jk.epipcg.pred.eb4_xgb_lm_10test_lg_seeds_gc_cg_noX_final','EB4_cnt'),
                intervals = over_doms,iterator = 20,colnames = c('pred_brz','truth_brz','pred_sns','truth_sns'))

    gext[is.na(gext)] = 0

over_doms_df = as.data.frame(tgs_matrix_tapply(t(as.matrix(gext[,c('pred_brz','truth_brz','pred_sns','truth_sns')])),gext$intervalID,mean))

over_doms_df  %>% pivot_longer(-c(3,4)) %>% ggplot(aes(x = name, y = value)) + geom_boxplot()+theme_bw()+
   coord_cartesian(ylim=c(4.5,9))
    }

plt_underdoms_in_xgb = function(){
#over_doms = (mod$gw$over_doms)
test_chroms = c('chr4','chr10','chr14','chr15')
f_test3 = mod$gw$under_doms$chrom %in% test_chroms
print(table(f_test3))
over_doms = (mod$gw$under_doms[f_test3,])
gext = gextract(c('seq.IQ.pcg.flashzoi.mm10.rf524k_EB4_cnt','EB4_cnt',
                  'jk.epipcg.pred.eb4_xgb_lm_10test_lg_seeds_gc_cg_noX_final','EB4_cnt'),
                intervals = over_doms,iterator = 20,colnames = c('pred_brz','truth_brz','pred_sns','truth_sns'))

    gext[is.na(gext)] = 0

over_doms_df = as.data.frame(tgs_matrix_tapply(t(as.matrix(gext[,c('pred_brz','truth_brz','pred_sns','truth_sns')])),gext$intervalID,mean))

over_doms_df  %>% pivot_longer(-c(3,4)) %>% ggplot(aes(x = name, y = value)) + geom_boxplot()+theme_bw()+
   coord_cartesian(ylim=c(4.5,9))
    }

plt_underdoms_in_borzoi = function(){
over_doms = (mod$gw$under_doms_bor)

gext = gextract(c('seq.IQ.pcg.flashzoi.mm10.rf524k_EB4_cnt','EB4_cnt',
                  'jk.epipcg.pred.eb4_xgb_lm_10test_lg_seeds_gc_cg_noX_final','EB4_cnt'),
                intervals = over_doms,iterator = 20,colnames = c('pred_brz','truth_brz','pred_sns','truth_sns'))

    gext[is.na(gext)] = 0

over_doms_df = as.data.frame(tgs_matrix_tapply(t(as.matrix(gext[,c('pred_brz','truth_brz','pred_sns','truth_sns')])),gext$intervalID,mean))

over_doms_df  %>% pivot_longer(-c(1,2)) %>% ggplot(aes(x = name, y = value)) + geom_boxplot()+theme_bw()+
    coord_cartesian(ylim=c(4.5,9))
}


plt_overdoms_in_borzoi = function(){
over_doms = (mod$gw$over_doms_bor)

gext = gextract(c('seq.IQ.pcg.flashzoi.mm10.rf524k_EB4_cnt','EB4_cnt',
                  'jk.epipcg.pred.eb4_xgb_lm_10test_lg_seeds_gc_cg_noX_final','EB4_cnt'),
                intervals = over_doms,iterator = 20,colnames = c('pred_brz','truth_brz','pred_sns','truth_sns'))

    gext[is.na(gext)] = 0

over_doms_df = as.data.frame(tgs_matrix_tapply(t(as.matrix(gext[,c('pred_brz','truth_brz','pred_sns','truth_sns')])),gext$intervalID,mean))

over_doms_df  %>% pivot_longer(-c(1,2)) %>% ggplot(aes(x = name, y = value)) + geom_boxplot()+theme_bw()+
    coord_cartesian(ylim=c(4.5,9))
}


###Fig4EDF

gen_regional_rna = function (mod) 
{
    #gen_track_vext()
    #gen_k27_vt(mod)
    #gen_k4_vt(mod)
    k4_tns = gen_k4_vt(mod)
    tss = mod$tss
    rownames(tss) = tss$geneSymbol
    tss = tss[tss$chrom != "chrM", ]
    eb_rna = mod$eb_legc_wt
    epi_rna = mod$emb_type_legc[, "Epiblast"]
    amn_rna = mod$emb_type_legc[, "Amnion/Chorion"]
    tss$eb_rna = rep(log2(1e-05), nrow(tss))
    tss$epi_rna = rep(log2(1e-05), nrow(tss))
    tss$amn_rna = rep(log2(1e-05), nrow(tss))
    gnms = intersect(rownames(tss), names(epi_rna))
    tss[gnms, "eb_rna"] = eb_rna[gnms]
    tss[gnms, "epi_rna"] = epi_rna[gnms]
    tss[gnms, "amn_rna"] = amn_rna[gnms]
    tss_ext = tss
    tss_ext$start = tss_ext$start - 1000
    tss_ext$end = tss_ext$end + 1000
    #prof = gextract(c("EB4_cnt", "EB4_cnt_k4", "ES_ch_k4", "atac_ext", 
    #    "seq.CG_500_mean_new", "seq.GC500_bin20"), iterator = 20, 
    #    intervals = tss_ext)
    #prof[is.na(prof)] = 0
    #T_k27 = mod$k27_track_thresh["EB4_cnt", "X0.99"]
    #T_k4 = mod$k4_track_thresh["EB4_cnt", "X0.99"]
    #T_atac = 497
    #tss$n2k_eb4_k27 = 20 * tapply(prof[, "EB4_cnt"] > T_k27, 
    #    prof$intervalID, sum)
    #tss$n2k_eb4_k4 = 20 * tapply(prof[, "EB4_cnt_k4"] > T_k4, 
    #    prof$intervalID, sum)
    #tss$eb4_k27_max = tapply(prof[, "EB4_cnt"], prof$intervalID, 
    #    max)
    #tss$eb4_k4_max = tapply(prof[, "EB4_cnt_k4"], prof$intervalID, 
    #    max)
    #tss$n2k_atac = 20 * tapply(prof[, "atac_ext"] > T_atac, prof$intervalID, 
    #    sum)
    #tss$atac_max = tapply(prof[, "atac_ext"], prof$intervalID, 
    #    max)
    #tss$cg_max = tapply(prof[, "seq.CG_500_mean_new"], prof$intervalID, 
    #    max)
    #tss$gc_max = tapply(prof[, "seq.GC500_bin20"], prof$intervalID, 
    #    max)
    #tss$gc_mean = tapply(prof[, "seq.GC500_bin20"], prof$intervalID, 
    #    mean)
    mod$epi_tss_high_cg_int = tss
    doms = mod$cgdom_ann[, c("chrom", "start", "end")]
    doms$id = 1:nrow(doms)
    dom_tss = gintervals.neighbors(doms, tss[, c("chrom", "start", 
        "end", "strand", "epi_rna", "eb_rna")], maxneighbors = 100, 
        maxdist = 1e+05, mindist = -1e+05)
    dom_region_rna = data.frame(id = doms$id)
    dom_region_rna = cbind(dom_region_rna, matrix(NA, nrow = length(doms$id), 
        ncol = 6))
    rownames(dom_region_rna) = doms$id
    ci = 2
    for (horiz in c(1000, 2000, 4000, 10000, 50000, 1e+05)) {
        dh = dom_tss[abs(dom_tss$dist) <= horiz, ]
        gmax = tapply(dh$eb_rna, dh$id, max)
        dom_region_rna[names(gmax), ci] = gmax
        ci = ci + 1
    }
    colnames(dom_region_rna) = c("ID", "rna1k", "rna2k", "rna4k", 
        "rna10k", "rna50k", "rna100k")
    mod$cg_high_int_rna = dom_region_rna
    mod$cgdom_ann_rna = merge(mod$cg_high_int_rna, mod$cgdom_ann, 
        by = "ID")
    return(mod)
}
pcg_test_miss_over_chroms_fcs_no_pred = function (mod, fn = NA,T_pos_t = 6.48, T_pos_p = 6.22, k_reg = 0.002, pred_max = 0.2, 
    miss_smooth_win = 40, smooth_win = 8, focus_chrom = NULL) 
{
    cg_trace = mod$gw$cg_trace
    chrom_y = 1:19
    names(chrom_y) = c("chr1", "chr2", "chr3", "chr4", "chr5", 
        "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
        "chr19")
    a = cg_trace %>% mutate(chrom = chrom_y[as.character(chrom)]) %>% 
        mutate(sbin = floor(cg_trace$start/50000)) %>% group_by(chrom, 
        sbin) %>% summarize(tpred = mean(pred_sns_lm > T_pos_p), tobs = mean(lk27_1k > 
        T_pos_t)) #T_pos_t T_pos_p .98 quantile
    a = a %>% mutate(lf = log2((k_reg + tobs)/(k_reg + tpred)))
    
    # Filter for specific chromosome if requested
    if (!is.null(focus_chrom)) {
        if (focus_chrom %in% names(chrom_y)) {
            focus_chrom_num = chrom_y[focus_chrom]
            a = a %>% filter(chrom == focus_chrom_num)
            cat("Focusing on", focus_chrom, "(chromosome", focus_chrom_num, ")\n")
        } else {
            stop("Invalid chromosome name. Use: ", paste(names(chrom_y), collapse = ", "))
        }
    }
    apred = a %>% reshape2::dcast(chrom ~ sbin, value.var = "tpred")
    apred = apred[, -1]
    apreds = t(apply(apred, 1, function(x) zoo::rollmean(x, k = smooth_win, 
        na.rm = T, f = "e")))
    aobs = a %>% reshape2::dcast(chrom ~ sbin, value.var = "tobs")
    aobs = aobs[, -1]
    aobss = t(apply(aobs, 1, function(x) zoo::rollmean(x, k = smooth_win, 
        na.rm = T, f = "e")))
    am = log2((k_reg + aobss)/(k_reg + apreds))
    ams = t(apply(am, 1, function(x) zoo::rollmean(x, k = miss_smooth_win, 
        na.rm = T, f = "e")))
    
    # Color palette for log-fold change only
    shades = colorRampPalette(c("cyan", "green", "lightgreen", 
        "white", "white", "white", "white", "gray", "pink", "darkred", 
        "black"))(1000)
    
    nbins = ncol(apreds)
    nrows = ifelse(!is.null(focus_chrom), 1, 19)  # Single row if focusing on one chromosome
    
    # Create nested matrix based on chromosome focus
    if (!is.null(focus_chrom)) {
        # Only log-fold change row for specific chromosome
        nested = matrix(as.vector(t(pmin(t(ams), 3))), ncol = 1)
    } else {
        # Log-fold change with NA spacers for all chromosomes
        nested = matrix(as.vector(rbind(pmin(t(ams), 3), 
            matrix(rep(NA, nbins * nrows), nrow = nbins, ncol = nrows))), 
            ncol = nrows * 2)
    }
    
    # Create x-axis labels for genomic coordinates (every 10 million base pairs)
    # Each bin represents 50kb, so 200 bins = 10Mb
    x_labels = rep("", nbins)
    x_at = seq(1, nbins, by = 200)  # Every 200th bin (10Mb intervals)
    x_labels[x_at] = paste0(((x_at - 1) * 50 / 1000), "M")  # Convert to Mb coordinates
    
    # Create plot title
    plot_title = ifelse(!is.null(focus_chrom), paste("Chromosome", focus_chrom), "All Chromosomes")
    
    if (!is.na(fn)) {
        # Main heatmap
        plot_height = ifelse(!is.null(focus_chrom), 4, 14)
        pheatmap::pheatmap(t(nested), col = shades, 
            breaks = seq(-3, 2.99, l = 1001), 
            cluster_rows = F, cluster_cols = F, 
            filename = fn, width = 18, height = plot_height,
                           
            #na_col = rgb(0, 0, 0, 0), 
                           legend = FALSE, main = plot_title,
            labels_col = x_labels)
        
        # Legend: Log-fold change only
        legend_data = matrix(seq(-3, 3, length.out = 100), ncol = 1)
        pheatmap::pheatmap(legend_data, col = shades, 
            breaks = seq(-3, 2.99, l = 1001), 
            cluster_rows = F, cluster_cols = F,
            filename = paste0(gsub("\\.[^.]*$", "", fn), "_legend_logfold.png"), 
            width = 2, height = 6,
            main = "Log2(Observed/Predicted)")
    }
    else {
        # For interactive display
        library(gridExtra)
        library(grid)
        
        # Main heatmap
        p_main = pheatmap::pheatmap(t(nested), col = shades, 
            breaks = seq(-3, 2.99, l = 1001), 
            cluster_rows = F, cluster_cols = F,
            legend = FALSE, silent = TRUE, main = plot_title,
            labels_col = x_labels)
        
        # Legend
        legend_data = matrix(seq(-3, 3, length.out = 100), ncol = 1)
        p_legend = pheatmap::pheatmap(legend_data, col = shades, 
            breaks = seq(-3, 2.99, l = 1001), 
            cluster_rows = F, cluster_cols = F, silent = TRUE,
            main = "Log2(Obs/Pred)")
        
        # Arrange plots
        grid.arrange(p_main[[4]], p_legend[[4]], 
                    ncol = 2, widths = c(6, 1))
    }
    
    return(t(nested))
}
###Fig3


pcg_test_miss_over_chroms_fcs_no_pred_hg = function (mod, fn = NA, T_pos = 4.5, k_reg = 0.002, pred_max = 0.2, 
    miss_smooth_win = 40, smooth_win = 8, focus_chrom = NULL) 
{
    cg_trace = mod$gw$cg_trace
    chrom_y = 1:23
    names(chrom_y) = c("chr1", "chr2", "chr3", "chr4", "chr5", 
        "chr6", "chr7", "chr8", "chr9", "chr10", "chr11", "chr12", 
        "chr13", "chr14", "chr15", "chr16", "chr17", "chr18", 
        "chr19", "chr20", "chr21", "chr22", "chrX")
    a = cg_trace %>% mutate(chrom = chrom_y[as.character(chrom)]) %>% 
        mutate(sbin = floor(cg_trace$start/50000)) %>% group_by(chrom, 
        sbin) %>% summarize(tpred = mean(pred > T_pos), tobs = mean(lk27_1k > 
        T_pos))
    a = a %>% mutate(lf = log2((k_reg + tobs)/(k_reg + tpred)))
    
    # Filter for specific chromosome if requested
    if (!is.null(focus_chrom)) {
        if (focus_chrom %in% names(chrom_y)) {
            focus_chrom_num = chrom_y[focus_chrom]
            a = a %>% filter(chrom == focus_chrom_num)
            cat("Focusing on", focus_chrom, "(chromosome", focus_chrom_num, ")\n")
        } else {
            stop("Invalid chromosome name. Use: ", paste(names(chrom_y), collapse = ", "))
        }
    }
    apred = a %>% reshape2::dcast(chrom ~ sbin, value.var = "tpred")
    apred = apred[, -1]
    apreds = t(apply(apred, 1, function(x) zoo::rollmean(x, k = smooth_win, 
        na.rm = T, f = "e")))
    aobs = a %>% reshape2::dcast(chrom ~ sbin, value.var = "tobs")
    aobs = aobs[, -1]
    aobss = t(apply(aobs, 1, function(x) zoo::rollmean(x, k = smooth_win, 
        na.rm = T, f = "e")))
    am = log2((k_reg + aobss)/(k_reg + apreds))
    ams = t(apply(am, 1, function(x) zoo::rollmean(x, k = miss_smooth_win, 
        na.rm = T, f = "e")))
    
    # Color palette for log-fold change only
    shades = colorRampPalette(c("cyan", "green", "lightgreen", 
        "white", "white", "white", "white", "gray", "pink", "darkred", 
        "black"))(1000)
    
    nbins = ncol(apreds)
    nrows = ifelse(!is.null(focus_chrom), 1, 23)  # Single row if focusing on one chromosome
    
    # Create nested matrix based on chromosome focus
    if (!is.null(focus_chrom)) {
        # Only log-fold change row for specific chromosome
        nested = matrix(as.vector(t(pmin(t(ams), 3))), ncol = 1)
    } else {
        # Log-fold change with NA spacers for all chromosomes
        nested = matrix(as.vector(rbind(pmin(t(ams), 3), 
            matrix(rep(NA, nbins * nrows), nrow = nbins, ncol = nrows))), 
            ncol = nrows * 2)
    }
    
    # Create x-axis labels for genomic coordinates (every 10 million base pairs)
    # Each bin represents 50kb, so 200 bins = 10Mb
    x_labels = rep("", nbins)
    x_at = seq(1, nbins, by = 200)  # Every 200th bin (10Mb intervals)
    x_labels[x_at] = paste0(((x_at - 1) * 50 / 1000), "M")  # Convert to Mb coordinates
    
    # Create plot title
    plot_title = ifelse(!is.null(focus_chrom), paste("Chromosome", focus_chrom), "All Chromosomes")
    
    if (!is.na(fn)) {
        # Main heatmap
        plot_height = ifelse(!is.null(focus_chrom), 4, 14)
        pheatmap::pheatmap(t(nested), col = shades, 
            breaks = seq(-3, 2.99, l = 1001), 
            cluster_rows = F, cluster_cols = F, 
            filename = fn, width = 18, height = plot_height, 
            #na_col = rgb(0, 0, 0, 0),
                           legend = FALSE, main = plot_title,
            labels_col = x_labels)
        
        # Legend: Log-fold change only
        legend_data = matrix(seq(-3, 3, length.out = 100), ncol = 1)
        pheatmap::pheatmap(legend_data, col = shades, 
            breaks = seq(-3, 2.99, l = 1001), 
            cluster_rows = F, cluster_cols = F,
            filename = paste0(gsub("\\.[^.]*$", "", fn), "_legend_logfold.png"), 
            width = 2, height = 6,
            main = "Log2(Observed/Predicted)")
    }
    else {
        # For interactive display
        library(gridExtra)
        library(grid)
        
        # Main heatmap
        p_main = pheatmap::pheatmap(t(nested), col = shades, 
            breaks = seq(-3, 2.99, l = 1001), 
            cluster_rows = F, cluster_cols = F,
            legend = FALSE, silent = TRUE, main = plot_title,
            labels_col = x_labels)
        
        # Legend
        legend_data = matrix(seq(-3, 3, length.out = 100), ncol = 1)
        p_legend = pheatmap::pheatmap(legend_data, col = shades, 
            breaks = seq(-3, 2.99, l = 1001), 
            cluster_rows = F, cluster_cols = F, silent = TRUE,
            main = "Log2(Obs/Pred)")
        
        # Arrange plots
        grid.arrange(p_main[[4]], p_legend[[4]], 
                    ncol = 2, widths = c(6, 1))
    }
    
    return(t(nested))
}


plt_genome_pred_CG_ppt_st_end = function (mod, g = NA, off5 = NA, off3 = NA, mark_reg = NA, chrom = NULL, 
    locus = NULL,start,end, win_smooth = 1, show_center = T, plot_pred_doms = F, 
    label_tss = F, pred_lwd = 1, plot_base_pred = F, fn = NULL, 
    fn_w = 15, fn_h = 10, more_tracks = c()) 
{
    cg_trace = mod$gw$cg_trace
    if (is.null(chrom)) {
        f = mod$tss$geneSymbol == g
        if (sum(f) == 0) {
            return
        }
        hits = mod$tss[f, ]
        locus = mean(mod$tss$start[f])
        chrom = hits[1, "chrom"]
    }
    f = as.character(cg_trace$chrom) == chrom & (cg_trace$start > 
        start + off5 & cg_trace$start < end + off3)
    tss = mod$epi_tss[mod$epi_tss$chrom == chrom & (mod$epi_tss$start > 
        start + off5 & mod$epi_tss$start < end + off3), ]
    add_n = length(more_tracks)
    layout(matrix(1:(2 + add_n), ncol = 1), h = c(1,1, 1.2))
    par(mar = c(0, 3, 2, 3))
    maxy = max(max(cg_trace$lk27_1k[f]), 7)
    pred_sm = zoo::rollmean(cg_trace$pred[f], win_smooth, f = "e")
    pred_segs = which(pred_sm[-1] > 4 & pred_sm[-length(pred_sm)] < 
        4)
    pred_segs = c(pred_segs, which(pred_sm[-1] < 4 & pred_sm[-length(pred_sm)] > 
        4))
    pred_segs_xs = cg_trace$start[f][pred_segs]
    obs = cg_trace$lk27_1k[f]
    r2 = cor(pmax(obs, 2), pmax(pred_sm, 2))^2
    obs_given_pred = round(sum(obs > 4 & pred_sm > 4)/(1 + sum(pred_sm > 
        4)), 3)
    pred_given_obs = round(sum(obs > 4 & pred_sm > 4)/(1 + sum(obs > 
        4)), 3)
    plot(cg_trace$start[f], cg_trace$lk27_1k[f], t = "l", ylim = c(1, 
        maxy), lwd = 1, col = "blue", xaxt = "n", main = sprintf("EB4 locus r2 = %s, obs|pred %s, pred|obs %s", 
        round(r2, 3), obs_given_pred, pred_given_obs), xlab = sprintf("chr %s, %s", 
        cg_trace$chrom[f][1], g))
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(cg_trace$lk27_1k[f], 
        rep(1, length(cg_trace$lk27_1k[f]))), col = "blue", border = NA)
    abline(h = log2(19), lty = 2)
    if (plot_pred_doms) {
        segments(x0 = pred_segs_xs, x1 = pred_segs_xs, y0 = rep(0, 
            length(pred_segs_xs)), y1 = rep(maxy, length(pred_segs_xs)), 
            lty = 1)
    }
    if (show_center) {
        abline(v = locus, col = "black", lwd = 1)
    }
    if (!is.na(mark_reg[1])) {
        for (mark in 1:length(mark_reg)) {
            abline(v = mark_reg[mark], col = "black", lwd = 1)
        }
    }
    par(mar = c(0, 3, 2, 3))
    plot(cg_trace$start[f], pred_sm, t = "l", col = "darkblue", 
        ylim = c(1, maxy), xaxt = "n", lwd = pred_lwd, main = "prediction")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])), c(pred_sm, 
        rep(1, length(pred_sm))), col = "darkblue", border = NA)
    if (!is.na(mark_reg[1])) {
        for (mark in 1:length(mark_reg)) {
            abline(v = mark_reg[mark], col = "black", lwd = 1)
        }
    }
    if (plot_base_pred) {
        pred_sm_b = zoo::rollmean(cg_trace$pred_base[f], win_smooth, 
            f = "e")
        lines(cg_trace$start[f], pred_sm_b, t = "l", col = "cyan", 
            ylim = c(1, maxy), lwd = pred_lwd, xlab = cg_trace$chrom[f][1])
    }
    if (label_tss) {
        if (nrow(tss) > 0) {
            par(srt = -45)
            text(x = tss$start, y = rep(maxy - 2, t = nrow(tss)), 
                labels = tss$geneSymbol, cex = 1.5)
        }
    }
    if (plot_pred_doms) {
        segments(x0 = pred_segs_xs, x1 = pred_segs_xs, y0 = rep(0, 
            length(pred_segs_xs)), y1 = rep(maxy, length(pred_segs_xs)), 
            lty = 1)
    }
    abline(h = log2(19), lty = 2)
    if (show_center) {
        abline(v = locus, col = "black", lwd = 2)
    }
    if (length(more_tracks) != 0) {
        par(mar = c(4, 3, 2, 3))
        profs = gextract(more_tracks, intervals = cg_trace[f, 
            1:3], iterator = cg_trace[f, 1:3])
        for (i in 1:length(more_tracks)) {
            sm_prf = profs[, more_tracks[i]]
            plot(cg_trace$start[f], sm_prf, t = "l", col = "black", 
                ylim = c(0, 0.1), lwd = pred_lwd, main = more_tracks[[i]])
        }
    }
    if (!is.null(fn)) {
    }
}

plt_genome_pred_CG_ppt_st_end_borzoi = function (
    mod, g = NA, off5 = NA, off3 = NA, mark_reg = NA,
    chrom = NULL, locus = NULL, start, end,
    win_smooth = 1, show_center = TRUE, plot_pred_doms = FALSE,
    label_tss = FALSE, pred_lwd = 1, plot_base_pred = FALSE,
    fn = NULL, fn_w = 15, fn_h = 10, more_tracks = c()
) {

    cg_trace = mod$gw$cg_trace

    # Identify locus / chromosome
    if (is.null(chrom)) {
        f = mod$tss$geneSymbol == g
        if (sum(f) == 0) return()
        hits = mod$tss[f, ]
        locus = mean(mod$tss$start[f])
        chrom = hits[1, "chrom"]
    }

    # Window filter
    if (is.null(locus)){
    f = as.character(cg_trace$chrom) == chrom &
        (cg_trace$start > start + off5 & cg_trace$start < end + off3)
        tss = mod$epi_tss[
        mod$epi_tss$chrom == chrom &
        (mod$epi_tss$start > start + off5 & mod$epi_tss$start < end + off3),
    ]
    }else{
      f = as.character(cg_trace$chrom) == chrom & (cg_trace$start > locus + off5 & cg_trace$start < locus + off3) 
     tss = mod$epi_tss[mod$epi_tss$chrom == chrom & (mod$epi_tss$start > locus + off5 & mod$epi_tss$start < locus + off3),]   
    }


    add_n = length(more_tracks)
    n_panels = 5 + add_n

    # --- TIGHT LAYOUT, Y LABELS OK ---
    layout(
        matrix(seq_len(n_panels), ncol = 1),
        heights = rep(1, n_panels)
    )

    # PRECOMPUTE
    maxy = max(max(cg_trace$lk27_1k[f]), 7)
    pred_sm = cg_trace$pred_sns_brz2k[f]
	pred_sm_lm = cg_trace$pred_sns_lm[f]
    pred_sm_borzoi = cg_trace$pred_borzoi[f]
	pred_sm_silicus = cg_trace$pred_sil[f]
    pred_segs = which(pred_sm[-1] > 4 & pred_sm[-length(pred_sm)] < 4)
    pred_segs = c(pred_segs,
                  which(pred_sm[-1] < 4 & pred_sm[-length(pred_sm)] > 4))
    pred_segs_xs = cg_trace$start[f][pred_segs]

    ### -------------------------
    ### PANEL 1: OBSERVED
    ### -------------------------
    par(mar=c(1,4,0,2), mgp=c(2,0.6,0))
    plot(cg_trace$start[f], cg_trace$lk27_1k[f],
         type="l", ylim=c(4.75,maxy), xaxt="n",
         lwd=1, col="blue", xlab="", ylab="obs")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])),
            c(cg_trace$lk27_1k[f], rep(1, sum(f))),
            col="blue", border=NA)
	abline(h=6.49, lty=2)#.98 quant
    if (plot_pred_doms)
        segments(pred_segs_xs, 0, pred_segs_xs, maxy)

    if (show_center) abline(v=locus)
    if (!is.na(mark_reg[1]))
        for (mark in mark_reg) abline(v=mark)

    ### -------------------------
    ### PANEL 2: MODEL PRED
    ### -------------------------
    par(mar=c(1,4,0,2))
    plot(cg_trace$start[f], pred_sm,
         type="l", ylim=c(4.8,maxy), xaxt="n",
         lwd=pred_lwd, col="darkblue", xlab="", ylab="pred_sns_brz2k")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])),
            c(pred_sm, rep(1, sum(f))),
            col="darkblue", border=NA)
	abline(h=6.26, lty=2)#.98 quant
    if (!is.na(mark_reg[1]))
        for (mark in mark_reg) abline(v=mark)
    ### -------------------------
    ### PANEL 3: MODEL PRED lm sns
    ### -------------------------
    par(mar=c(1,4,0,2))
    plot(cg_trace$start[f], pred_sm_lm,
         type="l", ylim=c(4.83,maxy), xaxt="n",
         lwd=pred_lwd, col="darkblue", xlab="", ylab="pred_sn_lm")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])),
            c(pred_sm_lm, rep(1, sum(f))),
            col="darkblue", border=NA)
	abline(h=6.22, lty=2)#.98 quant
    if (!is.na(mark_reg[1]))
        for (mark in mark_reg) abline(v=mark)
    ### -------------------------
    ### PANEL 3: BORZOI
    ### -------------------------
    par(mar=c(1,4,0,2))
    plot(cg_trace$start[f], pred_sm_borzoi,
         type="l", ylim=c(5,maxy), xaxt="n",
         lwd=pred_lwd, col="darkgreen", xlab="", ylab="borzoi")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])),
            c(pred_sm_borzoi, rep(1, sum(f))),
            col="darkgreen", border=NA)
	abline(h=6.5, lty=2)#.98 quant
    if (plot_pred_doms)
        segments(pred_segs_xs, 0, pred_segs_xs, maxy)
    if (show_center) abline(v=locus)
    ### -------------------------
    ### PANEL 3: SILICUS
    ### -------------------------
    par(mar=c(1,4,0,2))
    plot(cg_trace$start[f], pred_sm_silicus,
         type="l", ylim=c(5.1,maxy), xaxt="n",
         lwd=pred_lwd, col="darkviolet", xlab="", ylab="silicus")
    polygon(c(cg_trace$start[f], rev(cg_trace$start[f])),
            c(pred_sm_silicus, rep(1, sum(f))),
            col="darkviolet", border=NA)
	abline(h=6.27, lty=2)#.98 quant
    ### -------------------------
    ### PANELS 4+: EXTRA TRACKS
    ### -------------------------
    if (add_n > 0) {
        profs = gextract(more_tracks,
                         intervals = cg_trace[f,1:3],
                         iterator  = cg_trace[f,1:3])

        for (i in seq_along(more_tracks)) {
            # Last panel should have x-axis shown
            is_last = (i == add_n)
            par(mar=c(ifelse(is_last, 2, 1), 4, 0, 2))

            plot(cg_trace$start[f], profs[, more_tracks[i]],
                 type="l",
				 ylim=c(0,0.1),
                 xaxt=ifelse(is_last, "s", "n"),
                 lwd=pred_lwd, col="black",
                 xlab=ifelse(is_last, "position", ""),
                 ylab=more_tracks[i])
        }
    }

    if (!is.null(fn)) {
        # file saving logic here
    }
}




plt_regionAT_ds_cg_cons = function (mod, chr, st, end,k, ds_k27 = TRUE, marg5 = 2000, marg3 = 2000, 
    mark = c(), wk4 = F, k4_max_fact = 1) 
{
    tracks_k27 = c("jk.epipcg.pcg.CRJK_0363_k27me3_eb_j1_d3_a_norm", 
        "jk.epipcg.pcg.CRJK_0365_k27me3_eb_dko_d3_a_norm")

    short_namesk27 = c('wt_k27','dko_k27')
    for (i in 1:length(tracks_k27)) {
        gvtrack.create(short_namesk27[i], tracks_k27[i], "avg")
        #gvtrack.iterator(short_namesk27[i], sshift = -140, eshift = 140)
    }
    
    prof = gextract(short_namesk27, intervals = gintervals(chr, 
        st - marg5, end + marg3), iterator = 20)
    prof[is.na(prof)] = 0
ylim <- ( max(apply(prof[, short_namesk27], 2, function(x) {
  max(zoo::rollmean(x, k = k, fill = 0))
})))
	
ndx = data.frame(rep(1,length(short_namesk27)))
    #gvtrack.create("WT_Epi_rep1", "jk.epipcg.meth.meissN20.WT_Epi_rep1", 
    #    "avg")
    #gvtrack.iterator("WT_Epi_rep1", sshift = -100, eshift = 100)
    #c = (c("seq.CG_500_mean_new", NA, NA, NA, NA, NA, "seq.CG_500_mean_new", 
    #    NA, 1, 1, 1, 1, 1, 1))
    #ndx = rbind(ndx, c)
    #prof2 = gextract(c("WT_Epi_rep1", "seq.CG_500_mean_new"), 
      #  intervals = gintervals(chr, st - marg5, end + marg3), 
      #  iterator = 20, colnames = c("meth", "cg"))
    #prof2[is.na(prof2)] = 0
    #prof$cg = prof2$cg
    #prof$intervalID = NULL
    layout(matrix(1:nrow(ndx), ncol = 1), h = c(1.4, rep(1, nrow(ndx) - 
        2), 1.4))
    for (i in 1:nrow(ndx)) {
        if (i == nrow(ndx)) {
            par(mar = c(2, 4, 0.8, 4))
        }
        else if (i == 1) {
            par(mar = c(0, 4, 2, 4))
        }
        else {
            par(mar = c(0, 4, 0.8, 4))
        }
        if (i < nrow(ndx)) {
            c((prof[, 3 + i]), rep(0, length(prof$start)))
            plot(prof$start, (zoo::rollmean(prof[, 3 + i],k = k,fill = 0)), pch = 19, type = "l", 
                lwd = 3, ylim = c(4.78, ylim), xaxt = ifelse(i == 
                  nrow(ndx), "s", "n"), main = short_namesk27[i], col = "darkblue", 
                ylab = NA)
            polygon(c(prof$start, rev(prof$start)), c((zoo::rollmean(prof[, 
                3 + i],k = k,fill = 0)), rep(0, length(prof$start))), col = "darkblue", 
                border = NA)
        }
        else {
            plot(prof$start, (zoo::rollmean(prof[, 3 + i],k = k,fill = 0)), pch = 19, type = "l", 
                lwd = 3, ylim = c(4.69, ylim), xaxt = ifelse(i == 
                  nrow(ndx), "s", "n"), main = short_namesk27[i], col = "turquoise", 
                ylab = NA)
                polygon(c(prof$start, rev(prof$start)), c((zoo::rollmean(prof[, 
                3 + i],k = k,fill = 0)), rep(0, length(prof$start))), col = "turquoise", 
                border = NA)
        }
        abline(v = st)
        abline(v = end)
        if (!is.null(mark) & length(mark) > 0) {
            for (x in mark) {
                abline(v = x, col = "blue", lwd = 2)
            }
        }
    }
}
ebpcg_diff_doms = function(mod)
{
	gen_k27_vt(mod)
	a_doms = list()
	for(sname in rownames(mod$k27_track_thresh)) {
		message("screen ", sname)
		thresh = mod$k27_track_thresh[sname,"X0.99"]
		d = gscreen(sprintf("%s > %s", sname, thresh))
		d$l = d$end - d$start
		a_doms[[sname]] = d
	}
	a_doms_e = list()
	for(sname in names(a_doms)) {
		doms = a_doms[[sname]]
		doms = doms[!doms$chrom %in% c("chrY","chrM"),]
		doms$start = doms$start - 1000
		doms$end = doms$end + 1000
		doms_e = gintervals.canonic(doms)
		doms_e$start = doms_e$start + 1000
		doms_e$end = doms_e$end - 1000
		doms_e$l = doms_e$end - doms_e$start
		a_doms_e[[sname]] = doms_e
	}
	eb_seeds = a_doms_e[["EB4_cnt"]]
	eb_seeds = eb_seeds[eb_seeds$l > 300,]
	colnames(eb_seeds)[4] = "base_l"
	
	for(sname in names(a_doms)) {
		doms_e = a_doms_e[[sname]]
		doms_e$ID = 1:nrow(doms_e)
#first find longes overlapping base seed
		seed_d0 = gintervals.neighbors(doms_e, eb_seeds, maxdist=0, maxneighbors=1000)

		seed_d0 = seed_d0[order(-seed_d0$base_l),]
		seed_d0 = seed_d0[!duplicated(seed_d0$ID),]
		seed_d = gintervals.neighbors(doms_e[!doms_e$ID %in% seed_d0$ID,], eb_seeds)
		seed_d = rbind(seed_d0, seed_d)
		seed_d = seed_d[order(seed_d$ID),]
		doms_e$d_seed = seed_d$dist
		doms_e$seed_size = seed_d$base_l
		a_doms_e[[sname]] = doms_e
	}
	D = 5000
#	layout(matrix(c(1:14,-1),nrow=3))
	layout(matrix(c(1:2),nrow=1))
	par(mar=c(2,2,2,2))
	x_color = NA
	for(sname in c("e75_ecto_cnt","e75_emeso_cnt")) {
	#for(sname in names(a_doms)) {
		d = a_doms_e[[sname]]
		f_long = d$l > D
		plot(y=pmin(log2(d$l[f_long]),17), x=pmin(log2(256+d$seed_size[f_long]),17), 
						main=sname, ylim=c(log2(D),17), xlim=c(8,17), pch=19, 
						col=ifelse(d$chrom[f_long]=="chrX", x_color, ifelse(d$d_seed[f_long]>0,"red","black")))
		abline(a=0,b=1)
	}

	mod$a_doms_e = a_doms_e

	return(mod)
}

test_overlap_graph = function(mod, nm1, nm2, l_thresh = 20e+3)
{
#	d_s1 = mod$a_doms_e[["EB4_cnt"]]
#	d_s2 = mod$a_doms_e[["e75_s2_cnt"]]
	d_s1 = mod$a_doms_e[[nm1]]
	d_s2 = mod$a_doms_e[[nm2]]

	d_s2_long = d_s2[d_s2$l > l_thresh,]
	d_s2_long = d_s2_long[order(d_s2_long$l),]
	d_s2_long$s2_id = 1:nrow(d_s2_long)
	d_s1_s2_long = gintervals.neighbors(d_s1, d_s2_long)

	d_s2$strand = 1
	d_s2_s2_pos = gintervals.neighbors(d_s2, d_s2, mindist = 1)
	d_s2_s2_neg = gintervals.neighbors(d_s2, d_s2, maxdist = -1)
	
	s2_marg_3 = d_s2_s2_pos$dist/2
	names(s2_marg_3) = d_s2_s2_pos$s2_id
	s2_marg_5 = -d_s2_s2_neg$dist/2
	names(s2_marg_5) = d_s2_s2_neg$s2_id

	d_s1_s2_long_in = d_s1_s2_long[d_s1_s2_long$dist == 0,]

	d_s1_s2_long_in$mid = (d_s1_s2_long_in$start1 + d_s1_s2_long_in$end1)/2
	over_count = table(d_s1_s2_long_in$s2_id)
	d_s1_s2_long_in$s2_dup = over_count[d_s1_s2_long_in$s2_id]

	a = d_s1_s2_long_in

	a$start_clip = pmax(a$start, a$start1 - s2_marg_5[a$s2_id])
	a$end_clip = pmin(a$end, a$end1 + s2_marg_3[a$s2_id])

	n_tot = nrow(d_s2_long)
	n_miss = n_tot-length(unique(a$s2_id))
	f_exp1 = a$start1 - a$start_clip > 5000
	f_exp2 = a$end_clip - a$end1 > 5000
	tot_expand = sum(f_exp1 | f_exp2)

	ids = sort(unique(a$s2_id))
	compact_ids = 1:length(ids)
	names(compact_ids) = ids
	a$s2_idc = compact_ids[as.character(a$s2_id)]

	plot(c(a$start1 - a$mid, a$end1 - a$mid), c(a$s2_idc, a$s2_idc), 
			pch=19, cex=0.4, 
			main = sprintf("ord %s dom > %s, exp in %s = %s ", nm2, l_thresh, nm2, tot_expand),
			xlab = sprintf("N = %s, miss = %s", n_tot, n_miss))
	segments(x0=a$start1 - a$mid, x1 = a$end1 - a$mid, y0 = a$s2_idc, y1 = a$s2_idc, lwd=1, col="darkviolet")
#	segments(x0=a$start - a$mid, x1 = a$end - a$mid, y0 = a$s2_idc, y1 = a$s2_idc, lwd=0.5, col="blue")
	segments(x0=a$start_clip - a$mid, x1 = a$end_clip - a$mid, y0 = a$s2_idc, y1 = a$s2_idc, lwd=.8, col="darkgreen")
}


gen_diff_mat = function(mod){
cgd = mod$cgdom_ann
	all_seed_k27 = gextract(mod$epi_tracks$short_name, mod$cgdom_ann, iterator=mod$cgdom_ann)
	len = all_seed_k27$end - all_seed_k27$start
	len_norm = 300/len

	eb_tname = mod$epi_tracks$track_k27[mod$epi_tracks$short_name == "EB4_cnt"]
	emeso_tname = mod$epi_tracks$track_k27[mod$epi_tracks$short_name == 'e75_e_meso_cnt']	
	ecto_tname = mod$epi_tracks$track_k27[mod$epi_tracks$short_name == "e75_ecto_cnt"]	
	cov_eb = gsummary(eb_tname)[5]
	cov_emeso = gsummary(emeso_tname)[5]
	cov_ecto = gsummary(ecto_tname)[5]

	norm_emeso  = cov_eb / cov_emeso
	norm_ecto = cov_eb / cov_ecto
	all_seed_k27$e75_emeso_cnt = all_seed_k27$e75_e_meso_cnt * norm_emeso
	all_seed_k27$e75_ecto_cnt = all_seed_k27$e75_ecto_cnt * norm_ecto

	mat_e75= cbind( log2(4+all_seed_k27$EB4_cnt * len_norm), 
						 log2(4+all_seed_k27$e75_e_meso_cnt * len_norm), 
						 log2(4+all_seed_k27$e75_ecto_cnt * len_norm))
						 
	colnames(mat_e75) = c('eb4','e75_emeso','e75_ecto')

	mat_e75 = as.data.frame(mat_e75)


    return(mat_e75)
    }
	
test_diff_seeds = function(mod)
{

	cgd = mod$cgdom_ann
	all_seed_k27 = gextract(mod$epi_tracks$short_name, mod$cgdom_ann, iterator=mod$cgdom_ann)
	len = all_seed_k27$end - all_seed_k27$start
	len_norm = 300/len

	eb_tname = mod$epi_tracks$track_k27[mod$epi_tracks$short_name == "EB4_cnt"]
	emeso_tname = mod$epi_tracks$track_k27[mod$epi_tracks$short_name == "e75_emeso_cnt"]	
	ecto_tname = mod$epi_tracks$track_k27[mod$epi_tracks$short_name == "e75_ecto_cnt"]	
	cov_eb = gsummary(eb_tname)[5]
	cov_emeso = gsummary(emeso_tname)[5]
	cov_ecto = gsummary(ecto_tname)[5]

	norm_emeso  = cov_eb / cov_emeso
	norm_ecto = cov_eb / cov_ecto
	all_seed_k27$e75_emeso_cnt = all_seed_k27$e75_emeso_cnt * norm_emeso
	all_seed_k27$e75_ecto_cnt = all_seed_k27$e75_ecto_cnt * norm_ecto

	mat_e75= cbind( log2(4+all_seed_k27$EB4_cnt * len_norm), 
						 log2(4+all_seed_k27$e75_emeso_cnt * len_norm), 
						 log2(4+all_seed_k27$e75_ecto_cnt * len_norm))

	emeso_cons = mat_e75[,1] > 6 & mat_e75[,2] > 6
	emeso_loss = mat_e75[,1] > 6 & mat_e75[,2] < 4

	ecto_cons = mat_e75[,1] > 6 & mat_e75[,3] > 6
	ecto_loss = mat_e75[,1] > 6 & mat_e75[,3] < 4

	png("figs/seed_eb_emeso.png", w=500, h =500)
	plot(mat_e75[,1], mat_e75[,2], pch=19, cex=0.5, col=ifelse(emeso_loss,"gold","black"))
	dev.off()
	png("figs/seed_eb_ecto.png", w=500, h =500)
	plot(mat_e75[,1], mat_e75[,3], pch=19, cex=0.5, col=ifelse(ecto_loss,"gold","black"))
	dev.off()
	


	gvtrack.create("d_seed_cons_emeso", cgd[emeso_cons,], "distance")
	gvtrack.create("d_seed_loss_emeso", cgd[emeso_loss,], "distance")

	gvtrack.create("d_seed_cons_ecto", cgd[ecto_cons,], "distance")
	gvtrack.create("d_seed_loss_ecto", cgd[ecto_loss,], "distance")

	trace = gextract(c("seq.CG_500_mean_new", "EB4_cnt", "norm_emeso*e75_emeso_cnt", "norm_ecto*e75_ecto_cnt", "d_seed_cons_emeso", "d_seed_loss_emeso", "d_seed_cons_ecto", "d_seed_loss_ecto"), 
				intervals = gintervals.all(), 
				iterator = 100,
				colnames=c("cg", "eb4", "emeso", "ecto", "d_cons_emeso", "d_loss_emeso", "d_cons_ecto", "d_loss_ecto"))

	f_cons_emeso = trace$d_cons_emeso < trace$d_loss_emeso
	f_cons_ecto = trace$d_cons_ecto < trace$d_loss_ecto

	cg_rng = as.numeric(cut(trace$cg, c(seq(0,0.04,0.01),1)))

	eb_cons_emeso_d_cg = list()
	eb_cons_ecto_d_cg = list()
	emeso_cons_d_cg = list()
	ecto_cons_d_cg = list()

	eb_div_emeso_d_cg = list()
	eb_div_ecto_d_cg = list()
	emeso_div_d_cg = list()
	ecto_div_d_cg = list()

	for(cg in 1:5) {	
		f = f_cons_emeso & cg_rng==cg
		eb_cons_emeso_d_cg[[cg]] = tapply(log2(1+trace$eb4[f]),
								 floor(2*log2(trace$d_cons_emeso[f])), mean)
		emeso_cons_d_cg[[cg]] = tapply(log2(1+trace$emeso[f]),
								 floor(2*log2(trace$d_cons_emeso[f])), mean)


		f = !f_cons_emeso & cg_rng==cg
		eb_div_emeso_d_cg[[cg]] = tapply(log2(1+trace$eb4[f]),
								 floor(2*log2(trace$d_loss_emeso[f])), mean)
		emeso_div_d_cg[[cg]] = tapply(log2(1+trace$emeso[f]),
								 floor(2*log2(trace$d_loss_emeso[f])), mean)

		f = f_cons_ecto & cg_rng==cg
		eb_cons_ecto_d_cg[[cg]] = tapply(log2(1+trace$eb4[f]),
								 floor(2*log2(trace$d_cons_ecto[f])), mean)
		ecto_cons_d_cg[[cg]] = tapply(log2(1+trace$ecto[f]),
								 floor(2*log2(trace$d_cons_ecto[f])), mean)

		f = !f_cons_ecto & cg_rng==cg
		eb_div_ecto_d_cg[[cg]] = tapply(log2(1+trace$eb4[f]),
								 floor(2*log2(trace$d_loss_ecto[f])), mean)
		ecto_div_d_cg[[cg]] = tapply(log2(1+trace$ecto[f]),
								 floor(2*log2(trace$d_loss_ecto[f])), mean)
	}


	png("figs/decay_eb_emeso.png", w=1000, h =500)
	layout(matrix(1:2, nrow=1))
	plot(as.numeric(names(eb_cons_emeso_d_cg[[4]]))/2, eb_cons_emeso_d_cg[[4]], type="l", lwd=2, col="black", ylim=c(0,7), xlim=c(8,18))
	for(cg in 1:4) {
		lines(as.numeric(names(eb_cons_emeso_d_cg[[cg]]))/2, eb_cons_emeso_d_cg[[cg]], type="l", lwd=2)
		lines(as.numeric(names(emeso_cons_d_cg[[cg]]))/2, emeso_cons_d_cg[[cg]], type="l", lwd=2, col="gold")
	}
	plot(as.numeric(names(eb_div_emeso_d_cg[[4]]))/2, eb_div_emeso_d_cg[[4]], type="l", lwd=2, col="black", ylim=c(0,7), xlim=c(8,18))
	for(cg in 1:4) {
		lines(as.numeric(names(eb_div_emeso_d_cg[[cg]]))/2, eb_div_emeso_d_cg[[cg]], type="l", lwd=2)
		lines(as.numeric(names(emeso_div_d_cg[[cg]]))/2, emeso_div_d_cg[[cg]], type="l", lwd=2, col="gold")
	}
	dev.off()

	png("figs/decay_eb_ecto.png", w=1000, h =500)
	layout(matrix(1:2, nrow=1))
	plot(as.numeric(names(eb_cons_ecto_d_cg[[4]]))/2, log2(1+eb_cons_ecto_d_cg[[4]]), type="l", lwd=2, col="black", ylim=c(0,3.5), xlim=c(8,18))
	for(cg in 1:4) {
		lines(as.numeric(names(eb_cons_ecto_d_cg[[cg]]))/2, log2(1+eb_cons_ecto_d_cg[[cg]]), type="l", lwd=2)
		lines(as.numeric(names(ecto_cons_d_cg[[cg]]))/2, log2(1+ecto_cons_d_cg[[cg]]), type="l", lwd=2, col="green")
	}
	plot(as.numeric(names(eb_div_ecto_d_cg[[4]]))/2, log2(1+eb_div_ecto_d_cg[[4]]), type="l", lwd=2, col="black", ylim=c(0,3.5), xlim=c(8,18))
	for(cg in 1:4) {
		lines(as.numeric(names(eb_div_ecto_d_cg[[cg]]))/2, log2(1+eb_div_ecto_d_cg[[cg]]), type="l", lwd=2)
		lines(as.numeric(names(ecto_div_d_cg[[cg]]))/2, log2(1+ecto_div_d_cg[[cg]]), type="l", lwd=2, col="green")
	}
	dev.off()

	eb_cons_emeso_d_cg_q = list()
	emeso_cons_d_cg_q = list()
	for(cg in 1:5) {	
		f = f_cons_emeso & cg_rng==cg
		eb_cons_emeso_d_cg_q[[cg]] = tapply(log2(1+trace$eb4[f]),
								 floor(2*log2(trace$d_cons_emeso[f])), quantile, c(0.5,0.75,0.9,0.95))
		emeso_cons_d_cg_q[[cg]] = tapply(log2(1+trace$emeso[f]),
								 floor(2*log2(trace$d_cons_emeso[f])), quantile, c(0.5,0.75,0.9,0.95))
	}

	png("figs/cons_emeso_cg_qs.png", w=1200, h =300)
	layout(matrix(1:4,nrow=1))
	for(cg in 1:4) {
		mat = do.call('rbind', eb_cons_emeso_d_cg_q[[cg]])
		mate = do.call('rbind', emeso_cons_d_cg_q[[cg]])
		plot(as.numeric(rownames(mat))/2, mat[,1], 
					type="l", lwd=2, col="black", ylim=c(0,8), xlim=c(8,18))
		lines(as.numeric(rownames(mat))/2, mat[,2], type="l", lwd=2, col="gray")
		lines(as.numeric(rownames(mat))/2, mat[,3], type="l", 
												lwd=2, col="gray", lty=2) 
		lines(as.numeric(rownames(mate))/2, mate[,1], type="l", lwd=2, col="yellow")
		lines(as.numeric(rownames(mate))/2, mate[,2], type="l", lwd=2, col="gold")
		lines(as.numeric(rownames(mate))/2, mate[,3], type="l", lwd=2, col="gold", lty=2)
	}
	dev.off()

}


cres_near_pcg_seed = function(seeds,nm1,nm2,more,delta,atac_normed,nm1_a,nm2_a){
atac = atac_normed
#atac = as.data.frame(fread('./data/atac_ap.csv'))
#atac = atac[ !atac$chr %in% c('chrY','chrM','chrX'),]
#seeds$type_diff = ifelse(as.numeric(seeds[,nm1])> 6 & as.numeric(seeds[,nm2])>6 ,'cons','else')
#seeds$type_diff = ifelse(as.numeric(seeds[,nm1])> more & #(as.numeric(seeds[,nm1])-as.numeric(seeds[,nm2]))>delta ,paste0(nm2,'_loss'),seeds$type_diff)
#seeds$type_diff = ifelse(as.numeric(seeds[,nm2])> 5.5 & as.numeric(seeds[,nm1])<4.5 #,paste0(nm1,'_loss'),seeds$type_diff)
seeds$type_diff = seeds$col2  
print(table(seeds$type_diff))
seeds_cons = seeds[ seeds$type_diff=='cons',] 
seeds_emeso_loss = seeds[ seeds$type_diff==paste0(nm2,'_loss'),] 
seeds_epi_loss = seeds[ seeds$type_diff==paste0(nm1,'_loss'),] 




neigh_cons = atac %>% gintervals.neighbors(seeds_cons) %>%
        mutate(dist = cut(abs(dist), breaks = c(0,100, 1e3,  4e3, 1e4, 2e4, 5e4, 1e5, 2e5,5e5, 5e9), include.lowest = TRUE)) 


neigh_emeso_loss = atac %>% gintervals.neighbors(seeds_emeso_loss) %>%
        mutate(dist = cut(abs(dist), breaks = c(0,100, 1e3,  4e3, 1e4, 2e4, 5e4, 1e5, 2e5,5e5, 5e9), include.lowest = TRUE)) 



neigh_epi_loss = atac %>% gintervals.neighbors(seeds_epi_loss) %>%
        mutate(dist = cut(abs(dist), breaks = c(0,100, 1e3,  4e3, 1e4, 2e4, 5e4, 1e5, 2e5,5e5, 5e9), include.lowest = TRUE)) 





neigh_cons$g = 'cons'
neigh_emeso_loss$g = paste0(nm2,'_loss')

neigh_comb = rbind(neigh_cons,neigh_emeso_loss)

neigh_comb_f = neigh_comb[ !neigh_comb$dist %in% c('[0,100]','(5e+05,5e+09]'),]
print(table(neigh_comb$dist))
options( repr.plot.width=14,repr.plot.height=10)
gg=neigh_comb_f %>%
  ggplot(aes(x = dist, y = track2_a - track1_a ,fill=g)) +
  geom_boxplot(outlier.shape = NA) +
    labs(
   
    
    y = paste0(nm2_a, ' - ', nm1_a)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme_bw() + coord_cartesian(ylim = c(-1.5, 1.5))
    print(gg)
    return(list(gg=gg,neigh_comb_f=neigh_comb_f))
}

cres_near_pcg_seed_iq_pred = function(seeds,nm1,nm2,more,delta){
atac = as.data.frame(fread('./data/atac_ap_pred.csv'))
atac = atac[ !atac$chr %in% c('chrY','chrM','chrX'),]
seeds$type_diff = ifelse(as.numeric(seeds[,nm1])> 6 & as.numeric(seeds[,nm2])>6 ,'cons','else')
seeds$type_diff = ifelse(as.numeric(seeds[,nm1])> more & (as.numeric(seeds[,nm1])-as.numeric(seeds[,nm2]))>delta ,paste0(nm2,'_loss'),seeds$type_diff)
seeds$type_diff = ifelse(as.numeric(seeds[,nm2])> 5.5 & as.numeric(seeds[,nm1])<4.5 ,paste0(nm1,'_loss'),seeds$type_diff)
print(table(seeds$type_diff))
seeds_cons = seeds[ seeds$type_diff=='cons',] 
seeds_emeso_loss = seeds[ seeds$type_diff==paste0(nm2,'_loss'),] 
seeds_epi_loss = seeds[ seeds$type_diff==paste0(nm1,'_loss'),] 




neigh_cons = atac %>% gintervals.neighbors(seeds_cons) %>%
        mutate(dist = cut(abs(dist), breaks = c(0,100, 1e3,  4e3, 1e4, 2e4, 5e4, 1e5, 2e5,5e5, 5e9), include.lowest = TRUE)) 


neigh_emeso_loss = atac %>% gintervals.neighbors(seeds_emeso_loss) %>%
        mutate(dist = cut(abs(dist), breaks = c(0,100, 1e3,  4e3, 1e4, 2e4, 5e4, 1e5, 2e5,5e5, 5e9), include.lowest = TRUE)) 



neigh_epi_loss = atac %>% gintervals.neighbors(seeds_epi_loss) %>%
        mutate(dist = cut(abs(dist), breaks = c(0,100, 1e3,  4e3, 1e4, 2e4, 5e4, 1e5, 2e5,5e5, 5e9), include.lowest = TRUE)) 





neigh_cons$g = 'cons'
neigh_emeso_loss$g = paste0(nm2,'_loss')

neigh_comb = rbind(neigh_cons,neigh_emeso_loss)

neigh_comb_f = neigh_comb[ neigh_comb$dist != '(5e+05,5e+09]',]

options( repr.plot.width=14,repr.plot.height=10)
gg=neigh_comb_f %>%
  ggplot(aes(x = dist, y = pred,fill=g)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+ theme_bw() + coord_cartesian(ylim = c(-.25, .25))
    print(gg)
    return(gg)
}




plt_pwm_in_a_peaks = function(motif_db,seeds,nm1,nm2,more,delta,atac_normed,nm1_a,nm2_a){
atac = atac_normed
#atac = as.data.frame(fread('./data/atac_ap.csv'))
#atac = atac[ !atac$chr %in% c('chrY','chrM','chrX'),]
mots = unique(motif_db$motif)
    for (mot in mots) {
        gvtrack.create(sprintf("PWM%s", mot), func = "pwm", params = list(pssm = motif_db[motif_db$motif == 
            mot, ]))
    }
    mots_vt = paste("PWM", mots, sep = "")




	seeds$type_diff = ifelse(as.numeric(seeds[,nm1])> 6 & as.numeric(seeds[,nm2])>6 ,'cons','else')
	seeds$type_diff = ifelse(as.numeric(seeds[,nm1])> more & (as.numeric(seeds[,nm1])-as.numeric(seeds[,nm2]))>delta ,paste0(nm2,'_loss'),seeds$type_diff)
	seeds$type_diff = ifelse(as.numeric(seeds[,nm2])> 5.5 & as.numeric(seeds[,nm1])<4.5 ,paste0(nm1,'_loss'),seeds$type_diff)
	print(table(seeds$type_diff))
	seeds_cons = seeds[ seeds$type_diff=='cons',] 
	seeds_emeso_loss = seeds[ seeds$type_diff==paste0(nm2,'_loss'),] 
	seeds_epi_loss = seeds[ seeds$type_diff==paste0(nm1,'_loss'),] 
			
	
    atac_tss = gintervals.neighbors1(atac, mod$tss)
    atac_f = atac_tss[atac_tss$dist == 0, c("chrom", "start", 
        "end", "epi", "exe_meso")]
    neigh_cons = atac %>% gintervals.neighbors(seeds_cons) %>% 
        mutate(dist = cut(abs(dist), breaks = c(0, 100, 1000, 
            2000, 4000, 10000, 20000, 50000, 1e+05, 2e+05, 5e+05, 
            5e+09), include.lowest = TRUE))
    neigh_emeso_loss = atac %>% gintervals.neighbors(seeds_emeso_loss) %>% 
        mutate(dist = cut(abs(dist), breaks = c(0, 100, 1000, 
            2000, 4000, 10000, 20000, 50000, 1e+05, 2e+05, 5e+05, 
            5e+09), include.lowest = TRUE))
    neigh_epi_loss = atac %>% gintervals.neighbors(seeds_epi_loss) %>% 
        mutate(dist = cut(abs(dist), breaks = c(0, 100, 1000, 
            2000, 4000, 10000, 20000, 50000, 1e+05, 2e+05, 5e+05, 
            5e+09), include.lowest = TRUE))
    neigh_cons$g = "cons"
    neigh_emeso_loss$g = "emeso_loss"
    neigh_comb = rbind(neigh_cons, neigh_emeso_loss)
	neigh_comb_f = neigh_comb[ !neigh_comb$dist %in% c('[0,100]','(5e+05,5e+09]'),]

    neigh_comb_f = neigh_comb_f %>% arrange(chrom,start)

gext = gextract(mots_vt,intervals = neigh_comb_f,iterator = neigh_comb_f)

a_energ = cbind(neigh_comb_f,gext)

for (motif in mots_vt){

a_energ_f = a_energ[, c('dist','g',motif)]

colnames(a_energ_f) = c('dist','g','motif')

gg = a_energ_f %>%
  ggplot(aes(x = dist, y = motif, fill = g)) +
  geom_boxplot(outlier.shape = NA) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme_bw() +
  #coord_cartesian(ylim = c(-0.25, 0.25)) +
  labs(y = motif)
print(gg)
save_baseR_to_ppt(plot(1,1),link_ppt = paste0('./figs/epi_to_meso_mtfs_bxplts/',motif,'_bxplt.pptx'))
save_gg_to_ppt(gg = gg,link = paste0('./figs/epi_to_meso_mtfs_bxplts/',motif,'_bxplt.pptx'))
    }

return(a_energ)
}


ebpcg_cg_bias_test = function(mod)
{
	gen_k27_vt(mod)
   gvtrack.create("line_d", "intervs.global.rmsk_line", "distance")
   gvtrack.create("ltr_d", "intervs.global.rmsk_ltr", "distance")
   gvtrack.create("simp_d", "intervs.global.rmsk_simple_repeat", "distance")
   gvtrack.create("sine_d","intervs.global.rmsk_sine", "distance")
   gvtrack.create("lowcomplex_d", "intervs.global.rmsk_low_complexity", "distance")
	autosomes = gintervals.all()[1:19,]
	auto_nonrpt = gscreen("abs(ltr_d) > 0 & abs(line_d) > 0 & abs(sine_d) > 0", intervals=autosomes, iterator = 20)
	cg_dists = NULL
	for(sname in rownames(mod$k27_track_thresh)) {
		thresh = gquantiles(sname, intervals=auto_nonrpt, 0.985)
		message("screen ", sname, " autosome thresh ", thresh)
		d = gscreen(sprintf("%s > %s", sname, thresh), intervals = auto_nonrpt)
		cg_dist = gextract(c("seq.CG_500_mean_new"), iterator=20, intervals=d, colnames=c("CG"))

		if(is.null(cg_dists)) {
			cg_dists = tabulate(pmin(floor(cg_dist$CG/0.01), 7))
		} else {
			cg_dists = rbind(cg_dists, tabulate(pmin(floor(cg_dist$CG/0.01), 7)))
		}
	}
	rownames(cg_dists) = rownames(mod$k27_track_thresh)

	pdf("figs/cg_density_bars.pdf", w = 800, h = 400)
	select_nms = c("EB3_cnt", "EB4_cnt", "e75_ecto_cnt", "e75_emeso_cnt", "e105_limb")
	cg_dists_p = cg_dists / rowSums(cg_dists)
	barplot(cg_dists_p[select_nms,], beside=T, las=2)	
	dev.off()
}




plt_genome_pred_png_gg = function (mod, g = NA, off5 = NA, off3 = NA, mark_reg = NA, chrom = NULL, 
    locus = NULL, win_smooth = 1, show_center = F, plot_pred_doms = F, 
    label_tss = F, pred_lwd = 1, plot_base_pred = F, fn = NULL, 
    fn_w = 15, fn_h = 10, more_tracks = c()) 
{
    cg_trace = mod$gw$cg_trace
    if (is.null(chrom)) {
        f = mod$tss$geneSymbol == g
        if (sum(f) == 0) {
            return
        }
        hits = mod$tss[f, ]
        locus = mean(mod$tss$start[f])
        chrom = hits[1, "chrom"]
    }
    f = as.character(cg_trace$chrom) == chrom & (cg_trace$start > 
        locus + off5 & cg_trace$start < locus + off3)
    tss = mod$epi_tss[mod$epi_tss$chrom == chrom & (mod$epi_tss$start > 
        locus + off5 & mod$epi_tss$start < locus + off3), ]
    add_n = length(more_tracks)

    #maxy = max(max(cg_trace$lk27_1k[f]), 7)
    pred_sm = cg_trace$pred[f]
    maxy = max(max(cg_trace$lk27_1k[f]), pred_sm)
    pred_segs = which(pred_sm[-1] > 4 & pred_sm[-length(pred_sm)] < 
        4)
    pred_segs = c(pred_segs, which(pred_sm[-1] < 4 & pred_sm[-length(pred_sm)] > 
        4))
    pred_segs_xs = cg_trace$start[f][pred_segs]
    obs = cg_trace$lk27_1k[f]
    r2 = cor(pmax(obs, 2), pmax(pred_sm, 2))^2
    obs_given_pred = round(sum(obs > 6.5 & pred_sm > 5.88)/(1 + 
        sum(pred_sm > 5.88)), 3)
    pred_given_obs = round(sum(obs > 6.5 & pred_sm > 5.88)/(1 + 
        sum(obs > 6.5)), 3)

    df = data.frame(
        x = cg_trace$start[f],
        obs = cg_trace$lk27_1k[f],
        pred = pred_sm,
        chrom = cg_trace$chrom[f]
    )

    vlines = c()
    if (show_center) vlines = c(vlines, locus)
    if (!is.na(mark_reg[1])) vlines = c(vlines, mark_reg)
ymin_obs <- 4.74
obs_clip <- pmax(df$obs, ymin_obs)


p1 = ggplot2::ggplot(df, ggplot2::aes(x = x, y = obs)) +
  ggplot2::geom_line(linewidth = 0.4, color = "blue") +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin_obs, ymax = obs_clip),
                       fill = "blue", alpha = 1) +
  ggplot2::geom_hline(yintercept = 6.85, linetype = 2, linewidth = .7) +
  ggplot2::coord_cartesian(ylim = c(ymin_obs, maxy), expand = FALSE) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.text.x  = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin  = ggplot2::margin(t = 6, r = 6, b = 0, l = 6)
        )

    if (plot_pred_doms && length(pred_segs_xs) > 0) {
        p1 = p1 + ggplot2::geom_vline(xintercept = pred_segs_xs, linetype = 1)
    }
    if (length(vlines) > 0) {
        p1 = p1 + ggplot2::geom_vline(xintercept = vlines, linewidth = 0.4, color = "black")
    }
ymin_pred <- 4.83
pred_clip <- pmax(df$pred, ymin_pred)
p2 = ggplot2::ggplot(df, ggplot2::aes(x = x, y = pred)) +
  ggplot2::geom_line(linewidth = 0.4, color = "darkblue") +
  ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin_pred, ymax = pred_clip),
                       fill = "darkblue", alpha = 1) +
  ggplot2::geom_hline(yintercept = 6.22, linetype = 2, linewidth = .7) +
  ggplot2::coord_cartesian(ylim = c(ymin_pred, maxy), expand = FALSE) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(t = 0, r = 6, b = 6, l = 6)
        )

    if (length(vlines) > 0) {
        p2 = p2 + ggplot2::geom_vline(xintercept = vlines, linewidth = 0.4, color = "black")
    }

    if (plot_base_pred) {
        pred_sm_b = zoo::rollmean(cg_trace$pred_base[f], win_smooth, 
            f = "e")
        df$pred_base = pred_sm_b
        p2 = p2 + ggplot2::geom_line(ggplot2::aes(y = pred_base), linewidth = pred_lwd * 0.4, color = "cyan")
    }

    if (label_tss) {
        if (nrow(tss) > 0) {
            df_tss = data.frame(x = tss$start, y = rep(maxy - 2, times = nrow(tss)), lab = tss$geneSymbol)
            p2 = p2 + ggplot2::geom_text(data = df_tss, ggplot2::aes(x = x, y = y, label = lab),
                                         angle = -45, size = 5, inherit.aes = FALSE)
        }
    }

    if (plot_pred_doms && length(pred_segs_xs) > 0) {
        p2 = p2 + ggplot2::geom_vline(xintercept = pred_segs_xs, linetype = 1)
    }
    if (show_center) {
        p2 = p2 + ggplot2::geom_vline(xintercept = locus, linewidth = 1, color = "black")
    }

    extra_plots = list()
    if (length(more_tracks) != 0) {
        profs = gextract(more_tracks, intervals = cg_trace[f, 
            1:3], iterator = cg_trace[f, 1:3])
        for (i in 1:length(more_tracks)) {
            sm_prf = zoo::rollmean(profs[, more_tracks[i]], 5, 
                f = "e")
            df_tr = data.frame(
                x = cg_trace$start[f],
                y = log2(2 + sm_prf)
            )
            p_tr = ggplot2::ggplot(df_tr, ggplot2::aes(x = x, y = y)) +
                ggplot2::geom_line(linewidth = pred_lwd * 0.4, color = "gray") +
                ggplot2::coord_cartesian(ylim = c(1, maxy), expand = FALSE) +
                ggplot2::labs(title = more_tracks[[i]]) +
                ggplot2::theme_classic() +
                ggplot2::theme(
                    plot.margin = ggplot2::margin(t = 0, r = 6, b = 6, l = 6)
                )
            extra_plots[[length(extra_plots) + 1]] = p_tr
        }
    }

    if (!requireNamespace("patchwork", quietly = TRUE)) {
        print(p1)
        print(p2)
        if (length(extra_plots) != 0) {
            for (p in extra_plots) print(p)
        }
    } else {
        p_all = p1 / p2
        if (length(extra_plots) != 0) {
            for (p in extra_plots) p_all = p_all / p
        }
        return(p_all)
    }

    if (!is.null(fn)) {
    }
}


# Save a single ggplot (or patchwork) object to a PPTX file using officer + rvg
save_gg_to_pptx <- function(gg, path) {
  stopifnot(requireNamespace("officer", quietly = TRUE))
  stopifnot(requireNamespace("rvg", quietly = TRUE))

  doc <- officer::read_pptx()
  doc <- officer::add_slide(doc, layout = "Blank", master = "Office Theme")
  doc <- officer::ph_with(
    doc,
    value = rvg::dml(ggobj = gg),
    location = officer::ph_location_fullsize()
  )
  print(doc, target = path)
  message("Saved PPTX to: ", normalizePath(path, winslash = "/", mustWork = FALSE))
  invisible(path)
}

plot_dense_scatter_ylim = function (data, x, y, lim = c(8.5)) 
{
    options(repr.plot.width = 10, repr.plot.height = 10)
    temp = data[, c(x, y)]
    colnames(temp) = c("x", "y")
    p_coldens = densCols(x = (temp[, "x"]), y = (temp[, "y"]), 
        colramp = colorRampPalette(c("darkgray", "blue3", "red", 
            "yellow")))
    temp$col = p_coldens
    gg = temp %>% ggplot(aes(x = (x), y = y, col = col)) + geom_point(alpha = 1, 
        size = 0.8) + labs(x = x, y = y) + theme(axis.title.y = element_text(angle = 0, 
        vjust = 0.5), axis.text.x = element_text(angle = -45, 
        hjust = 0)) + theme_bw() + scale_color_identity() + xlim(c(3.3, 
        8.5))+ ylim(c(3.3,8.5))
    print(gg)
}




plot_dense_scatter_ylim_legend = function(data, x, y, lim_d = 3.3,lim_u=8.5) {

    options(repr.plot.width = 10, repr.plot.height = 10)
my_pal <- colorRampPalette(c("darkgray", "blue3", "red", "yellow"))
    temp = data[, c(x, y)]
    colnames(temp) = c("x", "y")

    # --- Compute 2D density using kde2d ---
    kd = MASS::kde2d(temp$x, temp$y, n = 200)

    # --- Convert each point to nearest kde2d grid cell ---
    ix = findInterval(temp$x, kd$x)
    iy = findInterval(temp$y, kd$y)

    # Clip to valid grid (avoid 0 or >length)
    ix[ix < 1] = 1;  ix[ix > length(kd$x)] = length(kd$x)
    iy[iy < 1] = 1;  iy[iy > length(kd$y)] = length(kd$y)

    # Extract density for each point
    temp$density = kd$z[cbind(ix, iy)]

    # --- Plot with ggplot2 and density legend ---
    gg = ggplot(temp, aes(x = x, y = y, color = density)) +
        geom_point(size = 0.8, alpha = 1) +
        scale_color_gradientn(
            colours = my_pal(200),
            name = "Density"
        ) +
        labs(x = x, y = y) +
        theme_bw() +
        theme(
            #axis.title.y = element_text(angle = 0, vjust = 0.5),
            axis.text.x = element_text(angle = -45, hjust = 0)
        ) +
      xlim(c(2.5, 8.1)) +
        ylim(c(3.4, 8.8))

    print(gg)
}
plt_regionAtac_k27_k4_qlim = function(mod, chr, st, end, marg5=2e+3, marg3=2e+3, mark=c(),  k4_max_fact=1){
#'jk.epipcg.multieb24.mcEpiblast'
 #   'jk.epipcg.atac.CRJK_0389_atac_wt_to_wt_eb_d3'
#tracks = c(
#    'jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm',
##'jk.epipcg.pcg.CRJK_0403_k27me3_wt_to_wt_eb_d3_a70ls',
#'jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3')
short_nms = c('k27','k4')
#	for(i in tracks) {
#		gvtrack.create(gsub('jk.epipcg.','',i), i, "sum") 
#		gvtrack.iterator(gsub('jk.epipcg.','',i), sshift = -140, eshift= 140)
#	}
 		gvtrack.create('k4', 'jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3', "sum") 
		gvtrack.iterator('k4', sshift = -140, eshift= 140)   
    
     		gvtrack.create('k27', 'jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm', "avg") 
		#gvtrack.iterator('k27', sshift = -500, eshift= 500)
    
    
prof = gextract(short_nms, intervals=gintervals(chr, st-marg5, end+marg3), iterator=20,colnames=short_nms)
prof[is.na(prof)] = 0 

#ds_rat = 
#k27 cov6487273.4
#k4 coverage5481174

#k27_vs = (2^as.numeric(prof[, 'k27']))-6
k27_vs = as.numeric(prof[, 'k27'])   
k4_vs = log2(6+as.numeric(prof[, 'k4']))# 
ymax_k27 =  (6.05)+3 #.99 #48.70000
ymax_k4 = log2(49.45+6)+3#.99
tmax = max(as.numeric(c(k27_vs,k4_vs)))
#tmax =  max(prof[, 'k27'])
# Set up the plot area with equal-sized panels and no spacing between them
# Important: call this BEFORE any plots are created
par(mfrow=c(2,1), oma=c(4,4,2,4), mar=c(0,0,0,0), mgp=c(2,1,0))



# Plot 2: K27 signal (scaled)



#message("scaled k27 max: ", tmax_k27, " (original max: ", max(prof[, 'k27']), ")")
plot(prof$start, pmin(k27_vs, ymax_k27),
     type="l", lwd=1, col="darkblue",
     ylim=c(4.9, ymax_k27),
     xaxt='n',  # No x-axis
     yaxt='s',  # Show y-axis
     main="", xlab="", ylab="")
polygon(c(prof$start, rev(prof$start)), c(pmin(k27_vs, ymax_k27), rep(1, length(pmin(k27_vs, ymax_k27)))), 
        col="darkblue", border=NA)
mtext("K27 Signal", side=2, line=2, cex=0.9)
		abline(v=st)
		abline(v=end)
    abline(h=6.05)
		if(!is.null(mark) & length(mark)>0) {
			for(x in mark) {
				abline(v=x, col="blue", lwd=2)
			}
		}
# Plot 1: ATAC signal
#tmax_atac = max(prof[, 'atac'])
#atac_vs = prof[, 'atac'] * tmax/tmax_atac
#message("scaled atac max: ", max(atac_vs), " (original max: ", max(prof[, 'atac']), ")")
#plot(prof$start, pmin(atac_vs, tmax), 
#     type="l", col='darkred', lwd=3,
#     ylim=c(0, tmax), 
#     #xaxt='n',  # Hide x-axis
#     yaxt='s',  # Show y-axis
#     main="",   # No individual plot title
#     xlab="", ylab="")
#mtext("Genomic Position", side=1, line=2.5, outer=TRUE)
#mtext("ATAC Signal", side=2, line=2, cex=0.9)
#		abline(v=st)
#		abline(v=end)
#		if(!is.null(mark) & length(mark)>0) {
#			for(x in mark) {
#				abline(v=x, col="blue", lwd=2)
#			}
#		}
# Plot 3: K4 signal (scaled)


#message("scaled k4 max: ", max(k4_vs), " (original max: ", max(prof[, 'k4']), ")")
plot(prof$start, pmin(k4_vs, ymax_k4),
     type="l", lwd=1, col="darkred",
     ylim=c(2.5, ymax_k4),
     yaxt='s',  # Show y-axis
     main="", xlab="", ylab="")
polygon(c(prof$start, rev(prof$start)), c(pmin(k4_vs, ymax_k4), rep(1, length(pmin(k4_vs, ymax_k4)))), 
        col="darkred", border=NA)
mtext("K4 Signal", side=2, line=2, cex=0.9)


# Add the x-axis label to the bottom of the entire figure
mtext("Genomic Position", side=1, line=2.5, outer=TRUE)
		abline(v=st)
		abline(v=end)
    abline(h=log2(6+49.45))
		if(!is.null(mark) & length(mark)>0) {
			for(x in mark) {
				abline(v=x, col="blue", lwd=2)
			}
		}
   }
pcg_report_rsqr_auc_cgdds_pptx_nofX = function(df, f_test,resp10,pred_k27,pred_k4,fit_type='IQ')
{

    cgd = df
	#f_500 = cgd$l>500
	lk27 = cgd$eb4_k27_mean
	lk4 = cgd$l6_eb4_k4_mean
	#k4_cov_1 = cgd$eb4_k4_cov_1
	#k27_cov_1 = cgd$eb4_k27_cov_1
	#f_out = lk27<6 & lk4<6
    #f_out = iq$f_out
	#f_ambig = k4_cov_1 > 0 & k27_cov_1 > 0
	f_train = !f_test
	f_test = f_test
	#f_tss = cgd$tss_dist == 0
    l10diff = lk4- lk27
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)

    pred_diff = pred_k4 - pred_k27
#	to_fit = mod$locmod[[fit_type]]$to_fit
	#f_X = cgd$chr == 'chrX'
	to_fit = resp10

	#png(sprintf("figs/model_%s_preds.png", fit_type), w=1000,h=800)
	plt1 = function(){
    layout(matrix(1:4,nrow=2))


		f1 =  f_train #& !f_X 
		f2 =  f_test#& !f_X
		#f3 = filt & f_500 & f_ambig & !f_X
		#f4 = filt & f_500 & f_out& !f_X
	for(i in 1:2) {
		if(i == 1) {
			stat = lk27
		} else {
			stat = l10diff
		}
		tr_cor = round(cor(stat[f1], pred_diff[f1],m="p")^2,3)
		te_cor = round(cor(stat[f2], pred_diff[f2], m="p")^2,3)
		
		message("comp ranges")
		ranges = c(min(pred_diff)-1,seq(ceiling(quantile(pred_diff,0.05)*100)/100,
                                        floor(quantile(pred_diff,0.95)*100)/100, 1), max(pred_diff)+1)
        #ranges = c(3.5,4.5,5,5.5,6,6.5,7,7.5,8.5)#seq(3.5, 8.5, 0.5)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred_diff[f1],ranges)), boxwex=0.15, col="blue",outline=FALSE, 
								main=sprintf("%s, rsqr train = %s,rsqr test = %s", fit_type, tr_cor, te_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred_diff[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col="gold", add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}
    save_baseR_to_ppt(plt1(),paste0('./figs/bxplts_cgdds_',fit_type,'.pptx'))
	print(plt1())
	#dev.off()

	#pdf(sprintf("figs/model_roc_%s.pdf", fit_type), w=6, h=6)
    plt_auc = function(){
	f1 = f_train
	f2 = f_test
	#f3 = !f_tss & f_500 & !f_ambig & !f_out & f_train
	#f4 = !f_tss& f_500 & !f_ambig & !f_out & f_test
	roc1 = comp_auc(pred_diff, to_fit, f1)
	roc2 = comp_auc(pred_diff, to_fit, f2)
	#roc3 = comp_auc(pred, to_fit, f3)
	#roc4 = comp_auc(pred, to_fit, f4)
        to_fit2 = to_fit
to_fit2[is.na(to_fit2)]=3
	plot(roc1$fn, roc1$tp, t="l", lwd=2, col="blue", xlab="fn", ylab="tp",
					 main=sprintf("auc train =  %s,auc test = %s, test n = %s, n1 = %s, n0 = %s", 
										round(roc1$auc,3), round(roc2$auc,3),
										sum(f2 & !is.na(to_fit) & to_fit2!=3),
										sum(f2 & !is.na(to_fit) & to_fit2==1),
										sum(f2 & !is.na(to_fit) & to_fit2==0)))
	lines(roc2$fn, roc2$tp, t="l", lwd=2, col="gold")
	#lines(roc3$fn, roc3$tp, t="l", lwd=2, "blue", lty=2)
	#lines(roc4$fn, roc4$tp, t="l", lwd=2, col="gold", lty=2)
	abline(a=0,b=1)
    }
    save_baseR_to_ppt(plt_auc(),paste0('./figs/auc_cgdds_',fit_type,'.pptx'))
	print(plt_auc())
	#dev.off()

}


pcg_report_rsqr_auc_cgdds_pptx = function(df, f_test,resp10,pred_k27,pred_k4,fit_type='IQ')
{

    cgd = df
	#f_500 = cgd$l>500
	lk27 = cgd$eb4_k27_mean
	lk4 = cgd$l6_eb4_k4_mean
	#k4_cov_1 = cgd$eb4_k4_cov_1
	#k27_cov_1 = cgd$eb4_k27_cov_1
	#f_out = lk27<6 & lk4<6
    #f_out = iq$f_out
	#f_ambig = k4_cov_1 > 0 & k27_cov_1 > 0
	f_train = !f_test
	f_test = f_test
	#f_tss = cgd$tss_dist == 0
    l10diff = lk4- lk27
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)

    pred_diff = pred_k4 - pred_k27
#	to_fit = mod$locmod[[fit_type]]$to_fit
	f_X = cgd$chr == 'chrX'
	to_fit = resp10

	#png(sprintf("figs/model_%s_preds.png", fit_type), w=1000,h=800)
	plt1 = function(){
    layout(matrix(1:4,nrow=2))


		f1 =  f_train & !f_X 
		f2 =  f_test& !f_X
		#f3 = filt & f_500 & f_ambig & !f_X
		#f4 = filt & f_500 & f_out& !f_X
	for(i in 1:2) {
		if(i == 1) {
			stat = lk27
		} else {
			stat = l10diff
		}
		tr_cor = round(cor(stat[f1], pred_diff[f1],m="p")^2,3)
		te_cor = round(cor(stat[f2], pred_diff[f2], m="p")^2,3)
		
		message("comp ranges")
		ranges = c(min(pred_diff)-1,seq(ceiling(quantile(pred_diff,0.05)*100)/100,
                                        floor(quantile(pred_diff,0.95)*100)/100, 1), max(pred_diff)+1)
        #ranges = c(3.5,4.5,5,5.5,6,6.5,7,7.5,8.5)#seq(3.5, 8.5, 0.5)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred_diff[f1],ranges)), boxwex=0.15, col="blue",outline=FALSE, 
								main=sprintf("%s, rsqr train = %s,rsqr test = %s", fit_type, tr_cor, te_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred_diff[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col="gold", add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}
    save_baseR_to_ppt(plt1(),paste0('./figs/bxplts_cgdds_',fit_type,'.pptx'))
	print(plt1())
	#dev.off()

	#pdf(sprintf("figs/model_roc_%s.pdf", fit_type), w=6, h=6)
    plt_auc = function(){
	f1 = f_train
	f2 = f_test
	#f3 = !f_tss & f_500 & !f_ambig & !f_out & f_train
	#f4 = !f_tss& f_500 & !f_ambig & !f_out & f_test
	roc1 = comp_auc(pred_diff, to_fit, f1)
	roc2 = comp_auc(pred_diff, to_fit, f2)
	#roc3 = comp_auc(pred, to_fit, f3)
	#roc4 = comp_auc(pred, to_fit, f4)
        to_fit2 = to_fit
to_fit2[is.na(to_fit2)]=3
	plot(roc1$fn, roc1$tp, t="l", lwd=2, col="blue", xlab="fn", ylab="tp",
					 main=sprintf("auc train =  %s,auc test = %s, test n = %s, n1 = %s, n0 = %s", 
										round(roc1$auc,3), round(roc2$auc,3),
										sum(f2 & !is.na(to_fit) & to_fit2!=3),
										sum(f2 & !is.na(to_fit) & to_fit2==1),
										sum(f2 & !is.na(to_fit) & to_fit2==0)))
	lines(roc2$fn, roc2$tp, t="l", lwd=2, col="gold")
	#lines(roc3$fn, roc3$tp, t="l", lwd=2, "blue", lty=2)
	#lines(roc4$fn, roc4$tp, t="l", lwd=2, col="gold", lty=2)
	abline(a=0,b=1)
    }
    save_baseR_to_ppt(plt_auc(),paste0('./figs/auc_cgdds_',fit_type,'.pptx'))
	print(plt_auc())
	#dev.off()

}

plt_region_cnt_vs_chip_gg = function(mod, chr, st, end, marg5=2e+3, marg3=2e+3, mark=c(),  k4_max_fact=1){
    gvtrack.create('ES_chip','jk.epipcg.lit.encode.es.ENCFF595SIA_H3K27me3_es', "sum") 
    gvtrack.iterator('ES_chip', sshift = -140, eshift= 140)

    gvtrack.create('ES_cnt_n','jk.epipcg.pcg.CRJK_0211_k27me3_es_50k_norm', "avg") 
    gvtrack.iterator('ES_cnt_n', sshift = -500, eshift= 500)

    short_nms = c('k27','k4')
    prof = gextract(c('ES_cnt_n','ES_chip'),
                    intervals=gintervals(chr, st-marg5, end+marg3),
                    iterator=20,
                    colnames=short_nms)
    prof[is.na(prof)] = 0 

    #ds_rat = 
    #k27 cov6487273.4
    #k4 coverage5481174

    k27_vs = as.matrix(prof[, 'k27']) # * tmax/tmax_k27
    #message("into dsamp k27",  " targ ", floor(sum(k27_vs)*1))
    #post_ds = downsample_matrix(k27_vs,  target_n = floor(sum(k27_vs)*1), seed=19)
    #k27_vs = as.numeric(post_ds)
    
    k4_vs = as.matrix(prof[, 'k4']) # * tmax/tmax_k4
    #message("into dsamp k4",  " targ ", floor(sum(k4_vs)*0.1))
    #post_ds = downsample_matrix(k4_vs,  target_n = floor(sum(k4_vs)*0.1), seed=19)
    #k4_vs = as.numeric(post_ds)

    k27_vs = zoo::rollmean(k27_vs,k=50,f='e')
    k4_vs  = zoo::rollmean(k4_vs, k=50,f='e')	
    k27_vs[is.na(k27_vs)] = 0
    k4_vs[is.na(k4_vs)]   = 0	
    #k27_vs = (2^(k27_vs))-6
    k4_vs = log2(6+k4_vs)

    tmax = max(as.numeric(c(k27_vs,k4_vs)))
    #tmax =  max(prof[, 'k27'])

    #tmax_k27 = max(prof[, 'k27']) # * k4_max_fact
    #tmax_k4  = max(as.numeric(prof[, 'k4']))

    # Build common data frame
    df <- data.frame(
        x   = prof$start,
        k27 = pmin(as.numeric(k27_vs), tmax),
        k4  = pmin(as.numeric(k4_vs),  tmax)
    )

    ## Panel 1: CnT Signal (k27)
    df_k27 <- data.frame(
        x = df$x,
        y = df$k27
    )
    poly_k27 <- data.frame(
        x = c(df_k27$x, rev(df_k27$x)),
        y = c(df_k27$y, rep(5, length(df_k27$y)))   # << was 1, now 5
    )

    p1 <- ggplot2::ggplot(df_k27, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line(colour = "darkviolet", linewidth = 0.3) +
        ggplot2::geom_polygon(
            data = poly_k27,
            ggplot2::aes(x = x, y = y),
            inherit.aes = FALSE,
            fill  = "darkviolet",
            colour = NA
        ) +
        ggplot2::coord_cartesian(ylim = c(5, 8.1), expand = FALSE) +
        ggplot2::labs(x = NULL, y = "CnT Signal", title = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            axis.title.x = ggplot2::element_blank(),
            axis.text.x  = ggplot2::element_blank(),
            axis.ticks.x = ggplot2::element_blank(),
            plot.margin  = ggplot2::margin(t = 4, r = 4, b = 0, l = 4)
        ) +
        ggplot2::geom_vline(xintercept = st) +
        ggplot2::geom_vline(xintercept = end)

    if (!is.null(mark) & length(mark) > 0) {
        for (xv in mark) {
            p1 <- p1 + ggplot2::geom_vline(xintercept = xv, colour = "blue", linewidth = 2)
        }
    }

    ## Panel 2: ChIP Signal (k4)
    df_k4 <- data.frame(
        x = df$x,
        y = df$k4
    )
    poly_k4 <- data.frame(
        x = c(df_k4$x, rev(df_k4$x)),
        y = c(df_k4$y, rep(2.5, length(df_k4$y)))  # << was 1, now 2.5
    )

    p2 <- ggplot2::ggplot(df_k4, ggplot2::aes(x = x, y = y)) +
        ggplot2::geom_line(colour = "darkgreen", linewidth = 0.3) +
        ggplot2::geom_polygon(
            data = poly_k4,
            ggplot2::aes(x = x, y = y),
            inherit.aes = FALSE,
            fill  = "darkgreen",
            colour = NA
        ) +
        ggplot2::coord_cartesian(ylim = c(2.5, 10), expand = FALSE) +
        ggplot2::labs(x = "Genomic Position", y = "ChIP Signal", title = NULL) +
        ggplot2::theme_classic() +
        ggplot2::theme(
            plot.margin = ggplot2::margin(t = 0, r = 4, b = 4, l = 4)
        ) +
        ggplot2::geom_vline(xintercept = st) +
        ggplot2::geom_vline(xintercept = end)

    if (!is.null(mark) & length(mark) > 0) {
        for (xv in mark) {
            p2 <- p2 + ggplot2::geom_vline(xintercept = xv, colour = "blue", linewidth = 2)
        }
    }

    # Stack panels like mfrow=c(2,1)
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        print(p1)
        print(p2)
    } else {
        p_all <- p1 / p2
        return(p_all)
    }
}

pcg_report_cooperativity_bxplts = function(df, resp10,pred_k27,pred_k4,fit_type='IQ',
                                              t_k4=6,t_k27=6.4,t_k4_on_k27 = 5.3)
{

    cgd = df
	#f_500 = cgd$l>500
	lk27 = cgd$eb4_k27_mean
	lk4 = cgd$l6_eb4_k4_mean
    k4 = cgd[cgd$pred_seed_k4 > t_k4 & cgd$pred_seed_k27 < t_k4_on_k27, 
        ]
    k27 = cgd[cgd$pred_seed_k27 > t_k27, ]
    cgd$d_k4 = misha.ext::gintervals.neighbors1(cgd, k4, mindist = 1)$dist
    cgd$d_k27 = misha.ext::gintervals.neighbors1(cgd, k27, mindist = 1)$dist
    #train near k4
    f_train = cgd$d_k4 < 10000 & 
        cgd$d_k27 > 10000
    #test near k27
    f_test = cgd$d_k4 > 10000 & 
        cgd$d_k27 < 10000
	#f_tss = cgd$tss_dist == 0
    l10diff = lk4- lk27
	#li = which(mod$locmod[[fit_type]]$lmod$lambda.1se==mod$locmod[[fit_type]]$lmod$lambda)

    pred_diff = pred_k4 - pred_k27
#	to_fit = mod$locmod[[fit_type]]$to_fit
	#f_X = cgd$chr == 'chrX'
	to_fit = resp10

	#png(sprintf("figs/model_%s_preds.png", fit_type), w=1000,h=800)
	plt1 = function(){
    layout(matrix(1:4,nrow=2))


		f1 =  f_train #& !f_X 
		f2 =  f_test#& !f_X
		#f3 = filt & f_500 & f_ambig & !f_X
		#f4 = filt & f_500 & f_out& !f_X
	for(i in 1:2) {
		if(i == 1) {
			stat = lk27
		} else {
			stat = l10diff
		}
		tr_cor = round(cor(stat[f1], pred_diff[f1],m="p")^2,3)
		te_cor = round(cor(stat[f2], pred_diff[f2], m="p")^2,3)
		
		message("comp ranges")
		ranges = c(min(pred_diff)-1,seq(ceiling(quantile(pred_diff,0.05)*100)/100,
                                        floor(quantile(pred_diff,0.95)*100)/100, 1), max(pred_diff)+1)
        #ranges = c(3.5,4.5,5,5.5,6,6.5,7,7.5,8.5)#seq(3.5, 8.5, 0.5)
		message("done ranges")
		nr = length(ranges)-1
		boxplot(split(stat[f1], cut(pred_diff[f1],ranges)), boxwex=0.15, col='#D60093',outline=FALSE, 
								main=sprintf("%s, rsqr nearK4 = %s,rsqr nearK27 = %s", fit_type, tr_cor, te_cor))
		message("bplt")
		boxplot(split(stat[f2], cut(pred_diff[f2],ranges)), at=0.2+1:nr,outline=FALSE,
												boxwex=0.15, col='#6600FF', add=T,xaxt='n')
		#message("bplt")
		#boxplot(split(stat[f3], cut(pred[f3],ranges)), at=0.5+1:nr,
		#										boxwex=0.15, col="darkgray", add=T, xaxt='n')
		#boxplot(split(stat[f4], cut(pred[f4],ranges)), at=0.7+1:nr,
		#										boxwex=0.15, col="lightgray", add=T, xaxt='n')
	}
	}
    save_baseR_to_ppt(plt1(),paste0('./figs/bxplts_cgdds_cooperativity',fit_type,'.pptx'))
	print(plt1())
	#dev.off()

	#pdf(sprintf("figs/model_roc_%s.pdf", fit_type), w=6, h=6)
    plt_auc = function(){
	f1 = f_train
	f2 = f_test
	#f3 = !f_tss & f_500 & !f_ambig & !f_out & f_train
	#f4 = !f_tss& f_500 & !f_ambig & !f_out & f_test
	roc1 = comp_auc(pred_diff, to_fit, f1)
	roc2 = comp_auc(pred_diff, to_fit, f2)
	#roc3 = comp_auc(pred, to_fit, f3)
	#roc4 = comp_auc(pred, to_fit, f4)
        to_fit2 = to_fit
to_fit2[is.na(to_fit2)]=3
	plot(roc1$fn, roc1$tp, t="l", lwd=2, col='#D60093', xlab="fn", ylab="tp",
					 main=sprintf("auc nearK4 =  %s,auc nearK27 = %s, nearK27 n = %s, n1 = %s, n0 = %s", 
										round(roc1$auc,3), round(roc2$auc,3),
										sum(f2 & !is.na(to_fit) & to_fit2!=3),
										sum(f2 & !is.na(to_fit) & to_fit2==1),
										sum(f2 & !is.na(to_fit) & to_fit2==0)))
	lines(roc2$fn, roc2$tp, t="l", lwd=2, col='#6600FF')
	#lines(roc3$fn, roc3$tp, t="l", lwd=2, "blue", lty=2)
	#lines(roc4$fn, roc4$tp, t="l", lwd=2, col="gold", lty=2)
	abline(a=0,b=1)
    }
    save_baseR_to_ppt(plt_auc(),paste0('./figs/auc_cgdds_cooperativity',fit_type,'.pptx'))
	print(plt_auc())
	#dev.off()

}
pcg_report_cooperativity_bxplts_hg = function (df, resp10,  pred, fit_type = "IQ", t_k4 = 3.8, t_k27 = 5.4) 
{
    cgd = df
    lk27 =log2(1+ cgd$h1_k27_sum)
    lk4 = log2(1+cgd$h1_k4_sum)
    k4 = cgd[cgd$pred5mc > 0.5 & cgd$pred < t_k4, ]
    k27 = cgd[cgd$pred5mc > 0.5 & cgd$pred > t_k27, ]
    cgd$d_k4 = misha.ext::gintervals.neighbors1(cgd, k4, mindist = 1)$dist
    cgd$d_k27 = misha.ext::gintervals.neighbors1(cgd, k27, mindist = 1)$dist
    f_train  =  cgd$d_k4 < 10000 & cgd$d_k27 > 
        10000
     f_test = cgd$d_k4 > 10000 & cgd$d_k27 < 
        10000
    l10diff = lk4 - lk27
    pred_diff = pred
    to_fit = resp10
    plt1 = function() {
        layout(matrix(1:4, nrow = 2))
        f1 = f_train
        f2 = f_test
        for (i in 1:2) {
            if (i == 1) {
                stat = lk4
            }
            else {
                stat = l10diff
            }
            tr_cor = round(cor(stat[f1], pred_diff[f1], m = "p")^2, 
                3)
            te_cor = round(cor(stat[f2], pred_diff[f2], m = "p")^2, 
                3)
            message("comp ranges")
            ranges = c(-1,4,5,5.5,6,6.5,9)
            message("done ranges")
            nr = length(ranges) - 1
            boxplot(split(stat[f1], cut(pred_diff[f1], ranges)), 
                boxwex = 0.15, col = "#D60093", outline = FALSE, 
                main = sprintf("%s, rsqr nearK4 = %s,rsqr nearK27 = %s", 
                  fit_type, tr_cor, te_cor))
            message("bplt")
            boxplot(split(stat[f2], cut(pred_diff[f2], ranges)), 
                at = 0.2 + 1:nr, outline = FALSE, boxwex = 0.15, 
                col = "#6600FF", add = T, xaxt = "n")
        }
    }
    save_baseR_to_ppt(plt1(), paste0("./figs/bxplts_cgdds_cooperativity_hg", 
        fit_type, ".pptx"))
    print(plt1())
    plt_auc = function() {
        f1 = f_train
        f2 = f_test
        roc1 = comp_auc(-pred_diff, to_fit, f1)
        roc2 = comp_auc(-pred_diff, to_fit, f2)
        to_fit2 = to_fit
        to_fit2[is.na(to_fit2)] = 3
        plot(roc1$fn, roc1$tp, t = "l", lwd = 2, col = "#D60093", 
            xlab = "fn", ylab = "tp", main = sprintf("auc nearK4 =  %s,auc nearK27 = %s, nearK27 n = %s, n1 = %s, n0 = %s", 
                round(roc1$auc, 3), round(roc2$auc, 3), sum(f2 & 
                  !is.na(to_fit) & to_fit2 != 3), sum(f2 & !is.na(to_fit) & 
                  to_fit2 == 1), sum(f2 & !is.na(to_fit) & to_fit2 == 
                  0)))
        lines(roc2$fn, roc2$tp, t = "l", lwd = 2, col = "#6600FF")
        abline(a = 0, b = 1)
    }
    save_baseR_to_ppt(plt_auc(), paste0("./figs/auc_cgdds_cooperativity_hg", 
        fit_type, ".pptx"))
    print(plt_auc())
}


plt_genome_pred_CG_ppt_st_end_borzoi_forA_lm_cg_gg = function (mod, g = NA, off5 = NA, off3 = NA, mark_reg = NA, chrom = NULL, 
    locus = NULL, start, end, win_smooth = 1, show_center = TRUE, 
    plot_pred_doms = FALSE, label_tss = FALSE, pred_lwd = 1, 
    plot_base_pred = FALSE, fn = NULL, fn_w = 15, fn_h = 10, 
    more_tracks = c()) 
{
    cg_trace = mod$gw$cg_trace
    if (is.null(chrom)) {
        f = mod$tss$geneSymbol == g
        if (sum(f) == 0) 
            return()
        hits = mod$tss[f, ]
        locus = mean(mod$tss$start[f])
        chrom = hits[1, "chrom"]
    }
    if (is.null(locus)) {
        f = as.character(cg_trace$chrom) == chrom & (cg_trace$start > 
            start + off5 & cg_trace$start < end + off3)
        tss = mod$epi_tss[mod$epi_tss$chrom == chrom & (mod$epi_tss$start > 
            start + off5 & mod$epi_tss$start < end + off3), ]
    }
    else {
        f = as.character(cg_trace$chrom) == chrom & (cg_trace$start > 
            locus + off5 & cg_trace$start < locus + off3)
        tss = mod$epi_tss[mod$epi_tss$chrom == chrom & (mod$epi_tss$start > 
            locus + off5 & mod$epi_tss$start < locus + off3), 
            ]
    }
    add_n = length(more_tracks)
    n_panels = 2 + add_n
    maxy = max(max(cg_trace$lk27_1k[f]), 7)
    pred_sm = cg_trace$pred_sns_brz2k[f]
    pred_sm_lm = cg_trace$pred_sns_lm[f]
    pred_sm_borzoi = cg_trace$pred_borzoi[f]
    pred_sm_silicus = cg_trace$pred_sil[f]
    pred_segs = which(pred_sm[-1] > 4 & pred_sm[-length(pred_sm)] < 
        4)
    pred_segs = c(pred_segs, which(pred_sm[-1] < 4 & pred_sm[-length(pred_sm)] > 
        4))
    pred_segs_xs = cg_trace$start[f][pred_segs]
    df = data.frame(x = cg_trace$start[f], obs = cg_trace$lk27_1k[f], 
        pred_sns_lm = pred_sm_lm, chrom = cg_trace$chrom[f])
    vlines = c()
    if (show_center) 
        vlines = c(vlines, locus)
    if (!is.na(mark_reg[1])) 
        vlines = c(vlines, mark_reg)
    ymin_obs <- 4.74
    obs_clip <- pmax(df$obs, ymin_obs)
    p1 = ggplot2::ggplot(df, ggplot2::aes(x = x, y = obs)) + 
        ggplot2::geom_line(linewidth = 0.4, color = "blue") + 
        ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin_obs, ymax = obs_clip), 
            fill = "blue", alpha = 1) + ggplot2::geom_hline(yintercept = 6.85, 
        linetype = 2) + ggplot2::coord_cartesian(ylim = c(ymin_obs, 
        maxy), expand = FALSE) + ggplot2::labs(x = "", y = "obs") + 
        ggplot2::theme_classic() + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
        plot.margin = ggplot2::margin(t = 0, r = 6, b = 0, l = 6))
    if (plot_pred_doms && length(pred_segs_xs) > 0) {
        p1 = p1 + ggplot2::geom_vline(xintercept = pred_segs_xs)
    }
    if (length(vlines) > 0) {
        p1 = p1 + ggplot2::geom_vline(xintercept = vlines, linewidth = 0.4, 
            color = "black")
    }
    ymin_pred <- 4.83
    pred_clip <- pmax(df$pred_sns_lm, ymin_pred)
    p2 = ggplot2::ggplot(df, ggplot2::aes(x = x, y = pred_sns_lm)) + 
        ggplot2::geom_line(linewidth = pred_lwd * 0.4, color = "darkblue") + 
        ggplot2::geom_ribbon(ggplot2::aes(ymin = ymin_pred, ymax = pred_clip), 
            fill = "darkblue", alpha = 1) + ggplot2::geom_hline(yintercept = 6.22, 
        linetype = 2) + ggplot2::coord_cartesian(ylim = c(ymin_pred, 
        maxy), expand = FALSE) + ggplot2::labs(x = "", y = "pred_sns_lm") + 
        ggplot2::theme_classic() + ggplot2::theme(axis.title.x = ggplot2::element_blank(), 
        axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), 
        plot.margin = ggplot2::margin(t = 0, r = 6, b = 0, l = 6))
    if (!is.na(mark_reg[1])) {
        p2 = p2 + ggplot2::geom_vline(xintercept = mark_reg, 
            linewidth = 0.4, color = "black")
    }
    if (length(vlines) > 0) {
        p2 = p2 + ggplot2::geom_vline(xintercept = vlines, linewidth = 0.4, 
            color = "black")
    }
    extra_plots = list()
    if (add_n > 0) {
        profs = gextract(more_tracks, intervals = cg_trace[f, 
            1:3], iterator = cg_trace[f, 1:3])
        for (i in seq_along(more_tracks)) {
            is_last = (i == add_n)
            df_tr = data.frame(x = cg_trace$start[f], y = profs[, 
                more_tracks[i]])
            p_tr = ggplot2::ggplot(df_tr, ggplot2::aes(x = x, y = y)) +
                ggplot2::geom_line(linewidth = pred_lwd * 0.4, color = "black") +
                ggplot2::coord_cartesian(ylim = c(0, 0.1), expand = FALSE) +
                ggplot2::labs(x = ifelse(is_last, "position", ""), y = more_tracks[i]) +
                ggplot2::theme_classic() +
                ggplot2::theme(
                    plot.margin = ggplot2::margin(t = 0, r = 6, b = ifelse(is_last, 6, 0), l = 6)
                )
            if (!is_last) {
                p_tr = p_tr + ggplot2::theme(
                    axis.title.x = ggplot2::element_blank(),
                    axis.text.x  = ggplot2::element_text(size = 0),
                    axis.ticks.x = ggplot2::element_line(linewidth = 0)
                )
            }
            extra_plots[[length(extra_plots) + 1]] = p_tr
        }
    }
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        print(p1)
        print(p2)
        if (length(extra_plots) != 0) {
            for (p in extra_plots) print(p)
        }
    }
    else {
        p_all = p1/p2
        if (length(extra_plots) != 0) {
            for (p in extra_plots) p_all = p_all/p
        }
        return(p_all)
    }
    if (!is.null(fn)) {
    }
}
plt_regionAT_ds_cg_meth_gg = function (mod, chr, st, end, k, ds_k27 = TRUE, marg5 = 2000, marg3 = 2000, 
    mark = c(), wk4 = F, k4_max_fact = 1) 
{
    tracks_k27 = c("jk.epipcg.pcg.CRJK_0363_k27me3_eb_j1_d3_a_norm", 
        "jk.epipcg.pcg.CRJK_0365_k27me3_eb_dko_d3_a_norm")

    short_namesk27 = c('wt_k27','dko_k27')
    for (i in 1:length(tracks_k27)) {
        gvtrack.create(short_namesk27[i], tracks_k27[i], "avg")
        #gvtrack.iterator(short_namesk27[i], sshift = -140, eshift = 140)
    }
    
    prof = gextract(short_namesk27, intervals = gintervals(chr, 
        st - marg5, end + marg3), iterator = 20)
    prof[is.na(prof)] = 0
    ylim <- ( max(apply(prof[, short_namesk27], 2, function(x) {
      max(zoo::rollmean(x, k = k, fill = 0))
    })))
	
    ndx = data.frame(rep(1,length(short_namesk27)))
    layout(matrix(1:nrow(ndx), ncol = 1), h = c(1.4, rep(1, nrow(ndx) - 
        2), 1.4))

    ## --- BEGIN ggplot replacement (no other changes) ---
    p_list <- list()
    for (i in seq_len(nrow(ndx))) {
        track_name <- short_namesk27[i]
        yvals <- zoo::rollmean(prof[[track_name]], k = k, fill = 0)
        df_plot <- data.frame(x = prof$start, y = yvals)

        is_last <- (i == nrow(ndx))
        # choose colors and y-limits exactly as in original
        col_line <- if (is_last) "turquoise" else "darkblue"
        y_min <- if (is_last) 4.69 else 4.78
        y_lim <- c(y_min, ylim)

        ## CLAMP the ribbon top so it never extends below the panel bottom
        df_plot$y_clip <- pmax(df_plot$y, y_min)

        theme_base <- ggplot2::theme_classic() +
            ggplot2::theme(
                plot.margin = ggplot2::margin(t = ifelse(i == 1, 8, 0), r = 6, b = ifelse(is_last, 8, 0), l = 6)
            )
# build polygon for fill (mimics base polygon)
df_poly <- data.frame(x = c(df_plot$x, rev(df_plot$x)),
                      y = c(df_plot$y_clip, rep(y_min, length(df_plot$y_clip))))

        
        
        p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_line(size = 1, colour = col_line) +
            #ggplot2::geom_ribbon(ggplot2::aes(ymin = y_min, ymax = y_clip), fill = col_line, colour = NA) +
            ggplot2::coord_cartesian(ylim = y_lim, expand = FALSE) +
            ggplot2::labs(x = ifelse(i == nrow(ndx), "position", ""), y = NA, title = track_name) +
            theme_base
p <- p + ggplot2::geom_polygon(data = df_poly, ggplot2::aes(x = x, y = y), inherit.aes = FALSE, fill = col_line, colour = NA)

        # add axis elements conditionally (avoid ifelse inside theme)
        if (is_last) {
            p <- p + ggplot2::theme(
                axis.title.x = ggplot2::element_text(),
                axis.text.x  = ggplot2::element_text(),
                axis.ticks.x = ggplot2::element_line()
            )
        } else {
            p <- p + ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x  = ggplot2::element_text(size = 0),
                axis.ticks.x = ggplot2::element_line(linewidth = 0)
            )
        }

        # add vertical lines corresponding to st, end and mark points (as in original)
        p <- p + ggplot2::geom_vline(xintercept = st)
        p <- p + ggplot2::geom_vline(xintercept = end)
        if (!is.null(mark) & length(mark) > 0) {
            for (x in mark) {
                p <- p + ggplot2::geom_vline(xintercept = x, colour = "blue", size = 1)
            }
        }

        p_list[[length(p_list) + 1]] <- p
    }

    # --- ADDITIONAL PLOT: use optional `add_track` provided by caller (no signature change) ---
    # If the caller provided a variable named `add_track` in their environment (string of a single track name),
    # extract it and append a black line plot (ymin = 0, ymax = value) with coord_cartesian(ylim = c(0,0.1))
    add_track <- 'seq.CG_500_mean_new'
    if (!is.null(add_track) && length(add_track) > 0) {
        # extract the add_track similarly to earlier gextract usage
        prof_add <- gextract(add_track, intervals = gintervals(chr, st - marg5, end + marg3), iterator = 20)
        prof_add[is.na(prof_add)] <- 0
        # assume single track name — take first if vector
        at_name <- add_track[1]
        if (at_name %in% colnames(prof_add)) {
            yvals_add <- prof_add[[at_name]]
            df_add <- data.frame(x = prof_add$start, y = yvals_add)
            p_add <- ggplot2::ggplot(df_add, ggplot2::aes(x = x, y = y)) +
                ggplot2::geom_line(size = 1, colour = "black") +
                ggplot2::coord_cartesian(ylim = c(0, 0.06), expand = FALSE) +
                ggplot2::labs(x = "position", y = NA, title = at_name) +
                ggplot2::theme_classic() +
                ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 6, b = 8, l = 6))
               p_add <- p_add + ggplot2::geom_vline(xintercept = st)
        p_add<- p_add + ggplot2::geom_vline(xintercept = end)
            p_list[[length(p_list) + 1]] <- p_add
        }
    }
    # --- END additional plot addition ---

    # print or return stacked plots like the original layout behavior
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        for (p in p_list) print(p)
    } else {
        p_all <- p_list[[1]]
        if (length(p_list) > 1) {
            for (j in 2:length(p_list)) p_all <- p_all / p_list[[j]]
        }
        return(p_all)
    }
    ## --- END ggplot replacement ---
}

plt_regionAT_ds_cg_meth_gg_nlg = function (mod, chr, st, end, k, ds_k27 = TRUE, marg5 = 2000, marg3 = 2000, 
    mark = c(), wk4 = F, k4_max_fact = 1) 
{
    tracks_k27 = c("jk.epipcg.pcg.CRJK_0363_k27me3_eb_j1_d3_a_norm", 
        "jk.epipcg.pcg.CRJK_0365_k27me3_eb_dko_d3_a_norm")

    short_namesk27 = c('wt_k27','dko_k27')
    for (i in 1:length(tracks_k27)) {
        gvtrack.create(short_namesk27[i], tracks_k27[i], "avg")
        #gvtrack.iterator(short_namesk27[i], sshift = -140, eshift = 140)
    }
    
    prof = gextract(c('2^wt_k27','2^dko_k27'), intervals = gintervals(chr, 
        st - marg5, end + marg3), iterator = 20,colnames = c(short_namesk27))
    prof[is.na(prof)] = 0
    ylim <- ( max(apply(prof[, short_namesk27], 2, function(x) {
      max(zoo::rollmean(x, k = k, fill = 0))
    })))
	
    ndx = data.frame(rep(1,length(short_namesk27)))
    layout(matrix(1:nrow(ndx), ncol = 1), h = c(1.4, rep(1, nrow(ndx) - 
        2), 1.4))

    ## --- BEGIN ggplot replacement (no other changes) ---
    p_list <- list()
    for (i in seq_len(nrow(ndx))) {
        track_name <- short_namesk27[i]
        yvals <- zoo::rollmean((prof[[track_name]]), k = k, fill = 0)
        df_plot <- data.frame(x = prof$start, y = yvals)

        is_last <- (i == nrow(ndx))
        # choose colors and y-limits exactly as in original
        col_line <- if (is_last) "turquoise" else "darkblue"
        y_min <- if (is_last) 25.81 else 27.47
        y_lim <- c(y_min, ylim)

        ## CLAMP the ribbon top so it never extends below the panel bottom
        df_plot$y_clip <- pmax(df_plot$y, y_min)

        theme_base <- ggplot2::theme_classic() +
            ggplot2::theme(
                plot.margin = ggplot2::margin(t = ifelse(i == 1, 8, 0), r = 6, b = ifelse(is_last, 8, 0), l = 6)
            )
# build polygon for fill (mimics base polygon)
df_poly <- data.frame(x = c(df_plot$x, rev(df_plot$x)),
                      y = c(df_plot$y_clip, rep(y_min, length(df_plot$y_clip))))

        
        
        p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_line(size = 1, colour = col_line) +
            #ggplot2::geom_ribbon(ggplot2::aes(ymin = y_min, ymax = y_clip), fill = col_line, colour = NA) +
            ggplot2::coord_cartesian(ylim = y_lim, expand = FALSE) +
            ggplot2::labs(x = ifelse(i == nrow(ndx), "position", ""), y = NA, title = track_name) +
            theme_base
p <- p + ggplot2::geom_polygon(data = df_poly, ggplot2::aes(x = x, y = y), inherit.aes = FALSE, fill = col_line, colour = NA)

        # add axis elements conditionally (avoid ifelse inside theme)
        if (is_last) {
            p <- p + ggplot2::theme(
                axis.title.x = ggplot2::element_text(),
                axis.text.x  = ggplot2::element_text(),
                axis.ticks.x = ggplot2::element_line()
            )
        } else {
            p <- p + ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x  = ggplot2::element_text(size = 0),
                axis.ticks.x = ggplot2::element_line(linewidth = 0)
            )
        }

        # add vertical lines corresponding to st, end and mark points (as in original)
        p <- p + ggplot2::geom_vline(xintercept = st)
        p <- p + ggplot2::geom_vline(xintercept = end)
        if (!is.null(mark) & length(mark) > 0) {
            for (x in mark) {
                p <- p + ggplot2::geom_vline(xintercept = x, colour = "blue", size = 1)
            }
        }

        p_list[[length(p_list) + 1]] <- p
    }

    # --- ADDITIONAL PLOT: use optional `add_track` provided by caller (no signature change) ---
    # If the caller provided a variable named `add_track` in their environment (string of a single track name),
    # extract it and append a black line plot (ymin = 0, ymax = value) with coord_cartesian(ylim = c(0,0.1))
    add_track <- 'seq.CG_500_mean_new'
    if (!is.null(add_track) && length(add_track) > 0) {
        # extract the add_track similarly to earlier gextract usage
        prof_add <- gextract(add_track, intervals = gintervals(chr, st - marg5, end + marg3), iterator = 20)
        prof_add[is.na(prof_add)] <- 0
        # assume single track name — take first if vector
        at_name <- add_track[1]
        if (at_name %in% colnames(prof_add)) {
            yvals_add <- prof_add[[at_name]]
            df_add <- data.frame(x = prof_add$start, y = yvals_add)
            p_add <- ggplot2::ggplot(df_add, ggplot2::aes(x = x, y = y)) +
                ggplot2::geom_line(size = 1, colour = "black") +
                ggplot2::coord_cartesian(ylim = c(0, 0.06), expand = FALSE) +
                ggplot2::labs(x = "position", y = NA, title = at_name) +
                ggplot2::theme_classic() +
                ggplot2::theme(plot.margin = ggplot2::margin(t = 0, r = 6, b = 8, l = 6))
               p_add <- p_add + ggplot2::geom_vline(xintercept = st)
        p_add<- p_add + ggplot2::geom_vline(xintercept = end)
            p_list[[length(p_list) + 1]] <- p_add
        }
    }
    # --- END additional plot addition ---

    # print or return stacked plots like the original layout behavior
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        for (p in p_list) print(p)
    } else {
        p_all <- p_list[[1]]
        if (length(p_list) > 1) {
            for (j in 2:length(p_list)) p_all <- p_all / p_list[[j]]
        }
        return(p_all)
    }
    ## --- END ggplot replacement ---
}

plt_regionAT_ds_cg_cons_gg = function (mod, chr, st, end, k, ds_k27 = TRUE, marg5 = 2000, marg3 = 2000, 
    mark = c(), wk4 = F, k4_max_fact = 1) 
{
    tracks_k27 = c("jk.epipcg.pcg.CRJK_0363_k27me3_eb_j1_d3_a_norm", 
        "jk.epipcg.pcg.CRJK_0365_k27me3_eb_dko_d3_a_norm")

    short_namesk27 = c('wt_k27','dko_k27')
    for (i in 1:length(tracks_k27)) {
        gvtrack.create(short_namesk27[i], tracks_k27[i], "avg")
        #gvtrack.iterator(short_namesk27[i], sshift = -140, eshift = 140)
    }
    
    prof = gextract(short_namesk27, intervals = gintervals(chr, 
        st - marg5, end + marg3), iterator = 20)
    prof[is.na(prof)] = 0
    ylim <- ( max(apply(prof[, short_namesk27], 2, function(x) {
      max(zoo::rollmean(x, k = k, fill = 0))
    })))
	
    ndx = data.frame(rep(1,length(short_namesk27)))
    layout(matrix(1:nrow(ndx), ncol = 1), h = c(1.4, rep(1, nrow(ndx) - 
        2), 1.4))

    ## --- BEGIN ggplot replacement (no other changes) ---
    p_list <- list()
    for (i in seq_len(nrow(ndx))) {
        track_name <- short_namesk27[i]
        yvals <- zoo::rollmean(prof[[track_name]], k = k, fill = 0)
        df_plot <- data.frame(x = prof$start, y = yvals)

        is_last <- (i == nrow(ndx))
        # choose colors and y-limits exactly as in original
        col_line <- if (is_last) "turquoise" else "darkblue"
        y_min <- if (is_last) 4.69 else 4.78
        y_lim <- c(y_min, ylim)

        theme_base <- ggplot2::theme_classic() +
            ggplot2::theme(
                plot.margin = ggplot2::margin(t = ifelse(i == 1, 8, 0), r = 6, b = ifelse(is_last, 8, 0), l = 6)
            )

        ## use polygon instead of ribbon (minimal change)
        df_poly <- data.frame(
            x = c(df_plot$x, rev(df_plot$x)),
            y = c(df_plot$y, rep(0, length(df_plot$y)))
        )

        p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_line(size = 1, colour = col_line) +
            ggplot2::geom_polygon(
                data = df_poly,
                ggplot2::aes(x = x, y = y),
                inherit.aes = FALSE,
                fill = col_line,
                colour = NA
            ) +
            ggplot2::coord_cartesian(ylim = y_lim, expand = FALSE) +
            ggplot2::labs(x = ifelse(i == nrow(ndx), "position", ""), y = NA, title = track_name) +
            theme_base
   
        # add axis elements conditionally (avoid ifelse inside theme)
        if (is_last) {
            p <- p + ggplot2::theme(
                axis.title.x = ggplot2::element_text(),
                axis.text.x  = ggplot2::element_text(),
                axis.ticks.x = ggplot2::element_line()
            )
        } else {
            p <- p + ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x  = ggplot2::element_text(size = 0),
                axis.ticks.x = ggplot2::element_line(linewidth = 0)
            )
        }

        # add vertical lines corresponding to st, end and mark points (as in original)
        p <- p + ggplot2::geom_vline(xintercept = st)
        p <- p + ggplot2::geom_vline(xintercept = end)
        if (!is.null(mark) & length(mark) > 0) {
            for (x in mark) {
                p <- p + ggplot2::geom_vline(xintercept = x, colour = "blue", size = 1)
            }
        }

        p_list[[length(p_list) + 1]] <- p
    }

    # print or return stacked plots like the original layout behavior
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        for (p in p_list) print(p)
    } else {
        p_all <- p_list[[1]]
        if (length(p_list) > 1) {
            for (j in 2:length(p_list)) p_all <- p_all / p_list[[j]]
        }
        return(p_all)
    }
    ## --- END ggplot replacement ---
}


plt_regionAT_ds_cg_cons_gg_nlg = function (mod, chr, st, end, k, ds_k27 = TRUE, marg5 = 2000, marg3 = 2000, 
    mark = c(), wk4 = F, k4_max_fact = 1) 
{
    tracks_k27 = c("jk.epipcg.pcg.CRJK_0363_k27me3_eb_j1_d3_a_norm", 
        "jk.epipcg.pcg.CRJK_0365_k27me3_eb_dko_d3_a_norm")

    short_namesk27 = c('wt_k27','dko_k27')
    for (i in 1:length(tracks_k27)) {
        gvtrack.create(short_namesk27[i], tracks_k27[i], "avg")
        #gvtrack.iterator(short_namesk27[i], sshift = -140, eshift = 140)
    }
    
    prof = gextract(c('2^wt_k27','2^dko_k27'), intervals = gintervals(chr, 
        st - marg5, end + marg3), iterator = 20,colnames = short_namesk27 )
    prof[is.na(prof)] = 0
    ylim <- ( max(apply(prof[, short_namesk27], 2, function(x) {
      max(zoo::rollmean(x, k = k, fill = 0))
    })))
	
    ndx = data.frame(rep(1,length(short_namesk27)))
    layout(matrix(1:nrow(ndx), ncol = 1), h = c(1.4, rep(1, nrow(ndx) - 
        2), 1.4))

    ## --- BEGIN ggplot replacement (no other changes) ---
    p_list <- list()
    for (i in seq_len(nrow(ndx))) {
        track_name <- short_namesk27[i]
        yvals <- zoo::rollmean(prof[[track_name]], k = k, fill = 0)
        df_plot <- data.frame(x = prof$start, y = yvals)

        is_last <- (i == nrow(ndx))
        # choose colors and y-limits exactly as in original
        col_line <- if (is_last) "turquoise" else "darkblue"
        y_min <- if (is_last) 2^4.69 else 2^4.78
        y_lim <- c(y_min, ylim)

        theme_base <- ggplot2::theme_classic() +
            ggplot2::theme(
                plot.margin = ggplot2::margin(t = ifelse(i == 1, 8, 0), r = 6, b = ifelse(is_last, 8, 0), l = 6)
            )

        ## use polygon instead of ribbon (minimal change)
        df_poly <- data.frame(
            x = c(df_plot$x, rev(df_plot$x)),
            y = c(df_plot$y, rep(0, length(df_plot$y)))
        )

        p <- ggplot2::ggplot(df_plot, ggplot2::aes(x = x, y = y)) +
            ggplot2::geom_line(size = 1, colour = col_line) +
            ggplot2::geom_polygon(
                data = df_poly,
                ggplot2::aes(x = x, y = y),
                inherit.aes = FALSE,
                fill = col_line,
                colour = NA
            ) +
            ggplot2::coord_cartesian(ylim = y_lim, expand = FALSE) +
            ggplot2::labs(x = ifelse(i == nrow(ndx), "position", ""), y = NA, title = track_name) +
            theme_base
   
        # add axis elements conditionally (avoid ifelse inside theme)
        if (is_last) {
            p <- p + ggplot2::theme(
                axis.title.x = ggplot2::element_text(),
                axis.text.x  = ggplot2::element_text(),
                axis.ticks.x = ggplot2::element_line()
            )
        } else {
            p <- p + ggplot2::theme(
                axis.title.x = ggplot2::element_blank(),
                axis.text.x  = ggplot2::element_text(size = 0),
                axis.ticks.x = ggplot2::element_line(linewidth = 0)
            )
        }

        # add vertical lines corresponding to st, end and mark points (as in original)
        p <- p + ggplot2::geom_vline(xintercept = st)
        p <- p + ggplot2::geom_vline(xintercept = end)
        if (!is.null(mark) & length(mark) > 0) {
            for (x in mark) {
                p <- p + ggplot2::geom_vline(xintercept = x, colour = "blue", size = 1)
            }
        }

        p_list[[length(p_list) + 1]] <- p
    }

    # print or return stacked plots like the original layout behavior
    if (!requireNamespace("patchwork", quietly = TRUE)) {
        for (p in p_list) print(p)
    } else {
        p_all <- p_list[[1]]
        if (length(p_list) > 1) {
            for (j in 2:length(p_list)) p_all <- p_all / p_list[[j]]
        }
        return(p_all)
    }
    ## --- END ggplot replacement ---
}
