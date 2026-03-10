suppressWarnings(
  suppressPackageStartupMessages({
    library(glmnet)
    library(doParallel)
  })
)
gen_tableS1 = function(){
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

tss = mod$tss
cgd = mod$cgdom_ann
rna = as.data.frame(mod$eb_legc)

df = gintervals.neighbors(tss,cgd,maxdist = 0,mindist = 0,na.if.notfound = F,use_intervals1_strand = TRUE)



#gen_k27_vt

cgtrace = mod$gw$cg_trace[f_norp,c('chrom','start','end','pred_sns_lm','lk27_1k')]

df2 = gintervals.neighbors(df,cgtrace,maxdist = 0,mindist = 0,na.if.notfound = F,use_intervals1_strand = TRUE)

df3 = df2[,c('chrom','start','end','strand','pred_sns_lm','lk27_1k','geneSymbol','chrom1','start1','end1')]
colnames(df3) = c('chrom','start','end','strand','pred','lk27_1k','geneSymbol','chrom_cgdd','start_cgdd','end_cgdd')



colnames(rna) = c('rna')

rna = rna %>% rownames_to_column(var = 'geneSymbol')

df4 = merge(df3,rna,by='geneSymbol')



df5 = df4[!(grepl('^Gm',df4$geneSymbol)|grepl('Rik$',df4$geneSymbol)),]


df5$delta = df5$pred - df5$lk27_1k

range(pmin(df5$rna,-12))

df5$rna_fl = pmin(df5$rna,-12)

df5$rna_bin = cut(df5$rna,c(-18,-16,-15.5,-15,-14,-13,-6))
df_d = df5 %>% filter(abs(delta)>2)#3
dim(df_d)

df_d$type = ifelse(df_d$delta >0 , 'over','under')
return(df_d)
    }
	
annotate_misses = function(){
miss = mod$gw$over_doms

#write.csv(miss3,'./data/gw_01back_misses.csv')

gintervals.ls("intervs.global.rmsk_simple_repeat")

gext = gextract(c("seq.CG_500_mean_new","seq.GC500_bin20" ),intervals = miss,iterator = 20,colnames=c('cg','gc'))

gext[is.na(gext)] = 0

miss = cbind(miss,tgs_matrix_tapply(t(as.matrix(gext[,c('cg','gc')])), gext$intervalID,max))

miss2 = gintervals.neighbors1(miss,mod$cgd_rpt[, c('chrom','start','end','l',
                                            'exon','line','sine',
                                                   'ltr','simp','low_complex','tot_rpt')],
                             na.if.notfound = F,maxdist = 1e9)

miss2 = miss2 %>% arrange(chrom,start)
miss2 = as.data.frame(miss2)

    gvtrack.create("d_ltr", "intervs.global.rmsk_ltr", "distance")
    gvtrack.create("d_line", "intervs.global.rmsk_line", "distance")
    gvtrack.create("d_sine", "intervs.global.rmsk_sine", "distance")
	gvtrack.create("simp_d", "intervs.global.rmsk_simple_repeat", "distance")
	gvtrack.create("lowcomplex_d", "intervs.global.rmsk_low_complexity", "distance")
rpt_trace = gextract(c("d_line", "d_ltr", "d_sine", 
       "simp_d","lowcomplex_d"), 
        intervals = miss2, iterator = miss2, colnames = c("d_line", 
            "d_ltr", "d_sine","simp_d","lowcomplex_d"))
    miss2$d_ltr = rpt_trace$d_ltr
    miss2$d_line = rpt_trace$d_line
    miss2$d_sine = rpt_trace$d_sine

	miss2$simp_d = rpt_trace$simp_d
	miss2$lowcomplex_d = rpt_trace$lowcomplex_d	


miss3 = gintervals.neighbors1(miss2,mod$cgdom_ann[,c('chrom','start','end')])

miss3$type = ifelse( abs(miss3$d_ltr)==0 |  abs(miss3$dist)==0 | abs(miss3$d_line)==0 ,'rpt','else')

miss3$type = ifelse( miss3$type == 'else' & miss3$cg <= .022  ,'low_cg',miss3$type)
miss3$type = ifelse( miss3$type == 'else' & miss3$dist1<1e4  ,'cgdd',miss3$type)
miss3$type = ifelse( miss3$type == 'else' & miss3$cg>.02& miss3$cg <= .04  ,'int_cg',miss3$type)

miss3$type2 = ifelse( miss3$type == 'low_cg' & miss3$dist > .5e5& miss3$dist1 > .5e5 ,'lowcg_far',miss3$type)
miss3$type2 = ifelse( miss3$type == 'low_cg' & (miss3$dist <= .5e5 | miss3$dist1 <= .5e5)  ,'lowcg_close',miss3$type2)

miss3$type2 = ifelse( miss3$type == 'int_cg' & miss3$dist > .5e5& miss3$dist1 > .5e5 ,'intcg_far',miss3$type2)
miss3$type2 = ifelse( miss3$type == 'int_cg' & (miss3$dist <= .5e5 | miss3$dist1 <= .5e5)  ,'intcg_close',miss3$type2)

miss3$type2 = ifelse( miss3$type == 'else'  ,'cgdd',miss3$type2)

table(miss3$type2)
return(miss3)
}
# ---- Leave-one-chromosome-out CV wrapper ----

loco_cv <- function(y, x, chrom=as.character((mod$cgdom_ann$chrom))) {
  if (length(y) != nrow(x)) {
    stop("y and x must have the same number of rows / observations.")
  }
  if (length(chrom) != length(y)) {
    stop("chrom must have the same length as y.")
  }
  
  chrom <- as.factor(chrom)
  chr_levels <- levels(chrom)
  
  # store predictions for each row
  y_pred <- rep(NA_real_, length(y))
  
  # per-chromosome metrics container
  per_chr <- data.frame(
    chrom = chr_levels,
    cor   = NA_real_,
    R2    = NA_real_,
    stringsAsFactors = FALSE
  )
  
  for (i in seq_along(chr_levels)) {
    chr <- chr_levels[i]
    message("LOCO CV: holding out ", chr, " ...")
    
    f_test <- (chrom == chr)
    
    # train on all other chromosomes
    glmod <- glm_train(f_test = f_test, y = y, x = x)
    
    # predict on held-out chromosome
    preds_chr <- glm_predict(glmod, x[f_test, , drop = FALSE])
    preds_chr <- as.numeric(preds_chr)
    
    y_pred[f_test] <- preds_chr
    
    # per-chromosome metrics
    y_chr <- y[f_test]

    # Pearson correlation
    cor_chr <- suppressWarnings(
      cor(y_chr, preds_chr, method = "pearson")
    )
    
    # R^2 as squared Pearson correlation
    R2_chr <- cor_chr^2
    
    per_chr$cor[i] <- cor_chr
    per_chr$R2[i]  <- R2_chr
  }
  
  # overall metrics (again via Pearson correlation)
  cor_overall <- suppressWarnings(
    cor(y, y_pred, method = "pearson")
  )
  R2_overall <- cor_overall^2
  
  list(
    y_true    = y,
    y_pred    = y_pred,
    overall   = list(
      cor = cor_overall,
      R2  = R2_overall
    ),
    per_chrom = per_chr
  )
}
generate_loco_models = function(m_en1_nr=m_en1_nr,feats_en1_nr=feats_en1_nr,cgd = cgd){
mod$seqmod_loc_feats_qpmax = m_en1_nr

mod$seqmod_loc_feat = feats_en1_nr

lk4 = log2(6 + cgd$eb4_k4_mean)

lk27 = cgd$eb4_k27_mean

resp = as.vector(lk4 - lk27)

test_lm = cgd$chrom %in% borzoi_chrs

fmat_all = as.matrix(cbind(cgd$l,
                           log2(cgd$l),
                           cgd$cg_max - 0.04,
                           cgd$gc_max,
                           mod$seqmod_loc_feats,
                           mod$seqmod_loc_feats_qpmax))

y = lk27

x = fmat_all

chrom = as.character(mod$cgdom_ann$chrom)

cv_res <- loco_cv(y = y, x = x, chrom = chrom)

cv_res$overall

cv_res$per_chrom

# for TrxG

y = lk4

x = fmat_all

chrom = as.character(mod$cgdom_ann$chrom)

cv_res_k4 <- loco_cv(y = y, x = x, chrom = chrom)

cv_res_k4$overall

cv_res_k4$per_chrom
if (!file.exists("./data/lm_cv_loco_k4_fig1.rds")) {
  saveRDS(cv_res_k4, "./data/lm_cv_loco_k4_fig1.rds")
}

if (!file.exists("./data/lm_cv_loco_fig1.rds")) {
  saveRDS(cv_res, "./data/lm_cv_loco_fig1.rds")
}
return(list(loco_pcg = cv_res, loco_trxg = cv_res_k4))
}

glm_train = function(f_test,y , x ){
    set.seed(42)
    # choose cores
    ncores <- max(1, parallel::detectCores() - 10)
    cl <- makeCluster(ncores)
    doParallel::registerDoParallel(cl)

    f = !f_test
    #f = train_bool
    #k27
   
    glmod = cv.glmnet(
      x = x[f,,drop=FALSE],
      y = as.vector(y)[f],
      parallel = TRUE
    )
    stopCluster(cl)
       
    
    return(glmod) 
    
}

glm_predict = function(glmod,x){
    glpred = predict.glmnet(glmod$glmnet.fit, x)
    res = glpred[, which(glmod$lambda==glmod$lambda.1se)]
    return(res)
    }
cluster_pwm_hits <- function(raw_pwm, q_energ = 0.95, k = 20) {

  pwm_mat = raw_pwm[, grepl('_center$',colnames(raw_pwm)) | grepl('_left1$',colnames(raw_pwm)) | grepl('_right1$',colnames(raw_pwm))]

  pwm_mat1 = raw_pwm[, grepl('_center$',colnames(raw_pwm))]
  pwm_mat2 = raw_pwm[, grepl('_left1$',colnames(raw_pwm))]
  pwm_mat3 = raw_pwm[, grepl('_right1$',colnames(raw_pwm))]
  colnames(pwm_mat1) = gsub('_center','',colnames(pwm_mat1))
  colnames(pwm_mat2) = gsub('_left1','',colnames(pwm_mat2))
  colnames(pwm_mat3) = gsub('_right1','',colnames(pwm_mat3))

  a = rbind(pwm_mat1,pwm_mat2,pwm_mat3)
  a <- apply(a, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
  })

  energy_95 = apply(a, 2, quantile, q_energ)

  ehit = t(t(a) > energy_95)

  n = nrow(ehit)

  d2 = t(ehit) %*% (ehit)

  nhit = n * (1 - q_energ)
  hc = hclust(as.dist(nhit - d2), "complete")

  #cl = cutree(hc, h=nhit*0.6)
  cl = cutree(hc, k = k)

  return(list(hc = hc, cl = cl, d2=d2,ehit=ehit,nhit=nhit))
}
save_rds_if_missing <- function(x, path) {
  if (!is.null(path) && !file.exists(path)) {
    saveRDS(x, path)
  }
}

prepare_model_features <- function(
    mod,
    dinucs_path = "./data/seqmod_loc_feats_dinucs_fig1_lm.rds",
    motifs20_path = "./data/motifs20_fig1.rds",
    motif_mode = c("all", "motif20"),
    feature_mode = c("all_bins", "center_only")
) {
  motif_mode <- match.arg(motif_mode)
  feature_mode <- match.arg(feature_mode)

  dinucs <- readRDS(dinucs_path)
  mod$seqmod_loc_feats_raw <- readRDS("./data/seqmod_loc_feats_raw_fig1_lm.rds")
  cgd <- mod$cgdom_ann
  f_x <- cgd$chrom == "chrX"

  raw_feats <- mod$seqmod_loc_feats_raw

  if (motif_mode == "motif20") {
    motifs20 <- readRDS(motifs20_path)
    cn <- colnames(raw_feats)
    pat <- paste0("^(", paste(motifs20, collapse = "|"), ")_(center|left1|left2|right1|right2)$")
    motif20_cols <- grep(pat, cn, value = TRUE)
    raw_feats <- raw_feats[, motif20_cols, drop = FALSE]
  }

  energ_norm <- iceqream::norm_energy_matrix(
    x = raw_feats,
    dataset_x = raw_feats,
    min_energy = -7,
    q = 0.995,
    norm_energy_max = 10
  )
  logist_fmat <- iceqream::create_logist_features(energ_norm)

  length_bin <- cut(cgd$l, c(-1, 400, 600, 1000, 2000, 4000, 1e6))
  feats_len <- data.frame(id = seq_along(length_bin))
  for (b in unique(length_bin)) {
    feats_len <- cbind(feats_len, as.numeric(length_bin == b))
  }
  feats_len <- feats_len[, -1, drop = FALSE]
  colnames(feats_len) <- c("l400", "l600", "l1k", "l2k", "l4k", "llong")

  if (feature_mode == "center_only") {
    logist_fmat <- logist_fmat[, grepl("_center", colnames(logist_fmat)), drop = FALSE]
    dinucs <- dinucs[, grepl("_center$", colnames(dinucs)), drop = FALSE]
  }

  list(
    cgd = cgd,
    f_x = f_x,
    logist_fmat = logist_fmat,
    dinucs = dinucs,
    feats_len = feats_len
  )
}

fit_model_and_plot <- function(
    prep,
    use_dinucs_logist = FALSE,
    use_dinucs_sq = FALSE,
    use_feats_len = FALSE,
    seed = 42,
    test_frac = 0.1,
    save_path = NULL
) {
  cgd <- prep$cgd
  f_x <- prep$f_x
  cgd_noX <- cgd[!f_x, , drop = FALSE]

  set.seed(seed)
  test <- rep(FALSE, nrow(cgd_noX))
  test[sample(nrow(cgd_noX), floor(test_frac * nrow(cgd_noX)))] <- TRUE

  lk4 <- log2(6 + cgd_noX$eb4_k4_mean)
  lk27 <- cgd_noX$eb4_k27_mean
  resp <- as.vector(lk4 - lk27)
  test_lm <- test

  feature_blocks <- list(
    cbind(
      cgd_noX$l,
      log2(cgd_noX$l),
      cgd_noX$cg_max - 0.04,
      cgd_noX$gc_max
    ),
    prep$logist_fmat[!f_x, , drop = FALSE],
    prep$dinucs[!f_x, , drop = FALSE]
  )

  if (use_dinucs_logist) {
    dinucs_norm <- iceqream::norm01(prep$dinucs) * 10
    dinucs_logist <- iceqream::create_logist_features(dinucs_norm)
    feature_blocks <- c(feature_blocks, list(dinucs_logist[!f_x, , drop = FALSE]))
  }

  if (use_dinucs_sq) {
    feature_blocks <- c(feature_blocks, list((prep$dinucs[!f_x, , drop = FALSE])^2))
  }

  if (use_feats_len) {
    feature_blocks <- c(feature_blocks, list(prep$feats_len[!f_x, , drop = FALSE]))
  }

  fmat_all <- as.matrix(do.call(cbind, feature_blocks))

  mod_fmat_all_k27 <- glm_train(f_test = test_lm, y = lk27, x = fmat_all)
  p_fmat_all_k27 <- glm_predict(glmod = mod_fmat_all_k27, x = fmat_all)

  mod_fmat_all_k4 <- glm_train(f_test = test_lm, y = lk4, x = fmat_all)
  p_fmat_all_k4 <- glm_predict(glmod = mod_fmat_all_k4, x = fmat_all)

  pred <- p_fmat_all_k4 - p_fmat_all_k27

  message(cor(pred[test_lm], resp[test_lm], m = "s")^2)
  message(cor(pred[!test_lm], resp[!test_lm], m = "s")^2)
  message(cor(pred[test_lm], resp[test_lm], m = "p")^2)
  message(cor(pred[!test_lm], resp[!test_lm], m = "p")^2)

  cgd_noX$k27_pred <- p_fmat_all_k27
  cgd_noX$k4_pred <- p_fmat_all_k4

  cgd_noX$modality <- ifelse(cgd_noX$k27_pred > 6.4, "pcg", "mix")
  cgd_noX$modality <- ifelse(
    cgd_noX$k4_pred >= 6 & cgd_noX$k27_pred < 5.3,
    "txg",
    cgd_noX$modality
  )

  print(table(cgd_noX$modality))

  sps(8, 8)
  p1 <- cgd_noX %>%
    ggplot(aes(x = k4_pred, y = k27_pred)) +
    geom_dense_scatter() +
    coord_cartesian(xlim = c(3, 8.5), ylim = c(4, 7)) +
    theme_bw() +
    geom_hline(yintercept = 6) +
    geom_vline(xintercept = 6) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 28, face = "bold"),
      axis.text.y = element_text(size = 28, face = "bold"),
      axis.title.x = element_text(size = 28, face = "bold"),
      axis.title.y = element_text(size = 28, face = "bold")
    )
  print(p1)

  p2 <- cgd_noX %>%
    ggplot(aes(x = k4_pred, y = k27_pred, col = modality)) +
    geom_point(size = .5, alpha = .7) +
    coord_cartesian(xlim = c(3, 8.5), ylim = c(4, 7)) +
    theme_bw() +
    geom_hline(yintercept = 6) +
    geom_vline(xintercept = 6) +
    scale_color_manual(values = c("black", "darkblue", "darkred")) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 28, face = "bold"),
      axis.text.y = element_text(size = 28, face = "bold"),
      axis.title.x = element_text(size = 28, face = "bold"),
      axis.title.y = element_text(size = 28, face = "bold")
    )
  print(p2)

  cgd_noX$f_test <- test

  save_rds_if_missing(cgd_noX, save_path)

  list(
    cgd_noX = cgd_noX,
    fmat_all = fmat_all,
    mod_k27 = mod_fmat_all_k27,
    mod_k4 = mod_fmat_all_k4,
    pred = pred,
    plot_dense = p1,
    plot_modality = p2
  )
}
get_pmw_min_pv_per_cl = function(hc_list){
set.seed(42)
##load the pwms and select the most enriched, 458 motifs
	mod$mot_stats[["cgd_inner"]] = readRDS("data/mot_stats_cgd_inner.RDS")
	mstm = mod$mot_stats[["cgd_inner"]]
	mstm$dist_bin = mstm$dist_bin + 20 - 11
	mstm = mstm[!is.na(mstm$pv),]
	mstm = mstm[mstm$pv < 1e-2,]
	mstm_dom = mstm[mstm$dist_bin >= -5 & mstm$dist_bin <= 5,]
	mstm_3 = mstm[mstm$dist_bin > 5,] 
	mstm_5 = mstm[mstm$dist_bin < -5,] 
	mot_n = table(mstm_dom$motif)
	mot_n3 = table(mstm_3$motif)
	mot_n5 = table(mstm_5$motif)
	tfs = names(mot_n[mot_n>3])
	tfs5 = names(mot_n5[mot_n5>3])
	tfs3 = names(mot_n3[mot_n3>3])
	tfs_a = unique(c(tfs,tfs5,tfs3))
	message('n tfs = ', length(tfs_a))


  ann <- data.frame(cluster = factor(hc_list$cl))
  rownames(ann) <- names(hc_list$cl)
  cl_order = unique(hc_list$cl[hc_list$hc$order])
min_pv_per_motif <- mstm %>% filter(motif %in% tfs_a) %>%
  group_by(motif) %>%
  slice_min(order_by = pv, n = 1, with_ties = FALSE) %>%
  ungroup()


length(intersect(rownames(ann),min_pv_per_motif$motif))

min_pv_per_motif = min_pv_per_motif %>%
  mutate(
    motif = gsub("\\.", "_", motif),
    motif = gsub("/", "_", motif),
    motif = gsub("-", "_", motif),
    motif = gsub("::", "_", motif)
  )


pv_cl_df = ann %>% rownames_to_column(var ='motif') %>% left_join(.,min_pv_per_motif,by = 'motif')

#library(dplyr)

min_pv_per_cluster <-pv_cl_df %>%
  group_by(cluster) %>%
  slice_min(order_by = pv, n = 1, with_ties = FALSE) %>%
  ungroup()
min_pv_per_cluster = min_pv_per_cluster %>%
  mutate(cluster = factor(cluster, levels = cl_order)) %>%
  arrange(cluster)
return(min_pv_per_cluster)
}

library("xgboost")
norm_CnT_loess = function (track, atac_nm, winsize, out_control_track = NULL) 
{
    set.seed(42)
    gvtrack.create("tr_win", track, "avg")
    gvtrack.iterator("tr_win", sshift = -winsize, eshift = winsize)
    gvtrack.create("atac_win", atac_nm, "avg")
    gvtrack.iterator("atac_win", sshift = -winsize, eshift = winsize)
    gvtrack.create("tr_2000", track, "avg")
    gvtrack.iterator("tr_2000", sshift = -1000, eshift = 1000)
    quant = gquantiles("tr_2000", c(0.98, 0.985))
    qk27 = as.numeric(quant[1])
    k27d = gscreen(paste("tr_2000 > ", qk27, sep = ""))
    pcg_d = gextract(c("tr_2000", "atac_win"), iterator = k27d, 
        intervals = k27d)
    pcg_d$l = pcg_d$end - pcg_d$start
    qk27_2 = as.numeric(quant[2])
    strict_pcg_d = pcg_d[pcg_d$tr_2000 > qk27_2, ]
    print("yalla")

	atac_q9995 = gquantiles("atac_win", 0.9995)
	max_atac = floor(2*atac_q9995)	

    a = gextract(c("seq.GC500_bin20"), iterator = 1e+06, intervals = gintervals.all())
    samp = a[a$seq.GC500_bin20 > 0.1 & runif(n = nrow(a)) < 0.05, ]
    profs = gextract(c("tr_win", "atac_win"), iterator = 20, intervals = samp)
    profs$chrom = as.character(profs$chrom)
    profs = profs[!profs$chrom %in% c("chrY", "chrM"), ]
    profs[is.na(profs)] = 0
    profs = profs %>% arrange(chrom, start)
    profs_d = gintervals.neighbors(profs, strict_pcg_d)
    profs_d = profs_d[!is.na(profs_d$dist), ]
    f_out = abs(profs_d$dist) > 1000
	
    tst_eb = tapply(profs$tr_win[f_out], pmin(floor(profs$atac_win[f_out] * 2), max_atac), mean)
    tst_all = tapply(profs$tr_win[f_out], floor(profs$atac_win[f_out] * 2), mean)
    ebin1 = min(floor(quantile(profs$atac_win, 0.9992) * 2), max_atac)
    ebin = min(floor(quantile(profs$atac_win, 0.999) * 2), max_atac)
    sbin = min(floor(quantile(profs$atac_win, 0.998) * 2), max_atac)
    f = loess(y ~ x, data = data.frame(x = 1:ebin1, y = tst_eb[1:ebin1]), 
        span = 0.8, control = loess.control(surface = "direct"))
    max_bin = 1000
    slope = (f$fitted[ebin] - f$fitted[sbin])/(ebin - sbin)
    norm_trend = c(f$fitted[1:ebin], f$fitted[ebin] + ((1 + (ebin:max_bin)) - 
        ebin) * slope)
    names(norm_trend) = as.character(1:length(norm_trend))
    plot(norm_trend)
    points(tst_eb, col = "blue")
    plot(tst_eb, col = "blue")
    points(norm_trend)
    plot(tst_all, col = "blue")
    points(norm_trend)
    plot(f$fitted)


    #norm_exp_simple = "pmax(ifelse(seq.GC500_bin20==0, 0, tr_win - norm_trend[1+pmin(floor(atac_win*2),max_bin)]),0)"


    #norm_exp = "ifelse(seq.GC500_bin20==0, 0, ifelse(tr_win - as.numeric(norm_trend[names(norm_trend)%in% as.character(1+pmin(floor(ifelse(is.na(atac_win),0,atac_win))*2,max_bin))])<0,0,tr_win - #as.numeric(norm_trend[names(norm_trend)%in% as.character(1+pmin(floor(ifelse(is.na(atac_win),0,atac_win))*2,max_bin))])) )"
	
	norm_exp_simple = "pmax(
	ifelse(
		seq.GC500_bin20 == 0,
		0,
		tr_win - norm_trend[1 + pmin(floor(ifelse(is.na(atac_win), 0, atac_win) * 2), max_bin)]
	),
	0
	)"
    if (!is.null(out_control_track)) {
        gtrack.create(track = out_control_track, description = "normalized using ATAC table", 
            expr = norm_exp_simple, iterator = 20)
    }
    return(list(tst_eb = tst_eb, norm_trend = norm_trend, f = f, 
        tst_all = tst_all))
}

norm_CnT_loess_k4 = function(track, atac_nm, winsize, out_control_track = NULL){
    
	set.seed(42)
	gvtrack.create("tr_win", track, "avg")
	gvtrack.iterator("tr_win", sshift=-winsize, eshift=winsize)
	gvtrack.create("atac_win", atac_nm, "avg")
	gvtrack.iterator("atac_win", sshift=-winsize, eshift=winsize)
	gvtrack.create("tr_2000", track, "avg")
	gvtrack.iterator("tr_2000", sshift=-1000, eshift=1000)		
	tss = readRDS('./data/tss_ref.rds')
	strict_txg_d = tss
	
	print("yalla")
	atac_q9995 = gquantiles("atac_win", 0.9995)
	max_atac = floor(2*atac_q9995)
	
	a = gextract(c("seq.GC500_bin20"), iterator=1e+6, intervals=gintervals.all())
	samp = a[a$seq.GC500_bin20 > 0.1 & runif(n=nrow(a))<0.05,]
	profs = gextract(c("tr_win", "atac_win"), iterator=20, intervals=samp)   
	profs$chrom = as.character(profs$chrom)
	profs = profs[ !profs$chrom%in% c('chrY','chrM'),]##
    profs[is.na(profs)] = 0
	profs = profs %>% arrange(chrom,start)	 
    profs_d = gintervals.neighbors(profs, strict_txg_d)
	profs_d = profs_d[!is.na(profs_d$dist),]
	f_out = abs(profs_d$dist) > 1000
	
	tst_eb = tapply(profs$tr_win[f_out], pmin(floor(profs$atac_win[f_out] * 2), max_atac), mean)
    tst_all = tapply(profs$tr_win[f_out], floor(profs$atac_win[f_out] * 2), mean)
    ebin1 = min(floor(quantile(profs$atac_win, 0.9992) * 2), max_atac)
    ebin = min(floor(quantile(profs$atac_win, 0.999) * 2), max_atac)
    sbin = min(floor(quantile(profs$atac_win, 0.998) * 2), max_atac)
    f = loess(y ~ x, data = data.frame(x = 1:ebin1, y = tst_eb[1:ebin1]), 
        span = 0.8, control = loess.control(surface = "direct"))
    max_bin = 1000
    slope = (f$fitted[ebin] - f$fitted[sbin])/(ebin - sbin)
    norm_trend = c(f$fitted[1:ebin], f$fitted[ebin] + ((1 + (ebin:max_bin)) - 
        ebin) * slope)
    names(norm_trend) = as.character(1:length(norm_trend))
    plot(norm_trend)
    points(tst_eb, col = "blue")
    plot(tst_eb, col = "blue")
    points(norm_trend)
    plot(tst_all, col = "blue")
    points(norm_trend)
    plot(f$fitted)
	norm_exp_simple = "pmax(
	ifelse(
		seq.GC500_bin20 == 0,
		0,
		tr_win - norm_trend[1 + pmin(floor(ifelse(is.na(atac_win), 0, atac_win) * 2), max_bin)]
	),
	0
	)"
	if(!is.null(out_control_track)) {
		gtrack.create(track=out_control_track,
					description="normalized using ATAC table", 
					expr=norm_exp_simple, iterator=20)}
	return(list(tst_eb=tst_eb,norm_trend=norm_trend,f=f,tst_all = tst_all))
         }
init_pipe = function(override = F)
{	
	if(!override & file.exists("data/cache_mod.RDS")) {
	mod = readRDS("data/cache_mod.RDS")
	compute_track_quantiles(mod = mod)
	message("done track quant")
	compute_cnt_quantiles(mod = mod)
	message("done cnt track quant")
	mod = pcg_init_epi_track_lib(mod) 
	message("done lib")
	
	return(mod)
	}else{ 
	message("data/cache_mod.RDS doesn't exist. Add it to the directory or use init_pipe_new()")
	
	}
}

init_pipe_new = function(override = F)
{	
	if(!override & file.exists("./data/cache_mod_new.RDS")) {
		return(readRDS("./data/cache_mod_new.RDS"))
	}else{ 
	mod = pcg_init() 
	message("done base")
	compute_track_quantiles(mod = mod)
	message("done track quant")
	compute_cnt_quantiles(mod = mod)
	message("done cnt track quant")
	mod = pcg_init_epi_track_lib(mod) 
	message("done lib")
	#mod = pcg_update_track_q_thresh(mod)
	#mod = pcg_compute_k27_doms(mod) 
	#message("done comp k27 dom")
	#mod = pcg_gen_eb4_doms(mod) 
	#message("done comp eb4 doms")
	#mod = pcg_gen_scale50k_k27(mod=mod) 
	#message("done comp scale50")
	#mod = pcg_gen_doms_interv_tiling(mod) 
	#message("done gen dom tilings")
	#mod = pcg_gen_tss_rna_k4k27_annots(mod)
	#message("done gen tss rna annots")
	#mod = pcg_load_mot_screens(mod)
	mod = pcg_gen_cg_h_doms(mod)
	message("generated CGDDs")
	mod = pcg_gen_cgd_tiling(mod)
	message("done CGDD tilling")
	mod$cgdom_ann$ID = 1:nrow(mod$cgdom_ann)
	mod$test_chroms = c("chr2", "chr8", "chr12", "chr18")
	message("saving")
	saveRDS(mod, file="./data/cache_mod_new.RDS")
	return(mod)
	
	}
}


pcg_init = function()
{
	mod = list()

	mod$tss_jk = readRDS("data/tss_jk.RDS")
	#mod$tss = gintervals.load("intervs.global.tss")
	mod$tss = gintervals.load("intervs.global.tss_10x")
	mod$exon = gintervals.load("intervs.global.exon")
	options(gmax.data.size=1e+9)
	options(scipen=999)

	if(!file.exists("data/multi/npeaks.RDS")) {
		mod$raw_npeaks = join_norm_peaks()
		mod$npeaks = mod$raw_npeaks
		gen_track_vext()
		for(i in 1:4) {	
			mod$npeaks = center_and_union_intervs(mod$npeaks, 
								"atac_ext", 140,min_space=100)
		}
		mod = add_peak_cggc(mod)
		saveRDS(mod$raw_npeaks, file="data/multi/npeaks.RDS")
		saveRDS(mod$npeaks, file="data/multi/npeaks_canonic.RDS")
	} else {
		mod$raw_npeaks = readRDS(file="data/multi/npeaks.RDS")
		mod$npeaks = readRDS(file="data/multi/npeaks_canonic.RDS")
	}
	mod$npeaks_tss = mod$npeaks[mod$npeaks$d_tss==0,]
	if(!file.exists("data/multi/multi_type_cov.RDS")) {
		mod  = collect_cov_r(mod)
	} else {
		mod$cov_stat = readRDS(file="data/multi/multi_type_cov.RDS")
		rownames(mod$cov_stat) = mod$cov_stat$nm
		mod$cov_peak_chr = readRDS(file="data/multi/cov_peak_chr.RDS")
		mod$cov_chr = readRDS(file="data/multi/cov_chr.RDS")
	}

	#mod = gen_cgi(mod)
	#mod = gen_ctcf(mod)
	#mod = init_exp(mod)

	mod$emb_type_legc =  readRDS(file="data/spatwt_type_legc.RDS")
 	mod$emb_mc_legc = readRDS(file="data/spatwt_mc_legc.RDS")
	mod$mat_eb_rna = readRDS("data/L67_bulk_per_cond.RDS")
	eb_umi = mod$mat_eb_rna[,"wt_t_wt_d3"]
	eb_umi_n = eb_umi/sum(eb_umi)
	mod$eb_legc = log2(1e-5+eb_umi_n)

	# message("building ehs modes")
	#if(!file.exists("data/epimod_ehs.RDS")) {
	#	mod = build_epi_hotspots_classes(mod)
	#	saveRDS(mod$ehs, "data/epimod_ehs.RDS")
	#} else {
	#	mod$ehs = readRDS("data/epimod_ehs.RDS")
	#}

#	if(file.exists("data/peak_stat.RDS")) {
#		mod$peak_stat = readRDS(file="data/peak_stat.RDS")
#	} else {
#		peak_stat = gen_peaks_z_mat(mod, npeaks)
#		mod$peak_stat = peak_stat
#		saveRDS(peak_stat, file="data/peak_stat.RDS")
#	}
	return(mod)
}


gen_track_vext = function(tname = "marginal", w_ext = 140)
{
	message("set vt for ", tname)
	nm = sprintf("jk.epipcg.multieb.%s", tname)
	gvtrack.create("atac_ext", nm, "sum")
	gvtrack.iterator("atac_ext", sshift = -w_ext, eshift = w_ext)
}
compute_track_quantiles = function(mod)
{
	ndx = as.data.frame(fread("./data/index_tracks.txt"))	
	
	for(i in 1:nrow(ndx)) {
#		message("running ", ndx$short_name[i])
		gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "avg") 
		gvtrack.iterator(ndx$short_name[i], sshift = -500, eshift= 500)
	}
	if(file.exists("./data/track_thresh.txt")) {
		mod$k27_track_thresh_cnt = read.table("./data/track_thresh.txt", sep="\t", stringsAsFactors=F,header=T)
		th_mat = mod$k27_track_thresh_cnt
		mod$k4_track_thresh_cnt = read.table("./data/track_thresh_k4.txt", sep="\t", stringsAsFactors=F,header=T)
	} else {
		thresh = list()
		for(i in 1:nrow(ndx)) {
			thresh[[i]] = gquantiles(ndx$short_name[i], c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
		}
		th_mat = do.call('rbind',thresh)
		rownames(th_mat) = ndx$short_name
		write.table(th_mat, file="./data/track_thresh.txt", quote=F, sep="\t")
		mod$k27_track_thresh_cnt = th_mat
		thresh_k4 = list()
		tnm_k4 = gen_k4_vt_cnt(mod)
		for(tnm in tnm_k4) {
			thresh_k4[[tnm]] = gquantiles(tnm, c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
		}
		th_mat_k4 = do.call('rbind',thresh_k4)
		rownames(th_mat_k4) = tnm_k4
		write.table(th_mat_k4, file="./data/track_thresh_k4.txt", quote=F, sep="\t")
		mod$k4_track_thresh_cnt = th_mat_k4
	}
	#return(mod)
}

compute_cnt_quantiles = function(mod)
{
	ndx = as.data.frame(fread("./data/index_tracks_cnt.txt"))	
	
	for(i in 1:nrow(ndx)) {
#		message("running ", ndx$short_name[i])
		gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "avg") 
		gvtrack.iterator(ndx$short_name[i], sshift = -500, eshift= 500)
	}
	if(file.exists("./data/track_thresh_cnt.txt")) {
		mod$k27_track_thresh_cnt = read.table("./data/track_thresh_cnt.txt", sep="\t", stringsAsFactors=F,header=T)
		th_mat = mod$k27_track_thresh_cnt
		mod$k4_track_thresh_cnt = read.table("./data/track_thresh_k4_cnt.txt", sep="\t", stringsAsFactors=F,header=T)
	} else {
		thresh = list()
		for(i in 1:nrow(ndx)) {
			thresh[[i]] = gquantiles(ndx$short_name[i], c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
		}
		th_mat = do.call('rbind',thresh)
		rownames(th_mat) = ndx$short_name
		write.table(th_mat, file="./data/track_thresh_cnt.txt", quote=F, sep="\t")
		mod$k27_track_thresh_cnt = th_mat
		thresh_k4 = list()
		tnm_k4 = gen_k4_vt_cnt(mod)
		for(tnm in tnm_k4) {
			thresh_k4[[tnm]] = gquantiles(tnm, c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
		}
		th_mat_k4 = do.call('rbind',thresh_k4)
		rownames(th_mat_k4) = tnm_k4
		write.table(th_mat_k4, file="./data/track_thresh_k4_cnt.txt", quote=F, sep="\t")
		mod$k4_track_thresh_cnt = th_mat_k4
	}
	#return(mod)
}
gen_k4_vt_cnt = function(mod) 
{
	ndx = as.data.frame(fread("./data/index_tracks_cnt.txt"))
	for(i in 1:nrow(ndx)) {
		if(!is.na(ndx$track_k4[i])) {
			gvtrack.create(ndx$short_name_k4[i], ndx$track_k4[i], "sum") 
			gvtrack.iterator(ndx$short_name_k4[i], sshift = -140, eshift= 140)
		}
	}
	k4_tns =  ndx$short_name_k4[!is.na(ndx$track_k4)]
	
	
	return(k4_tns)
}
compute_tot_cnt_cov_vs_hcg = function(mod){
ndx = as.data.frame(fread("./data/index_tracks_cnt.txt"))		
	for(i in 1:nrow(ndx)) {
#		message("running ", ndx$short_name[i])
		gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "sum") 
		
        }

	if(file.exists("./data/track_cov_cnt.txt")) {
		mod$cnt_track_cov = read.table("data/track_cov_cnt.txt", sep="\t", stringsAsFactors=F,header=T)
		

	} else {
		thresh = list()
		for(i in 1:nrow(ndx)) {
			thresh[[i]] = gsummary(ndx$short_name[i])
		}
		th_mat = do.call('rbind',thresh)
		rownames(th_mat) = ndx$short_name
        cg = mod$cgdom_ann

        cg = cg[ cg$type%in% c('cl1','cl2'),]



        gext = gextract(ndx$short_name,intervals = cg,iterator = cg,colnames = ndx$short_name)

        gext[is.na(gext)] = 0


        
        th_mat = cbind(th_mat,as.numeric(colSums(gext[,-c(1,2,3,ncol(gext))])))
		write.table(th_mat, file="./data/track_cov_cnt.txt", quote=F, sep="\t")
		mod$cnt_track_cov = th_mat       
	}
 return(mod)
}
save_gw_to_misha = function(mod,track_nm = "jk.epipcg.pred.eb4_xgb_lm_10test_lg_seeds_gc_cg_noX_final"){
##save the genome wide prediction as misha track
mod$gw$cg_trace$pred = mod$gw$xg_pred$base2iqdn

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
cg_trace_f = cg_trace[ f_norp,]

data_tr = cg_trace_f

data_tr[is.na(data_tr)] = 0


data_tr$end = data_tr$end - 1

gtrack.create_sparse(track = track_nm,
                     description ='eb4_xgb_lm_seeds_on_normed_cnt_gc_cg_feats_10test_noX' ,
                     intervals = data_tr,values = data_tr$pred)
}
compute_tot_cnt_k4_cov_vs_hcg = function(mod){
ndx = as.data.frame(fread("./data/index_tracks_cnt.txt"))
ndx = ndx[!is.na(ndx$track_k4),]	
	for(i in 1:nrow(ndx)) {
#		message("running ", ndx$short_name[i])
		gvtrack.create(ndx$short_name_k4[i], ndx$track_k4[i], "sum") 
		
        }

	if(file.exists("./data/track_cov_cnt_k4.txt")) {
		mod$cnt_track_cov_k4 = read.table("data/track_cov_cnt_k4.txt", sep="\t", stringsAsFactors=F,header=T,row.names = 1)###be aware colnames shifted
		colnames(mod$cnt_track_cov_k4) = c('Total.intervals','NaN.intervals','Min','Max','Sum','Mean','Std.dev','X')

	} else {
		thresh = list()
		for(i in 1:nrow(ndx)) {
			thresh[[i]] = gsummary(ndx$short_name_k4[i])
		}
		th_mat = do.call('rbind',thresh)
		rownames(th_mat) = ndx$short_name_k4
        cg = mod$cgdom_ann

        cg = cg[ cg$type%in% c('cl4','cl3'),]



        gext = gextract(ndx$short_name_k4,intervals = cg,iterator = cg,colnames = ndx$short_name_k4)

        gext[is.na(gext)] = 0


        
        th_mat = cbind(th_mat,as.numeric(colSums(gext[,-c(1,2,3,ncol(gext))])))
		write.table(th_mat, file="./data/track_cov_cnt_k4.txt", quote=F, sep="\t")
		mod$cnt_track_cov_k4 = th_mat       
	}
 return(mod)
}

compute_tot_atac_cov_vs_hcg = function(mod){
ndx = mod$track_atac_ribo	
	for(i in 1:nrow(ndx)) {
#		message("running ", ndx$short_name[i])
		gvtrack.create(ndx$short_name[i], ndx$tr_atac[i], "sum") 
		
        }

	if(file.exists("./data/track_cov_atac.txt")) {
		mod$atac_track_cov = read.table("./data/track_cov_atac.txt", sep="\t", stringsAsFactors=F,header=T)
		

	} else {
		thresh = list()
		for(i in 1:nrow(ndx)) {
			thresh[[i]] = gsummary(ndx$short_name[i])
		}
		th_mat = do.call('rbind',thresh)
		rownames(th_mat) = ndx$short_name
        cg = mod$cg_high_int

        cg = cg[ cg$type%in% c('cl4','cl3'),]



        gext = gextract(ndx$short_name,intervals = cg,iterator = cg,colnames = ndx$short_name)

        gext[is.na(gext)] = 0


        
        th_mat = cbind(th_mat,as.numeric(colSums(gext[,-c(1,2,3,ncol(gext))])))
		write.table(th_mat, file="./data/track_cov_atac.txt", quote=F, sep="\t")
		mod$atac_track_cov = th_mat       
	}
 return(mod)
}




pcg_gen_cgd_tiling = function(mod)
{
	cgdoms = mod$cgdom_ann[,c("chrom","start","end","l","tss_start", "tss_strand")]
	cgdoms$id = 1:nrow(cgdoms)
	cgdoms$strand = cgdoms$tss_strand
	cgdoms300 = cgdoms
	cgdoms300$start = cgdoms$start - 600
	cgdoms300$end = cgdoms$end + 600
	cgdom_i = giterator.intervals(expr = c("1"), iterator=100, intervals=cgdoms300)
	cgdom_i = gintervals.neighbors(cgdom_i, cgdoms)
	near_tss = gintervals(cgdoms$chrom, cgdoms$tss_start, cgdoms$tss_start+20, cgdoms$tss_strand)
	gvtrack.create("dom_d", mod$cgdom_ann, "distance")
	gvtrack.create("tss_d", near_tss, "distance")
	cgdom_ds = gextract(c("tss_d", "dom_d", "seq.GC20_bin20", "seq.AA_bin20", "seq.TT_bin20"), iterator=cgdom_i, intervals=cgdom_i, colnames=c("tss_d", "dom_d", "GC", "AA", "TT"))

	cgdom_tile = cbind(cgdom_i[,c("chrom","start","end","strand","l","id")], cgdom_ds[,c("tss_d","dom_d","GC", "AA","TT")])

#this should have be done by the vtrack - didn't it?
	cgdom_tile$tss_d = ifelse(cgdom_tile$strand==1, cgdom_tile$tss_d, -cgdom_tile$tss_d)

	gvtrack.create("inner_d", cgdoms, "distance", 1000)
	tile_in_dist = gextract(c("inner_d"), iterator=cgdom_tile, intervals=cgdom_tile)
	
	cgdom_tile$inner_dom_d = round(tile_in_dist$inner_d)
	
	mod$cgd_tile = cgdom_tile
	return(mod)
}

pcg_screen_motifs = function(mod, doms, foc_tiles, dist_bins, label)
{
	set.seed(42)
	load("data/motif_db.rda")
	motif_db_s <- motif_db %>% filter(dataset %in% c("HOMER", "JASPAR", "JOLMA"))

	mot_hits = c()

#work in groups of k mots?
	tf_mots = unique(motif_db_s$motif)

	n_gcbin = 9
	foc_tiles$GC[is.na(foc_tiles$GC)]=0
	foc_tiles$gcb = pmin(pmax(floor((foc_tiles$GC/20-0.35)/0.05),0),n_gcbin-1)
	bins = foc_tiles$gcb + n_gcbin * dist_bins
	
	k27= doms[foc_tiles$id,"eb4_k27_mean"]
	k4= log2(6+doms[foc_tiles$id,"eb4_k4_mean"])
	resp = k4 - k27

	max_sbin = max(bins)

	all_stats = NULL

	for(i in seq(1, length(tf_mots), 100)) {
		message("tf group ", i)
		mot_set = tf_mots[i:min(i+99, length(tf_mots))]
		en = prego::gextract_pwm(intervals=foc_tiles, dataset=motif_db_s[motif_db_s$motif %in% mot_set,])
		eq95 = apply(en[,-(1:13)], 2, quantile, 0.95)

		for(sbin in 1:max_sbin) {
			if(sbin/40 == round(sbin/40)) {
				message("sbin ", sbin)
			}
			fbin = (bins == sbin)
			resp_b = resp[fbin]
			hits = t(t(en[fbin,-(1:13)])>eq95)
			stats = apply(hits, 2, function(h) { 
					if(sum(h) > 20 & sum(!h) > 20) {
						ks = ks.test(resp_b[h], resp_b[!h])
						return(c(sum(h), ks$p.value, ks$statistic, mean(resp_b[h]), mean(resp_b[!h])))
					} else {
						return(c(sum(h), NA, NA, NA, NA))
					}
			})
			
			if(is.null(all_stats)) {
				all_stats = cbind(colnames(stats), rep(sbin,ncol(stats)), t(stats))
			} else {
				all_stats = rbind(all_stats, cbind(colnames(stats), rep(sbin,ncol(stats)), t(stats)))
			}
		}
		saveRDS(all_stats,"data/mot_stats.RDS")
	}
	mot_stats = data.frame(motif = all_stats[,1], 
							  bin = as.numeric(all_stats[,2]),
							  gc_bin = as.numeric(all_stats[,2])%%n_gcbin,
							  dist_bin = floor(as.numeric(all_stats[,2])/n_gcbin)-20,
							  n = as.numeric(all_stats[,3]),
							  pv = as.numeric(all_stats[,4]),
							  D = as.numeric(all_stats[,5]),
							  e_hit = as.numeric(all_stats[,6]),
							  e_nhit = as.numeric(all_stats[,7]))

	mot_stats$delta = mot_stats$e_hit-mot_stats$e_nhit
	if(!is.list(mod$mot_stat)) {
		mod$mot_stat = list()
	}
	mod$mot_stats[[label]] = mot_stats
	saveRDS(mot_stats,sprintf("data/mot_stats_%s.RDS", label))
	return(mod)
}

pcg_screen_cgd_dom_mots = function(mod, dist_mod = "tss", test_chroms=mod$test_chroms)
{
	set.seed(42)
	foc_tiles = mod$cgd_tile
	doms = mod$cgdom_ann
	if(!is.null(test_chroms)) {
		foc_tiles = foc_tiles[!as.character(foc_tiles$chrom) %in% test_chroms,]
	}
	if(dist_mod == "tss") {
		foc_tiles = foc_tiles[foc_tiles$l > 500 & doms[foc_tiles$id,"tss_dist"]==0,]
		dist_bins = pmin(pmax(floor(foc_tiles$tss_d/100)+20,0),39)
	} else if(dist_mod == "inner") {
		foc_tiles = foc_tiles[foc_tiles$l > 500 & doms[foc_tiles$id,"tss_dist"]==0,]
		dist_bins = pmin(pmax(floor(foc_tiles$inner_dom_d/100)+11,0),21)

	}
	mod = pcg_screen_motifs(mod, doms, foc_tiles, dist_bins,
                        label= sprintf("cgd_%s", dist_mod))

}

pcg_build_local_seq_feats = function(mod, add_dinucs=F)
{
	set.seed(42)
	mod$mot_stats[["cgd_inner"]] = readRDS("data/mot_stats_cgd_inner.RDS")
	mstm = mod$mot_stats[["cgd_inner"]]
	mstm$dist_bin = mstm$dist_bin + 20 - 11
	mstm = mstm[!is.na(mstm$pv),]
	mstm = mstm[mstm$pv < 1e-2,]

	mstm_dom = mstm[mstm$dist_bin >= -5 & mstm$dist_bin <= 5,]
	mstm_3 = mstm[mstm$dist_bin > 5,] 
	mstm_5 = mstm[mstm$dist_bin < -5,] 

	mot_n = table(mstm_dom$motif)
	mot_n3 = table(mstm_3$motif)
	mot_n5 = table(mstm_5$motif)
	tfs = names(mot_n[mot_n>3])
	tfs5 = names(mot_n[mot_n5>3])
	tfs3 = names(mot_n[mot_n3>3])
	tfs_a = unique(c(tfs,tfs5,tfs3))
	message('n tfs = ', length(tfs_a))
	cgd = mod$cgdom_ann

	cgd = mod$cgdom_ann
	cgd$center = cgd$start + (cgd$end-cgd$start)/2
	cgd1 = cgd[,c('chrom','start','end','center')]
	cgd1$start = cgd$center - 200
	cgd1$end = cgd$center + 200
	cgd2 = cgd[,c('chrom','start','end','center')]
	cgd2$start  = cgd1$start - 401
	cgd2$end  = cgd1$start - 1
	cgd3 = cgd[,c('chrom','start','end','center')]
	cgd3$start  = cgd2$start - 401
	cgd3$end  = cgd2$start - 1
	
	
	## 1st bin to the right of the center bin (cgd1)
	cgd4 = cgd[,c('chrom','start','end','center')]
	cgd4$start  = cgd1$end + 1
	cgd4$end    = cgd1$end + 401
	
	## 2nd bin to the right
	cgd5 = cgd[,c('chrom','start','end','center')]
	cgd5$start  = cgd4$end + 1
	cgd5$end    = cgd4$end + 401

	##check that bining is okay
	n = 10
	plot(c(cgd1[n,]$start,cgd1[n,]$end),c(1,1),xlim = c(cgd1[n,]$center-1200,cgd1[n,]$center+1200))
	abline(v = cgd1[n,]$center)
	lines(c(cgd2[n,]$start,cgd2[n,]$end),c(1,1),col='darkblue')
	lines(c(cgd3[n,]$start,cgd3[n,]$end),c(1.1,1.1),col='blue')
	lines(c(cgd4[n,]$start,cgd4[n,]$end),c(1,1),col='darkred')
	lines(c(cgd5[n,]$start,cgd5[n,]$end),c(1.1,1.1),col='red')
	
	message('calculating pwm cgd1')
	## cgd1 (center)
	en1 = prego::gextract_pwm(intervals = cgd1, motifs = tfs_a)
	m_en1 = as.matrix(en1[,-(1:4)])
    m_en1 <- apply(m_en1, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })
	m_en1_n = t(t(m_en1) - apply(m_en1, 2, quantile, 0.995))
	m_en1_nr = apply(m_en1_n, 2, function(x) pmax(pmin(x, 0), -12))
	
	message('calculating pwm cgd2')	
	## cgd2 (left 1)
	en2 = prego::gextract_pwm(intervals = cgd2, motifs = tfs_a)
	m_en2 = as.matrix(en2[,-(1:4)])
    m_en2 <- apply(m_en2, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })
	m_en2_n = t(t(m_en2) - apply(m_en2, 2, quantile, 0.995))
	m_en2_nr = apply(m_en2_n, 2, function(x) pmax(pmin(x, 0), -12))
	
	message('calculating pwm cgd3')	
	## cgd3 (left 2)
	en3 = prego::gextract_pwm(intervals = cgd3, motifs = tfs_a)
	m_en3 = as.matrix(en3[,-(1:4)])
    m_en3 <- apply(m_en3, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })
	m_en3_n = t(t(m_en3) - apply(m_en3, 2, quantile, 0.995))
	m_en3_nr = apply(m_en3_n, 2, function(x) pmax(pmin(x, 0), -12))
	
	message('calculating pwm cgd4')	
	## cgd4 (right 1)
	en4 = prego::gextract_pwm(intervals = cgd4, motifs = tfs_a)
	m_en4 = as.matrix(en4[,-(1:4)])
    m_en4 <- apply(m_en4, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })
	m_en4_n = t(t(m_en4) - apply(m_en4, 2, quantile, 0.995))
	m_en4_nr = apply(m_en4_n, 2, function(x) pmax(pmin(x, 0), -12))
	
	message('calculating pwm cgd5')	
	## cgd5 (right 2)
	en5 = prego::gextract_pwm(intervals = cgd5, motifs = tfs_a)
	m_en5 = as.matrix(en5[,-(1:4)])
    m_en5 <- apply(m_en5, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })
	m_en5_n = t(t(m_en5) - apply(m_en5, 2, quantile, 0.995))
	m_en5_nr = apply(m_en5_n, 2, function(x) pmax(pmin(x, 0), -12))
	


	# center
	colnames(m_en1_nr) = paste(colnames(m_en1_nr), "_center", sep = "")
	
	# left 1
	colnames(m_en2_nr) = paste(colnames(m_en2_nr), "_left1", sep = "")
	
	# left 2
	colnames(m_en3_nr) = paste(colnames(m_en3_nr), "_left2", sep = "")
	
	# right 1
	colnames(m_en4_nr) = paste(colnames(m_en4_nr), "_right1", sep = "")
	
	# right 2
	colnames(m_en5_nr) = paste(colnames(m_en5_nr), "_right2", sep = "")


	feats = cbind(m_en1_nr, m_en2_nr,m_en3_nr, m_en4_nr,m_en5_nr)
	#feats3 = cbind(feats2, m_en1_nr)
	#feats = cbind(feats3, m_en1e_nr)

	if(add_dinucs) {
	message('calculating dinucs')
		adins = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
		for(din in adins) {
			gvtrack.create(din, sprintf("seq.%s", din), "sum")
		}
		din1 = gextract(adins, intervals= cgd1[,c(1:3)], iterator=cgd1[,c(1:3)], colnames=paste(adins,"center", sep="_"))
		din2 = gextract(adins, intervals= cgd2[,c(1:3)], iterator=cgd2[,c(1:3)], colnames=paste(adins,"left1", sep="_"))
		din3 = gextract(adins, intervals= cgd3[,c(1:3)], iterator=cgd3[,c(1:3)], colnames=paste(adins,"left2", sep="_"))
		din4 = gextract(adins, intervals= cgd4[,c(1:3)], iterator=cgd4[,c(1:3)], colnames=paste(adins,"right1", sep="_"))
		din5 = gextract(adins, intervals= cgd5[,c(1:3)], iterator=cgd5[,c(1:3)], colnames=paste(adins,"right2", sep="_"))		
		din_feats = cbind(din1[,4:19], din2[,4:19], din3[,4:19], din4[,4:19], din5[,4:19])
		din_feats[is.na(din_feats)] = 0
		feats = cbind(feats, din_feats)
	}

	mod$seqmod_loc_feats = feats

	return(mod)
}
pcg_build_local_seq_feats_raw = function(mod,mot_db, add_dinucs=F)
{
	set.seed(42)
	mod$mot_stats[["cgd_inner"]] = readRDS("data/mot_stats_cgd_inner.RDS")
	mstm = mod$mot_stats[["cgd_inner"]]
	mstm$dist_bin = mstm$dist_bin + 20 - 11
	mstm = mstm[!is.na(mstm$pv),]
	mstm = mstm[mstm$pv < 1e-2,]

	mstm_dom = mstm[mstm$dist_bin >= -5 & mstm$dist_bin <= 5,]
	mstm_3 = mstm[mstm$dist_bin > 5,] 
	mstm_5 = mstm[mstm$dist_bin < -5,] 

	mot_n = table(mstm_dom$motif)
	mot_n3 = table(mstm_3$motif)
	mot_n5 = table(mstm_5$motif)
	tfs = names(mot_n[mot_n>3])
	tfs5 = names(mot_n5[mot_n5>3])
	tfs3 = names(mot_n3[mot_n3>3])
	tfs_a = unique(c(tfs,tfs5,tfs3))
	message('n tfs = ', length(tfs_a))
	cgd = mod$cgdom_ann

	cgd$center = floor(cgd$start + (cgd$end-cgd$start)/2)
	cgd1 = cgd[,c('chrom','start','end','center')]
	cgd1$start = cgd$center - 200
	cgd1$end = cgd$center + 200
	cgd2 = cgd[,c('chrom','start','end','center')]
	cgd2$start  = cgd1$start - 401
	cgd2$end  = cgd1$start - 1
	cgd3 = cgd[,c('chrom','start','end','center')]
	cgd3$start  = cgd2$start - 401
	cgd3$end  = cgd2$start - 1
	
	
	## 1st bin to the right of the center bin (cgd1)
	cgd4 = cgd[,c('chrom','start','end','center')]
	cgd4$start  = cgd1$end + 1
	cgd4$end    = cgd1$end + 401
	
	## 2nd bin to the right
	cgd5 = cgd[,c('chrom','start','end','center')]
	cgd5$start  = cgd4$end + 1
	cgd5$end    = cgd4$end + 401

	##check that bining is okay
	n = 10
	plot(c(cgd1[n,]$start,cgd1[n,]$end),c(1,1),xlim = c(cgd1[n,]$center-1200,cgd1[n,]$center+1200))
	abline(v = cgd1[n,]$center)
	lines(c(cgd2[n,]$start,cgd2[n,]$end),c(1,1),col='darkblue')
	lines(c(cgd3[n,]$start,cgd3[n,]$end),c(1.1,1.1),col='blue')
	lines(c(cgd4[n,]$start,cgd4[n,]$end),c(1,1),col='darkred')
	lines(c(cgd5[n,]$start,cgd5[n,]$end),c(1.1,1.1),col='red')
	
	message('calculating pwm cgd1')
	## cgd1 (center)
	mot_db_f = motif_db %>% filter(motif%in% tfs_a)
	mot_db_f <- mot_db_f %>% 
    mutate(motif = gsub("\\.", "_", motif), motif = gsub("/", "_", motif), motif = gsub("-", "_", motif), motif = gsub("::", "_", motif))

	
mots = unique(mot_db_f$motif)
    for (mot in mots) {
        gvtrack.create(sprintf("PWM%s", mot), func = "pwm", params = list(pssm = mot_db_f[ mot_db_f$motif == 
            mot, ]))
    }
    mots_vt = paste("PWM", mots, sep = "")
	
	
	en1 = gextract(mots_vt, intervals = cgd1, iterator = cgd1,colnames = mots)
	m_en1 = as.matrix(en1[,mots])
    #m_en1 <- apply(m_en1, 2, function(x) {
    #  x[is.na(x)] <- min(x, na.rm = TRUE)
    #  x
    #})
	#m_en1_n = t(t(m_en1) - apply(m_en1, 2, quantile, 0.995))
	#m_en1_nr = apply(m_en1_n, 2, function(x) pmax(pmin(x, 0), -12))
	#
	message('calculating pwm cgd2')	
	## cgd2 (left 1)
	en2 = gextract(mots_vt, intervals = cgd2, iterator = cgd2,colnames = mots)
	m_en2 = as.matrix(en2[,mots])
    #m_en2 <- apply(m_en2, 2, function(x) {
    #  x[is.na(x)] <- min(x, na.rm = TRUE)
    #  x
    #})
	#m_en2_n = t(t(m_en2) - apply(m_en2, 2, quantile, 0.995))
	#m_en2_nr = apply(m_en2_n, 2, function(x) pmax(pmin(x, 0), -12))
	#
	message('calculating pwm cgd3')	
	## cgd3 (left 2)
	en3 = gextract(mots_vt, intervals = cgd3, iterator = cgd3,colnames = mots)
	m_en3 = as.matrix(en3[,mots])
    #m_en3 <- apply(m_en3, 2, function(x) {
    #  x[is.na(x)] <- min(x, na.rm = TRUE)
    #  x
    #})
	#m_en3_n = t(t(m_en3) - apply(m_en3, 2, quantile, 0.995))
	#m_en3_nr = apply(m_en3_n, 2, function(x) pmax(pmin(x, 0), -12))
	
	message('calculating pwm cgd4')	
	## cgd4 (right 1)
	en4 = gextract(mots_vt, intervals = cgd4, iterator = cgd4,colnames = mots)
	m_en4 = as.matrix(en4[,mots])
    #m_en4 <- apply(m_en4, 2, function(x) {
    #  x[is.na(x)] <- min(x, na.rm = TRUE)
    #  x
    #})
	#m_en4_n = t(t(m_en4) - apply(m_en4, 2, quantile, 0.995))
	#m_en4_nr = apply(m_en4_n, 2, function(x) pmax(pmin(x, 0), -12))
	
	message('calculating pwm cgd5')	
	## cgd5 (right 2)
	en5 = gextract(mots_vt, intervals = cgd5, iterator = cgd5,colnames = mots)
	m_en5 = as.matrix(en5[,mots])
    #m_en5 <- apply(m_en5, 2, function(x) {
    #  x[is.na(x)] <- min(x, na.rm = TRUE)
    #  x
    #})
	#m_en5_n = t(t(m_en5) - apply(m_en5, 2, quantile, 0.995))
	#m_en5_nr = apply(m_en5_n, 2, function(x) pmax(pmin(x, 0), -12))
	


	# center
	colnames(m_en1) = paste(colnames(m_en1), "_center", sep = "")
	
	# left 1
	colnames(m_en2) = paste(colnames(m_en2), "_left1", sep = "")
	
	# left 2
	colnames(m_en3) = paste(colnames(m_en3), "_left2", sep = "")
	
	# right 1
	colnames(m_en4) = paste(colnames(m_en4), "_right1", sep = "")
	
	# right 2
	colnames(m_en5) = paste(colnames(m_en5), "_right2", sep = "")


	feats = cbind(m_en1, m_en2,m_en3, m_en4,m_en5)
	#feats3 = cbind(feats2, m_en1_nr)
	#feats = cbind(feats3, m_en1e_nr)

	if(add_dinucs) {
	message('calculating dinucs')
		adins = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
		for(din in adins) {
			gvtrack.create(din, sprintf("seq.%s", din), "sum")
		}
		din1 = gextract(adins, intervals= cgd1[,c(1:3)], iterator=cgd1[,c(1:3)], colnames=paste(adins,"center", sep="_"))
		din2 = gextract(adins, intervals= cgd2[,c(1:3)], iterator=cgd2[,c(1:3)], colnames=paste(adins,"left1", sep="_"))
		din3 = gextract(adins, intervals= cgd3[,c(1:3)], iterator=cgd3[,c(1:3)], colnames=paste(adins,"left2", sep="_"))
		din4 = gextract(adins, intervals= cgd4[,c(1:3)], iterator=cgd4[,c(1:3)], colnames=paste(adins,"right1", sep="_"))
		din5 = gextract(adins, intervals= cgd5[,c(1:3)], iterator=cgd5[,c(1:3)], colnames=paste(adins,"right2", sep="_"))		
		din_feats = cbind(din1[,paste(adins,"center", sep="_")], 
		din2[,paste(adins,"left1", sep="_")],
		din3[,paste(adins,"left2", sep="_")],
		din4[,paste(adins,"right1", sep="_")],
		din5[,paste(adins,"right2", sep="_")])
		din_feats[is.na(din_feats)] = 0
		feats = cbind(feats, din_feats)
	}

	mod$seqmod_loc_feats_raw = feats

	return(mod)
}
pcg_build_local_seq_feats_q_pmax = function(mod, add_dinucs=F)
{
	set.seed(42)
	mod$mot_stats[["cgd_inner"]] = readRDS("data/mot_stats_cgd_inner.RDS")
	mstm = mod$mot_stats[["cgd_inner"]]
	mstm$dist_bin = mstm$dist_bin + 20 - 11
	mstm = mstm[!is.na(mstm$pv),]
	mstm = mstm[mstm$pv < 1e-2,]

	mstm_dom = mstm[mstm$dist_bin >= -5 & mstm$dist_bin <= 5,]
	mstm_3 = mstm[mstm$dist_bin > 5,] 
	mstm_5 = mstm[mstm$dist_bin < -5,] 

	mot_n = table(mstm_dom$motif)
	mot_n3 = table(mstm_3$motif)
	mot_n5 = table(mstm_5$motif)
	tfs = names(mot_n[mot_n>3])
	tfs5 = names(mot_n[mot_n5>3])
	tfs3 = names(mot_n[mot_n3>3])
	tfs_a = unique(c(tfs,tfs5,tfs3))
	message('n tfs = ', length(tfs_a))
	cgd = mod$cgdom_ann

	cgd = mod$cgdom_ann
	cgd$center = floor(cgd$start + (cgd$end-cgd$start)/2)
	cgd1 = cgd[,c('chrom','start','end','center')]
	cgd1$start = cgd$center - 200
	cgd1$end = cgd$center + 200
	cgd2 = cgd[,c('chrom','start','end','center')]
	cgd2$start  = cgd1$start - 401
	cgd2$end  = cgd1$start - 1
	cgd3 = cgd[,c('chrom','start','end','center')]
	cgd3$start  = cgd2$start - 401
	cgd3$end  = cgd2$start - 1
	
	
	## 1st bin to the right of the center bin (cgd1)
	cgd4 = cgd[,c('chrom','start','end','center')]
	cgd4$start  = cgd1$end + 1
	cgd4$end    = cgd1$end + 401
	
	## 2nd bin to the right
	cgd5 = cgd[,c('chrom','start','end','center')]
	cgd5$start  = cgd4$end + 1
	cgd5$end    = cgd4$end + 401

	##check that bining is okay
	n = 10
	plot(c(cgd1[n,]$start,cgd1[n,]$end),c(1,1),xlim = c(cgd1[n,]$center-1200,cgd1[n,]$center+1200))
	abline(v = cgd1[n,]$center)
	lines(c(cgd2[n,]$start,cgd2[n,]$end),c(1,1),col='darkblue')
	lines(c(cgd3[n,]$start,cgd3[n,]$end),c(1.1,1.1),col='blue')
	lines(c(cgd4[n,]$start,cgd4[n,]$end),c(1,1),col='darkred')
	lines(c(cgd5[n,]$start,cgd5[n,]$end),c(1.1,1.1),col='red')
	
	message('calculating pwm cgd1')
	## cgd1 (center)
	en1 = prego::gextract_pwm(intervals = cgd1, motifs = tfs_a)
	m_en1 = as.matrix(en1[,-(1:4)])
    m_en1 <- apply(m_en1, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })

	m_en1_n = t(t(m_en1) - apply(m_en1, 2, quantile, 0.98))
	m_en1_nr = apply(m_en1_n, 2, function(x) pmax(x,0))
	
	message('calculating pwm cgd2')	
	## cgd2 (left 1)
	en2 = prego::gextract_pwm(intervals = cgd2, motifs = tfs_a)
	m_en2 = as.matrix(en2[,-(1:4)])
    m_en2 <- apply(m_en2, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })
	m_en2_n = t(t(m_en2) - apply(m_en2, 2, quantile, 0.98))
	m_en2_nr = apply(m_en2_n, 2, function(x) pmax(x,0))
	
	message('calculating pwm cgd3')	
	## cgd3 (left 2)
	en3 = prego::gextract_pwm(intervals = cgd3, motifs = tfs_a)
	m_en3 = as.matrix(en3[,-(1:4)])
    m_en3 <- apply(m_en3, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })
	m_en3_n = t(t(m_en3) - apply(m_en3, 2, quantile, 0.98))
	m_en3_nr = apply(m_en3_n, 2, function(x) pmax(x,0))
	
	message('calculating pwm cgd4')	
	## cgd4 (right 1)
	en4 = prego::gextract_pwm(intervals = cgd4, motifs = tfs_a)
	m_en4 = as.matrix(en4[,-(1:4)])
    m_en4 <- apply(m_en4, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })
	m_en4_n = t(t(m_en4) - apply(m_en4, 2, quantile, 0.98))
	m_en4_nr = apply(m_en4_n, 2, function(x) pmax(x,0))
	
	message('calculating pwm cgd5')	
	## cgd5 (right 2)
	en5 = prego::gextract_pwm(intervals = cgd5, motifs = tfs_a)
	m_en5 = as.matrix(en5[,-(1:4)])
    m_en5 <- apply(m_en5, 2, function(x) {
      x[is.na(x)] <- min(x, na.rm = TRUE)
      x
    })
	m_en5_n = t(t(m_en5) - apply(m_en5, 2, quantile, 0.98))
	m_en5_nr = apply(m_en5_n, 2, function(x) pmax(x,0))
	


	# center
	colnames(m_en1_nr) = paste(colnames(m_en1_nr), "_center", sep = "")
	
	# left 1
	colnames(m_en2_nr) = paste(colnames(m_en2_nr), "_left1", sep = "")
	
	# left 2
	colnames(m_en3_nr) = paste(colnames(m_en3_nr), "_left2", sep = "")
	
	# right 1
	colnames(m_en4_nr) = paste(colnames(m_en4_nr), "_right1", sep = "")
	
	# right 2
	colnames(m_en5_nr) = paste(colnames(m_en5_nr), "_right2", sep = "")


	feats = cbind(m_en1_nr, m_en2_nr,m_en3_nr, m_en4_nr,m_en5_nr)
	#feats3 = cbind(feats2, m_en1_nr)
	#feats = cbind(feats3, m_en1e_nr)

	if(add_dinucs) {
	message('calculating dinucs')
		adins = c("AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT")
		for(din in adins) {
			gvtrack.create(din, sprintf("seq.%s", din), "sum")
		}
		din1 = gextract(adins, intervals= cgd1[,c(1:3)], iterator=cgd1[,c(1:3)], colnames=paste(adins,"center", sep="_"))
		din2 = gextract(adins, intervals= cgd2[,c(1:3)], iterator=cgd2[,c(1:3)], colnames=paste(adins,"left1", sep="_"))
		din3 = gextract(adins, intervals= cgd3[,c(1:3)], iterator=cgd3[,c(1:3)], colnames=paste(adins,"left2", sep="_"))
		din4 = gextract(adins, intervals= cgd4[,c(1:3)], iterator=cgd4[,c(1:3)], colnames=paste(adins,"right1", sep="_"))
		din5 = gextract(adins, intervals= cgd5[,c(1:3)], iterator=cgd5[,c(1:3)], colnames=paste(adins,"right2", sep="_"))		
		din_feats = cbind(din1[,4:19], din2[,4:19], din3[,4:19], din4[,4:19], din5[,4:19])
		din_feats[is.na(din_feats)] = 0
		feats = cbind(feats, din_feats)
	}

	mod$seqmod_loc_feats_qpmax = feats

	return(mod)
}



pcg_init_epi_track_lib = function(mod)
{
	mod$epi_tracks = as.data.frame(fread("data/index_tracks.txt"))
    mod$k27_track_thresh = read.table("data/track_thresh.txt",row.names=1, 
        sep = "\t", stringsAsFactors = F, header = T)
	mod$epi_tracks$short_name_k4 = paste(mod$epi_tracks$short_name,"k4", sep="_")
	mod$epi_tracks_all = mod$epi_tracks
	mod$epi_tracks = mod$epi_tracks[mod$epi_tracks$use_comb==1,]
	mod$epi_tracks_cnt = as.data.frame(fread("data/index_tracks_cnt.txt"))
	mod$epi_tracks_atac = as.data.frame(fread('data/index_tracks_atac.csv'))
	mod$atac_track_thresh = read.table("data/track_thresh_atac.txt", sep="\t", stringsAsFactors=F,header=T)
	mod$track_atac_ribo = fread("data/index_tracks_atac.csv", stringsAsFactors=F)
	mod$k27_track_thresh_cnt = read.table("data/track_thresh_cnt.txt", sep="\t", stringsAsFactors=F,header=T)
	mod$k4_track_thresh_cnt = read.table("data/track_thresh_k4_cnt.txt", sep="\t", stringsAsFactors=F,header=T,row.names = 1)
	mod$k4_track_thresh = read.table("data/track_thresh_k4_cnt.txt", sep="\t", stringsAsFactors=F,header=T,row.names = 1)
	#mod$cnt_track_cov = read.table("data/track_cov_cnt.txt", sep="\t", stringsAsFactors=F,header=T)
		#mod$cnt_track_cov_k4 = read.table("data/track_cov_cnt_k4.txt", sep="\t", stringsAsFactors=F,header=T,row.names = 1)###be aware colnames shifted
		#colnames(mod$cnt_track_cov_k4) = c('Total.intervals','NaN.intervals','Min','Max','Sum','Mean','Std.dev','X')
mod$atac_track_cov = read.table("data/track_cov_atac.txt", sep="\t", stringsAsFactors=F,header=T)
				


	return(mod)
}

compute_specific_atac_quantiles_ribo = function(mod, tracks,sh_nms,file_name,th_percent="X0.99",calc_ribo=TRUE){
    ndx = data.frame(tr_atac = tracks, short_name = sh_nms)
	for(i in 1:nrow(ndx)) {
#		message("running ", ndx$short_name[i])
		gvtrack.create(ndx$short_name[i], ndx$tr_atac[i], "sum") 
		gvtrack.iterator(ndx$short_name[i], sshift = -140, eshift= 140)
        }

		#mod$track_atac_ribo = fread(paste0("data/index_tracks_atac_",file_name,'.csv'), stringsAsFactors=F)

	
	tss = mod$tss
    tss = tss[ grepl('Rp[l/s]',tss$geneSymbol),]
    tss$start = tss$start - 400
    tss$end = tss$end + 400
    tss = tss %>% arrange(chrom,start)
    message('calcultating ribo coverage')
    gext = gextract(ndx$short_name,intervals = tss,iterator = tss,colnames = ndx$short_name)
    gext[is.na(gext)] = 0
    csum = colSums(gext[,-c(1,2,3,ncol(gext))])
    ndx$ribo_c = as.numeric(csum)
    message('calcultating total coverage')
    gext = gextract(ndx$tr_atac,intervals = gintervals.all(),iterator = 100,colnames = ndx$short_name)
    gext[is.na(gext)] = 0
    csum = colSums(gext[,-c(1,2,3,ncol(gext))])
    ndx$cov = as.numeric(csum)
    #mod$track_atac_ribo = ndx
    fwrite(ndx,paste0("data/index_tracks_atac_",file_name,'.csv'))
		 message('calcultating quantiles')
		thresh = list()
		for(i in 1:nrow(ndx)) {
			thresh[[i]] = gquantiles(ndx$short_name[i], c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
		}
		th_mat = do.call('rbind',thresh)
		rownames(th_mat) = ndx$short_name
		write.table(th_mat, file=paste0("data/track_thresh_atac_",file_name,'.txt'), quote=F, sep="\t")
		#mod$atac_track_thresh = th_mat       
	
 return(list(ndx = ndx,th_mat = th_mat))
}

pcg_update_track_q_thresh = function(mod)
{
	ndx = mod$epi_tracks
	thresh = list()
	for(i in 1:nrow(ndx)) {
		thresh[[i]] = gquantiles(ndx$short_name[i], c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
	}
	th_mat = do.call('rbind',thresh)
	rownames(th_mat) = ndx$short_name
	write.table(th_mat, file="data/track_thresh.txt", quote=F, sep="\t")

	mod$k27_track_thresh = th_mat

	thresh_k4 = list()
	tnm_k4 = gen_k4_vt(mod)
	for(tnm in tnm_k4) {
		thresh_k4[[tnm]] = gquantiles(tnm, c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
	}
	th_mat_k4 = do.call('rbind',thresh_k4)
	rownames(th_mat_k4) = tnm_k4
	write.table(th_mat_k4, file="data/track_thresh_k4.txt", quote=F, sep="\t")
	mod$k4_track_thresh = th_mat_k4
	return(mod)
}

pcg_compute_k27_doms = function(mod, th_percent="X0.99")
{
	ndx = mod$epi_tracks
	for(i in 1:nrow(ndx)) {
#		message("running ", ndx$short_name[i])
		gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "sum") 
		gvtrack.iterator(ndx$short_name[i], sshift = -140, eshift= 140)
	}
	if(file.exists("data/track_thresh.txt")) {
		mod$k27_track_thresh = read.table("data/track_thresh.txt", sep="\t", stringsAsFactors=F)
		th_mat = mod$k27_track_thresh
		mod$k4_track_thresh = read.table("data/track_thresh_k4.txt", sep="\t", stringsAsFactors=F)
	} else {
		thresh = list()
		for(i in 1:nrow(ndx)) {
			thresh[[i]] = gquantiles(ndx$short_name[i], c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
		}
		th_mat = do.call('rbind',thresh)
		rownames(th_mat) = ndx$short_name
		write.table(th_mat, file="data/track_thresh.txt", quote=F, sep="\t")

		thresh_k4 = list()
		tnm_k4 = gen_k4_vt(mod)
		for(tnm in tnm_k4) {
			thresh_k4[[tnm]] = gquantiles(tnm, c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
		}
		th_mat_k4 = do.call('rbind',thresh_k4)
		rownames(th_mat_k4) = tnm_k4
		write.table(th_mat_k4, file="data/track_thresh_k4.txt", quote=F, sep="\t")
		mod$k4_track_thresh = th_mat_k4
	}

	gvtrack.create("cg_max", "seq.CG_500_mean_new", "max")
	gvtrack.create("gc_max", "seq.GC500_bin20", "max")

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

	mat = gextract(ndx$short_name, intervals=doms1k, iterator=doms1k)
	mat_n = t(t(mat[,4:(ncol(mat)-1)])/apply(mat[,4:(ncol(mat)-1)], 2, median))
	mat_n2 = mat_n / rowMeans(mat_n)
	mod$k27_mat = mat
	mod$k27_mat_n2 = mat_n2

	return(mod)
}

gen_k27_vt = function(mod) 
{
	ndx = mod$epi_tracks
	ndx = mod$epi_tracks
	ndx = ndx%>% filter(short_name != 'EB4_cnt')
	for(i in 1:nrow(ndx)) {
		if(!is.na(ndx$track_k27[i])) {
			gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "avg") 
			gvtrack.iterator(ndx$short_name[i], sshift = -500, eshift= 500)
		}
	}
		gvtrack.create('EB4_cnt', 'jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm', "avg") 
		gvtrack.iterator('EB4_cnt', sshift = -500, eshift= 500)
}

gen_k4_vt = function(mod) 
{
	ndx = mod$epi_tracks
	for(i in 1:nrow(ndx)) {
		if(!is.na(ndx$track_k4[i])) {
			gvtrack.create(ndx$short_name_k4[i], ndx$track_k4[i], "sum") 
			gvtrack.iterator(ndx$short_name_k4[i], sshift = -140, eshift= 140)
		}
	}
	k4_tns =  ndx$short_name_k4[!is.na(ndx$track_k4)]
	return(k4_tns)
}

add_dom_k4_mat = function(mod)
{
	k4_tns = gen_k4_vt(mod)
	a = gextract(k4_tns, iterator=20, intervals=mod$k27_doms)
	b = tgs_matrix_tapply(t(a[,4:(ncol(a)-1)]), a$intervalID, max)
	mod$k27_mat_k4 = b
	return(mod)
}

pcg_cluster_k27_doms = function(mod, K=60)
{
	km = tglkmeans::TGL_kmeans(log2(mod$k27_mat_n2+0.2), k=K, id_column=F, seed=19)

	cls_cg = tapply(mod$k27_doms$CG_max, km$cluster, median)
	cls_gc = tapply(mod$k27_doms$GC, km$cluster, median)
#	ordmat = log2(0.2+mat_n2[order(cls_cg[km$cluster]),])
	ordmat = log2(0.2+mod$k27_mat_n2[order(km$cluster),])
	ordmat_s = apply(ordmat, 2, rollmean, 20, f='e')
	
	mod$k27_km = km

	shades = colorRampPalette(c("white","gray","darkblue","purple"))

	pdf("figs/k27_doms_cmp.pdf", h=16,w=6)
	image(t(pmin(pmax(ordmat_s,-2),2)), col=shades(1000), breaks = seq(-2,2,l=1001), xaxt='n', yaxt='n')
	dev.off()

}

gen_cg_corplot = function(mod)
{
	ndx = mod$epi_tracks

	gvtrack.create("CG2k", "seq.CG_500_mean_new", "avg")
	gvtrack.iterator("CG2k", sshift=-1000,eshift=1000)
	gvtrack.create("CG1k", "seq.CG_500_mean_new", "avg")
	gvtrack.iterator("CG1k", sshift=-500,eshift=500)

	samp_d = mod$k27_doms[sample(1:nrow(mod$k27_doms), 1000),]
	samp_d$start = samp_d$start - 3000
	samp_d$end = samp_d$end + 3000
#	prf = gextract(c(ndx$short_name,"seq.CG_500_mean_new", "CG1k", "CG2k"), iterator=20, intervals=mod$k27_doms)
	prf = gextract(c(ndx$short_name,"seq.CG_500_mean_new", "CG1k", "CG2k"), iterator=20, intervals=samp_d)

	prf_v  = prf[,ndx$short_name]
	lprf_v = log2(0.05+t(t(prf_v)/colMeans(prf_v)))
	c500 = apply(prf_v, 2,cor,prf$seq.CG_500_mean_new,m="s")
	c1k = apply(prf_v, 2,cor,prf$CG1k,m="s")
	c2k = apply(prf_v, 2,cor,prf$CG2k,m="s")

	pdf("figs/CG2k_dom_strat.pdf", h=3, w=24)
	layout(matrix(1:ncol(prf_v),nrow=1))
	f_hcg = prf$CG2k > 0.05
	f_lcg = prf$CG2k < 0.02
	f_icg = !f_hcg & !f_lcg
	for(i in 1:ncol(prf_v)) {
		par(mar=c(2,1,2,1))
		plot(density(lprf_v[f_hcg,i]), col="red", yaxt='n', main=colnames(prf_v)[i])
		lines(density(lprf_v[f_lcg,i]), col="black")
		lines(density(lprf_v[f_icg,i]), col="pink")
	}
	dev.off()
	pdf("figs/CG_dom_strat.pdf", h=3, w=24)
	layout(matrix(1:ncol(prf_v),nrow=1))
	f_hcg = prf$seq.CG_500_mean_new > 0.05
	f_lcg = prf$seq.CG_500_mean_new < 0.02
	f_icg = !f_hcg & !f_lcg
	for(i in 1:ncol(prf_v)) {
		par(mar=c(2,1,2,1))
		plot(density(lprf_v[f_hcg,i]), col="red", yaxt='n', main=colnames(prf_v)[i])
		lines(density(lprf_v[f_lcg,i]), col="black")
		lines(density(lprf_v[f_icg,i]), col="pink")
	}
	dev.off()

}

pcg_gen_cnt_comparisons = function(mod)
{
	tr_es = mod$epi_tracks_all[grepl("ES_", mod$epi_tracks_all$short_name),]
	samp_d = mod$k27_doms[sample(1:nrow(mod$k27_doms), 1000),]
	samp_d$start = samp_d$start - 3000
	samp_d$end = samp_d$end + 3000

	horiz = 140
	for(i in 1:nrow(tr_es)) {
		gvtrack.create(tr_es$short_name[i], tr_es$track_k27[i], "sum") 
		gvtrack.iterator(tr_es$short_name[i], sshift = -horiz, eshift= horiz)
	}
	#cnt_doms_a = gscreen("EB4_cnt > 18.615989")
	cnt_doms_a = gscreen("ES_cnt > 5.802099")
	cnt_doms_a$l = cnt_doms_a$end - cnt_doms_a$start
	cnt_doms = cnt_doms_a[cnt_doms_a$l>20,]
	cnt_bord5 = cnt_doms
	cnt_bord5$end = cnt_bord5$start+1
	cnt_bord5$strand = 1
	cnt_bord3 = cnt_doms
	cnt_bord3$start = cnt_bord3$end-1
	cnt_bord3$strand = -1
	cnt_bord = rbind(cnt_bord3, cnt_bord5)
	gvtrack.create("cnt_dom_dist", cnt_bord, "distance")

	chip_doms_a = gscreen("ES_ch > 117.201529")
	chip_doms_a$l = chip_doms_a$end - chip_doms_a$start
	chip_doms = chip_doms_a[chip_doms_a$l>20,]
	chip_bord5 = chip_doms
	chip_bord5$end = chip_bord5$start+1
	chip_bord5$strand = 1
	chip_bord3 = chip_doms
	chip_bord3$start = chip_bord3$end-1
	chip_bord3$strand = -1
	chip_bord = rbind(chip_bord3, chip_bord5)
	gvtrack.create("chip_dom_dist", chip_bord, "distance")


	prof = gextract(c(tr_es$short_name,"cnt_dom_dist", "chip_dom_dist", "seq.CG_500_mean_new", "seq.GC500_bin20"), intervals=samp_d, iterator=20)
	prof = prof[rowSums(is.na(prof))==0,]
	prf_v  = prof[,tr_es$short_name]
	lprf_v = log2(0.05+t(t(prf_v)/colMeans(prf_v)))

	cnt_dst_bin = 20*pmin(pmax(floor(prof$cnt_dom_dist/20),-250),250)
	chip_dst_bin = 20*pmin(pmax(floor(prof$chip_dom_dist/20),-250),250)
	
	cnt_cnt_trend = tapply(lprf_v[,"ES_cnt"], cnt_dst_bin, mean)
	chip_cnt_trend = tapply(lprf_v[,"ES_ch"], cnt_dst_bin, mean)

	cnt_chip_trend = tapply(lprf_v[,"ES_cnt"], chip_dst_bin, mean)
	chip_chip_trend = tapply(lprf_v[,"ES_ch"], chip_dst_bin, mean)

	pdf("figs/cnt_chip_borders.pdf", w=12,h=6)
	layout(matrix(1:2,nrow=1))
	plot(names(cnt_cnt_trend), cnt_cnt_trend, xlim=c(-4000,2000), type="l", lwd=3, main="at CnT domain borders", ylim=c(-3.7,3.2), xlab="dist from border")
	lines(names(chip_cnt_trend), chip_cnt_trend, lwd=3, col="blue")
	grid()
	abline(v=0)
	plot(names(cnt_chip_trend), cnt_chip_trend, xlim=c(-4000,2000), type="l", lwd=3, main="at ChIP domain borders", ylim=c(-3.7,3.2), xlab="dist from border")
	lines(names(chip_chip_trend), chip_chip_trend, lwd=3, col="blue")
	grid()
	abline(v=0)
	dev.off()

}

pcg_gen_eb4_doms = function(mod)
{
	gen_track_vext()
	k4_tns = gen_k4_vt(mod)

	t27 = mod$k27_track_thresh["EB4_cnt",6]
	t4 = mod$k4_track_thresh["EB4_cnt_k4",6]
   doms = gscreen(sprintf("EB4_cnt > %s | EB4_cnt_k4 > %s", t27, t4))
	doms = doms[doms$chrom != "chrM" & doms$chrom != "chrY",]

	doms_ext = doms
	doms_ext$start = doms_ext$start - 100
	doms_ext$end = doms_ext$end + 100
	doms_ext = gintervals.canonic(doms_ext)
	doms_ext$start = doms_ext$start + 100
	doms_ext$end = doms_ext$end - 100
	doms = doms_ext

   doms$l = doms$end - doms$start

	a = gextract(c("EB4_cnt", "EB4_cnt_k4", "ES_ch_k4","atac_ext"), iterator=20, intervals=doms)
	a[is.na(a)]= 0
	b = tgs_matrix_tapply(t(a[,4:(ncol(a)-1)]), a$intervalID, max)
	doms$eb4_k27 = b[,"EB4_cnt"]
	doms$eb4_k4 = b[,"EB4_cnt_k4"]
	doms$es_k4 = b[,"ES_ch_k4"]
	doms$eb_atac = b[,"atac_ext"]

#fix this 
	max_loc = a[a$atac_ext == doms[a$intervalID,"eb_atac"], c("start","intervalID")]
	doms$loc_max_atac = tapply(max_loc$start, max_loc$intervalID, mean)

# 	doms$loc_max_atac = a[a$atac_ext == doms[a$intervalID,"atac_ext"], "start"]
# 	doms$loc_max_atac = a[a$atac_ext == doms[a$intervalID,"atac_ext"], "start"]

	gvtrack.create("cg_max", "seq.CG_500_mean_new", "max")
	gvtrack.create("gc_max", "seq.GC500_bin20", "max")

   doms_cg = gextract(c("seq.GC500_bin20", "gc_max", "cg_max"), iterator=doms, intervals=doms, colnames=c("GC","GC_max", "CG_max"))
   doms_tss = gintervals.neighbors(doms, mod$tss[,c("chrom","start","end","strand","geneSymbol")])
	doms_tad = gintervals.neighbors(doms,"intervs.global.tad_names")
	base_coli = 10
	colnames(doms_tad)[base_coli + 4] = "tad_index"
	doms_tad$tad_5_dist = doms_tad$start - doms_tad[, base_coli+1]
	doms_tad$tad_3_dist = doms_tad[,base_coli + 2] - doms_tad$end

   doms = cbind(doms, doms_cg[,c("GC", "GC_max", "CG_max")])
   doms = cbind(doms, doms_tss[,c("dist","strand","geneSymbol")])
	colnames(doms)[base_coli+3] = "tss_dist"
   doms = cbind(doms, doms_tad[,c("tad_name", "tad_index", "tad_5_dist", "tad_3_dist")])

	sines = gintervals.load("intervs.global.rmsk_sine")
	gvtrack.create("line_d", "intervs.global.rmsk_line", "distance")
	gvtrack.create("ltr_d", "intervs.global.rmsk_ltr", "distance")
	gvtrack.create("sine_d",sines, "distance")
	gvtrack.create("lowcomplex_d", "intervs.global.rmsk_low_complexity", "distance")
	dom_rpts_20 = gextract(c("ifelse(line_d==0,1,0)", "ifelse(ltr_d==0,1,0)", "ifelse(sine_d==0,1,0)", "ifelse(lowcomplex_d==0,1,0)"), iterator=20, intervals=doms, colnames=c("line","ltr","sine","low_complex"))
	dom_rpts = data.frame(
						line = tapply(dom_rpts_20$line, dom_rpts_20$intervalID, mean),
						sine = tapply(dom_rpts_20$sine, dom_rpts_20$intervalID, mean),
						ltr = tapply(dom_rpts_20$ltr, dom_rpts_20$intervalID, mean),
						low_complex = tapply(dom_rpts_20$low_complex, dom_rpts_20$intervalID, mean))
	doms = cbind(doms, dom_rpts)

	cgbin = as.numeric(cut(doms$CG_max, c(-1,0.02,0.04,1)))
	doms$cgtss_mode = ifelse(doms$tss_dist == 0, 0, 1)*3 + cgbin
	mod$eb4_biv_doms = doms

	return(mod)
}

pcg_gen_tss_rna_k4k27_annots = function(mod)
{
	gen_track_vext()
	k4_tns = gen_k4_vt(mod)

	tss = mod$tss
	rownames(tss) = tss$geneSymbol
	tss = tss[tss$chrom != "chrM",]

	eb_rna = mod$eb_legc
	epi_rna = mod$emb_type_legc[,"Epiblast"]
	amn_rna = mod$emb_type_legc[,"Amnion/Chorion"]

	tss$epi_rna = rep(log2(1e-5), nrow(tss))
	tss$amn_rna = rep(log2(1e-5), nrow(tss))
	gnms = intersect(rownames(tss), names(epi_rna))
	tss[gnms,"eb_rna"] = eb_rna[gnms]
	tss[gnms,"epi_rna"] = epi_rna[gnms]
	tss[gnms,"amn_rna"] = amn_rna[gnms]

	tss_ext = tss
	tss_ext$start = tss_ext$start - 1000
	tss_ext$end = tss_ext$end + 1000
	prof = gextract(c("EB4_cnt", "EB4_cnt_k4", "ES_ch_k4","atac_ext","seq.CG_500_mean_new"), iterator=20, intervals=tss_ext)
	prof[is.na(prof)]= 0
	T_k27 = mod$k27_track_thresh["EB4_cnt","X0.99"]
	T_k4 = mod$k4_track_thresh["EB4_cnt","X0.99"]
	T_atac = 497
	tss$n2k_eb4_k27 = 20*tapply(prof[,"EB4_cnt"]>T_k27, prof$intervalID, sum)
	tss$n2k_eb4_k4 = 20*tapply(prof[,"EB4_cnt_k4"]>T_k4, prof$intervalID, sum)
	tss$eb4_k27_max = tapply(prof[,"EB4_cnt"], prof$intervalID, max)
	tss$eb4_k4_max = tapply(prof[,"EB4_cnt_k4"], prof$intervalID, max)
	tss$n2k_atac = 20*tapply(prof[,"atac_ext"]>T_atac, prof$intervalID, sum)
	tss$atac_max = tapply(prof[,"atac_ext"], prof$intervalID, max)
	tss$cg_max = tapply(prof[,"seq.CG_500_mean_new"], prof$intervalID, max)

	mod$epi_tss = tss

	doms = mod$eb4_biv_doms[,c("chrom","start","end")]
	doms$id = 1:nrow(doms)
	dom_tss = gintervals.neighbors(doms, tss[,c("chrom","start","end","strand", "epi_rna", "eb_rna")], maxneighbors=100, maxdist=1e+5, mindist=-1e+5)

	dom_region_rna = data.frame(id = doms$id)
	dom_region_rna = cbind(dom_region_rna, matrix(NA, nrow=length(doms$id), ncol=6))
	rownames(dom_region_rna) = doms$id
	ci = 2
	for(horiz in c(1000,2000,4000,10000,50000, 100000)) {
		dh = dom_tss[abs(dom_tss$dist) <= horiz,]
		gmax = tapply(dh$eb_rna, dh$id, max)
		dom_region_rna[names(gmax), ci] = gmax
		ci = ci + 1
	}

	colnames(dom_region_rna) = c("ID", "rna1k", "rna2k", "rna4k", "rna10k", "rna50k", "rna100k")
	mod$eb4_biv_rna = dom_region_rna
	return(mod)
}

pcg_gen_scale50k_k27 = function(mod)
{
	all_c = gintervals.all()
	all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]

	d = mod$eb4_biv_doms
	d = d[d$l > 300 & !is.na(d$line) & d$line < 0.5 & d$ltr < 0.5 & d$sine <0.5 & d$low_complex<0.5,]

	pall = gextract(c("EB4_cnt"), iterator=500, intervals = all_c, colnames=c("k27"))
	T_k27 = mod$k27_track_thresh["EB4_cnt","X0.99"]
	pall$dom = ifelse(pall$k27 > T_k27, 1, 0)
	pall$dom50 = rollmean(pall$dom,100,f='e')
	d_to50 = gintervals.neighbors(d, pall)

	pdf("figs/k27_cg_seeding.pdf", w=6, h=6)
	bin_cols = c("black","orange","orange","red","red","darkred")
	cg_b = as.numeric(cut(d_to50$CG_max, c(-1,0.02,0.03,0.04,0.05,0.06,1)))
	plot(density(log2(0.02+d_to50$dom50[cg_b==1])), col=bin_cols[1], xlab="log2(50kb density", xlim=c(-6,0), ylim=c(0,0.5), main="seeding density")
	for(i in 1:6) {
		lines(density(log2(0.02+d_to50$dom50[cg_b==i])), col=bin_cols[i], lwd=2)
	}
	dev.off()
	all_to50 = gintervals.neighbors(mod$eb4_biv_doms, pall)
	
	mod$eb4_biv_doms$scale50k = all_to50$dom50
	mod$eb4_k27_doms = mod$eb4_biv_doms[mod$eb4_biv_doms$eb4_k27 > 18,]
	return(mod)
}



####
pcg_gen_cg_h_doms = function(mod, force_update=F)
{
	if(!force_update & file.exists("data/cgdom_ann_new.RDS")) {
		mod$cgdom_ann = readRDS("data/cgdom_ann_new.RDS")
		mod$cgd_rpt = readRDS("data/cgd_rpt_new.RDS")
		return(mod)
	}
	options(gmax.data.size=1e+9)
	gen_k4_vt(mod)
	gen_k27_vt(mod)
	gen_track_vext()
	gvtrack.create("k27_raw", 'jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a', "sum")
	gvtrack.iterator("k27_raw", sshift = -140, eshif = 140)
	cgd = gscreen("seq.CG_500_mean_new>0.04")
	message("done gscreen")
	cgd = cgd[!cgd$chrom %in% c("chrM", "chrY"),]
	cgd$l = cgd$end-cgd$start
	cgd = cgd[cgd$l>300,]

	rownames(cgd) = 1:nrow(cgd)

#annotate: max CG, max k27, k4, GC
	prof = gextract(c("EB4_cnt", "EB3_cnt", "EB4_cnt_k4", "atac_ext","seq.CG_500_mean_new", "seq.GC500_bin20",'k27_raw'), iterator=20, intervals=cgd)
	prof[is.na(prof)]= 0
	message("done gextract")
	cgd$eb3_k27_max = tapply(prof[,"EB3_cnt"], prof$intervalID, max)
	cgd$eb4_k27_max = tapply(prof[,"EB4_cnt"], prof$intervalID, max)
	cgd$eb4_k27_mean = tapply(prof[,"EB4_cnt"], prof$intervalID, mean)
	cgd$eb4_k4_max = tapply(prof[,"EB4_cnt_k4"], prof$intervalID, max)
	cgd$eb4_k4_mean = tapply(prof[,"EB4_cnt_k4"], prof$intervalID, mean)
	cgd$l6_eb4_k4_mean = log2(6+cgd$eb4_k4_mean)
	cgd$l1_eb4_k4_mean = log2(1+cgd$eb4_k4_mean)
	cgd$k27_raw_max = tapply(prof[,"k27_raw"], prof$intervalID, max)
	T_k27_995 = mod$k27_track_thresh["EB4_cnt",7]
	T_k4_995 = mod$k4_track_thresh["EB3_wt_t_wt",7]
	T_k27_99 = mod$k27_track_thresh["EB4_cnt",6]
	T_k4_99 = mod$k4_track_thresh["EB3_wt_t_wt",6]
	up_k27 = ifelse(prof[,"EB4_cnt"]>T_k27_995,1,0)
	up_k4 = ifelse(prof[,"EB4_cnt_k4"]>T_k4_995,1,0)
	cgd$eb4_k27_cov = tapply(up_k27 * (1-up_k4), prof$intervalID, mean)
	cgd$eb4_k4_cov = tapply((1-up_k27) * up_k4, prof$intervalID, mean)
	cgd$eb4_biv_cov = tapply(up_k27 * up_k4, prof$intervalID, mean)

	up_k27_1 = ifelse(prof[,"EB4_cnt"]>T_k27_99,1,0)
	up_k4_1 = ifelse(prof[,"EB4_cnt_k4"]>T_k4_99,1,0)
	cgd$eb4_k27_cov_1 = tapply(up_k27_1 * (1-up_k4_1), prof$intervalID, mean)
	cgd$eb4_k4_cov_1 = tapply((1-up_k27_1) * up_k4_1, prof$intervalID, mean)
	cgd$eb4_biv_cov_1 = tapply(up_k27_1 * up_k4_1, prof$intervalID, mean)

	cgd$atac_max = tapply(prof[,"atac_ext"], prof$intervalID, max)
	cgd$cg_max = tapply(prof[,"seq.CG_500_mean_new"], prof$intervalID, max)
	cgd$gc_max = tapply(prof[,"seq.GC500_bin20"], prof$intervalID, max)

	#fcov = (prof$EB4_cnt_k4+prof$EB4_cnt)>20
	#prof[!fcov,"EB4_cnt_k4"] = NA
	#prof[!fcov,"EB4_cnt"] = NA
	#cgd$eb4_k4_rat = tapply(log2((10+prof[,"EB4_cnt_k4"])/(10+prof[,"EB4_cnt"])), #prof$intervalID, max, na.rm=T)
	#cgd$eb4_k27_rat = tapply(log2((10+prof[,"EB4_cnt"])/(10+prof[,"EB4_cnt_k4"])), #prof$intervalID, max, na.rm=T)

	sines = gintervals.load("intervs.global.rmsk_sine")
	gvtrack.create("exon_d", "intervs.global.exon", "distance")
	gvtrack.create("line_d", "intervs.global.rmsk_line", "distance")
	gvtrack.create("ltr_d", "intervs.global.rmsk_ltr", "distance")
	gvtrack.create("simp_d", "intervs.global.rmsk_simple_repeat", "distance")
	gvtrack.create("sine_d",sines, "distance")
	gvtrack.create("lowcomplex_d", "intervs.global.rmsk_low_complexity", "distance")

	cgd_rpts_20 = gextract(c("ifelse(exon_d==0,1,0)","ifelse(line_d==0,1,0)", "ifelse(ltr_d==0,1,0)", "ifelse(sine_d==0,1,0)", "ifelse(lowcomplex_d==0,1,0)", "ifelse(simp_d==0,1,0)"), iterator=20, intervals=cgd, colnames=c("exon","line","ltr","sine","low_complex", "simple"))
	message("done gextract repeats")
	cgd_rpts = data.frame(
						exon = tapply(cgd_rpts_20$exon, cgd_rpts_20$intervalID, mean),
						line = tapply(cgd_rpts_20$line, cgd_rpts_20$intervalID, mean),
						sine = tapply(cgd_rpts_20$sine, cgd_rpts_20$intervalID, mean),
						ltr = tapply(cgd_rpts_20$ltr, cgd_rpts_20$intervalID, mean),
						simp = tapply(cgd_rpts_20$simp, cgd_rpts_20$intervalID, mean),
						low_complex = tapply(cgd_rpts_20$low_complex, cgd_rpts_20$intervalID, mean))
	cgd_rpts$tot_rpt = rowSums(cgd_rpts[,c("line","sine","ltr","simp","low_complex")])
	
	cgd_ann = cbind(cgd, cgd_rpts)
	f_bad = cgd_ann$line > 0.5 | cgd_ann$ltr > 0.5 
	f_bad1 =  !f_bad & cgd_ann$tot_rpt > 0.5 
	f_exon = cgd_ann$exon > 0.8
	f_lowmapability = !f_bad & !f_bad1 & (cgd_ann$k27_raw_max+cgd_ann$eb4_k4_max) == 0
	f_mask = f_bad | f_bad1 | f_lowmapability
	cgd_ann_mask = cgd_ann[f_mask,]
	cgd_ann = cgd_ann[!f_mask,]

	cgd_tss = gintervals.neighbors(cgd_ann, mod$tss)
	cgd_ann$tss_dist = cgd_tss$dist
	cgd_ann$tss_strand = cgd_tss$strand
	cgd_ann$tss_start = cgd_tss$start
	cgd_ann$tss_gene = cgd_tss$geneSymbol

	gvtrack.create("d_cgd",cgd_ann[,1:3], "distance")
	gvtrack.create("d_cgdnotss",cgd_ann[abs(cgd_ann$tss_dist)>1000,1:3], "distance")

	#all_c = gintervals.all()
	#all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]
	#cg_trace = gextract(c("ifelse(seq.CG_500_mean_new > 0.04, 1, 0)"),
	#							 intervals=all_c, iterator = 200, 
	#								colnames("cg_all"))
#
	#cg_trace1 = gextract(c("ifelse(!is.na(d_cgd) & d_cgd==0, 1, 0)"),
	#							 intervals=all_c, iterator = 200, 
	#								colnames("cg_filt"))

	#cg_trace2 = gextract(c("ifelse(!is.na(d_cgdnotss) & d_cgdnotss==0, 1, 0)"),
	#							 intervals=all_c, iterator = 200, 
	#								colnames("cg_notss"))

	#cgd_regstat = data.frame( 
	#	chrom = cg_trace$chrom,
	#	start = cg_trace$start,
	#	end = cg_trace$end,
	#	cg_2k = rollmean(cg_trace[,4],10,f='e'),
	#	cg_4k = rollmean(cg_trace[,4],20,f='e'),
	#	cg_8k = rollmean(cg_trace[,4],40,f='e'),
	#	cg_16k = rollmean(cg_trace[,4],80,f='e'),
	#	cg_32k = rollmean(cg_trace[,4],160,f='e'),
	#	cgd_2k = rollmean(cg_trace1[,4],10,f='e'),
	#	cgd_4k = rollmean(cg_trace1[,4],20,f='e'),
	#	cgd_8k = rollmean(cg_trace1[,4],40,f='e'),
	#	cgd_16k = rollmean(cg_trace1[,4],80,f='e'),
	#	cgd_32k = rollmean(cg_trace1[,4],160,f='e'))
	#cgd_regstat$cgnotss_2k = rollmean(cg_trace2[,4],10,f='e')
	#cgd_regstat$cgnotss_4k = rollmean(cg_trace2[,4],20,f='e')
	#cgd_regstat$cgnotss_8k = rollmean(cg_trace2[,4],40,f='e')
	#cgd_regstat$cgnotss_16k = rollmean(cg_trace2[,4],80,f='e')
	#cgd_regstat$cgnotss_32k = rollmean(cg_trace2[,4],160,f='e')
#regional CG elements
	#cgd_cent = cgd_ann[,c("chrom","start","end")]
	#cgd_cent$start = (cgd_cent$start+cgd_cent$end)/2
	#cgd_cent$end = cgd_cent$start + 1
	#cgd_reg = gintervals.neighbors(cgd_cent, cgd_regstat)
	#cgd_ann  = cbind(cgd_ann, cgd_reg[,-c(1:6,ncol(cgd_reg))])
	####seed clust
	#all_c = gintervals.all()
	#all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]
	#pall = gextract(c("EB4_cnt"), iterator=500, intervals = all_c, colnames=c("k27"))
	#T_k27 = mod$k27_track_thresh["EB4_cnt","X0.99"]
	#pall$dom = ifelse(pall$k27 > T_k27, 1, 0)
	#pall$dom50 = rollmean(pall$dom,100,f='e')	
	#all_to50 = gintervals.neighbors(cgd_ann, pall)	
	#cgd_ann$scale50k = all_to50$dom50
	
	###divide into clusters
	d = cgd_ann
	d$type = ifelse(d$eb3_k27_max >= 6 & log2(1+d$eb4_k4_max)<6 ,'cl1','cl5')

	d$type = ifelse(d$eb3_k27_max >= 6 & log2(1+d$eb4_k4_max)>6 ,'cl2',d$type)

	d$type = ifelse(d$eb3_k27_max < 6 & d$eb3_k27_max > 4 & log2(1+d$eb4_k4_max)>5.5,'cl3',d$type)

	d$type = ifelse( d$eb3_k27_max < 5 & log2(1+d$eb4_k4_max)>5.5,'cl4',d$type)
	#d$type_50k = ifelse(d$scale50k<.06,'seed','clust')
	cgd_ann = d	
	
	mod$cgdom_ann = cgd_ann
	mod$cgd_rpt = cgd_ann_mask
	message("saving CGDDs")
	saveRDS(mod$cgdom_ann, "./data/cgdom_ann_new.RDS")
	saveRDS(mod$cgd_rpt, "./data/cgd_rpt_new.RDS")

	return(mod)
}
pcg_icg_gen = function(mod)
{

	pcg_gen_cg_vtracks(20)
	pcg_gen_rmsk_vtrakcs()
	gen_k27_vt(mod)
	gen_k4_vt(mod)

	gvtrack.create("eb4dko_k27","jk.epipcg.pcg.CRJK_0408_k27me3_dko_to_dko_eb_d4", "sum")
	gvtrack.create("eb3dko_k27","jk.epipcg.pcg.CRJK_0407_k27me3_dko_to_dko_eb_d3", "sum")
	gvtrack.iterator("eb4dko_k27", sshift=-140, eshift=140)
	gvtrack.iterator("eb3dko_k27", sshift=-140, eshift=140)

	icg = gscreen("cg500 > 10", iterator=20)
	icg$ID = 1:nrow(icg)

	icg_rmsk = gextract(c("line_d", "ltr_d", "sine_d", "lowcomplex_d"), interval=icg, iterator=icg)

	icg_stat = gextract(c("cg200","cg300","cg400","cg500","cg600","cg700","gc500", "EB4_cnt", "EB4_cnt_k4", "eb4dko_k27", "eb3dko_k27"), intervals=icg, iterator=20)
	icg_stat[is.na(icg_stat)] = 0 

	f= icg_rmsk$line_d != 0 & icg_rmsk$ltr_d !=0

	icg$l = icg$end - icg$start
	icg$cg200 = tapply(icg_stat$cg200/200, icg_stat$intervalID, max)
	icg$cg300 = tapply(icg_stat$cg300/300, icg_stat$intervalID, max)
	icg$cg400 = tapply(icg_stat$cg400/400, icg_stat$intervalID, max)
	icg$cg500 = tapply(icg_stat$cg500/500, icg_stat$intervalID, max)
	icg$cg600 = tapply(icg_stat$cg600/600, icg_stat$intervalID, max)
	icg$cg700 = tapply(icg_stat$cg700/700, icg_stat$intervalID, max)
	icg$gc500 = tapply(icg_stat$gc500/500, icg_stat$intervalID, max)
	icg$max_eb4_k27 = tapply(icg_stat$EB4_cnt, icg_stat$intervalID, max)
	icg$max_eb4_k4 = tapply(icg_stat$EB4_cnt_k4, icg_stat$intervalID, max)
	icg$max_eb4_k27_dko = tapply(icg_stat$eb4dko_k27, icg_stat$intervalID, max)
	icg$max_eb3_k27_dko = tapply(icg_stat$eb3dko_k27, icg_stat$intervalID, max)

	icg_stat$maxdom_cg700 = icg[icg_stat$intervalID, "cg700"]
	icg_stat$dom_l = icg[icg_stat$intervalID, "l"]
	icg_stat_r = icg_stat[f[icg_stat$intervalID],]

	mod$meth = list()
	mod$meth$icg = icg
	mod$meth$icg_stat = icg_stat
	mod$meth$icg_stat_r = icg_stat_r
	return(mod)

}

pcg_gen_cg_vtracks = function(win_len=20)
{
	wl2 = win_len/2
	gvtrack.create("cg200", "seq.CG", "sum")
	gvtrack.iterator("cg200", sshift=-100+wl2, eshift=100-wl2)
	gvtrack.create("gc200", "seq.G_or_C", "sum")
	gvtrack.iterator("gc200", sshift=-100+wl2, eshift=100-wl2)
	gvtrack.create("cg300", "seq.CG", "sum")
	gvtrack.iterator("cg300", sshift=-150+wl2, eshift=150-wl2)
	gvtrack.create("cg400", "seq.CG", "sum")
	gvtrack.iterator("cg400", sshift=-200+wl2, eshift=200-wl2)
	gvtrack.create("gc500", "seq.G_or_C", "sum")
	gvtrack.iterator("gc500", sshift=-250+wl2, eshift=250-wl2)
	gvtrack.create("cg500", "seq.CG", "sum")
	gvtrack.iterator("cg500", sshift=-250+wl2, eshift=250-wl2)
	gvtrack.create("cg600", "seq.CG", "sum")
	gvtrack.iterator("cg600", sshift=-300+wl2, eshift=300-wl2)
	gvtrack.create("cg700", "seq.CG", "sum")
	gvtrack.iterator("cg700", sshift=-350+wl2, eshift=350-wl2)
}
pcg_gen_rmsk_vtrakcs = function()
{
	gvtrack.create("line_d", "intervs.global.rmsk_line", "distance")
	gvtrack.create("ltr_d", "intervs.global.rmsk_ltr", "distance")
	gvtrack.create("sine_d","intervs.global.rmsk_sine", "distance")
	gvtrack.create("lowcomplex_d", "intervs.global.rmsk_low_complexity", "distance")
}

pcg_gen_dkocnt_vtracks = function(win_len=20)
{
	#wl2 = win_len/2
	#gvtrack.create("eb4_k27","jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm", "avg")
	#gvtrack.create("eb4dko_k27","jk.epipcg.pcg.CRJK_0408_k27me3_dko_to_dko_eb_d4_norm", "avg")
	#gvtrack.create("eb3dko_k27","jk.epipcg.pcg.CRJK_0407_k27me3_dko_to_dko_eb_d3", "sum")
	##gvtrack.iterator("eb4_k27", sshift=-500, eshift=500)
	##gvtrack.iterator("eb4dko_k27", sshift=-500, eshift=500)
	##gvtrack.iterator("eb3dko_k27", sshift=-500, eshift=500)
	#
	#gvtrack.create("eb3_k4","jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3", "sum")
	#gvtrack.create("eb4dko_k4","jk.epipcg.pcg.CRJK_0413_k4me3_dko_to_dko_eb_d3", "sum")
	#gvtrack.create("eb3dko_k4","jk.epipcg.pcg.CRJK_0414_k4me3_dko_to_dko_eb_d4", "sum")
	#gvtrack.iterator("eb3_k4", sshift=-150+wl2, eshift=150-wl2)
	#gvtrack.iterator("eb4dko_k4", sshift=-150+wl2, eshift=150-wl2)
	#gvtrack.iterator("eb3dko_k4", sshift=-150+wl2, eshift=150-wl2)
	#
}
pcg_gen_gw_distrib = function(mod, foc_chroms = c(1,2,3))
{
	
	win_len = 200
    wl2 = win_len/2
    gvtrack.create("eb4_k27", "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm", 
        "avg")
    gvtrack.create("eb4dko_k27", "jk.epipcg.pcg.CRJK_0408_k27me3_dko_to_dko_eb_d4_norm", 
        "avg")
    gvtrack.create("eb3dko_k27", "jk.epipcg.pcg.CRJK_0407_k27me3_dko_to_dko_eb_d3", 
        "sum")
    gvtrack.iterator("eb3dko_k27", sshift = -150 + wl2, eshift = 150 - 
        wl2)
    gvtrack.create("eb3_k4", "jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3", 
        "sum")
    gvtrack.create("eb4dko_k4", "jk.epipcg.pcg.CRJK_0413_k4me3_dko_to_dko_eb_d3", 
        "sum")
    gvtrack.create("eb3dko_k4", "jk.epipcg.pcg.CRJK_0414_k4me3_dko_to_dko_eb_d4", 
        "sum")
    gvtrack.iterator("eb3_k4", sshift = -150 + wl2, eshift = 150 - 
        wl2)
    gvtrack.iterator("eb4dko_k4", sshift = -150 + wl2, eshift = 150 - 
        wl2)
    gvtrack.iterator("eb3dko_k4", sshift = -150 + wl2, eshift = 150 - 
        wl2)
	
	gen_k27_vt(mod)
   gen_k4_vt(mod)
	pcg_gen_cg_vtracks(200)
	#pcg_gen_dkocnt_vtracks(200)
	pcg_gen_rmsk_vtrakcs()
	gvtrack.create("epi_meth", "jk.epipcg.meth.meissN20.WT_Epi_rep1", "avg")
	gvtrack.iterator("epi_meth", sshift=-150, eshift=150)

        message("lm 10% test noX final model")
        seeds = readRDS('./data/lm_10test_noX_seeds_logist_dinucs_logist_fig1.rds')
        seeds$pred_seed_k27 = seeds$k27_pred
        seeds$pred_seed_k4 = seeds$k4_pred
    
    f_pcg = seeds$modality == "pcg"
    f_txg = seeds$modality == "txg"
    f_mix = seeds$modality == "mix"
	cgdann = mod$cgdom_ann[,c("chrom","start","end", "tss_strand")]
	colnames(cgdann)[4] = "strand"
	cgd_pcg = cgdann[f_pcg,]
	cgd_txg = cgdann[f_txg,]
	cgd_mix = cgdann[f_mix,]

	gvtrack.create("d_pcg", cgd_pcg, "distance", 1000)
	gvtrack.create("d_txg", cgd_txg, "distance", 1000)
	gvtrack.create("d_mix", cgd_mix, "distance", 1000)

	dst = gextract(c("line_d", "ltr_d","cg200", "gc200", "cg500","gc500","cg700",
							"2^eb4_k27 - 27.85", "2^eb4dko_k27 - 26.17", "eb3dko_k27",
							"eb3_k4", "eb4dko_k4", "eb3dko_k4",
							"EB4_cnt", "EB4_cnt_k4", "epi_meth", 
							"d_pcg", "d_txg", "d_mix"), 
							intervals=gintervals.all()[foc_chroms,], iterator=200,
							colnames = c("line_d", "ltr_d","cg200", "gc200", "cg500","gc500","cg700",
							"eb4_k27", "eb4dko_k27", "eb3dko_k27",
							"eb3_k4", "eb4dko_k4", "eb3dko_k4",
							"EB4_cnt", "EB4_cnt_k4", "epi_meth", 
							"d_pcg", "d_txg", "d_mix"))
	dst[is.na(dst)] = 0

#	dst$gcb = floor(dst$gc500/10)
#	dst$gcb = pmin(pmax(dst$gcb, 15), 34)-15
	dst$gcb = floor(dst$gc500/20) #0:25
	dst$gcb = pmin(pmax(dst$gcb, 9), 18)-9
	dst$cgb = pmin(dst$cg500, 29)
	dst$methb = pmin(pmax(floor(dst$epi_meth*10),5),9)

	min_cgd_d = pmin( abs(dst$d_pcg) , abs(dst$d_mix) , abs(dst$d_txg))
	min_cgd_d = pmax( min_cgd_d - 500, 0)
	dst$farb = ifelse(min_cgd_d < 300, 0, ifelse(min_cgd_d < 2000, 1,2))

	return(dst)
}

plt_border_5mc_cg_pcg = function(dst, pdf_fn = "figs/5mc_on_dko_025_95.pdf")
{
	if(!is.na(pdf_fn)) {
		pdf(pdf_fn, w=12, h = 5)
	}

	f_pcg = abs(dst$d_pcg) < abs(dst$d_txg) & abs(dst$d_pcg) < abs(dst$d_mix)
	f_txg = abs(dst$d_txg) < abs(dst$d_pcg) & abs(dst$d_txg) < abs(dst$d_mix)
	f_mix = !f_pcg & !f_txg

	prox_pcg = cut(dst$d_pcg[f_pcg],c(-1e+5,seq(-2900, 2900, 200),1e+5))
	prox_txg = cut(dst$d_txg[f_txg],c(-1e+5,seq(-2900, 2900, 200),1e+5))
	prox_mix = cut(dst$d_mix[f_mix],c(-1e+5,seq(-2900, 2900, 200),1e+5))
	prox_pcg = as.numeric(prox_pcg)
	prox_mix = as.numeric(prox_mix)
	prox_txg = as.numeric(prox_txg)

	cg_trend_pcg = tapply(dst$cg200[f_pcg]/2, prox_pcg, mean)
	cg_trend_txg = tapply(dst$cg200[f_txg]/2, prox_txg, mean)
	cg_trend_mix = tapply(dst$cg200[f_mix]/2, prox_mix, mean)

	gc_trend_pcg = tapply(dst$gc200[f_pcg]/200, prox_pcg, mean)
	gc_trend_txg = tapply(dst$gc200[f_txg]/200, prox_txg, mean)
	gc_trend_mix = tapply(dst$gc200[f_mix]/200, prox_mix, mean)

	meth_trend_pcg = tapply(dst$epi_meth[f_pcg], prox_pcg, mean)
	meth_trend_txg = tapply(dst$epi_meth[f_txg], prox_txg, mean)
	meth_trend_mix = tapply(dst$epi_meth[f_mix], prox_mix, mean)

	k27_trend_pcg = tapply(dst$eb4_k27[f_pcg], prox_pcg, mean)
	k27_trend_txg = tapply(dst$eb4_k27[f_txg], prox_txg, mean)
	k27_trend_mix = tapply(dst$eb4_k27[f_mix], prox_mix, mean)

	r_cov = sum(dst$eb4_k27)/sum(dst$eb4dko_k27)
	dkok27_trend_pcg = tapply(dst$eb4dko_k27[f_pcg], prox_pcg, mean)*r_cov
	dkok27_trend_txg = tapply(dst$eb4dko_k27[f_txg], prox_txg, mean)*r_cov
	dkok27_trend_mix = tapply(dst$eb4dko_k27[f_mix], prox_mix, mean)*r_cov

	kpcor_trend_pcg = rep(NA, t=max(prox_pcg, na.rm=T))
	diffmcor_trend_pcg = rep(NA, t=max(prox_pcg, na.rm=T))
	kpcor_trend_mix = rep(NA, t=max(prox_pcg, na.rm=T))
	diffmcor_trend_mix = rep(NA, t=max(prox_pcg, na.rm=T))
	kpcor_trend_txg = rep(NA, t=max(prox_pcg, na.rm=T))

	diffmcor_trend_txg = rep(NA, t=max(prox_pcg, na.rm=T))
	for(i in 1:max(prox_pcg, na.rm=T)) {
		f = f_pcg 
		fb = prox_pcg == i
		kpcor_trend_pcg[i] = cor(dst$eb4_k27[f][fb], dst$epi_meth[f][fb], u="p")
		diffmcor_trend_pcg[i] = cor(-log2(1+dst$eb4_k27[f][fb])+log2(1+dst$eb4dko_k27[f][fb]), dst$epi_meth[f][fb], u="p")
		f = f_txg 
		fb = prox_txg == i
		kpcor_trend_txg[i] = cor(dst$eb4_k27[f][fb], dst$epi_meth[f][fb], u="p")
		diffmcor_trend_txg[i] = cor(-log2(1+dst$eb4_k27[f][fb])+log2(1+dst$eb4dko_k27[f][fb]), dst$epi_meth[f][fb], u="p")
		f = f_mix 
		fb = prox_mix == i
		kpcor_trend_mix[i] = cor(dst$eb4_k27[f][fb], dst$epi_meth[f][fb], u="p")
		diffmcor_trend_mix[i] = cor(-log2(1+dst$eb4_k27[f][fb])+log2(1+dst$eb4dko_k27[f][fb]), dst$epi_meth[f][fb], u="p")
	}

	layout(t(matrix(1:6, nrow=3)))
	plot(cg_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red", main="CpG content")
	lines(cg_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	lines(cg_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	abline(v=12.5)
	abline(v=17.5)

	plot(gc_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red", 
				ylim=c(0.4,0.7), main="GC content")
	lines(gc_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	lines(gc_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	abline(v=12.5)
	abline(v=17.5)

	plot(meth_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red", main="wt 5mc")
	lines(meth_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	lines(meth_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	abline(v=12.5)
	abline(v=17.5)

	plot(-kpcor_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red",ylim=c(-0.05,0.6), main="anti cor of wt 5mc and wt k27")
	lines(-kpcor_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	lines(-kpcor_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	abline(v=12.5)
	abline(v=17.5)

	plot(diffmcor_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red",ylim=c(0, 0.7), main="cor of wt 5mc and k27 dko-wt")
	lines(diffmcor_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	lines(diffmcor_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	abline(v=12.5)
	abline(v=17.5)

	plot(log2(0.5+k27_trend_txg[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="red", ylim=c(2,7), main="k27 in wt and dko (dashed)")
	lines(log2(0.5+k27_trend_pcg[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	lines(log2(0.5+k27_trend_mix[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="black")

	lines(log2(0.5+dkok27_trend_txg[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="red", lty=2)
	lines(log2(0.5+dkok27_trend_pcg[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="blue", lty=2)
	lines(log2(0.5+dkok27_trend_mix[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="black", lty=2)
	abline(v=12.5)
	abline(v=17.5)

	if(!is.na(pdf_fn)) {
		dev.off()
	}

}

plt_border_5mc_cg_pcg_meth_plot_only = function(dst, pdf_fn = "figs/5mc_on_dko_025_95.pdf")
{
	if(!is.na(pdf_fn)) {
		pdf(pdf_fn, w=12, h = 5)
	}

	f_pcg = abs(dst$d_pcg) < abs(dst$d_txg) & abs(dst$d_pcg) < abs(dst$d_mix)
	f_txg = abs(dst$d_txg) < abs(dst$d_pcg) & abs(dst$d_txg) < abs(dst$d_mix)
	f_mix = !f_pcg & !f_txg

	prox_pcg = cut(dst$d_pcg[f_pcg],c(-1e+5,seq(-2900, 2900, 200),1e+5))
	prox_txg = cut(dst$d_txg[f_txg],c(-1e+5,seq(-2900, 2900, 200),1e+5))
	prox_mix = cut(dst$d_mix[f_mix],c(-1e+5,seq(-2900, 2900, 200),1e+5))
	prox_pcg = as.numeric(prox_pcg)
	prox_mix = as.numeric(prox_mix)
	prox_txg = as.numeric(prox_txg)

	cg_trend_pcg = tapply(dst$cg200[f_pcg]/2, prox_pcg, mean)
	cg_trend_txg = tapply(dst$cg200[f_txg]/2, prox_txg, mean)
	cg_trend_mix = tapply(dst$cg200[f_mix]/2, prox_mix, mean)

	gc_trend_pcg = tapply(dst$gc200[f_pcg]/200, prox_pcg, mean)
	gc_trend_txg = tapply(dst$gc200[f_txg]/200, prox_txg, mean)
	gc_trend_mix = tapply(dst$gc200[f_mix]/200, prox_mix, mean)

	meth_trend_pcg = tapply(dst$epi_meth[f_pcg], prox_pcg, mean)
	meth_trend_txg = tapply(dst$epi_meth[f_txg], prox_txg, mean)
	meth_trend_mix = tapply(dst$epi_meth[f_mix], prox_mix, mean)

	#k27_trend_pcg = tapply(dst$eb4_k27[f_pcg], prox_pcg, mean)
	#k27_trend_txg = tapply(dst$eb4_k27[f_txg], prox_txg, mean)
	#k27_trend_mix = tapply(dst$eb4_k27[f_mix], prox_mix, mean)

	#r_cov = sum(dst$eb4_k27)/sum(dst$eb4dko_k27)
	#dkok27_trend_pcg = tapply(dst$eb4dko_k27[f_pcg], prox_pcg, mean)*r_cov
	#dkok27_trend_txg = tapply(dst$eb4dko_k27[f_txg], prox_txg, mean)*r_cov
	#dkok27_trend_mix = tapply(dst$eb4dko_k27[f_mix], prox_mix, mean)*r_cov

	#kpcor_trend_pcg = rep(NA, t=max(prox_pcg, na.rm=T))
	#diffmcor_trend_pcg = rep(NA, t=max(prox_pcg, na.rm=T))
	#kpcor_trend_mix = rep(NA, t=max(prox_pcg, na.rm=T))
	#diffmcor_trend_mix = rep(NA, t=max(prox_pcg, na.rm=T))
	#kpcor_trend_txg = rep(NA, t=max(prox_pcg, na.rm=T))

	#diffmcor_trend_txg = rep(NA, t=max(prox_pcg, na.rm=T))
	##for(i in 1:max(prox_pcg, na.rm=T)) {
	##	f = f_pcg 
	##	fb = prox_pcg == i
	##	kpcor_trend_pcg[i] = cor(dst$eb4_k27[f][fb], dst$epi_meth[f][fb], u="p")
	##	diffmcor_trend_pcg[i] = cor(-log2(1+dst$eb4_k27[f][fb])+log2(1+dst$eb4dko_k27[f][fb]), dst$epi_meth[f][fb], u="p")
	##	f = f_txg 
	##	fb = prox_txg == i
	##	kpcor_trend_txg[i] = cor(dst$eb4_k27[f][fb], dst$epi_meth[f][fb], u="p")
	##	diffmcor_trend_txg[i] = cor(-log2(1+dst$eb4_k27[f][fb])+log2(1+dst$eb4dko_k27[f][fb]), dst$epi_meth[f][fb], u="p")
	##	f = f_mix 
	##	fb = prox_mix == i
	##	kpcor_trend_mix[i] = cor(dst$eb4_k27[f][fb], dst$epi_meth[f][fb], u="p")
	##	diffmcor_trend_mix[i] = cor(-log2(1+dst$eb4_k27[f][fb])+log2(1+dst$eb4dko_k27[f][fb]), dst$epi_meth[f][fb], u="p")
	##}
######
	#layout(t(matrix(1, nrow=1)))
	#plot(cg_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red", main="CpG content")
	#lines(cg_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	#lines(cg_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	#abline(v=12.5)
	#abline(v=17.5)

	#plot(gc_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red", 
	#			ylim=c(0.4,0.7), main="GC content")
	#lines(gc_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	#lines(gc_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	#abline(v=12.5)
	#abline(v=17.5)
    plt = function(){
	plot(meth_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red", main="wt 5mc")
	lines(meth_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	lines(meth_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	abline(v=12.5)
	abline(v=17.5)
    }
    save_baseR_to_ppt(plot_func = plt(),link_ppt = './figs/meth_on_cgd.pptx')
	#plot(-kpcor_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red",ylim=c(-0.05,0.6), main="anti cor of wt 5mc and wt k27")
	#lines(-kpcor_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	#lines(-kpcor_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	#abline(v=12.5)
	#abline(v=17.5)

	#plot(diffmcor_trend_txg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="red",ylim=c(0, 0.7), main="cor of wt 5mc and k27 dko-wt")
	#lines(diffmcor_trend_pcg[-1], type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	#lines(diffmcor_trend_mix[-1], type="l", lwd = 2, cex=0.4, pch=19, col="black")
	#abline(v=12.5)
	#abline(v=17.5)

	#plot(log2(0.5+k27_trend_txg[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="red", ylim=c(2,7), main="k27 in wt and dko (dashed)")
	#lines(log2(0.5+k27_trend_pcg[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="blue")
	#lines(log2(0.5+k27_trend_mix[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="black")

	#lines(log2(0.5+dkok27_trend_txg[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="red", lty=2)
	#lines(log2(0.5+dkok27_trend_pcg[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="blue", lty=2)
	#lines(log2(0.5+dkok27_trend_mix[-1]), type="l", lwd = 2, cex=0.4, pch=19, col="black", lty=2)
	#abline(v=12.5)
	#abline(v=17.5)

	if(!is.na(pdf_fn)) {
		#dev.off()
	}

}



plt_far_5mc_cg_pcg = function(mod,  pdf_fn = "figs/far_5mc_pcg.pdf")
{
	#pcg_gen_dkocnt_vtracks(200)
	gvtrack.create("epi_meth", "jk.epipcg.meth.meissN20.WT_Epi_rep1", "avg")
	gvtrack.iterator("epi_meth", sshift=-150, eshift=150)
	win_len = 200
wl2 = win_len/2
	gvtrack.create("eb4_k27","jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm", "avg")
	gvtrack.create("eb4dko_k27","jk.epipcg.pcg.CRJK_0408_k27me3_dko_to_dko_eb_d4_norm", "avg")
	gvtrack.create("eb3dko_k27","jk.epipcg.pcg.CRJK_0407_k27me3_dko_to_dko_eb_d3", "sum")
	#gvtrack.iterator("eb4_k27", sshift=-150+wl2, eshift=150-wl2)
	#gvtrack.iterator("eb4dko_k27", sshift=-150+wl2, eshift=150-wl2)
	gvtrack.iterator("eb3dko_k27", sshift=-150+wl2, eshift=150-wl2)

	gvtrack.create("eb3_k4","jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3", "sum")
	gvtrack.create("eb4dko_k4","jk.epipcg.pcg.CRJK_0413_k4me3_dko_to_dko_eb_d3", "sum")
	gvtrack.create("eb3dko_k4","jk.epipcg.pcg.CRJK_0414_k4me3_dko_to_dko_eb_d4", "sum")
	gvtrack.iterator("eb3_k4", sshift=-150+wl2, eshift=150-wl2)
	gvtrack.iterator("eb4dko_k4", sshift=-150+wl2, eshift=150-wl2)
	gvtrack.iterator("eb3dko_k4", sshift=-150+wl2, eshift=150-wl2)
	cgt = readRDS('./data/cg_trace_mm10.rds')
	if (file.exists('./data/cgt_met_fig4.rds')) {
        cgt_met = readRDS('./data/cgt_met_fig4.rds')
       
    }else {
cgt_met = gextract(c("epi_meth", "2^eb4_k27 - 27.85", "eb3dko_k27", 
        "2^eb4dko_k27 - 26.17"), intervals = cgt, iterator = cgt, 
        colnames = c("epi_meth", "eb4_k27", "eb3dko_k27", "eb4dko_k27"))
saveRDS(cgt_met,'./data/cgt_met_fig4.rds')
        }
	
	#cgt_met = gextract(c("epi_meth", "2^eb4_k27 - 27.85","eb3dko_k27", "2^eb4dko_k27 - 26.17"), intervals = cgt, iterator = cgt,
	#colnames= c("epi_meth", "eb4_k27","eb3dko_k27", "eb4dko_k27"))
#cgt_met$eb4_k27     <- 2^cgt_met$eb4_k27 - 6
#cgt_met$eb4dko_k27  <- 2^cgt_met$eb4dko_k27 - 6

	#cgt_met$eb4_k27 = cgt_met$eb4_k27
	#cgt_met$eb4dko_k27 = cgt_met$eb4dko_k27
	f_far = abs(cgt$d_pcg)>3000 & abs(cgt$d_mix)>3000  & abs(cgt$d_txg)>3000
	#f_far = abs(cgt$d_pcg)>3000 #& abs(cgt$d_txg)>3000 & abs(cgt$d_pcg)>3000
	cg_rng = c(seq(0,0.1,0.005),0.3)
	n_cg = tapply(cgt_met$eb4_k27, cut(cgt$CG, cg_rng), length)
	p_cg = n_cg/sum(n_cg)
	nk27_cg = tapply(cgt_met$eb4_k27 , cut(cgt$CG, cg_rng),sum)
	pk27_cg = nk27_cg/sum(nk27_cg)
	nk27dko_cg = tapply(cgt_met$eb4dko_k27 , cut(cgt$CG, cg_rng),sum)
	pk27dko_cg = nk27dko_cg/sum(nk27dko_cg)

	far_cg_rng = cut(cgt$CG[f_far], c(seq(0,0.04,0.005),0.3))
	farn_cg = tapply(cgt_met$eb4_k27[f_far], far_cg_rng, length)
	farp_cg = farn_cg/sum(farn_cg)
	farnk27_cg = tapply(cgt_met$eb4_k27[f_far], far_cg_rng, sum)
	farpk27_cg = farnk27_cg/sum(farnk27_cg)
	farnk27dko_cg = tapply(cgt_met$eb4dko_k27[f_far], far_cg_rng, sum)
	farpk27dko_cg = farnk27dko_cg/sum(farnk27dko_cg)

	f_near = !f_far
	near_cg_rng = cut(cgt$CG[f_near], c(seq(0,0.1,0.005),0.3))
	nearn_cg = tapply(cgt_met$eb4_k27[f_near], near_cg_rng, length)
	nearp_cg = nearn_cg/sum(nearn_cg)
	nearnk27_cg = tapply(cgt_met$eb4_k27[f_near], near_cg_rng, sum)
	nearpk27_cg = nearnk27_cg/sum(nearnk27_cg)
	nearnk27dko_cg = tapply(cgt_met$eb4dko_k27[f_near], near_cg_rng, sum)
	nearpk27dko_cg = nearnk27dko_cg/sum(nearnk27dko_cg)


	#pdf("figs/near_k27_dko_cg.pdf", h=5, w=5)
    plt1 = function(){
	plot(0.0025+cg_rng[1:21],nearnk27_cg/nearn_cg, lwd=1, type="l", ylim=c(0,70), col="darkblue", xlab="CG", ylab="mean H3k27me3")
	lines(0.0025+cg_rng[1:21], nearnk27dko_cg/nearn_cg, lwd=1, type="b", col="darkblue")
        }#1.4*
    save_baseR_to_ppt(plot_func =plt1() ,link_ppt ="figs/near_k27_dko_cg.pptx" )
	#dev.off()
	#pdf("figs/far_k27_dko_cg.pdf", h=5, w=5)
    plt2 = function(){
	plot(0.0025+cg_rng[1:9],farnk27_cg/farn_cg, lwd=1, type="l", xlim=c(0,0.035),ylim=c(0,35), xlab="CG", col="darkblue", ylab="mean H3k27me3")
	lines(0.0025+cg_rng[1:9],farnk27dko_cg/farn_cg, lwd=1, type="b", col="darkblue")}#1.4*
    save_baseR_to_ppt(plot_func =plt2() ,link_ppt ="figs/far_k27_dko_cg.pptx" )
        
    
	#dev.off()
	

#why does the k27 cor so poorly with meth while pred looks better?

	#boxplot(split(cgt_met$epi_meth[f], cut(cgt$lk27_1k[f], seq(-0.1,7.5,0.2))),add=F, col="blue",boxwex=0.4, pch=19, cex=0.3, las=2)
	#boxplot(split(cgt_met$epi_meth[f], cut(cgt$pred[f], seq(-0.1,7.5,0.2))),add=F, col="blue",boxwex=0.4, pch=19, cex=0.3, las=2)
}
pcg_gen_dkocnt_vtracks_2i = function (win_len = 20) 
{
   wl2 = win_len/2
   gvtrack.create("eb4_k27", "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_norm", 
       "avg")
   gvtrack.create("eb4dko_k27", "jk.epipcg.pcg.CRJK_0408_k27me3_dko_to_dko_eb_d4_norm", 
       "avg")
   gvtrack.create("eb3dko_k27", "jk.epipcg.pcg.CRJK_0407_k27me3_dko_to_dko_eb_d3_norm", 
       "avg")
   gvtrack.create("ES2i_cnt", "jk.epipcg.pcg.CRJK_0232_k27me3_es_2i_dmso_norm", 
       "avg")
   gvtrack.create("ES_cnt", "jk.epipcg.pcg.CRJK_0211_k27me3_es_50k_norm", 
       "avg")
   #gvtrack.iterator("eb4_k27", sshift = -150 + wl2, eshift = 150 - 
    #   wl2)
   #gvtrack.iterator("eb4dko_k27", sshift = -150 + wl2, eshift = 150 - 
   #    wl2)
   #gvtrack.iterator("eb3dko_k27", sshift = -150 + wl2, eshift = 150 - 
    #   wl2)
     #  gvtrack.iterator("ES2i_cnt", sshift = -150 + wl2, eshift = 150 - 
     #  wl2)
      # gvtrack.iterator("ES_cnt", sshift = -150 + wl2, eshift = 150 - 
     #  wl2)
   gvtrack.create("eb3_k4", "jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3", 
       "sum")
   gvtrack.create("eb3dko_k4", "jk.epipcg.pcg.CRJK_0413_k4me3_dko_to_dko_eb_d3", 
       "sum")
   gvtrack.create("eb4dko_k4", "jk.epipcg.pcg.CRJK_0414_k4me3_dko_to_dko_eb_d4", 
       "sum")
   gvtrack.iterator("eb3_k4", sshift = -150 + wl2, eshift = 150 - 
       wl2)
   gvtrack.iterator("eb4dko_k4", sshift = -150 + wl2, eshift = 150 - 
       wl2)
    gvtrack.iterator("eb3dko_k4", sshift = -150 + wl2, eshift = 150 - 
        wl2)
}
plt_far_5mc_cg_pcg_k4 = function (mod, pdf_fn = "figs/far_5mc_pcg_2i.pdf") 
{
    pcg_gen_dkocnt_vtracks_2i(200)
    gvtrack.create("epi_meth", "jk.epipcg.meth.meissN20.WT_Epi_rep1", 
        "avg")
    gvtrack.iterator("epi_meth", sshift = -150, eshift = 150)
    cgt = readRDS('./data/cg_trace_mm10.rds')
	
	
if (file.exists('./data/cgt_met_k4_fig4.rds')) {
        cgt_met = readRDS('./data/cgt_met_k4_fig4.rds')
       
    }else {
cgt_met = gextract(c("epi_meth", "eb3_k4", "eb3dko_k4", "ES_cnt"), 
        intervals = cgt, iterator = cgt)
saveRDS(cgt_met,'./data/cgt_met_k4_fig4.rds')
        }
	
    #cgt_met = gextract(c("epi_meth", "eb3_k4", "eb3dko_k4", 
    #    "ES_cnt"), intervals = cgt, iterator = cgt)
    f_far = abs(cgt$d_pcg) > 3000 & abs(cgt$d_txg) > 3000 & abs(cgt$d_pcg) > 
        3000
    cg_rng = c(seq(0, 0.1, 0.005), 0.3)
    n_cg = tapply(cgt_met$eb3_k4, cut(cgt$CG, cg_rng), length)
    p_cg = n_cg/sum(n_cg)
    nk27_cg = tapply(cgt_met$eb3_k4, cut(cgt$CG, cg_rng), sum)
    pk27_cg = nk27_cg/sum(nk27_cg)
    nk27dko_cg = tapply(cgt_met$eb3dko_k4, cut(cgt$CG, cg_rng), 
        sum)
    pk27dko_cg = nk27dko_cg/sum(nk27dko_cg)
    far_cg_rng = cut(cgt$CG[f_far], c(seq(0, 0.04, 0.005), 0.3))
    farn_cg = tapply(cgt_met$eb3_k4[f_far], far_cg_rng, length)
    farp_cg = farn_cg/sum(farn_cg)
    farnk27_cg = tapply(cgt_met$eb3_k4[f_far], far_cg_rng, sum)
    farpk27_cg = farnk27_cg/sum(farnk27_cg)
    farnk27dko_cg = tapply(cgt_met$eb3dko_k4[f_far], far_cg_rng, 
        sum)
    farpk27dko_cg = farnk27dko_cg/sum(farnk27dko_cg)
    f_near = !f_far
    near_cg_rng = cut(cgt$CG[f_near], c(seq(0, 0.1, 0.005), 0.3))
    nearn_cg = tapply(cgt_met$eb3_k4[f_near], near_cg_rng, length)
    nearp_cg = nearn_cg/sum(nearn_cg)
    nearnk27_cg = tapply(cgt_met$eb3_k4[f_near], near_cg_rng, 
        sum)
    nearpk27_cg = nearnk27_cg/sum(nearnk27_cg)
    nearnk27dko_cg = tapply(cgt_met$eb3dko_k4[f_near], near_cg_rng, 
        sum)
    nearpk27dko_cg = nearnk27dko_cg/sum(nearnk27dko_cg)
    plt1 = function() {
        plot(0.0025 + cg_rng[1:21], nearnk27_cg/nearn_cg, lwd = 1,col = "darkred", 
            type = "l", ylim = c(0, 200), xlab = "CG", ylab = "mean H3k27me3")
        lines(0.0025 + cg_rng[1:21],  nearnk27dko_cg/nearn_cg, 
            lwd = 1, type = "b", col = "darkred")
    }
    save_baseR_to_ppt(plot_func = plt1(), link_ppt = "figs/near_k27_dko_cg_k4.pptx")
    plt2 = function() {
        plot(0.0025 + cg_rng[1:9], farnk27_cg/farn_cg, lwd = 1,col = "darkred", 
            type = "l", ylim = c(0, 85), xlab = "CG", ylab = "mean H3k27me3")
        lines(0.0025 + cg_rng[1:9],  farnk27dko_cg/farn_cg, 
            lwd = 1, type = "b", col = "darkred")
    }
    save_baseR_to_ppt(plot_func = plt2(), link_ppt = "figs/far_k27_dko_cg_k4.pptx")
}
plt_regionAT_ds_cg_meth = function(mod, chr , st , end, ds_k27 =TRUE, marg5=2e+3, marg3=2e+3, mark=c(), wk4=F, k4_max_fact=1)
{
    tracks_k27 = c('jk.epipcg.pcg.CRJK_0363_k27me3_eb_j1_d3_a',
 
  'jk.epipcg.pcg.CRJK_0365_k27me3_eb_dko_d3_a'
)
	ndx = mod$epi_tracks_cnt
    ndx = ndx[ ndx$track_k27%in% tracks_k27,]
    ndx = ndx[ match(tracks_k27,ndx$track_k27),]
    short_namesk27 = ndx$short_name
	for(i in 1:nrow(ndx)) {
		gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "sum") 
		gvtrack.iterator(ndx$short_name[i], sshift = -140, eshift= 140)
	}
    cov = mod$cnt_track_cov
    cov = cov[ short_namesk27,]
    ###hox
    hoxs = grep('Hox',mod$tss$geneSymbol,v=T)

    hox_int = mod$tss[ mod$tss$geneSymbol %in% hoxs,]
    
    hox_int$start = hox_int$start - 500
    hox_int$end = hox_int$end + 500

    hox_int = (gintervals.canonic(hox_int))

    h = gextract(c(short_namesk27),intervals = hox_int,iterator =hox_int)

    hs = colSums(h[ ,short_namesk27])
    hs_df = as.data.frame(hs)

    cov$hox = hs_df$hs
    ####
    
	prof = gextract(ndx$short_name, intervals=gintervals(chr, st-marg5, end+marg3), iterator=20)
    
    prof[is.na(prof)] = 0
    ylim = max(apply((prof[, short_namesk27]),2,max))
    ####ds
    ribo = cov$hox

    tr_dsamp_r = min(ribo)/ribo
    names(tr_dsamp_r) = short_namesk27
    prof_rel = prof
    if (ds_k27 ==T){
    for (nm in short_namesk27) {
      if(!is.null(tr_dsamp_r[nm])) {
        to_ds = matrix(prof_rel[,nm],ncol=1)
        r_ds = tr_dsamp_r[nm]
        message("into dsamp, ", nm, " targ ", floor(sum(to_ds)*r_ds))
        post_ds = downsample_matrix(to_ds,  target_n = floor(sum(to_ds)*r_ds), seed=19)
        prof_rel[,nm] = as.numeric(post_ds)
        }
    }
    } else {
        for (nm in short_namesk27) {
      
        to_ds = matrix(prof_rel[,nm],ncol=1)
        
        
        post_ds = (to_ds/sum(to_ds))*1e6
        prof_rel[,nm] = post_ds
        
    }
        
    }
    
    prof_rel = as.data.frame(prof_rel)
    #####
    prof = prof_rel
    
  
    gvtrack.create('WT_Epi_rep1', 'jk.epipcg.meth.meissN20.WT_Epi_rep1', "avg")
    gvtrack.iterator('WT_Epi_rep1', sshift=-100, eshift=100)
 
   

    #m = (c("WT_Epi_rep1",NA,NA,NA,NA,NA,"WT_Epi_rep1",NA,1,1,1,1,1,1))
    #ndx = rbind(ndx,m)
    c = (c("seq.CG_500_mean_new",NA,NA,NA,NA,NA,"seq.CG_500_mean_new",NA,1,1,1,1,1,1))
    ndx = rbind(ndx,c)
    prof2 = gextract(c("WT_Epi_rep1","seq.CG_500_mean_new"), intervals=gintervals(chr, st-marg5, end+marg3), iterator=20,
                    colnames= c('meth','cg'))
    prof2[is.na(prof2)] = 0
    #prof$meth = prof2$meth
    prof$cg = prof2$cg
    prof$intervalID = NULL
#	if(wk4) {
#        ndxK4 = ndx[ !is.na(ndx$track_k4),]
#        	for(i in 1:nrow(ndxK4)) {
#		gvtrack.create(ndxK4$short_name_k4[i], ndxK4$track_k4[i], "sum") 
#		gvtrack.iterator(ndxK4$short_name_k4[i], sshift = -140, eshift= 140)
#	}
		
#        k4_tns = ndxK4$short_name_k4
#		prof_k4 = gextract(k4_tns, intervals=gintervals(chr, st-marg5, end+marg3), iterator=20)
#	}
	layout(matrix(1:nrow(ndx),ncol=1),h=c(1.4, rep(1,nrow(ndx)-2), 1.4))
	for(i in 1:nrow(ndx)) {
		if(i == nrow(ndx)) {
			par(mar=c(2,4,0.8,4))
		} else if(i==1) {
			par(mar=c(0,4,2,4))
		} else {
			par(mar=c(0,4,0.8,4))
		}
		if ( i < nrow(ndx)){c((prof[,3+i]), rep(0, length(prof$start)))
                plot(prof$start, (prof[,3+i]), pch=19, type="l", 
                    lwd=3, ylim = c(0,ylim), xaxt=ifelse(i==nrow(ndx), 's', 'n'), 
                    main = ndx$short_name[i], ylab=NA)   
                polygon(c(prof$start, rev(prof$start)), c((prof[,3+i]), rep(0, length(prof$start))), 
                    col="black", border=NA)
        } #else if(i==nrow(ndx)-1){
           #     plot(prof$start, (prof[,3+i]), pch=19, type="l", 
           #         lwd=3, ylim = c(0,1), xaxt=ifelse(i==nrow(ndx), 's', 'n'), 
            #        main = ndx$short_name[i], ylab=NA)    
                #polygon(c(prof$start, rev(prof$start)), c(log2(1+prof[,3+i]), rep(0, length(prof$start))), 
                #    col="black", border=NA)
         else {
                plot(prof$start, (prof[,3+i]), pch=19, type="l", 
                    lwd=3, ylim = c(0,.06), xaxt=ifelse(i==nrow(ndx), 's', 'n'), 
                    main = ndx$short_name[i], ylab=NA)  
               # polygon(c(prof$start, rev(prof$start)), c(log2(1+prof[,3+i]), rep(0, length(prof$start))), 
                #    col="black", border=NA)
        }

#		if(wk4 & !is.na(ndx$track_k4[i])) {
#			tmax_k4 = mod$k4_track_thresh_cnt[ mod$k4_track_thresh_cnt$`X0.5`==ndx$short_name_k4[i],8]* k4_max_fact
#			message("added ", k4_max_fact, " to k4 thresh", " ", tmax_k4, " ", tmax)
#			k4_vs = prof_k4[,ndx$short_name_k4[i]]*tmax/tmax_k4
#			lines(prof_k4$start, pmin(k4_vs,tmax), pch=19, type="l", lwd=2, col="red")
#		}
        #cg = gextract("seq.CG_500_mean_new",intervals = gintervals(chr, st-marg5, end+marg3),iterator = 20,colnames = 'cg')
        
        
		abline(v=st)
		abline(v=end)
#        abline(h=0.02)
#        abline(h=0.04)
		if(!is.null(mark) & length(mark)>0) {
			for(x in mark) {
				abline(v=x, col="blue", lwd=2)
			}
		}
	}
  
}
plt_cg_h_distribs = function(mod)
{
	cgd = mod$cgdom_ann
	cgdr = mod$cgd_rpt

	pdf("figs/cgd_length.pdf", w=6, h=4)
	plot(density(log2(cgd$l)), lwd=3, ylim=c(0,1.4), xlim=c(8,12), col="darkgreen")
	lines(density(log2(cgdr$l[cgdr$line>0.5])), lwd=2, col="darkgray")
	lines(density(log2(cgdr$l[cgdr$ltr>0.5])), lwd=2, col="darkgray",lty=2)
	dev.off()

	lbin = as.numeric(cut(cgd$l,c(0,seq(500,3000,250))))
#plot various features vs length
	pdf("figs/cgd_k27.pdf", w=6,h=4)
	boxplot(split(log2(10+cgd$eb4_k27_max), lbin), col="darkblue")
	dev.off()
	pdf("figs/cgd_k4.pdf", w=6,h=4)
	boxplot(split(log2(10+cgd$eb4_k4_max), lbin), col="darkred")
	dev.off()
	pdf("figs/cgd_atac.pdf", w=6,h=4)
	boxplot(split(log2(10+cgd$atac_max), lbin), col="darkgray")
	dev.off()
	pdf("figs/cgd_cg.pdf", w=6,h=4)
	boxplot(split(cgd$cg_max, lbin), col="gray")
	dev.off()
	pdf("figs/cgd_gc.pdf", w=6,h=4)
	boxplot(split(cgd$gc_max, lbin), col="gray")
	dev.off()
	pdf("figs/cgd_ratio_cg.pdf", w=6,h=4)
	boxplot(split(cgd$cg_max*4/(cgd$gc_max**2), lbin), col="gray")
	dev.off()
	pdf("figs/cgd_simp.pdf", w=6,h=4)
	boxplot(split(cgd$simp+cgd$low_complex, lbin), col="gold")
	dev.off()
	pdf("figs/cgd_tss.pdf", w=6,h=4)
	boxplot(split(log2(100+abs(cgd$tss_dist)), lbin), col="gold")
	dev.off()
	pdf("figs/cgd_exon.pdf", w=6,h=4)
	boxplot(split(cgd$exon, lbin), col="gold")
	dev.off()

	gen_trend = function(strat, w=401) 
	{
		x = rollmean(sort(strat), w, f='e')
		y = rollmean(lk27[order(strat)], w, f='e')
		return(list(x=x, y=y))
	}

	lk27 = log2(1+cgd$eb4_k27_max)
	lk4 = log2(1+cgd$eb4_k4_max)
	k27_rat = cgd$eb4_k27_rat
	k4_rat = cgd$eb4_k4_rat
	l_class = as.numeric(cut(cgd$l, c(-1,500,750,1000,1e+6)))
	pdf("figs/cgd_k4k27_strat.pdf", w=12.5,h=10)
	strat = lk4
	xlim = c(0,8)
	#strat = cgd$gc_max
	#xlim = c(0.4,0.9)
	#strat = log2(1+cgd$atac)
	#xlim = c(5, 14)
	tss_class = as.numeric(cut(abs(cgd$tss_dist),c(-1,0,500,2000,8000,1e+8)))
	layout(t(matrix(1:20,ncol=4)))
	for(li in 1:4) {
		par(mar=c(2,1,1,1))
		f = l_class==li
		for(i in 1:5) { 
			n0 = sum(f & tss_class==i & lk27<6 & lk4 < 6)
			n1 = sum(f & tss_class==i & lk27<4 & lk4>6)
			n2 = sum(f & tss_class==i & lk27<6 & lk4>6) - n1
			n3 = sum(f & tss_class==i) - n1 - n2 - n0
			n = n1+n2+n3+n0
			f1 = round(n1/n,2)
			f2 = round(n2/n,2)
			f3 = round(n3/n,2)
			f0 = round(n0/n,2)
			plot(strat[f], lk27[f], pch=19, cex=0.3, col="gray", 
					xlim = xlim, ylim=c(0,8),
					main=sprintf("%s (pcg %s, %s, k4 %s, no %s)", n, f3,f2,f1,f0),
					xaxt='n',yaxt='n') 

			points(strat[f & tss_class==i],lk27[f & tss_class==i], pch=19,cex=0.7,col="red") 
			abline(h=6, lwd=2)
			segments(x0=c(6,6),x1=c(9,6), y0=c(4,6), y1=c(4,-2), lwd=2)
		}
	}
	dev.off()

	f_out = lk27 < 6 & lk4 < 6

	png("figs/k4cov_k27max.png", w=400,h=1200)
	layout(matrix(1:3,nrow=3))
	for(li in 2:4) {
		par(mar = c(2,2,1,1))
		f = l_class==li & !f_out
		n_ambig = sum(f & lk27>6 & cgd$eb4_k4_cov_1>0.1)
		n_pcg = sum(f & lk27>6 & cgd$eb4_k4_cov_1==0)
		n_wambig = sum(f & lk27>6 & cgd$eb4_k4_cov_1<0.1) - n_pcg
		plot(cgd$eb4_k4_cov_1[f], lk27[f], xaxt='n', yaxt='n', xlim=c(0,1),ylim=c(0,8), pch=19, cex=0.5, main=sprintf("n_ambig = %s,%s,%s", n_ambig, n_wambig, n_pcg))
	}
	dev.off()

	for(use_1_percentile in c(F,T)) {

		if(use_1_percentile) {
			pdf("figs/cgd_biv_dist_1perc.pdf", w=2.5,h=10)
			k27cov = cgd$eb4_k27_cov_1
			k4cov = cgd$eb4_k4_cov_1
			bivcov = cgd$eb4_biv_cov_1
			dst = cgd[,c("eb4_k27_cov_1", "eb4_biv_cov_1", "eb4_k4_cov_1")]
		} else {
			pdf("figs/cgd_biv_dist_0.5perc.pdf", w=2.5,h=10)
			k27cov = cgd$eb4_k27_cov
			k4cov = cgd$eb4_k4_cov
			bivcov = cgd$eb4_biv_cov
			dst = cgd[,c("eb4_k27_cov", "eb4_biv_cov", "eb4_k4_cov")]
		}
		layout(matrix(1:4, ncol=1))
		f_out = lk27 < 6 & lk4 < 6
		dst_n = dst/rowSums(dst)
		for(li in 1:4) {
			par(mar = c(2,2,1,1))
			f = l_class==li & !f_out

			dst_i = dst_n[f,]

			n = nrow(dst_i)
			ord = order((1e-3+dst_i[,1]+dst_i[,2])/(1e-3+dst_i[,3]+dst_i[,2]))
			dst_i = dst_i[ord,]
			plot(-0,-0, xlim=c(0,1), ylim=c(1,nrow(dst_i)), xaxt='n',yaxt='n')
			polygon(x = c(dst_i[,3], rep(0,t=n)), y=c(1:n, n:1), col="darkred", border=NA)
			polygon(x = c(dst_i[,3], rev(dst_i[,2]+ dst_i[,3])), y=c(1:n, n:1), col="black", border=NA)
			polygon(x = c(rep(1,t=n), rev(dst_i[,3]+dst_i[,2])), y=c(1:n, n:1), col="darkblue", NA)
		}
		dev.off()

		if(use_1_percentile) {
			pdf("figs/cgd_biv_dens_1perc.pdf", w=2.5,h=10)
		} else {
			pdf("figs/cgd_biv_dens_0.5perc.pdf", w=2.5,h=10)
		}
		layout(t(matrix(1:8, nrow=2)))
		l_class = as.numeric(cut(cgd$l, c(-1,500,750,1000,1e+6)))
		for(li in 1:4) {
			par(mar = c(2,2,1,1))
			f = l_class==li & !f_out

			dst_i = dst[f,]
			hist(dst_i[,3]+dst_i[,2], freq=F, col="darkred", ylim=c(0,5), main=NA,breaks=seq(-0.1,1,0.1))
			hist(dst_i[,1]+dst_i[,2], freq=F, col="darkblue", ylim=c(0,6), main=NA,breaks=seq(-0.1,1,0.1))
		}
		dev.off()
	}
}

plt_cg_interm_distribs = function(mod)
{
	cgd = mod$cgdom_ann_i
	cgdr = mod$cgd_rpt_i

	pdf("figs/cgd_length_i.pdf", w=6, h=4)
	plot(density(log2(cgd$l)), lwd=3, ylim=c(0,1.4), xlim=c(8,12), col="darkgreen")
	lines(density(log2(cgdr$l[cgdr$line>0.5])), lwd=2, col="darkgray")
	lines(density(log2(cgdr$l[cgdr$ltr>0.5])), lwd=2, col="darkgray",lty=2)
	dev.off()

	lbin = as.numeric(cut(cgd$l,c(0,seq(250,3000,250))))
#plot various features vs length
	pdf("figs/cgd_k27_i.pdf", w=6,h=4)
	boxplot(split(log2(10+cgd$eb4_k27_max), lbin), col="darkblue")
	dev.off()
	pdf("figs/cgd_k4_i.pdf", w=6,h=4)
	boxplot(split(log2(10+cgd$eb4_k4_max), lbin), col="darkred")
	dev.off()
	pdf("figs/cgd_atac_i.pdf", w=6,h=4)
	boxplot(split(log2(10+cgd$atac_max), lbin), col="darkgray")
	dev.off()
	pdf("figs/cgd_cg_i.pdf", w=6,h=4)
	boxplot(split(cgd$cg_max, lbin), col="gray")
	dev.off()
	pdf("figs/cgd_gc_i.pdf", w=6,h=4)
	boxplot(split(cgd$gc_max, lbin), col="gray")
	dev.off()
	pdf("figs/cgd_ratio_cg_i.pdf", w=6,h=4)
	boxplot(split(cgd$cg_max*4/(cgd$gc_max**2), lbin), col="gray")
	dev.off()
	pdf("figs/cgd_simp_i.pdf", w=6,h=4)
	boxplot(split(cgd$simp+cgd$low_complex, lbin), col="gold")
	dev.off()
	pdf("figs/cgd_tss_i.pdf", w=6,h=4)
	boxplot(split(log2(100+abs(cgd$tss_dist)), lbin), col="gold")
	dev.off()
	pdf("figs/cgd_exon_i.pdf", w=6,h=4)
	boxplot(split(cgd$exon, lbin), col="gold")
	dev.off()

	gen_trend = function(strat, w=401) 
	{
		x = rollmean(sort(strat), w, f='e')
		y = rollmean(lk27[order(strat)], w, f='e')
		return(list(x=x, y=y))
	}

	lk27 = log2(1+cgd$eb4_k27_max)
	lk4 = log2(1+cgd$eb4_k4_max)
	k27_rat = cgd$eb4_k27_rat
	k4_rat = cgd$eb4_k4_rat
	l_class = as.numeric(cut(cgd$l, c(-1,500,750,1000,1e+6)))
	pdf("figs/cgd_k4k27_strat_i.pdf", w=12.5,h=10)
	strat = lk4
	xlim = c(0,8)
	#strat = cgd$gc_max
	#xlim = c(0.4,0.9)
	#strat = log2(1+cgd$atac)
	#xlim = c(5, 14)
	tss_class = as.numeric(cut(abs(cgd$tss_dist),c(-1,0,500,2000,8000,1e+8)))
	layout(t(matrix(1:20,ncol=4)))
	for(li in 1:4) {
		par(mar=c(2,1,1,1))
		f = l_class==li
		for(i in 1:5) { 
			n0 = sum(f & tss_class==i & lk27<6 & lk4 < 6)
			n1 = sum(f & tss_class==i & lk27<4 & lk4>6)
			n2 = sum(f & tss_class==i & lk27<6 & lk4>6) - n1
			n3 = sum(f & tss_class==i) - n1 - n2 - n0
			n = n1+n2+n3+n0
			f1 = round(n1/n,2)
			f2 = round(n2/n,2)
			f3 = round(n3/n,2)
			f0 = round(n0/n,2)
			plot(strat[f], lk27[f], pch=19, cex=0.3, col="gray", 
					xlim = xlim, ylim=c(0,8),
					main=sprintf("%s (pcg %s, %s, k4 %s, no %s)", n, f3,f2,f1,f0),
					xaxt='n',yaxt='n') 

			points(strat[f & tss_class==i],lk27[f & tss_class==i], pch=19,cex=0.7,col="red") 
			abline(h=6, lwd=2)
			segments(x0=c(6,6),x1=c(9,6), y0=c(4,6), y1=c(4,-2), lwd=2)
		}
	}
	dev.off()

	f_out = lk27 < 6 & lk4 < 6

	png("figs/k4cov_k27max_i.png", w=400,h=1200)
	layout(matrix(1:3,nrow=3))
	for(li in 2:4) {
		par(mar = c(2,2,1,1))
		f = l_class==li & !f_out
		n_ambig = sum(f & lk27>6 & cgd$eb4_k4_cov_1>0.1)
		n_pcg = sum(f & lk27>6 & cgd$eb4_k4_cov_1==0)
		n_wambig = sum(f & lk27>6 & cgd$eb4_k4_cov_1<0.1) - n_pcg
		plot(cgd$eb4_k4_cov_1[f], lk27[f], xaxt='n', yaxt='n', xlim=c(0,1),ylim=c(0,8), pch=19, cex=0.5, main=sprintf("n_ambig = %s,%s,%s", n_ambig, n_wambig, n_pcg))
	}
	dev.off()

	for(use_1_percentile in c(F,T)) {

		if(use_1_percentile) {
			pdf("figs/cgd_biv_dist_1perc_i.pdf", w=2.5,h=10)
			k27cov = cgd$eb4_k27_cov_1
			k4cov = cgd$eb4_k4_cov_1
			bivcov = cgd$eb4_biv_cov_1
			dst = cgd[,c("eb4_k27_cov_1", "eb4_biv_cov_1", "eb4_k4_cov_1")]
		} else {
			pdf("figs/cgd_biv_dist_0.5perc_i.pdf", w=2.5,h=10)
			k27cov = cgd$eb4_k27_cov
			k4cov = cgd$eb4_k4_cov
			bivcov = cgd$eb4_biv_cov
			dst = cgd[,c("eb4_k27_cov", "eb4_biv_cov", "eb4_k4_cov")]
		}
		layout(matrix(1:4, ncol=1))
		f_out = lk27 < 6 & lk4 < 6
		dst_n = dst/rowSums(dst)
		for(li in 1:4) {
			par(mar = c(2,2,1,1))
			f = l_class==li & !f_out

			dst_i = dst_n[f,]

			n = nrow(dst_i)
			ord = order((1e-3+dst_i[,1]+dst_i[,2])/(1e-3+dst_i[,3]+dst_i[,2]))
			dst_i = dst_i[ord,]
			plot(-0,-0, xlim=c(0,1), ylim=c(1,nrow(dst_i)), xaxt='n',yaxt='n')
			polygon(x = c(dst_i[,3], rep(0,t=n)), y=c(1:n, n:1), col="darkred", border=NA)
			polygon(x = c(dst_i[,3], rev(dst_i[,2]+ dst_i[,3])), y=c(1:n, n:1), col="black", border=NA)
			polygon(x = c(rep(1,t=n), rev(dst_i[,3]+dst_i[,2])), y=c(1:n, n:1), col="darkblue", NA)
		}
		dev.off()

		if(use_1_percentile) {
			pdf("figs/cgd_biv_dens_1perc_i.pdf", w=2.5,h=10)
		} else {
			pdf("figs/cgd_biv_dens_0.5perc_i.pdf", w=2.5,h=10)
		}
		layout(t(matrix(1:8, nrow=2)))
		l_class = as.numeric(cut(cgd$l, c(-1,500,750,1000,1e+6)))
		for(li in 1:4) {
			par(mar = c(2,2,1,1))
			f = l_class==li & !f_out

			dst_i = dst[f,]
			hist(dst_i[,3]+dst_i[,2], freq=F, col="darkred", ylim=c(0,5), main=NA,breaks=seq(-0.1,1,0.1))
			hist(dst_i[,1]+dst_i[,2], freq=F, col="darkblue", ylim=c(0,6), main=NA,breaks=seq(-0.1,1,0.1))
		}
		dev.off()
	}
}

plt_domain_spatial_d_atac_tss = function(mod, doms_interv = NULL, prefix="pcg",prime5=T,colors = cols,plt_tss=F){
gen_track_vext()
	gen_k4_vt(mod)

	ylims = list()
	ylims[["gc"]] = c(0.3,0.8)
	ylims[["cg"]] = c(0,0.16)
	ylims[["k4"]] = c(-1,14)
	ylims[["k27"]] = c(-1,8)
	ylims[["atac"]] = c(0,12)
	gvtrack.create("tss_d", mod$tss, "distance")
	l_rng = c(500,800,1000,2000,4000,8000)
	l_rng = c(1000,3000)
	shades=colorRampPalette(c("lightblue", "darkblue", "black", "darkblue", "lightblue"))(33)
#doms_interv = doms_interv[ doms_interv$strand== -1,]
	doms_interv = doms_interv %>% arrange(chrom,start)
	a = mod$npeaks

	a$start = a$start + (a$end-a$start)/2
	a$end = a$start +1
	a$strand = a$gene_strand
	gvtrack.create("atac_d", a, "distance")
	
	for(i in 1:(length(l_rng)-1)) {
		l = l_rng[i+1]
		l_min = l_rng[i]
		if(sum(doms_interv$l > l_min  & doms_interv$l < l) > 10) {
		 doms_interv_f = doms_interv[doms_interv$l > l_min  & doms_interv$l < l,]
          doms_interv_fd = doms_interv_f  
            if(prime5==T){
         doms_interv_fd$start = ifelse(doms_interv_fd$strand==1,doms_interv_fd$start,doms_interv_fd$end) 
            
         doms_interv_fd$end = doms_interv_fd$start+1  
		  gvtrack.create("dom_intd", doms_interv_fd, "distance")
                }
            else {
         doms_interv_fd$start = ifelse(doms_interv_fd$strand==1,doms_interv_fd$end,doms_interv_fd$start) 
            
         doms_interv_fd$end = doms_interv_fd$start+1  
		 gvtrack.create("dom_intd", doms_interv_fd, "distance") 
            }
		 doms_interv_pad = doms_interv_f
		 #doms_interv_pad$start = ifelse(doms_interv_f$strand==1,doms_interv_f$start-200,doms_interv_f$start-2000)
		 #doms_interv_pad$end = ifelse(doms_interv_f$strand==1,doms_interv_f$end+2000,doms_interv_f$end+200)
         doms_interv_pad$start = doms_interv_pad$start#-2000
		 doms_interv_pad$end = doms_interv_pad$end#+2000
        #doms_interv_pad$strand =NULL 
           #doms_interv_f$strand = NULL 
       
		 prof = gextract(c("atac_ext","EB4_cnt", "ES_ch_k4", 
								"jk.epipcg.multieb.marginal_regcap",
								"seq.CG_500_mean_new", "seq.GC500_bin20", 
								"dom_intd", "tss_d",'atac_d'),
			 intervals=doms_interv_pad, iterator=20, 
			colnames=c("atac","k27","k4","reg_atac","cg", "gc", "dst", "tss_d",'atac_d'))
		 prof[is.na(prof)] = 0
		 prof$at_tss = ifelse(abs(prof$tss_d) < 10,1,0)
		 prof$at_atac = ifelse(abs(prof$atac_d) < 10,1,0)
         doms_interv_pad$intervalID = 1:nrow(doms_interv_pad)
         prof = merge(prof,doms_interv_pad[,c('intervalID','type')],by='intervalID')
		 #pdf(sprintf("figs/dom_spat/int_dom_tss_%s_%s.pdf", prefix, l),w=8,h=5)
          if(plt_tss==TRUE)  {
		 tr = tapply(prof[ prof$type=='cl4',]$at_tss, floor(prof[ prof$type=='cl4',]$dst/10), mean, na.rm=T)
		 tr = rollmean(tr, 14, f='e')
		 tr3 = tapply(prof[ prof$type=='cl3',]$at_tss, floor(prof[ prof$type=='cl3',]$dst/10), mean, na.rm=T)
		 tr3 = rollmean(tr3, 14, f='e')        
		 tr2 = tapply(prof[ prof$type=='cl2',]$at_tss, floor(prof[ prof$type=='cl2',]$dst/10), mean, na.rm=T)
		 tr2 = rollmean(tr2, 14, f='e')     
		 tr1 = tapply(prof[ prof$type=='cl1',]$at_tss, floor(prof[ prof$type=='cl1',]$dst/10), mean, na.rm=T)
		 tr1 = rollmean(tr1, 14, f='e') 
		 plot(as.numeric(names(tr))*10, tr, type="l", lwd=2, pch=19,cex=0.5,col=colors[1],xlim=c(-1200,0))
         lines(as.numeric(names(tr3))*10, tr3, type="l", lwd=2, pch=19,cex=0.5,col=colors[2])
         lines(as.numeric(names(tr2))*10, tr2, type="l", lwd=2, pch=19,cex=0.5,col=colors[3])
         lines(as.numeric(names(tr1))*10, tr1, type="l", lwd=2, pch=19,cex=0.5,col=colors[4])
        }else{
            
	 	 #dev.off()
		 #pdf(sprintf("figs/dom_spat/int_dom_dist_atac_%s_%s.pdf", prefix, l),w=8,h=5)
		 tr = tapply(prof[ prof$type=='cl4',]$at_atac, floor(prof[ prof$type=='cl4',]$dst/10), mean, na.rm=T)
		 tr = rollmean(tr, 20, f='e')
		 tr3 = tapply(prof[ prof$type=='cl3',]$at_atac, floor(prof[ prof$type=='cl3',]$dst/10), mean, na.rm=T)
		 tr3 = rollmean(tr3, 20, f='e')        
		 tr2 = tapply(prof[ prof$type=='cl2',]$at_atac, floor(prof[ prof$type=='cl2',]$dst/10), mean, na.rm=T)
		 tr2 = rollmean(tr2, 20, f='e')     
		 tr1 = tapply(prof[ prof$type=='cl1',]$at_atac, floor(prof[ prof$type=='cl1',]$dst/10), mean, na.rm=T)
		 tr1 = rollmean(tr1, 20, f='e') 
		 plot(as.numeric(names(tr))*10, tr, type="l", lwd=2, pch=19,cex=0.5,col=colors[1],ylim = c(0,0.05),xlim=c(-1200,0))
         lines(as.numeric(names(tr3))*10, tr3, type="l", lwd=2, pch=19,cex=0.5,col=colors[2])
         lines(as.numeric(names(tr2))*10, tr2, type="l", lwd=2, pch=19,cex=0.5,col=colors[3])
         lines(as.numeric(names(tr1))*10, tr1, type="l", lwd=2, pch=19,cex=0.5,col=colors[4])
	 	 #dev.off()
		 }
		 
		}}
    }
plt_eb4_k27_domainogram = function(mod)
{
	all_c = gintervals.all()
	all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]

	d = mod$eb4_k27_doms
	
	d = d[d$l > 300 & !is.na(d$line) & d$line < 0.5 & d$ltr < 0.5 & d$sine <0.5 & d$low_complex<0.5,]

	tads = gintervals.load("intervs.global.tad_names")
	tads$l = tads$end-tads$start
	tads = tads[tads$l > 1e+5,]

	shades = colorRampPalette(c("white","white","white","lightgray","lightblue", "blue","darkblue","gold"))(1000)

	for(ci in 1:nrow(all_c)) {
		prof = gextract(c("jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a_a70ls"), iterator=500, intervals = all_c[ci,], colnames=c("k27"))
		prof[is.na(prof)] = 0
		dgram = matrix(ifelse(prof$k27>18,1,0), ncol=1)
		for(si in seq(1,10,0.5)) {
			sc = floor(2**si)
			dgram = cbind(dgram, rollmean(prof$k27, sc, f='e'))
		}

		xmax = all_c$end[ci]
		chrom = as.character(all_c$chrom[ci])
		#pdf(sprintf("figs/dgram%s.pdf", all_c$chrom[ci]), h=5,w=floor(20*xmax/1e+8))
		png(sprintf("figs/dgram%s.png", all_c$chrom[ci]), h=500,w=floor(3000*xmax/1e+8))
		image(log2(1e-5+dgram), col=shades, xaxt='n',yaxt='n')
		c_tads = tads[tads$chrom == all_c$chrom[ci],]
		segments(x0=c_tads$start/xmax, x1=c_tads$start/xmax, 
					y0=rep(-0.5,nrow(c_tads)), y1 = rep(1.5, nrow(c_tads)), 
					lwd=1, col="darkgreen", lty=1)
		f_lcg = d$chrom==chrom & d$CG_max < 0.02
		f_icg = d$chrom==chrom & d$CG_max > 0.02 & d$CG_max < 0.04
		f_hcg = d$chrom==chrom & d$CG_max > 0.04
		points(d$start[f_hcg]/xmax, 0.1+rnorm(sum(f_hcg),0,0.03), col="red", pch=19, cex=1)
		points(d$start[f_icg]/xmax, 0.3+rnorm(sum(f_icg),0,0.03), col="orange", pch=19, cex=1)
		points(d$start[f_lcg]/xmax, 0.5+rnorm(sum(f_lcg), 0,0.03), col="black", pch=19, cex=1)
		dev.off()
	}
		
}

plt_eb4_dom_tad = function(mod)
{
	tads = gintervals.load("intervs.global.tad_names")
	tads$l = tads$end-tads$start

	d = mod$eb4_k27_doms
	d = d[d$l > 300 & d$line < 0.5 & d$ltr < 0.5 & d$sine <0.5 & d$low_complex<0.5,]
	cgmode = as.numeric(cut(d$CG_max, c(0,0.02,0.04,1)))
	n_tad1 = tabulate(d$tad_index[cgmode==1], nrow(tads))
	n_tad2 = tabulate(d$tad_index[cgmode==2], nrow(tads))
	n_tad3 = tabulate(d$tad_index[cgmode==3], nrow(tads))
	sz_tad1 = tapply(d$l[cgmode==1], d$tad_index[cgmode==1], sum)
	sz_tad2 = tapply(d$l[cgmode==2], d$tad_index[cgmode==2], sum)
	sz_tad3 = tapply(d$l[cgmode==3], d$tad_index[cgmode==3], sum)
	sz_tad1 = sz_tad1[as.character(1:nrow(tads))]
	sz_tad2 = sz_tad2[as.character(1:nrow(tads))]
	sz_tad3 = sz_tad3[as.character(1:nrow(tads))]
	names(sz_tad1) = 1:nrow(tads)
	sz_tad1[is.na(sz_tad1)] = 0
	names(sz_tad2) = 1:nrow(tads)
	sz_tad2[is.na(sz_tad2)] = 0
	names(sz_tad3) = 1:nrow(tads)
	sz_tad3[is.na(sz_tad3)] = 0

#1. for ranges of TAD size show scatter of cg ncg cov
#2. starta of TAD HCG cov - how much LCG, reverse: starta of LCG cov
#3. domainogram of a kind? 1,2,4,8,16,32,64,128k distrib for LCG/HCG

	f = tads$l < 10e+5 & tads$l > 5e+5
	smoothScatter(log2(2e-3+(sz_tad3[f])/tads$l[f]), 
							log2(2e-3+sz_tad1[f]/tads$l[f]),
			colramp=colorRampPalette(c("white","gray","lightblue","blue","darkblue","red","yellow","black")),
			xlim=c(-9,-2),ylim=c(-9,-2),nrpoints=500, pch=19, cex=0.5)

	pcgmod = mod$eb4_k27_doms$cgtss_mode
	tad_mods = table(mod$eb4_k27_doms1k$tad_index, pcgmod)
	pheatmap::pheatmap(log2(1+tad_mods),cluster_cols=F, 
					show_rownames=F,	
					treeheight_row = 0,
					col=colorRampPalette(c("white","darkblue","black"))(1000),
					filename = "figs/eb_pcg_dom_tad.pdf", w = 4, h =10)

	pdf("figs/eb_pcg_tad_sqrt_d.pdf", w=6,h=6)	
	d5 = pmax(mod$eb4_k27_doms1k$tad_5_dist,1)+mod$eb4_k27_doms1k$l/2
	d3 = pmax(mod$eb4_k27_doms1k$tad_3_dist,1)+mod$eb4_k27_doms1k$l/2

	plot(sqrt(d5), sqrt(d3), pch=19,cex=0.4, xlim=c(0,1e+3), ylim=c(0,1e+3), col=ifelse(mod$eb4_k27_doms$cgtss_mode==1, "red", "black"), xlab="sqrt -dist border", ylab="sqrt +dist border")
	dev.off()

	pdf("figs/eb_pcg_tad_mode4_sqrt_d.pdf", w=6,h=6)	
	plot(sqrt(d5), sqrt(d3), pch=19,cex=0.4, xlim=c(0,1e+3), ylim=c(0,1e+3), col=ifelse(mod$eb4_k27_doms$cgtss_mode==4, "orange", "black"), xlab="sqrt -dist border", ylab="sqrt +dist border")
	dev.off()

	pdf("figs/eb_pcg_tad_d.pdf", w=6,h=6)	
	f = (d3+d5)<1500000
	d5 = d5[f]
	d3 = d3[f]
	plot((-(d5+d3)/2+d3)[order(d3+d5)], 1:length(d3),
				pch=19,cex=0.4, 
				xlab="locationi in TAD", ylab="sqrt +dist border")
	segments(x0=(-(d5+d3)/2)[order(d3+d5)], x1=((d5+d3)/2)[order(d3+d5)], 
					y0 = 1:length(d3), y1=1:length(d3), col="lightgray")
	points((-(d5+d3)/2+d3)[order(d3+d5)], 1:length(d3),
				pch=19,cex=0.4, 
				col=ifelse(mod$eb4_k27_doms$cgtss_mode[f][order(d3+d5)]==1, "red", "black"))
	dev.off()
}
plt_cgd_in_context = function(mod)
{
	d = mod$cgdom_ann
	d$strand = d$tss_strand
	f_cl1tss = d$type%in%'cl1' & d$tss_dist==0
	f_cl2tss = d$type%in%'cl2' & d$tss_dist==0
	f_cl3tss = d$type%in%'cl3' & d$tss_dist==0
	f_cl4tss = d$type%in%'cl4' & d$tss_dist==0
	f_cl1n = d$type%in%'cl1' & d$tss_dist!=0
	f_cl2n = d$type%in%'cl2' & d$tss_dist!=0
	f_cl3n = d$type%in%'cl3' & d$tss_dist!=0
	f_cl4n = d$type%in%'cl4' & d$tss_dist!=0
	
	
	#f_biv_seed = (mod$cgdom_ann$type_50k=='seed' & mod$cgdom_ann$eb4_biv_cov_1>0)
	
	#plt_domain_spatial(mod, d[f_biv_seed,], prefix="biv_seed")
	plt_domain_spatial(mod, d[f_cl1tss,], prefix="cl1tss")
	plt_domain_spatial(mod, d[f_cl2tss,], prefix="cl2tss")
	plt_domain_spatial(mod, d[f_cl3tss,], prefix="cl3tss")
	plt_domain_spatial(mod, d[f_cl4tss,], prefix="cl4tss")
	plt_domain_spatial(mod, d[f_cl1n,], prefix="cl1n")
	plt_domain_spatial(mod, d[f_cl2n,], prefix="cl2n")
	plt_domain_spatial(mod, d[f_cl3n,], prefix="cl3n")
	plt_domain_spatial(mod, d[f_cl4n,], prefix="cl4n")
}

plt_domain_in_context = function(mod)
{
	d = mod$eb4_k27_doms
    
	f_hcg_seed = d$CG_max > 0.05 & d$scale50k < 0.05 & d$l>300
	f_lcg_seed = d$CG_max < 0.02 & d$scale50k < 0.05 & d$l>300

	f_hcg_clust = d$CG_max > 0.05 & d$scale50k > 0.05 & d$l>300
	f_lcg_clust = d$CG_max < 0.02 & d$scale50k > 0.05 & d$l>300

	
	plt_domain_spatial(mod, d[f_hcg_seed,], prefix="hcg_seed")
	plt_domain_spatial(mod, d[f_lcg_seed,], prefix="lcg_seed")
	plt_domain_spatial(mod, d[f_hcg_clust,], prefix="hcg_clust")
	plt_domain_spatial(mod, d[f_lcg_clust,], prefix="lcg_clust")

}




plt_domain_spatial = function(mod, doms_interv = NULL, prefix="pcg", l_pad=1000)
{
   gen_track_vext()
	gen_k4_vt(mod)

	ylims = list()
	ylims[["gc"]] = c(0.3,0.8)
	ylims[["cg"]] = c(0,0.16)
	ylims[["k4"]] = c(-1,14)
	ylims[["k27"]] = c(-1,8)
	ylims[["atac"]] = c(0,12)
	gvtrack.create("tss_d", mod$tss, "distance")
	l_rng = c(500,800,1000,2000,4000,8000)
	#l_rng = c(1000,2000)
	shades=colorRampPalette(c("lightblue", "darkblue", "black", "darkblue", "lightblue"))(33)
	doms_interv = doms_interv %>% arrange(chrom,start)
	a = mod$npeaks

	a$start = a$start + (a$end-a$start)/2
	a$end = a$start +1
	a$strand = a$gene_strand
	gvtrack.create("atac_d", a, "distance")
	
	for(i in 1:(length(l_rng)-1)) {
		l = l_rng[i+1]
		l_min = l_rng[i]
		if(sum(doms_interv$l > l_min  & doms_interv$l < l) > 10) {
		 doms_interv_f = doms_interv[doms_interv$l > l_min  & doms_interv$l < l,]
		 gvtrack.create("dom_intd", doms_interv_f, "distance", l)
		 doms_interv_pad = doms_interv_f
		 doms_interv_pad$start = doms_interv_f$start - l_pad
		 doms_interv_pad$end = doms_interv_f$end + l_pad

		 prof = gextract(c("atac_ext","EB4_cnt", "ES_ch_k4", 
								"jk.epipcg.multieb.marginal_regcap",
								"seq.CG_500_mean_new", "seq.GC500_bin20", 
								"dom_intd", "tss_d",'atac_d'),
			 intervals=doms_interv_pad, iterator=20, 
			colnames=c("atac","k27","k4","reg_atac","cg", "gc", "dst", "tss_d",'atac_d'))
		 prof[is.na(prof)] = 0
		 prof$at_tss = ifelse(abs(prof$tss_d) < 10,1,0)
		 prof$at_atac = ifelse(abs(prof$atac_d) < 10,1,0)
		 for(nm in c("atac","k27","k4","cg", "gc")) {
			trend = tapply(prof[,nm], floor(prof$dst/20), quantile, 
														seq(0.1,0.9,0.025))
			mtrend = do.call('rbind', trend)
			if(nm != "gc" & nm != "cg") {
				mtrend = log2(0.5+mtrend)
			}
			mtrend = apply(mtrend,2,rollmean, 5, f='e')
			mtrend = pmax(mtrend, ylims[[nm]][1])
			pdf(sprintf("figs/dom_spat/%s_int_dom_%s_%s.pdf", prefix, nm, l),w=8,h=5)
			maxx = nrow(mtrend)
			plot(-100,-100,col=NA, xlim=c(1,maxx), ylim=ylims[[nm]], 
						main = sprintf("%s l %s %s n=%s", 
											prefix, l, nm, nrow(doms_interv_f)))
			for(ri in 1:32) {
				polygon(x=c(1:maxx,maxx:1), y=c(mtrend[,ri], rev(mtrend[,ri+1])),
							border=NA, col=shades[ri])
			}
#			barplot(t(cbind(mtrend[,1],mtrend[,-1]-mtrend[,-33])),
#						ylim=ylims[[nm]], space=0, border=NA, 
#						col=shades)
			dev.off()
		 }
		 pdf(sprintf("figs/dom_spat/int_dom_tss_%s_%s.pdf", prefix, l),w=8,h=5)
		 tr = tapply(prof$at_tss, floor(prof$dst/10), mean, na.rm=T)
		 tr = rollmean(tr, 14, f='e')
		 plot(as.numeric(names(tr))*10, tr, type="l", lwd=2, pch=19,cex=0.5, ylim=c(0,0.05))
	 	 dev.off()
		 pdf(sprintf("figs/dom_spat/int_dom_dist_atac_%s_%s.pdf", prefix, l),w=8,h=5)
		 tr = tapply(prof$at_atac, floor(prof$dst/10), mean, na.rm=T)
		 tr = rollmean(tr, 14, f='e')
		 plot(as.numeric(names(tr))*10, tr, type="l", lwd=2, pch=19,cex=0.5, ylim=c(0,0.05))
	 	 dev.off()
		 
		 
		}
		
	}

}

plot_k27_k4_starta = function(mod)
{
	d = mod$eb4_biv_doms
	d = d[d$l > 300 & !is.na(d$line) & d$line < 0.5 & d$ltr < 0.5 & d$sine <0.5 & d$low_complex<0.5,]
	cgbin = as.numeric(cut(d$CG_max, c(-1,0.02,0.04,0.06,0.10,1)))
	gcbin = as.numeric(cut(d$GC, c(0,0.48,0.53,0.58,0.63,1)))
	gcmaxbin = as.numeric(cut(d$GC_max, c(0,0.5,0.6,0.65,0.7,1)))
	lbin = as.numeric(cut(d$l, c(0,750,1000,2000,4000,1e+6)))

	q_atac = gquantiles("atac_ext", c(0.975,0.99,0.995,0.999,1))

	message("ATAC bins ", paste(q_atac, collapse=","))
	atacbin = as.numeric(cut(d$eb_atac, c(0,q_atac)))
	q_scale50 = quantile(d$scale50k, c(0,0.2,0.4,0.6,0.8,1))

	message("scale50 bins ", paste(q_scale50, collapse=","))
	clsbin = as.numeric(cut(d$scale50k, q_scale50))

	dk27 = log2(10+d$eb4_k27)	
	dk4 = log2(10+d$eb4_k4)	
	png("figs/BIV_cg_strata.png", w = 1500,h=300)
	layout(matrix(1:5,nrow=1))
	for(i in 1:5) {
		plot(dk4, dk27, pch=19, cex=0.3, col="gray")
		f = cgbin==i
		points(dk4[f], dk27[f], pch=19, cex=0.3, col="darkblue")
	}
	dev.off()
	png("figs/BIV_gcmax_strata.png", w = 1500,h=300)
	layout(matrix(1:5,nrow=1))
	for(i in 1:5) {
		plot(dk4, dk27, pch=19, cex=0.3, col="gray")
		f = gcmaxbin==i
		points(dk4[f], dk27[f], pch=19, cex=0.3, col="darkblue")
	}
	dev.off()
	png("figs/BIV_gc_strata.png", w = 1500,h=300)
	layout(matrix(1:5,nrow=1))
	for(i in 1:5) {
		plot(dk4, dk27, pch=19, cex=0.3, col="gray")
		f = gcbin==i
		points(dk4[f], dk27[f], pch=19, cex=0.3, col="darkblue")
	}
	dev.off()
	png("figs/BIV_l_strata.png", w = 1500,h=300)
	layout(matrix(1:5,nrow=1))
	for(i in 1:5) {
		plot(dk4, dk27, pch=19, cex=0.3, col="gray", main=sprintf("s50 %s", q_scale50[i+1]))
		f = lbin==i
		points(dk4[f], dk27[f], pch=19, cex=0.3, col="darkblue")
	}
	dev.off()
	png("figs/BIV_atac_strata.png", w = 1500,h=300)
	layout(matrix(1:5,nrow=1))
	for(i in 1:5) {
		plot(dk4, dk27, pch=19, cex=0.3, col="gray" , main=sprintf("ATAC %s", q_atac[i+1]))
		f = atacbin==i
		points(dk4[f], dk27[f], pch=19, cex=0.3, col="darkblue")
	}
	dev.off()
	png("figs/BIV_scale50_strata.png", w = 1500,h=300)
	layout(matrix(1:5,nrow=1))
	for(i in 1:5) {
		plot(dk4, dk27, pch=19, cex=0.3, col="gray")
		f = clsbin==i
		points(dk4[f], dk27[f], pch=19, cex=0.3, col="darkblue")
	}


}

pcg_test_polarize_nucs = function(mod)
{
   gen_track_vext()
   gen_k4_vt(mod)

	prof = gextract(c("atac_ext"), iterator=20,  intervals=mod$eb4_biv_doms)
	atac_max = tapply(prof$atac_ext, prof$intervalID, max)
	prof$atac_max = atac_max[ prof$intervalID]
	prof_max = prof[prof$atac_ext == prof$atac_max,]; 
	prof_max = prof_max[!duplicated(prof_max$intervalID),]
	prof_max = prof_max[!is.na(prof_max$chrom),]
	prof_max$strand = mod$eb4_biv_doms$strand[prof_max$intervalID]
	max_intervs = gintervals(prof_max$chrom, prof_max$start, prof_max$end, prof_max$strand)
	gvtrack.create("atacmax_d", max_intervs, "distance")
	mnase_es = gtrack.ls("Sou.+mnase_mESC",ignore.case=T)

	d = mod$eb4_biv_doms
	dom_pad = d
	dom_pad$start = dom_pad$start - 1000
	dom_pad$end = dom_pad$end + 1000
	prof_rel = gextract(c("atac_ext", "EB4_cnt", "EB4_cnt_k4", mnase_es, "atacmax_d"), 
								iterator=20,  intervals=dom_pad)
	prof_rel[is.na(prof_rel)] = 0
	f_k27_hcg_seed = d$eb4_k27 > 18 & d$eb4_k4 < 64 & d$CG_max > 0.04 & d$l>500 & d$scale50k < 0.06
	f_k27b_hcg_seed = d$eb4_k27 > 18 & d$eb4_k4 > 64 & d$CG_max > 0.04 & d$l>500 & d$scale50k < 0.06
	f_k27_hcg_cls = d$eb4_k27 > 18 & d$CG_max > 0.04 & d$l>500 & d$scale50k > 0.06
	f_k4_hcg_seed = d$eb4_k27 < 6 & d$CG_max > 0.04 & d$l>500 & d$scale50k < 0.06
	f_k4_hcg_cls = d$eb4_k27 < 6 & d$CG_max > 0.04 & d$l>500 & d$scale50k > 0.06
	ff_k27_hcg_seed = f_k27_hcg_seed[prof_rel$intervalID]
	ff_k27b_hcg_seed = f_k27b_hcg_seed[prof_rel$intervalID]
	ff_k4_hcg_seed = f_k4_hcg_seed[prof_rel$intervalID]
	ff_k27_hcg_cls = f_k27_hcg_cls[prof_rel$intervalID]
	ff_k4_hcg_cls = f_k4_hcg_cls[prof_rel$intervalID]

	round_polar = function(x) { ifelse(x <= 0, floor(x), ceiling(x)) }
	trk27_seed_mn = matrix(tapply(prof_rel[ff_k27_hcg_seed,7], 
					round_polar(prof_rel$atacmax_d[ff_k27_hcg_seed]/20), length),ncol=1)
	trk27b_seed_mn = matrix(tapply(prof_rel[ff_k27b_hcg_seed,7], 
					round_polar(prof_rel$atacmax_d[ff_k27b_hcg_seed]/20), length),ncol=1)
	trk27_cls_mn = matrix(tapply(prof_rel[ff_k27_hcg_cls,7], 
						round_polar(prof_rel$atacmax_d[ff_k27_hcg_cls]/20), length),ncol=1)
	trk4_seed_mn = matrix(tapply(prof_rel[ff_k4_hcg_seed,7], 
						round_polar(prof_rel$atacmax_d[ff_k4_hcg_seed]/20), length),ncol=1)
	trk4_cls_mn = matrix(tapply(prof_rel[ff_k4_hcg_cls,7], 
						round_polar(prof_rel$atacmax_d[ff_k4_hcg_cls]/20), length),ncol=1)

	for(i in 7:10) {
		trk27_seed_mn= cbind(trk27_seed_mn, tapply(prof_rel[ff_k27_hcg_seed,i], 
						round_polar(prof_rel$atacmax_d[ff_k27_hcg_seed]/20), mean))
		trk27b_seed_mn= cbind(trk27b_seed_mn, tapply(prof_rel[ff_k27b_hcg_seed,i], 
						round_polar(prof_rel$atacmax_d[ff_k27b_hcg_seed]/20), mean))
		trk27_cls_mn = cbind(trk27_cls_mn, tapply(prof_rel[ff_k27_hcg_cls,i], 
						round_polar(prof_rel$atacmax_d[ff_k27_hcg_cls]/20), mean))
		trk4_seed_mn= cbind(trk4_seed_mn, tapply(prof_rel[ff_k4_hcg_seed,i], 
						round_polar(prof_rel$atacmax_d[ff_k4_hcg_seed]/20), mean))
		trk4_cls_mn = cbind(trk4_cls_mn, tapply(prof_rel[ff_k4_hcg_cls,i], 
						round_polar(prof_rel$atacmax_d[ff_k4_hcg_cls]/20), mean))
	}

	plt_nuc_tr = function(trace, pref, norm_type='none') {
		if(norm_type == 'region') {
			fcol = abs(as.numeric(rownames(trace))) < 200
			yl = c(0.5,2)
			trace[,2:5] = t(t(trace[,2:5])/colMeans(trace[fcol,2:5]))
		} else if(norm_type == 'nuc') {
			for(i in 2:5) {
				yl = c(0.7,1.5)
				trace[,i] = trace[,i] / rollmean(trace[,i],16, f='e')
			}
		}
		n = max(trace[,1])
		pdf(sprintf("figs/nuc_trace_%s.pdf", pref), h=5,w=6)
		plot(as.numeric(rownames(trace))*20, trace[,3], 
						xlim=c(-1000,1000), ylim=yl, type="l", col="gray",
						xlab = "atac offset", ylab="MNase",
						main = sprintf("%s n=%s", pref, n))
		
		lines(as.numeric(rownames(trace))*20, trace[,4], col="black")
		lines(as.numeric(rownames(trace))*20, trace[,2], col="pink")
		lines(as.numeric(rownames(trace))*20, trace[,5], col="red")
		for(x in seq(-1000,1000,200)) {
			abline(v=x, lty=2)
		}
		dev.off()
	}
	plt_nuc_tr(trk27_seed_mn, "k27_seed", "region")
	plt_nuc_tr(trk27b_seed_mn, "k27b_seed", "region")
	plt_nuc_tr(trk27_cls_mn, "k27_cls", "region")
	plt_nuc_tr(trk4_seed_mn, "k4_seed", "region")
	plt_nuc_tr(trk4_cls_mn, "k4_cls", "region")
	plt_nuc_tr(trk27_seed_mn, "norm_k27_seed", "nuc")
	plt_nuc_tr(trk27b_seed_mn, "norm_k27b_seed", "nuc")
	plt_nuc_tr(trk27_cls_mn, "norm_k27_cls", "nuc")
	plt_nuc_tr(trk4_seed_mn, "norm_k4_seed", "nuc")
	plt_nuc_tr(trk4_cls_mn, "norm_k4_cls", "nuc")

}

test_polarize_nucs_mnase_cl_officer = function(mod)
{	
	library(officer)
	library(rvg)
   gen_track_vext()
   gen_k4_vt(mod)
    cgd = mod$cgdom_ann
	cgd$strand = mod$cgdom_ann$tss_strand
	cgd = cgd %>% filter(type!='cl5')
	prof = gextract(c("atac_ext"), iterator=20,  intervals=cgd)
	atac_max = tapply(prof$atac_ext, prof$intervalID, max)
	prof$atac_max = atac_max[ prof$intervalID]
	prof_max = prof[prof$atac_ext == prof$atac_max,]; 
	prof_max = prof_max[!duplicated(prof_max$intervalID),]
	prof_max = prof_max[!is.na(prof_max$chrom),]
	prof_max$strand = cgd$strand[prof_max$intervalID]
	max_intervs = gintervals(prof_max$chrom, prof_max$start, prof_max$end, prof_max$strand)
	gvtrack.create("atacmax_d", max_intervs, "distance")
	mnase_es = gtrack.ls("Sou.+mnase_mESC",ignore.case=T)

	d = cgd
	dom_pad = d
	dom_pad$start = dom_pad$start - 1000
	dom_pad$end = dom_pad$end + 1000
	prof_rel = gextract(c("atac_ext", "EB4_cnt", "EB4_cnt_k4", mnase_es, "atacmax_d"), 
								iterator=20,  intervals=dom_pad,colnames = c("atac_ext", "EB4_cnt", "EB4_cnt_k4", mnase_es, "atacmax_d"))
	prof_rel[is.na(prof_rel)] = 0

	
	#f_cl1_seed = log2(1+d$eb4_k27) >= 6.5 & log2(1+d$eb4_k4)<5.5 & d$CG_max > 0.04 & d$scale50k < 0.06
	#f_cl1_clust = log2(1+d$eb4_k27) >= 6.5 & log2(1+d$eb4_k4)<5.5 & d$CG_max > 0.04 & d$scale50k > 0.06
	#f_cl2_seed = log2(1+d$eb4_k27) >= 6.5 & log2(1+d$eb4_k4)>5.5 & d$CG_max > 0.04 & d$scale50k < 0.06
	#f_cl2_clust = log2(1+d$eb4_k27) >= 6.5 & log2(1+d$eb4_k4)>5.5 & d$CG_max > 0.04 & d$scale50k > 0.06
	#f_cl3_seed = log2(1+d$eb4_k27) < 6.5 & log2(1+d$eb4_k27) > 5.6 & log2(1+d$eb4_k4)>5.5& d$CG_max > 0.04 & d$scale50k < 0.06
	#f_cl3_clust = log2(1+d$eb4_k27) < 6.5 & log2(1+d$eb4_k27) > 5.6 & log2(1+d$eb4_k4)>5.5& d$CG_max > 0.04 & d$scale50k > 0.06
	#f_cl4_seed = log2(1+d$eb4_k27) > 3.8 & log2(1+d$eb4_k27) < 5.6 & log2(1+d$eb4_k4)>5.5& d$CG_max > 0.04 & d$scale50k < 0.06
	#f_cl4_clust = log2(1+d$eb4_k27) > 3.8 & log2(1+d$eb4_k27) < 5.6 & log2(1+d$eb4_k4)>5.5& d$CG_max > 0.04 & d$scale50k > 0.06
	#f_cl5_seed = log2(1+d$eb4_k27) < 3.8  & log2(1+d$eb4_k4)>5.5& d$CG_max > 0.04 & d$scale50k < 0.06
	#f_cl5_clust = log2(1+d$eb4_k27) < 3.8  & log2(1+d$eb4_k4)>5.5& d$CG_max > 0.04 & d$scale50k > 0.06
	
	f_cl1_seed = d$type=='cl1' & d$tss_dist == 0
	f_cl1_clust = d$type=='cl1' & d$tss_dist != 0
	f_cl2_seed = d$type=='cl2' & d$tss_dist == 0
	f_cl2_clust = d$type=='cl2' & d$tss_dist != 0
	f_cl3_seed = d$type=='cl3' & d$tss_dist == 0
	f_cl3_clust = d$type=='cl3' & d$tss_dist != 0
	f_cl4_seed = d$type=='cl4' & d$tss_dist == 0
	f_cl4_clust = d$type=='cl4' & d$tss_dist != 0

	
	ff_cl1_seed = f_cl1_seed[prof_rel$intervalID]
	ff_cl1_clust = f_cl1_clust[prof_rel$intervalID]
	ff_cl2_seed = f_cl2_seed[prof_rel$intervalID]
	ff_cl2_clust = f_cl2_clust[prof_rel$intervalID]	
	ff_cl3_seed = f_cl3_seed[prof_rel$intervalID]
	ff_cl3_clust = f_cl3_clust[prof_rel$intervalID]	
	ff_cl4_seed = f_cl4_seed[prof_rel$intervalID]
	ff_cl4_clust = f_cl4_clust[prof_rel$intervalID]	
	#ff_cl5_seed = f_cl5_seed[prof_rel$intervalID]
	#ff_cl5_clust = f_cl5_clust[prof_rel$intervalID]	


	round_polar = function(x) { ifelse(x <= 0, floor(x), ceiling(x)) }

						
	trcl1_seed = matrix(tapply(prof_rel[ff_cl1_seed,4], 
					round_polar(prof_rel$atacmax_d[ff_cl1_seed]/20), length),ncol=1)
	trcl1_clust = matrix(tapply(prof_rel[ff_cl1_clust,4], 
					round_polar(prof_rel$atacmax_d[ff_cl1_clust]/20), length),ncol=1)
					
	trcl2_seed = matrix(tapply(prof_rel[ff_cl2_seed,4], 
					round_polar(prof_rel$atacmax_d[ff_cl2_seed]/20), length),ncol=1)
	trcl2_clust = matrix(tapply(prof_rel[ff_cl2_clust,4], 
					round_polar(prof_rel$atacmax_d[ff_cl2_clust]/20), length),ncol=1)
					
	trcl3_seed = matrix(tapply(prof_rel[ff_cl3_seed,4], 
					round_polar(prof_rel$atacmax_d[ff_cl3_seed]/20), length),ncol=1)
	trcl3_clust = matrix(tapply(prof_rel[ff_cl3_clust,4], 
					round_polar(prof_rel$atacmax_d[ff_cl3_clust]/20), length),ncol=1)

	trcl4_seed = matrix(tapply(prof_rel[ff_cl4_seed,4], 
					round_polar(prof_rel$atacmax_d[ff_cl4_seed]/20), length),ncol=1)
	trcl4_clust = matrix(tapply(prof_rel[ff_cl4_clust,4], 
					round_polar(prof_rel$atacmax_d[ff_cl4_clust]/20), length),ncol=1)
					
	#trcl5_seed = matrix(tapply(prof_rel[ff_cl5_seed,4], 
	#				round_polar(prof_rel$atacmax_d[ff_cl5_seed]/20), length),ncol=1)
	#trcl5_clust = matrix(tapply(prof_rel[ff_cl5_clust,4], 
	#				round_polar(prof_rel$atacmax_d[ff_cl5_clust]/20), length),ncol=1)


	for(i in 7:10) {
		trcl1_seed = cbind(trcl1_seed , tapply(prof_rel[ff_cl1_seed,i], 
						round_polar(prof_rel$atacmax_d[ff_cl1_seed]/20), mean))
		trcl1_clust = cbind(trcl1_clust , tapply(prof_rel[ff_cl1_clust,i], 
						round_polar(prof_rel$atacmax_d[ff_cl1_clust]/20), mean))

		trcl2_seed = cbind(trcl2_seed , tapply(prof_rel[ff_cl2_seed,i], 
						round_polar(prof_rel$atacmax_d[ff_cl2_seed]/20), mean))
		trcl2_clust = cbind(trcl2_clust , tapply(prof_rel[ff_cl2_clust,i], 
						round_polar(prof_rel$atacmax_d[ff_cl2_clust]/20), mean))

		trcl3_seed = cbind(trcl3_seed , tapply(prof_rel[ff_cl3_seed,i], 
						round_polar(prof_rel$atacmax_d[ff_cl3_seed]/20), mean))
		trcl3_clust = cbind(trcl3_clust , tapply(prof_rel[ff_cl3_clust,i], 
						round_polar(prof_rel$atacmax_d[ff_cl3_clust]/20), mean))	

		trcl4_seed = cbind(trcl4_seed , tapply(prof_rel[ff_cl4_seed,i], 
						round_polar(prof_rel$atacmax_d[ff_cl4_seed]/20), mean))
		trcl4_clust = cbind(trcl4_clust , tapply(prof_rel[ff_cl4_clust,i], 
						round_polar(prof_rel$atacmax_d[ff_cl4_clust]/20), mean))

		#trcl5_seed = cbind(trcl5_seed , tapply(prof_rel[ff_cl5_seed,i], 
		#				round_polar(prof_rel$atacmax_d[ff_cl5_seed]/20), mean))
		#trcl5_clust = cbind(trcl5_clust , tapply(prof_rel[ff_cl5_clust,i], 
		#				round_polar(prof_rel$atacmax_d[ff_cl5_clust]/20), mean))
						
}




  back = colMeans(prof_rel[abs(as.numeric(prof_rel$atacmax_d))>700,c(mnase_es)])


	plt_nuc_tr = function(trace, pref, norm_type='none') {
		if(norm_type == 'region') {
			fcol = abs(as.numeric(rownames(trace))) < 700
			yl = c(0.5,2.61)
			trace[,2:5] = t(t(trace[,2:5])/back)
		} else if(norm_type == 'nuc') {
			for(i in 2:5) {
				yl = c(0.7,1.35)
				trace[,i] = trace[,i] / rollmean(trace[,i],16, f='e')
			}
		}
		n = max(trace[,1])
		
		#pdf(sprintf("figs/nuc_trace_%s.pdf", pref), h=5,w=6)
		
		plot(as.numeric(rownames(trace))*20, trace[,3], 
						xlim=c(-1000,1000), ylim=yl, type="l", col='#FDAE6B',lwd = 2,
						xlab = "atac offset", ylab="MNase",
						main = sprintf("%s n=%s", pref, n))
####	"gray" 	"black" "pink" "red"
		lines(as.numeric(rownames(trace))*20, trace[,4], col='#F16913',lwd = 2 )
		lines(as.numeric(rownames(trace))*20, trace[,2], col='#D94801',lwd = 2)
		lines(as.numeric(rownames(trace))*20, trace[,5], col='#8C2D04',lwd = 2)
		for(x in seq(-1000,1000,200)) {
			abline(v=x, lty=2)
		}
		#dev.off()
	}


	
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl1_seed, "cl1_tss", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl1_tss.pptx'))
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl1_clust, "cl1_ntss", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl1_ntss.pptx'))	
	
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl2_seed, "cl2_tss", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl2_tss.pptx'))
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl2_clust, "cl2_ntss", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl2_ntss.pptx'))
	
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl3_seed, "cl3_tss", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl3_tss.pptx'))
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl3_clust, "cl3_ntss", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl3_ntss.pptx'))

	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl4_seed, "cl4_tss", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl4_tss.pptx'))
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl4_clust, "cl4_ntss", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl4_ntss.pptx'))

	#save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl5_seed, "cl5_seed", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl5_seed.pptx'))
	#save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl5_clust, "cl5_clust", "region") , link_ppt = paste0("./figs/mnase_nucs",'_cl5_clust.pptx'))
	
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl1_seed, "cl1_tss", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl1_tss.pptx'))
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl1_clust, "cl1_ntss", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl1_ntss.pptx'))	
	
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl2_seed, "cl2_tss", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl2_tss.pptx'))
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl2_clust, "cl2_ntss", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl2_ntss.pptx'))
	
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl3_seed, "cl3_tss", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl3_tss.pptx'))
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl3_clust, "cl3_ntss", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl3_ntss.pptx'))

	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl4_seed, "cl4_tss", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl4_tss.pptx'))
	save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl4_clust, "cl4_ntss", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl4_ntss.pptx'))

	#save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl5_seed, "cl5_seed", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl5_seed.pptx'))
	#save_baseR_to_ppt(plot_func =plt_nuc_tr(trcl5_clust, "cl5_ntss", "nuc") , link_ppt = paste0("./figs/norm_mnase_nucs",'_cl5_ntss.pptx'))
	
	
}





plt_reg_g = function(mod, g, marg5=2e+3, marg3=2e+3, add_mark=F, wk4=F)
{
	f = mod$tss$geneSymbol==g
	if(sum(f) !=0) {
		hits = mod$tss[f,]
		locus = mean(mod$tss$start[f])
		chrom = hits[1,"chrom"]
		if(add_mark) {
			x_tss = hits
		} else {
			x_tss = NULL
		}
		plt_region(mod, chr=hits$chrom, st = locus, end = locus+1,  marg5=marg5, marg3=marg3,
						mark=x_tss, wk4=wk4, title=g)
	} else {
		message("no hit")
	}
}
plt_region = function(mod, chr, st, end, marg5=2e+3, marg3=2e+3, mark=c(), wk4=F, k4_max_fact=1, title=NA)
{
	ndx = mod$epi_tracks
	for(i in 1:nrow(ndx)) {
		gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "sum") 
		gvtrack.iterator(ndx$short_name[i], sshift = -140, eshift= 140)
	}
	prof = gextract(ndx$short_name, intervals=gintervals(chr, st-marg5, end+marg3), iterator=20)

	if(wk4) {
		k4_tns = gen_k4_vt(mod)
		prof_k4 = gextract(k4_tns, intervals=gintervals(chr, st-marg5, end+marg3), iterator=20)
	}
	prof_cg = gextract(c("seq.CG_500_mean_new"), intervals=gintervals(chr, st-marg5, end+marg3), iterator=20, colnames=c("cg"))
	prof_cg$cg = pmin(prof_cg$cg, 0.1)-0.02

	layout(matrix(1:(1+nrow(ndx)),ncol=1),h=c(1.4,1.4, rep(1,nrow(ndx)-2), 1.4))
	
	par(mar=c(0,4,2,4))
	plot(prof$start, rep(0,length(prof$start)), col=NA, ylim=c(2,8), xaxt='n')
	cgd = mod$cgdom_ann
	fcgd = as.character(cgd$chrom)==chr & (fcgd = cgd$end > st | cgd$start < end)
	if(sum(fcgd) != 0) {
		cgd = cgd[fcgd,]
		pred = mod$locmod[["k27_eb4"]]$pred[fcgd]
		pred = pmin(pmax(5-pred,3),7)
		obs = pmin(pmax(log2(10+cgd$eb4_k27_max),3),7)
		segments(x0=cgd$start, x1=cgd$end, y0=obs, y1=obs,
					col="blue", lwd=3)
		segments(x0=cgd$start, x1=cgd$end, y0=pred, y1=pred,
					col="cyan", lwd=3)
	}
	abline(h=5)
	abline(h=6)
	abline(h=4)
	for(i in 1:nrow(ndx)) {
		if(i == nrow(ndx)) {
			par(mar=c(2,4,0.8,4))
		} else if(i==1) {
			par(mar=c(0,4,0.8,4))
		} else {
			par(mar=c(0,4,0.8,4))
		}
		tmax = mod$k27_track_thresh[i,7]	
		plot(prof$start, pmin(prof[,3+i],tmax), pch=19, type="l", 
				lwd=3, ylim=c(0,tmax), 
				xaxt=ifelse(i==nrow(ndx), 's', 'n'), 
				main = sprintf("%s : %s", title, ndx$short_name[i]), ylab=NA)
		if(wk4 & !is.na(ndx$track_k4[i])) {
			tmax_k4 = mod$k4_track_thresh[ndx$short_name_k4[i],7]	* k4_max_fact
			#message("added ", k4_max_fact, " to k4 thresh", " ", tmax_k4, " ", tmax)
			k4_vs = prof_k4[,ndx$short_name_k4[i]]*tmax/tmax_k4
			lines(prof_k4$start, pmin(k4_vs,tmax), pch=19, type="l", lwd=2, col="red")
		}
		points(prof_cg$start, prof_cg$cg * tmax/0.08, pch=19, cex=0.4, col="green")
		abline(v=st)
		abline(v=end)
		if(!is.null(mark) & length(mark)>0) {
			for(x in mark) {
				abline(v=x, col="blue", lwd=2)
			}
		}
	}
}



pcg_gen_cg_ih_doms = function(mod, force_update=F)
{
	if(!force_update & file.exists("data/cgdom_ann_i.RDS")) {
		mod$cgdom_ann_i = readRDS("data/cgdom_ann_i.RDS")
		mod$cgd_rpt_i = readRDS("data/cgd_rpt_i.RDS")
		return(mod)
	}else{
	options(gmax.data.size=1e+9)
	gen_k4_vt(mod)
	gen_k27_vt(mod)
	gen_track_vext()

	cgd = gscreen("seq.CG_500_mean_new>0.022 & seq.CG_500_mean_new <= 0.04")
	cgd = cgd[!cgd$chrom %in% c("chrM", "chrY"),]
	cgd$l = cgd$end-cgd$start
	cgd = cgd[cgd$l>300,]

	rownames(cgd) = 1:nrow(cgd)
	
#annotate: max CG, max k27, k4, GC
	prof = gextract(c("EB4_cnt", "EB3_cnt", "EB4_cnt_k4", "atac_ext","seq.CG_500_mean_new", "seq.GC500_bin20"), iterator=20, intervals=cgd)
	prof[is.na(prof)]= 0
	cgd$eb3_k27_max = tapply(prof[,"EB3_cnt"], prof$intervalID, max)
	cgd$eb4_k27_max = tapply(prof[,"EB4_cnt"], prof$intervalID, max)
	cgd$eb4_k4_max = tapply(prof[,"EB4_cnt_k4"], prof$intervalID, max)

	T_k27_995 = mod$k27_track_thresh["EB4_cnt",7]
	T_k4_995 = mod$k4_track_thresh["EB4_cnt",7]
	T_k27_99 = mod$k27_track_thresh["EB4_cnt",6]
	T_k4_99 = mod$k4_track_thresh["EB4_cnt",6]
	up_k27 = ifelse(prof[,"EB4_cnt"]>T_k27_995,1,0)
	up_k4 = ifelse(prof[,"EB4_cnt_k4"]>T_k4_995,1,0)
	cgd$eb4_k27_cov = tapply(up_k27 * (1-up_k4), prof$intervalID, mean)
	cgd$eb4_k4_cov = tapply((1-up_k27) * up_k4, prof$intervalID, mean)
	cgd$eb4_biv_cov = tapply(up_k27 * up_k4, prof$intervalID, mean)

	up_k27_1 = ifelse(prof[,"EB4_cnt"]>T_k27_99,1,0)
	up_k4_1 = ifelse(prof[,"EB4_cnt_k4"]>T_k4_99,1,0)
	cgd$eb4_k27_cov_1 = tapply(up_k27_1 * (1-up_k4_1), prof$intervalID, mean)
	cgd$eb4_k4_cov_1 = tapply((1-up_k27_1) * up_k4_1, prof$intervalID, mean)
	cgd$eb4_biv_cov_1 = tapply(up_k27_1 * up_k4_1, prof$intervalID, mean)

	cgd$atac_max = tapply(prof[,"atac_ext"], prof$intervalID, max)
	cgd$cg_max = tapply(prof[,"seq.CG_500_mean_new"], prof$intervalID, max)
	cgd$gc_max = tapply(prof[,"seq.GC500_bin20"], prof$intervalID, max)

	fcov = (prof$EB4_cnt_k4+prof$EB4_cnt)>20
	prof[!fcov,"EB4_cnt_k4"] = NA
	prof[!fcov,"EB4_cnt"] = NA
	cgd$eb4_k4_rat = tapply(log2((10+prof[,"EB4_cnt_k4"])/(10+prof[,"EB4_cnt"])), prof$intervalID, max, na.rm=T)
	cgd$eb4_k27_rat = tapply(log2((10+prof[,"EB4_cnt"])/(10+prof[,"EB4_cnt_k4"])), prof$intervalID, max, na.rm=T)

	sines = gintervals.load("intervs.global.rmsk_sine")
	gvtrack.create("exon_d", "intervs.global.exon", "distance")
	gvtrack.create("line_d", "intervs.global.rmsk_line", "distance")
	gvtrack.create("ltr_d", "intervs.global.rmsk_ltr", "distance")
	gvtrack.create("simp_d", "intervs.global.rmsk_simple_repeat", "distance")
	gvtrack.create("sine_d",sines, "distance")
	gvtrack.create("lowcomplex_d", "intervs.global.rmsk_low_complexity", "distance")

	cgd_rpts_20 = gextract(c("ifelse(exon_d==0,1,0)","ifelse(line_d==0,1,0)", "ifelse(ltr_d==0,1,0)", "ifelse(sine_d==0,1,0)", "ifelse(lowcomplex_d==0,1,0)", "ifelse(simp_d==0,1,0)"), iterator=20, intervals=cgd, colnames=c("exon","line","ltr","sine","low_complex", "simple"))
	cgd_rpts = data.frame(
						exon = tapply(cgd_rpts_20$exon, cgd_rpts_20$intervalID, mean),
						line = tapply(cgd_rpts_20$line, cgd_rpts_20$intervalID, mean),
						sine = tapply(cgd_rpts_20$sine, cgd_rpts_20$intervalID, mean),
						ltr = tapply(cgd_rpts_20$ltr, cgd_rpts_20$intervalID, mean),
						simp = tapply(cgd_rpts_20$simp, cgd_rpts_20$intervalID, mean),
						low_complex = tapply(cgd_rpts_20$low_complex, cgd_rpts_20$intervalID, mean))
	cgd_rpts$tot_rpt = rowSums(cgd_rpts[,c("line","sine","ltr","simp","low_complex")])
	cgd_ann = cbind(cgd, cgd_rpts)
	f_bad = cgd_ann$line > 0.5 | cgd_ann$ltr > 0.5 
	f_bad1 =  !f_bad & cgd_ann$tot_rpt > 0.5 
	f_exon = cgd_ann$exon > 0.8
	f_lowmapability = !f_bad & !f_bad1 & (cgd_ann$eb4_k27_max+cgd_ann$eb4_k4_max) == 0
	f_mask = f_bad | f_bad1 | f_lowmapability
	cgd_ann_mask = cgd_ann[f_mask,]
	cgd_ann = cgd_ann[!f_mask,]

	cgd_tss = gintervals.neighbors(cgd_ann, mod$tss)
	cgd_ann$tss_dist = cgd_tss$dist
	cgd_ann$tss_strand = cgd_tss$strand
	cgd_ann$tss_start = cgd_tss$start
	cgd_ann$tss_gene = cgd_tss$geneSymbol

	gvtrack.create("d_cgd",cgd_ann[,1:3], "distance")
	gvtrack.create("d_cgdnotss",cgd_ann[abs(cgd_ann$tss_dist)>1000,1:3], "distance")

	all_c = gintervals.all()
	all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]
	cg_trace = gextract(c("ifelse(seq.CG_500_mean_new > 0.04, 1, 0)"),
								 intervals=all_c, iterator = 200, 
									colnames("cg_all"))

	cg_trace1 = gextract(c("ifelse(!is.na(d_cgd) & d_cgd==0, 1, 0)"),
								 intervals=all_c, iterator = 200, 
									colnames("cg_filt"))

	cg_trace2 = gextract(c("ifelse(!is.na(d_cgdnotss) & d_cgdnotss==0, 1, 0)"),
								 intervals=all_c, iterator = 200, 
									colnames("cg_notss"))

	cgd_regstat = data.frame( 
		chrom = cg_trace$chrom,
		start = cg_trace$start,
		end = cg_trace$end,
		cg_2k = rollmean(cg_trace[,4],10,f='e'),
		cg_4k = rollmean(cg_trace[,4],20,f='e'),
		cg_8k = rollmean(cg_trace[,4],40,f='e'),
		cg_16k = rollmean(cg_trace[,4],80,f='e'),
		cg_32k = rollmean(cg_trace[,4],160,f='e'),
		cgd_2k = rollmean(cg_trace1[,4],10,f='e'),
		cgd_4k = rollmean(cg_trace1[,4],20,f='e'),
		cgd_8k = rollmean(cg_trace1[,4],40,f='e'),
		cgd_16k = rollmean(cg_trace1[,4],80,f='e'),
		cgd_32k = rollmean(cg_trace1[,4],160,f='e'))
	cgd_regstat$cgnotss_2k = rollmean(cg_trace2[,4],10,f='e')
	cgd_regstat$cgnotss_4k = rollmean(cg_trace2[,4],20,f='e')
	cgd_regstat$cgnotss_8k = rollmean(cg_trace2[,4],40,f='e')
	cgd_regstat$cgnotss_16k = rollmean(cg_trace2[,4],80,f='e')
	cgd_regstat$cgnotss_32k = rollmean(cg_trace2[,4],160,f='e')
# region al CG elements
	cgd_cent = cgd_ann[,c("chrom","start","end")]
	cgd_cent$start = (cgd_cent$start+cgd_cent$end)/2
	cgd_cent$end = cgd_cent$start + 1
	cgd_reg = gintervals.neighbors(cgd_cent, cgd_regstat)
	cgd_ann  = cbind(cgd_ann, cgd_reg[,-c(1:6,ncol(cgd_reg))])
	####seed clust
	all_c = gintervals.all()
	all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]
	pall = gextract(c("EB4_cnt"), iterator=500, intervals = all_c, colnames=c("k27"))
	T_k27 = mod$k27_track_thresh["EB4_cnt","X0.99"]
	pall$dom = ifelse(pall$k27 > T_k27, 1, 0)
	pall$dom50 = rollmean(pall$dom,100,f='e')	
	all_to50 = gintervals.neighbors(cgd_ann, pall)	
	cgd_ann$scale50k = all_to50$dom50
	
	###divide into clusters
	d = cgd_ann
	d$type = ifelse(log2(1+d$eb3_k27_max) >= 6 & log2(1+d$eb4_k4_max)<6 ,'cl1','cl5')

	d$type = ifelse(log2(1+d$eb3_k27_max) >= 6 & log2(1+d$eb4_k4_max)>6 ,'cl2',d$type)

	d$type = ifelse(log2(1+d$eb3_k27_max) < 6 & log2(1+d$eb3_k27_max) > 4 & log2(1+d$eb4_k4_max)>5.5,'cl3',d$type)

	d$type = ifelse( log2(1+d$eb3_k27_max) < 4 & log2(1+d$eb4_k4_max)>5.5,'cl4',d$type)
	d$type_50k = ifelse(d$scale50k<.06,'seed','clust')
	cgd_ann = d	
	
	mod$cgdom_ann_i = cgd_ann
	mod$cgd_rpt_i = cgd_ann_mask
	saveRDS(mod$cgdom_ann_i, "data/cgdom_ann_i.RDS")
	saveRDS(mod$cgd_rpt_i, "data/cgd_rpt_i.RDS")

	return(mod)
	}
}

plot_pcg_enr_wt_vs_inh = function(tracks,cl = c('cl1','cl2')) {


intervs = mod$cgdom_ann[, c('chrom','start','end','tss_strand','ID','type')]
colnames(intervs) = c('chrom','start','end','strand','intervalID','cluster')
expand = 20000#20000
intervs  = intervs %>% filter(cluster%in% c('cl4'))
intervs = intervs %>% arrange(chrom,start)
p1k=intervs
n_domains = dim(p1k)[1]
total_dom_length = sum(p1k$end-p1k$start)
p1k_pad = p1k
p1k_pad$start = p1k$start - expand
p1k_pad$end = p1k$end + expand
p1k_pad = gintervals.force_range(p1k_pad)
p1k_pad = gintervals.canonic(p1k_pad)
gvtrack.create(vtrack = 'dist_peak_enc_k27me3',src=p1k,func = 'distance' )
peak_distr = gextract(c(tracks,'dist_peak_enc_k27me3'),intervals = p1k_pad, iterator = 20)
colnames(peak_distr) = c('chrom','start','end',tracks,'dist','intervalID')
peak_distr$chrom = as.character(peak_distr$chrom)
temp = peak_distr[,c(tracks)]
temp[is.na(temp)]=0
temp[temp<0]=0
temp_n = temp
peak_distr = cbind(peak_distr[,c('chrom','start','end','dist','intervalID')],temp_n)
colnames(peak_distr) = c('chrom','start','end','dist','intervalID','cntr','inh')

peak_plot = peak_distr %>% group_by(dist)%>% summarise(m_1 = mean(cntr),m_2 = mean(inh))



back = c(mean(peak_plot[abs(peak_plot$dist) < 1000,]$m_1),mean(peak_plot[abs(peak_plot$dist) < 1000,]$m_2))
names(back) = c('cntr','inh')

intervs = mod$cgdom_ann[, c('chrom','start','end','tss_strand','ID','type')]
colnames(intervs) = c('chrom','start','end','strand','intervalID','cluster')
expand = 20000#20000
intervs  = intervs %>% filter(cluster%in% c(cl))
intervs = intervs %>% arrange(chrom,start)
p1k=intervs
n_domains = dim(p1k)[1]
total_dom_length = sum(p1k$end-p1k$start)
p1k_pad = p1k
p1k_pad$start = p1k$start - expand
p1k_pad$end = p1k$end + expand
p1k_pad = gintervals.force_range(p1k_pad)
p1k_pad = gintervals.canonic(p1k_pad)
gvtrack.create(vtrack = 'dist_peak_enc_k27me3',src=p1k,func = 'distance' )
peak_distr = gextract(c(tracks,'dist_peak_enc_k27me3'),intervals = p1k_pad, iterator = 20)
colnames(peak_distr) = c('chrom','start','end',tracks,'dist','intervalID')
peak_distr$chrom = as.character(peak_distr$chrom)
temp = peak_distr[,c(tracks)]
temp[is.na(temp)]=0
temp[temp<0]=0
temp_n = temp
peak_distr = cbind(peak_distr[,c('chrom','start','end','dist','intervalID')],temp_n)
colnames(peak_distr) = c('chrom','start','end','dist','intervalID','cntr','inh')

peak_plot = peak_distr %>% group_by(dist)%>% summarise(m_1 = mean(cntr),m_2 = mean(inh))



max=c(max(peak_plot$m_1),max(peak_plot$m_2))
names(max) = c('cntr','inh')



peak_plot$m_1 = peak_plot$m_1/as.numeric(back[1])

peak_plot$m_2 = peak_plot$m_2/as.numeric(back[2])


peak_plot$m_1 = rollmean(peak_plot$m_1, 30, na.rm=T, f='e')
peak_plot$m_2 = rollmean(peak_plot$m_2, 30, na.rm=T, f='e')

max/back
peak_plot_l = pivot_longer(peak_plot,-c(1))

options( repr.plot.width=10,repr.plot.height=10)
gg=peak_plot_l %>% ggplot(aes(x=dist,y=value,col=name))+geom_path(size=.5,alpha=1)+
scale_color_manual(values = c('darkgreen','darkviolet'))+
ggtitle(paste0('cl = ',cl,'  n = ', nrow(intervs), 'enr wt = ',as.numeric(max/back)[1], 'enr inh = ',as.numeric(max/back)[2]))+

theme_bw()+
      theme(plot.title = element_text(hjust = 0.5,size = 5),plot.subtitle = element_text(size = 5,hjust = 0.5),
              axis.text.y = element_text(size = 5),
            axis.text.x = element_text(size = 5,angle=90),axis.title = element_text(size = 5))

print(gg)
    return(gg)

}



#####jk







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
# endregion



pcg_build_gw_feats = function(mod, brz='lm')
{
	all_c = gintervals.all()
	all_c = all_c[!all_c$chrom %in% c("chrM", "chrY",'chrX'),]

	options(gmax.data.size=1e+9)
	gen_track_vext()
	gen_k27_vt(mod)
	if (brz == 'brz') {
	message('brz')
	seeds <- readRDS("brz2k_seeds_for_fig2.rds") %>% arrange(chrom, start)
	seeds$pred_seed_k27 <- seeds$brz2k_max
	seeds$pred_seed_k4  <- log2(6 + seeds$brz2k_k4_max)
	
	} else if (brz == 'lm_co') {
	message('lm_co')
	seeds <- readRDS("lm_co_seeds_for_fig2.rds")
	seeds$pred_seed_k27 <- seeds$pred_lm_k27_co35_d
	seeds$pred_seed_k4  <- seeds$pred_lm_k4_co35_d
	
	} else {
	message('lm 10% test noX final')
	seeds = readRDS('./data/lm_10test_noX_seeds_logist_dinucs_logist_fig1.rds')

	seeds$pred_seed_k27 = seeds$k27_pred
	seeds$pred_seed_k4 = seeds$k4_pred
		}
	
	f_pcg = seeds$modality=='pcg'
	f_txg = seeds$modality=='txg'
	f_mix = seeds$modality=='mix'
	
	cgd_pcg = seeds[f_pcg,1:3]
	cgd_txg = seeds[f_txg,1:3]
	cgd_mix = seeds[f_mix,1:3]
print('ok')
	gvtrack.create("d_cgd_pcg", cgd_pcg, "distance")
	gvtrack.create("d_cgd_txg", cgd_txg, "distance")
	gvtrack.create("d_cgd_mix", cgd_mix, "distance")
	cg_trace = gextract(c("d_cgd_pcg", "d_cgd_txg", "d_cgd_mix", "atac_ext", "EB4_cnt", "EB3_cnt", "ES_cnt+ES2i_cnt", "e75_ecto_cnt", "e75_emeso_cnt","ifelse(seq.CG_500_mean_new > 0.04, 1, 0)"), intervals=all_c, iterator = 200, colnames = c("d_pcg", "d_txg", "d_mix", "atac", "k27", "k27_eb3","k27_es", "k27_ecto","k27_emeso","cg_all"))

	cgd_pcg$strand = 1
	cgd_txg$strand = 1
	cgd_mix$strand = 1

	pcg_5 = gintervals.neighbors(cg_trace, cgd_pcg, mindist=0, na=T)
	pcg_3 = gintervals.neighbors(cg_trace, cgd_pcg, maxdist=0, na=T)
	mix_5 = gintervals.neighbors(cg_trace, cgd_mix, mindist=0, na=T)
	mix_3 = gintervals.neighbors(cg_trace, cgd_mix, maxdist=0, na=T)
	txg_5 = gintervals.neighbors(cg_trace, cgd_txg, mindist=0, na=T)
	txg_3 = gintervals.neighbors(cg_trace, cgd_txg, maxdist=0, na=T)

	cg_trace$d_pcg5 = ifelse(is.na(pcg_5$dist), 4e+6, pmin(pcg_5$dist, 4e+6))
	cg_trace$d_pcg3 = ifelse(is.na(pcg_3$dist), -4e+6, pmax(pcg_3$dist, -4e+6))
	cg_trace$d_mix5 = ifelse(is.na(mix_5$dist), 4e+6, pmin(mix_5$dist, 4e+6))
	cg_trace$d_mix3 = ifelse(is.na(mix_3$dist), -4e+6, pmax(mix_3$dist, -4e+6))
	cg_trace$d_txg5 = ifelse(is.na(txg_5$dist), 4e+6, pmin(txg_5$dist, 4e+6))
	cg_trace$d_txg3 = ifelse(is.na(txg_3$dist), -4e+6, pmax(txg_3$dist, -4e+6))

	gc_trace = gextract(c("seq.GC500_bin20", "seq.CG_500_mean_new",
								"seq.IQ.epiblast", "seq.IQ.deeptopic_cnn","seq.IQ.deeptopic_cnn"), 
						intervals=all_c, iterator = 200, 
						colnames = c("GC", "CG", "IQbase", "tn5bias","deep_base") )
						

	f = is.infinite(gc_trace$GC) | is.na(gc_trace$GC)
	gc_traceGC = 0.3

	gc_trace$deep_base[is.na(gc_trace$deep_base)] = 
									min(cg_trace$deep_base,na.rm=T)
	gc_trace$deep_base[is.infinite(gc_trace$deep_base)] = 
									max(cg_trace$deep_base,na.rm=T)
	gc_trace$GC[is.na(gc_trace$GC)] = 
									min(cg_trace$GC,na.rm=T)
	cg_trace$GC = gc_trace$GC
	cg_trace$CG = gc_trace$CG
	cg_trace$IQbase = gc_trace$IQbase
	cg_trace$deep_base = gc_trace$deep_base
	iq_T = quantile(cg_trace$IQbase, 0.98)
	dn_T = quantile(cg_trace$deep_base, 0.98)

	iqb_dens10 = rollmean(ifelse(cg_trace$IQbase>iq_T,1,0),10, f='e')
	iqb_dens100 = rollmean(ifelse(cg_trace$IQbase>iq_T,1,0),100, f='e')
	iqb_dens1k = rollmean(ifelse(cg_trace$IQbase>iq_T,1,0),1000, f='e')

	deepb_dens10 = rollmean(ifelse(cg_trace$deep_base > dn_T,1,0),10, f='e')
	deepb_dens100 = rollmean(ifelse(cg_trace$deep_base > dn_T,1,0),100, f='e')
	deepb_dens1k = rollmean(ifelse(cg_trace$deep_base > dn_T,1,0),1000, f='e')
	
	iqbase_feats = cbind(cg_trace$IQbase, iqb_dens10, iqb_dens100, iqb_dens1k)
	colnames(iqbase_feats) = c("IQ","IQ2k", "IQ20k", "IQ200k")
	deepbase_feats = cbind(cg_trace$deep_base, deepb_dens10, deepb_dens100, deepb_dens1k)
	colnames(deepbase_feats) = c("DeepTop","DeepTop2k", "DeepTop20k", "DeepTop200k")

	#cg_trace$tn5bias = gc_trace$tn5bias

	gvtrack.create("d_ltr", "intervs.global.rmsk_ltr", "distance")
	gvtrack.create("d_line", "intervs.global.rmsk_line", "distance")
	gvtrack.create("d_sine", "intervs.global.rmsk_sine", "distance")
	rpt_trace = gextract(c("d_line","d_ltr", "d_sine"), intervals=all_c,iterator=200)

	cg_trace$d_ltr = rpt_trace$d_ltr	
	cg_trace$d_line = rpt_trace$d_line
	cg_trace$d_sine = rpt_trace$d_sine

	cg_trace$k27[is.na(cg_trace$k27)] = 0
	cg_trace$atac[is.na(cg_trace$atac)] = 0

	is_pcg = ifelse(cg_trace$d_pcg==0, 1,0)
	is_txg = ifelse(cg_trace$d_txg==0, 1,0)
	is_mix = ifelse(cg_trace$d_mix==0, 1,0)

	#k27_1k = rollmean(cg_trace$k27, 5, f='e')
	k27_1k = cg_trace$k27# ##already loged smoothed -+500
	cg_trace$k27_1k = k27_1k
	#cg_trace$lk27_1k = log2(2+k27_1k)
	cg_trace$lk27_1k = k27_1k

	cg_sc = matrix(cg_trace$cg_all, ncol=1)
	cgpcg_sc = matrix(is_pcg, ncol=1)
	cgtxg_sc = matrix(is_txg, ncol=1)
	cgmix_sc = matrix(is_mix, ncol=1)
	#we rely on the facto that chromosomes are all padded by trivial MBs - so no wraparound problems
	for(scale in 2**(1:12)) {
		message("scale ", scale)
		cg_sc = cbind(cg_sc, rollmean(cg_trace$cg_all,1+scale,f='e'))
		cgpcg_sc = cbind(cgpcg_sc, rollmean(is_pcg, 1+scale, f='e'))
		cgtxg_sc = cbind(cgtxg_sc, rollmean(is_txg, 1+scale, f='e'))
		cgmix_sc = cbind(cgmix_sc, rollmean(is_mix, 1+scale, f='e'))
	}
	colnames(cgmix_sc) = paste(rep("cgmix",13), 1:13, sep="_")
	colnames(cgpcg_sc) = paste(rep("cgpcg",13), 1:13, sep="_")
	colnames(cgtxg_sc) = paste(rep("cgtxg",13), 1:13, sep="_")

	n = nrow(cg_trace)
	for(i in 6:12) {
		scale = 2**(i-1)
		half_win = scale/2
		pcg_pad3 = c(cgpcg_sc[-(1:half_win), i], rep(cgpcg_sc[n,i], half_win))	
		pcg_pad5 = c(rep(cgpcg_sc[1,i], half_win), cgpcg_sc[-((n-half_win+1):n), i])
		txg_pad3 = c(cgtxg_sc[-(1:half_win), i], rep(cgtxg_sc[n,i], half_win))	
		txg_pad5 = c(rep(cgtxg_sc[1,i], half_win), cgtxg_sc[-((n-half_win+1):n), i])
		mix_pad3 = c(cgmix_sc[-(1:half_win), i], rep(cgmix_sc[n,i], half_win))	
		mix_pad5 = c(rep(cgmix_sc[1,i], half_win), cgmix_sc[-((n-half_win+1):n), i])

		if(i == 6) {
			cgpcg_3 = matrix(pcg_pad3, ncol=1)
			cgpcg_5 = matrix(pcg_pad5, ncol=1)
			cgtxg_3 = matrix(txg_pad3, ncol=1)
			cgtxg_5 = matrix(txg_pad5, ncol=1)
			cgmix_3 = matrix(mix_pad3, ncol=1)
			cgmix_5 = matrix(mix_pad5, ncol=1)
		} else {
			cgpcg_3 = cbind(cgpcg_3, pcg_pad3)
			cgpcg_5 = cbind(cgpcg_5, pcg_pad5)
			cgtxg_3 = cbind(cgtxg_3, txg_pad3)
			cgtxg_5 = cbind(cgtxg_5, txg_pad5)
			cgmix_3 = cbind(cgmix_3, mix_pad3)
			cgmix_5 = cbind(cgmix_5, mix_pad5)
		}
	}
	colnames(cgmix_3) = paste(rep("cgmix3",7), 1:7, sep="_")
	colnames(cgpcg_3) = paste(rep("cgpcg3",7), 1:7, sep="_")
	colnames(cgtxg_3) = paste(rep("cgtxg3",7), 1:7, sep="_")

	colnames(cgmix_5) = paste(rep("cgmix5",7), 1:7, sep="_")
	colnames(cgpcg_5) = paste(rep("cgpcg5",7), 1:7, sep="_")
	colnames(cgtxg_5) = paste(rep("cgtxg5",7), 1:7, sep="_")

	inv_d_txg = 1/(cg_trace$d_txg+100)
	inv_d_mix = 1/(cg_trace$d_mix+100)
	inv_d_pcg = 1/(cg_trace$d_pcg+100)

	min_d_cg = pmax(pmax(inv_d_txg, inv_d_mix), inv_d_pcg)
	brz_score = misha.ext::gintervals.neighbors1(cg_trace[,c('chrom','start','end')],
	seeds[,c('chrom','start','end','pred_seed_k27','pred_seed_k4')],maxneighbors = 1)
	
	gvtrack.create('GC20_l1','seq.GC20_bin20','avg')
gvtrack.iterator('GC20_l1',sshift = -400,eshift = -100)
gvtrack.create('GC20_l2','seq.GC20_bin20','avg')
gvtrack.iterator('GC20_l2',sshift = -900,eshift = -600)
gvtrack.create('GC20_l3','seq.GC20_bin20','avg')
gvtrack.iterator('GC20_l3',sshift = -1400,eshift = -1100)
gvtrack.create('GC20_l4','seq.GC20_bin20','avg')
gvtrack.iterator('GC20_l4',sshift = -1900,eshift = -1600)

gvtrack.create('GC20_r1','seq.GC20_bin20','avg')
gvtrack.iterator('GC20_r1',sshift = 100,eshift = 400)
gvtrack.create('GC20_r2','seq.GC20_bin20','avg')
gvtrack.iterator('GC20_r2',sshift = 600,eshift = 900)
gvtrack.create('GC20_r3','seq.GC20_bin20','avg')
gvtrack.iterator('GC20_r3',sshift = 1100,eshift = 1400)
gvtrack.create('GC20_r4','seq.GC20_bin20','avg')
gvtrack.iterator('GC20_r4',sshift = 1600,eshift = 1900)

featsGC = gextract(c('GC20_l1/20','GC20_l2/20','GC20_l3/20','GC20_l4/20',
                    'GC20_r1/20','GC20_r2/20','GC20_r3/20','GC20_r4/20'),intervals = cg_trace,iterator = cg_trace,
                  colnames = c('GC20_l1','GC20_l2','GC20_l3','GC20_l4',
                    'GC20_r1','GC20_r2','GC20_r3','GC20_r4'))
featsGC = featsGC[,c('GC20_l1','GC20_l1','GC20_l2','GC20_l3','GC20_l4',
                    'GC20_r1','GC20_r2','GC20_r3','GC20_r4')]
featsGC[is.na(featsGC)] = 0


gvtrack.create('CG20_l1','seq.CG_500_mean_new','avg')
gvtrack.iterator('CG20_l1',sshift = -400,eshift = -100)
gvtrack.create('CG20_l2','seq.CG_500_mean_new','avg')
gvtrack.iterator('CG20_l2',sshift = -900,eshift = -600)
gvtrack.create('CG20_l3','seq.CG_500_mean_new','avg')
gvtrack.iterator('CG20_l3',sshift = -1400,eshift = -1100)
gvtrack.create('CG20_l4','seq.CG_500_mean_new','avg')
gvtrack.iterator('CG20_l4',sshift = -1900,eshift = -1600)

gvtrack.create('CG20_r1','seq.CG_500_mean_new','avg')
gvtrack.iterator('CG20_r1',sshift = 100,eshift = 400)
gvtrack.create('CG20_r2','seq.CG_500_mean_new','avg')
gvtrack.iterator('CG20_r2',sshift = 600,eshift = 900)
gvtrack.create('CG20_r3','seq.CG_500_mean_new','avg')
gvtrack.iterator('CG20_r3',sshift = 1100,eshift = 1400)
gvtrack.create('CG20_r4','seq.CG_500_mean_new','avg')
gvtrack.iterator('CG20_r4',sshift = 1600,eshift = 1900)
###dont divide cg by 20, but doesn matter just scaling
featsCG = gextract(c('CG20_l1/20','CG20_l2/20','CG20_l3/20','CG20_l4/20',
                    'CG20_r1/20','CG20_r2/20','CG20_r3/20','CG20_r4/20'),intervals = cg_trace,iterator = cg_trace,
                  colnames = c('CG20_l1','CG20_l2','CG20_l3','CG20_l4',
                    'CG20_r1','CG20_r2','CG20_r3','CG20_r4'))
featsCG = featsCG[,c('CG20_l1','CG20_l1','CG20_l2','CG20_l3','CG20_l4',
                    'CG20_r1','CG20_r2','CG20_r3','CG20_r4')]
featsCG[is.na(featsCG)] = 0

	
	
	feats = cbind(ld_mix=log2(inv_d_mix), 
							ld_pcg = log2(inv_d_pcg),
							ld_txg = log2(inv_d_txg),
							min_dcg = min_d_cg,
							cgmix_sc[,1:11],
							cgtxg_sc[,1:11], 
							cgpcg_sc[,1:11],
							GC = gc_trace$GC,
							CG = gc_trace$CG,
							#tn5b = gc_trace$tn5bias,
							dltr = abs(rpt_trace$d_ltr),
							dline = abs(rpt_trace$d_line),
							ldltr = log2(1+abs(rpt_trace$d_ltr)),
							ldline = log2(1+abs(rpt_trace$d_line)),
							pred_seed_k27 = as.numeric(brz_score$pred_seed_k27),
							pred_seed_k4 = as.numeric(brz_score$pred_seed_k4),
							featsGC,featsCG)

	feats35 = cbind(cgmix_3, cgmix_5,
							cgtxg_3, cgtxg_5,
							cgpcg_3, cgpcg_5)
	feats_iqdn = cbind(deepbase_feats,
							iqbase_feats)

	mod$gw = list()
	mod$gw$cg_trace = cg_trace
	mod$gw$cgpcg_sc = cgpcg_sc
	mod$gw$cgtxg_sc = cgtxg_sc
	mod$gw$cgmix_sc = cgmix_sc
	mod$gw$feats = feats
	mod$gw$feats35 = feats35
	mod$gw$feats_iqdn = feats_iqdn

	for(i in 1:ncol(mod$gw$feats_iqdn)) {
		minv = min(mod$gw$feats_iqdn[,i], na.rm=T)
		maxv = max(mod$gw$feats_iqdn[,i], na.rm=T)
		f = is.na(mod$gw$feats_iqdn[,i])
		f1 = is.infinite(mod$gw$feats_iqdn[,i])
		if(sum(f)>0) {
			message("replace ", sum(f), " nas with ", minv, " in ", colnames(mod$gw$feats_iqdn)[i])
			mod$gw$feats_iqdn[f,i] = minv
		}
		if(sum(f1)>0) {
			mod$gw$feats_iqdn[f1,i] = maxv
		}
	}
	f = is.infinite(mod$gw$feats[,"GC"]) | is.na(mod$gw$feats[,"GC"])
	mod$gw$feats[f,"GC"] = 0.3
	
	return(mod)
}
pcg_build_gw_feats_k4 = function(mod, quick_mode=F)
{
	all_c = gintervals.all()
	all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]

	options(gmax.data.size=1e+9)
	gen_track_vext()
	gen_k27_vt(mod)
	gen_k4_vt(mod)
	f_pcg = mod$seqmod_loc_pred > 6.8
	f_txg = mod$seqmod_loc_pred < 5.5
	f_mix = !f_pcg & !f_txg
	cgd_pcg = mod$cgdom_ann[f_pcg,1:3]
	cgd_txg = mod$cgdom_ann[f_txg,1:3]
	cgd_mix = mod$cgdom_ann[f_mix,1:3]
print('ok')
	gvtrack.create("d_cgd_pcg", cgd_pcg, "distance")
	gvtrack.create("d_cgd_txg", cgd_txg, "distance")
	gvtrack.create("d_cgd_mix", cgd_mix, "distance")
	cg_trace = gextract(c("d_cgd_pcg", "d_cgd_txg", "d_cgd_mix", "atac_ext", "EB4_cnt_k4", "EB3_cnt", "ES_cnt+ES2i_cnt", "e75_ecto_cnt", "e75_emeso_cnt","ifelse(seq.CG_500_mean_new > 0.04, 1, 0)"), intervals=all_c, iterator = 200, colnames = c("d_pcg", "d_txg", "d_mix", "atac", "k27", "k27_eb3","k27_es", "k27_ecto","k27_emeso","cg_all"))

	cgd_pcg$strand = 1
	cgd_txg$strand = 1
	cgd_mix$strand = 1

	pcg_5 = gintervals.neighbors(cg_trace, cgd_pcg, mindist=0, na=T)
	pcg_3 = gintervals.neighbors(cg_trace, cgd_pcg, maxdist=0, na=T)
	mix_5 = gintervals.neighbors(cg_trace, cgd_mix, mindist=0, na=T)
	mix_3 = gintervals.neighbors(cg_trace, cgd_mix, maxdist=0, na=T)
	txg_5 = gintervals.neighbors(cg_trace, cgd_txg, mindist=0, na=T)
	txg_3 = gintervals.neighbors(cg_trace, cgd_txg, maxdist=0, na=T)

	cg_trace$d_pcg5 = ifelse(is.na(pcg_5$dist), 4e+6, pmin(pcg_5$dist, 4e+6))
	cg_trace$d_pcg3 = ifelse(is.na(pcg_3$dist), -4e+6, pmax(pcg_3$dist, -4e+6))
	cg_trace$d_mix5 = ifelse(is.na(mix_5$dist), 4e+6, pmin(mix_5$dist, 4e+6))
	cg_trace$d_mix3 = ifelse(is.na(mix_3$dist), -4e+6, pmax(mix_3$dist, -4e+6))
	cg_trace$d_txg5 = ifelse(is.na(txg_5$dist), 4e+6, pmin(txg_5$dist, 4e+6))
	cg_trace$d_txg3 = ifelse(is.na(txg_3$dist), -4e+6, pmax(txg_3$dist, -4e+6))

	gc_trace = gextract(c("seq.GC500_bin20", "seq.CG_500_mean_new",
								"seq.IQ.epiblast", "seq.IQ.deeptopic_cnn","seq.IQ.deeptopic_cnn"), 
						intervals=all_c, iterator = 200, 
						colnames = c("GC", "CG", "IQbase", "tn5bias","deep_base") )
						

	f = is.infinite(gc_trace$GC) | is.na(gc_trace$GC)
	gc_traceGC = 0.3

	gc_trace$deep_base[is.na(gc_trace$deep_base)] = 
									min(cg_trace$deep_base,na.rm=T)
	gc_trace$deep_base[is.infinite(gc_trace$deep_base)] = 
									max(cg_trace$deep_base,na.rm=T)
	gc_trace$GC[is.na(gc_trace$GC)] = 
									min(cg_trace$GC,na.rm=T)
	cg_trace$GC = gc_trace$GC
	cg_trace$CG = gc_trace$CG
	cg_trace$IQbase = gc_trace$IQbase
	cg_trace$deep_base = gc_trace$deep_base
	iq_T = quantile(cg_trace$IQbase, 0.98)
	dn_T = quantile(cg_trace$deep_base, 0.98)

	iqb_dens10 = rollmean(ifelse(cg_trace$IQbase>iq_T,1,0),10, f='e')
	iqb_dens100 = rollmean(ifelse(cg_trace$IQbase>iq_T,1,0),100, f='e')
	iqb_dens1k = rollmean(ifelse(cg_trace$IQbase>iq_T,1,0),1000, f='e')

	deepb_dens10 = rollmean(ifelse(cg_trace$deep_base > dn_T,1,0),10, f='e')
	deepb_dens100 = rollmean(ifelse(cg_trace$deep_base > dn_T,1,0),100, f='e')
	deepb_dens1k = rollmean(ifelse(cg_trace$deep_base > dn_T,1,0),1000, f='e')
	
	iqbase_feats = cbind(cg_trace$IQbase, iqb_dens10, iqb_dens100, iqb_dens1k)
	colnames(iqbase_feats) = c("IQ","IQ2k", "IQ20k", "IQ200k")
	deepbase_feats = cbind(cg_trace$deep_base, deepb_dens10, deepb_dens100, deepb_dens1k)
	colnames(deepbase_feats) = c("DeepTop","DeepTop2k", "DeepTop20k", "DeepTop200k")

	#cg_trace$tn5bias = gc_trace$tn5bias

	gvtrack.create("d_ltr", "intervs.global.rmsk_ltr", "distance")
	gvtrack.create("d_line", "intervs.global.rmsk_line", "distance")
	gvtrack.create("d_sine", "intervs.global.rmsk_sine", "distance")
	rpt_trace = gextract(c("d_line","d_ltr", "d_sine"), intervals=all_c,iterator=200)

	cg_trace$d_ltr = rpt_trace$d_ltr	
	cg_trace$d_line = rpt_trace$d_line
	cg_trace$d_sine = rpt_trace$d_sine

	cg_trace$k27[is.na(cg_trace$k27)] = 0
	cg_trace$atac[is.na(cg_trace$atac)] = 0

	is_pcg = ifelse(cg_trace$d_pcg==0, 1,0)
	is_txg = ifelse(cg_trace$d_txg==0, 1,0)
	is_mix = ifelse(cg_trace$d_mix==0, 1,0)

	k27_1k = rollmean(cg_trace$k27, 5, f='e')
	cg_trace$k27_1k = k27_1k
	cg_trace$lk27_1k = log2(2+k27_1k)

	cg_sc = matrix(cg_trace$cg_all, ncol=1)
	cgpcg_sc = matrix(is_pcg, ncol=1)
	cgtxg_sc = matrix(is_txg, ncol=1)
	cgmix_sc = matrix(is_mix, ncol=1)
	#we rely on the facto that chromosomes are all padded by trivial MBs - so no wraparound problems
	for(scale in 2**(1:12)) {
		message("scale ", scale)
		cg_sc = cbind(cg_sc, rollmean(cg_trace$cg_all,1+scale,f='e'))
		cgpcg_sc = cbind(cgpcg_sc, rollmean(is_pcg, 1+scale, f='e'))
		cgtxg_sc = cbind(cgtxg_sc, rollmean(is_txg, 1+scale, f='e'))
		cgmix_sc = cbind(cgmix_sc, rollmean(is_mix, 1+scale, f='e'))
	}
	colnames(cgmix_sc) = paste(rep("cgmix",13), 1:13, sep="_")
	colnames(cgpcg_sc) = paste(rep("cgpcg",13), 1:13, sep="_")
	colnames(cgtxg_sc) = paste(rep("cgtxg",13), 1:13, sep="_")

	n = nrow(cg_trace)
	for(i in 6:12) {
		scale = 2**(i-1)
		half_win = scale/2
		pcg_pad3 = c(cgpcg_sc[-(1:half_win), i], rep(cgpcg_sc[n,i], half_win))	
		pcg_pad5 = c(rep(cgpcg_sc[1,i], half_win), cgpcg_sc[-((n-half_win+1):n), i])
		txg_pad3 = c(cgtxg_sc[-(1:half_win), i], rep(cgtxg_sc[n,i], half_win))	
		txg_pad5 = c(rep(cgtxg_sc[1,i], half_win), cgtxg_sc[-((n-half_win+1):n), i])
		mix_pad3 = c(cgmix_sc[-(1:half_win), i], rep(cgmix_sc[n,i], half_win))	
		mix_pad5 = c(rep(cgmix_sc[1,i], half_win), cgmix_sc[-((n-half_win+1):n), i])

		if(i == 6) {
			cgpcg_3 = matrix(pcg_pad3, ncol=1)
			cgpcg_5 = matrix(pcg_pad5, ncol=1)
			cgtxg_3 = matrix(txg_pad3, ncol=1)
			cgtxg_5 = matrix(txg_pad5, ncol=1)
			cgmix_3 = matrix(mix_pad3, ncol=1)
			cgmix_5 = matrix(mix_pad5, ncol=1)
		} else {
			cgpcg_3 = cbind(cgpcg_3, pcg_pad3)
			cgpcg_5 = cbind(cgpcg_5, pcg_pad5)
			cgtxg_3 = cbind(cgtxg_3, txg_pad3)
			cgtxg_5 = cbind(cgtxg_5, txg_pad5)
			cgmix_3 = cbind(cgmix_3, mix_pad3)
			cgmix_5 = cbind(cgmix_5, mix_pad5)
		}
	}
	colnames(cgmix_3) = paste(rep("cgmix3",7), 1:7, sep="_")
	colnames(cgpcg_3) = paste(rep("cgpcg3",7), 1:7, sep="_")
	colnames(cgtxg_3) = paste(rep("cgtxg3",7), 1:7, sep="_")

	colnames(cgmix_5) = paste(rep("cgmix5",7), 1:7, sep="_")
	colnames(cgpcg_5) = paste(rep("cgpcg5",7), 1:7, sep="_")
	colnames(cgtxg_5) = paste(rep("cgtxg5",7), 1:7, sep="_")

	inv_d_txg = 1/(cg_trace$d_txg+100)
	inv_d_mix = 1/(cg_trace$d_mix+100)
	inv_d_pcg = 1/(cg_trace$d_pcg+100)

	min_d_cg = pmax(pmax(inv_d_txg, inv_d_mix), inv_d_pcg)
	iq_score = misha.ext::gintervals.neighbors1(cg_trace[,c('chrom','start','end')],mod$cgdom_ann[,c('chrom','start','end','pred')],maxneighbors = 1)
	

	
	
	
	
	
	feats = cbind(ld_mix=log2(inv_d_mix), 
							ld_pcg = log2(inv_d_pcg),
							ld_txg = log2(inv_d_txg),
							min_dcg = min_d_cg,
							cgmix_sc[,1:11],
							cgtxg_sc[,1:11], 
							cgpcg_sc[,1:11],
							GC = gc_trace$GC,
							CG = gc_trace$CG,
							#tn5b = gc_trace$tn5bias,
							dltr = abs(rpt_trace$d_ltr),
							dline = abs(rpt_trace$d_line),
							ldltr = log2(1+abs(rpt_trace$d_ltr)),
							ldline = log2(1+abs(rpt_trace$d_line)),
							iq_score = as.numeric(iq_score$pred));

	feats35 = cbind(cgmix_3, cgmix_5,
							cgtxg_3, cgtxg_5,
							cgpcg_3, cgpcg_5)
	feats_iqdn = cbind(deepbase_feats,
							iqbase_feats)

	mod$gw = list()
	mod$gw$cg_trace = cg_trace
	mod$gw$cgpcg_sc = cgpcg_sc
	mod$gw$cgtxg_sc = cgtxg_sc
	mod$gw$cgmix_sc = cgmix_sc
	mod$gw$feats = feats
	mod$gw$feats35 = feats35
	mod$gw$feats_iqdn = feats_iqdn

	for(i in 1:ncol(mod$gw$feats_iqdn)) {
		minv = min(mod$gw$feats_iqdn[,i], na.rm=T)
		maxv = max(mod$gw$feats_iqdn[,i], na.rm=T)
		f = is.na(mod$gw$feats_iqdn[,i])
		f1 = is.infinite(mod$gw$feats_iqdn[,i])
		if(sum(f)>0) {
			message("replace ", sum(f), " nas with ", minv, " in ", colnames(mod$gw$feats_iqdn)[i])
			mod$gw$feats_iqdn[f,i] = minv
		}
		if(sum(f1)>0) {
			mod$gw$feats_iqdn[f1,i] = maxv
		}
	}
	f = is.infinite(mod$gw$feats[,"GC"]) | is.na(mod$gw$feats[,"GC"])
	mod$gw$feats[f,"GC"] = 0.3
	
	return(mod)
}
pcg_build_gw_feats_brz = function(mod, quick_mode=F)
{
	all_c = gintervals.all()
	all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]

	options(gmax.data.size=1e+9)
	gen_track_vext()
	gen_k27_vt(mod)
	f_pcg = mod$seqmod_loc_pred > 6.8
	f_txg = mod$seqmod_loc_pred < 5
	f_mix = !f_pcg & !f_txg
	cgd_pcg = mod$cgdom_ann[f_pcg,1:3]
	cgd_txg = mod$cgdom_ann[f_txg,1:3]
	cgd_mix = mod$cgdom_ann[f_mix,1:3]
print('ok')
	gvtrack.create("d_cgd_pcg", cgd_pcg, "distance")
	gvtrack.create("d_cgd_txg", cgd_txg, "distance")
	gvtrack.create("d_cgd_mix", cgd_mix, "distance")
	cg_trace = gextract(c("d_cgd_pcg", "d_cgd_txg", "d_cgd_mix", "atac_ext", "EB4_cnt", "EB3_cnt", "ES_cnt+ES2i_cnt", "e75_ecto_cnt", "e75_emeso_cnt","ifelse(seq.CG_500_mean_new > 0.04, 1, 0)"), intervals=all_c, iterator = 200, colnames = c("d_pcg", "d_txg", "d_mix", "atac", "k27", "k27_eb3","k27_es", "k27_ecto","k27_emeso","cg_all"))

	cgd_pcg$strand = 1
	cgd_txg$strand = 1
	cgd_mix$strand = 1

	pcg_5 = gintervals.neighbors(cg_trace, cgd_pcg, mindist=0, na=T)
	pcg_3 = gintervals.neighbors(cg_trace, cgd_pcg, maxdist=0, na=T)
	mix_5 = gintervals.neighbors(cg_trace, cgd_mix, mindist=0, na=T)
	mix_3 = gintervals.neighbors(cg_trace, cgd_mix, maxdist=0, na=T)
	txg_5 = gintervals.neighbors(cg_trace, cgd_txg, mindist=0, na=T)
	txg_3 = gintervals.neighbors(cg_trace, cgd_txg, maxdist=0, na=T)

	cg_trace$d_pcg5 = ifelse(is.na(pcg_5$dist), 4e+6, pmin(pcg_5$dist, 4e+6))
	cg_trace$d_pcg3 = ifelse(is.na(pcg_3$dist), -4e+6, pmax(pcg_3$dist, -4e+6))
	cg_trace$d_mix5 = ifelse(is.na(mix_5$dist), 4e+6, pmin(mix_5$dist, 4e+6))
	cg_trace$d_mix3 = ifelse(is.na(mix_3$dist), -4e+6, pmax(mix_3$dist, -4e+6))
	cg_trace$d_txg5 = ifelse(is.na(txg_5$dist), 4e+6, pmin(txg_5$dist, 4e+6))
	cg_trace$d_txg3 = ifelse(is.na(txg_3$dist), -4e+6, pmax(txg_3$dist, -4e+6))

	gc_trace = gextract(c("seq.GC500_bin20", "seq.CG_500_mean_new",
								"seq.IQ.epiblast", "seq.IQ.deeptopic_cnn","seq.IQ.deeptopic_cnn"), 
						intervals=all_c, iterator = 200, 
						colnames = c("GC", "CG", "IQbase", "tn5bias","deep_base") )
						

	f = is.infinite(gc_trace$GC) | is.na(gc_trace$GC)
	gc_traceGC = 0.3

	gc_trace$deep_base[is.na(gc_trace$deep_base)] = 
									min(cg_trace$deep_base,na.rm=T)
	gc_trace$deep_base[is.infinite(gc_trace$deep_base)] = 
									max(cg_trace$deep_base,na.rm=T)
	gc_trace$GC[is.na(gc_trace$GC)] = 
									min(cg_trace$GC,na.rm=T)
	cg_trace$GC = gc_trace$GC
	cg_trace$CG = gc_trace$CG
	cg_trace$IQbase = gc_trace$IQbase
	cg_trace$deep_base = gc_trace$deep_base
	iq_T = quantile(cg_trace$IQbase, 0.98)
	dn_T = quantile(cg_trace$deep_base, 0.98)

	iqb_dens10 = rollmean(ifelse(cg_trace$IQbase>iq_T,1,0),10, f='e')
	iqb_dens100 = rollmean(ifelse(cg_trace$IQbase>iq_T,1,0),100, f='e')
	iqb_dens1k = rollmean(ifelse(cg_trace$IQbase>iq_T,1,0),1000, f='e')

	deepb_dens10 = rollmean(ifelse(cg_trace$deep_base > dn_T,1,0),10, f='e')
	deepb_dens100 = rollmean(ifelse(cg_trace$deep_base > dn_T,1,0),100, f='e')
	deepb_dens1k = rollmean(ifelse(cg_trace$deep_base > dn_T,1,0),1000, f='e')
	
	iqbase_feats = cbind(cg_trace$IQbase, iqb_dens10, iqb_dens100, iqb_dens1k)
	colnames(iqbase_feats) = c("IQ","IQ2k", "IQ20k", "IQ200k")
	deepbase_feats = cbind(cg_trace$deep_base, deepb_dens10, deepb_dens100, deepb_dens1k)
	colnames(deepbase_feats) = c("DeepTop","DeepTop2k", "DeepTop20k", "DeepTop200k")

	#cg_trace$tn5bias = gc_trace$tn5bias

	gvtrack.create("d_ltr", "intervs.global.rmsk_ltr", "distance")
	gvtrack.create("d_line", "intervs.global.rmsk_line", "distance")
	gvtrack.create("d_sine", "intervs.global.rmsk_sine", "distance")
	rpt_trace = gextract(c("d_line","d_ltr", "d_sine"), intervals=all_c,iterator=200)

	cg_trace$d_ltr = rpt_trace$d_ltr	
	cg_trace$d_line = rpt_trace$d_line
	cg_trace$d_sine = rpt_trace$d_sine

	cg_trace$k27[is.na(cg_trace$k27)] = 0
	cg_trace$atac[is.na(cg_trace$atac)] = 0

	is_pcg = ifelse(cg_trace$d_pcg==0, 1,0)
	is_txg = ifelse(cg_trace$d_txg==0, 1,0)
	is_mix = ifelse(cg_trace$d_mix==0, 1,0)

	k27_1k = rollmean(cg_trace$k27, 5, f='e')
	cg_trace$k27_1k = k27_1k
	cg_trace$lk27_1k = log2(2+k27_1k)

	cg_sc = matrix(cg_trace$cg_all, ncol=1)
	cgpcg_sc = matrix(is_pcg, ncol=1)
	cgtxg_sc = matrix(is_txg, ncol=1)
	cgmix_sc = matrix(is_mix, ncol=1)
	#we rely on the facto that chromosomes are all padded by trivial MBs - so no wraparound problems
	for(scale in 2**(1:12)) {
		message("scale ", scale)
		cg_sc = cbind(cg_sc, rollmean(cg_trace$cg_all,1+scale,f='e'))
		cgpcg_sc = cbind(cgpcg_sc, rollmean(is_pcg, 1+scale, f='e'))
		cgtxg_sc = cbind(cgtxg_sc, rollmean(is_txg, 1+scale, f='e'))
		cgmix_sc = cbind(cgmix_sc, rollmean(is_mix, 1+scale, f='e'))
	}
	colnames(cgmix_sc) = paste(rep("cgmix",13), 1:13, sep="_")
	colnames(cgpcg_sc) = paste(rep("cgpcg",13), 1:13, sep="_")
	colnames(cgtxg_sc) = paste(rep("cgtxg",13), 1:13, sep="_")

	n = nrow(cg_trace)
	for(i in 6:12) {
		scale = 2**(i-1)
		half_win = scale/2
		pcg_pad3 = c(cgpcg_sc[-(1:half_win), i], rep(cgpcg_sc[n,i], half_win))	
		pcg_pad5 = c(rep(cgpcg_sc[1,i], half_win), cgpcg_sc[-((n-half_win+1):n), i])
		txg_pad3 = c(cgtxg_sc[-(1:half_win), i], rep(cgtxg_sc[n,i], half_win))	
		txg_pad5 = c(rep(cgtxg_sc[1,i], half_win), cgtxg_sc[-((n-half_win+1):n), i])
		mix_pad3 = c(cgmix_sc[-(1:half_win), i], rep(cgmix_sc[n,i], half_win))	
		mix_pad5 = c(rep(cgmix_sc[1,i], half_win), cgmix_sc[-((n-half_win+1):n), i])

		if(i == 6) {
			cgpcg_3 = matrix(pcg_pad3, ncol=1)
			cgpcg_5 = matrix(pcg_pad5, ncol=1)
			cgtxg_3 = matrix(txg_pad3, ncol=1)
			cgtxg_5 = matrix(txg_pad5, ncol=1)
			cgmix_3 = matrix(mix_pad3, ncol=1)
			cgmix_5 = matrix(mix_pad5, ncol=1)
		} else {
			cgpcg_3 = cbind(cgpcg_3, pcg_pad3)
			cgpcg_5 = cbind(cgpcg_5, pcg_pad5)
			cgtxg_3 = cbind(cgtxg_3, txg_pad3)
			cgtxg_5 = cbind(cgtxg_5, txg_pad5)
			cgmix_3 = cbind(cgmix_3, mix_pad3)
			cgmix_5 = cbind(cgmix_5, mix_pad5)
		}
	}
	colnames(cgmix_3) = paste(rep("cgmix3",7), 1:7, sep="_")
	colnames(cgpcg_3) = paste(rep("cgpcg3",7), 1:7, sep="_")
	colnames(cgtxg_3) = paste(rep("cgtxg3",7), 1:7, sep="_")

	colnames(cgmix_5) = paste(rep("cgmix5",7), 1:7, sep="_")
	colnames(cgpcg_5) = paste(rep("cgpcg5",7), 1:7, sep="_")
	colnames(cgtxg_5) = paste(rep("cgtxg5",7), 1:7, sep="_")

	inv_d_txg = 1/(cg_trace$d_txg+100)
	inv_d_mix = 1/(cg_trace$d_mix+100)
	inv_d_pcg = 1/(cg_trace$d_pcg+100)

	min_d_cg = pmax(pmax(inv_d_txg, inv_d_mix), inv_d_pcg)
	iq_score = misha.ext::gintervals.neighbors1(cg_trace[,c('chrom','start','end')],
	mod$cgdom_ann[,c('chrom','start','end','l10_borzoi_mean')],maxneighbors = 1)
	feats = cbind(ld_mix=log2(inv_d_mix), 
							ld_pcg = log2(inv_d_pcg),
							ld_txg = log2(inv_d_txg),
							min_dcg = min_d_cg,
							cgmix_sc[,1:11],
							cgtxg_sc[,1:11], 
							cgpcg_sc[,1:11],
							GC = gc_trace$GC,
							CG = gc_trace$CG,
							#tn5b = gc_trace$tn5bias,
							dltr = abs(rpt_trace$d_ltr),
							dline = abs(rpt_trace$d_line),
							ldltr = log2(1+abs(rpt_trace$d_ltr)),
							ldline = log2(1+abs(rpt_trace$d_line)),
							iq_score = as.numeric(iq_score$l10_borzoi_mean));

	feats35 = cbind(cgmix_3, cgmix_5,
							cgtxg_3, cgtxg_5,
							cgpcg_3, cgpcg_5)
	feats_iqdn = cbind(deepbase_feats,
							iqbase_feats)

	mod$gw = list()
	mod$gw$cg_trace = cg_trace
	mod$gw$cgpcg_sc = cgpcg_sc
	mod$gw$cgtxg_sc = cgtxg_sc
	mod$gw$cgmix_sc = cgmix_sc
	mod$gw$feats = feats
	mod$gw$feats35 = feats35
	mod$gw$feats_iqdn = feats_iqdn

	for(i in 1:ncol(mod$gw$feats_iqdn)) {
		minv = min(mod$gw$feats_iqdn[,i], na.rm=T)
		maxv = max(mod$gw$feats_iqdn[,i], na.rm=T)
		f = is.na(mod$gw$feats_iqdn[,i])
		f1 = is.infinite(mod$gw$feats_iqdn[,i])
		if(sum(f)>0) {
			message("replace ", sum(f), " nas with ", minv, " in ", colnames(mod$gw$feats_iqdn)[i])
			mod$gw$feats_iqdn[f,i] = minv
		}
		if(sum(f1)>0) {
			mod$gw$feats_iqdn[f1,i] = maxv
		}
	}
	f = is.infinite(mod$gw$feats[,"GC"]) | is.na(mod$gw$feats[,"GC"])
	mod$gw$feats[f,"GC"] = 0.3
	
	return(mod)
}


pcg_build_gw_pred = function(mod, rebuild_iqdn=T,perc_genome, test_borzoi=FALSE,rebuild_base=F, rebuild_atac=F)
{
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
	
	
	#f_norp = cg_trace$d_ltr !=0 & cg_trace$d_line!=0  & cg_trace$start > 3e+6
	p = runif(n=nrow(gw$cg_trace))
	f_sub = f_norp & (p >= perc_genome | cg_trace$CG>0.02) 
#& cg_trace$cg_all!=1
    message('f_norp is ', sum(f_norp))
	ptrain = runif(n=nrow(gw$cg_trace))
	if (test_borzoi==FALSE){
	f_test_chr = cg_trace$chrom %in% c("chr1", "chr4", "chr6","chr10","chr13", "chr16", "chr19")
	}else{
	borzoi_chrs = c('chr14','chr10','chr15','chr4')
	f_test_chr = cg_trace$chrom %in% borzoi_chrs
	
	}
	message('test chromosomes ',unique(cg_trace$chrom[f_test_chr]))
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
		formu2iq = paste("resp ~ ", paste(colnames(feats_iqdn),collapse="+"))
		formu2iq = paste(formu2iq," + GC*CG",sep="")
		mod = pcg_build_mod_from_feats(mod, resp = cg_trace$lk27_1k, 
								formu = formu2iq, feats= feats_iqdn, 
								f_train = f_train, f_test = f_test, 
								f_train_all = f_train_all, f_test_all = f_test_all, tag = "base2iqdn")
	}

	if(rebuild_base) {
		formu = paste("resp ~ ", paste(colnames(feats),collapse="+"))
		formu = paste(formu," + GC*CG",sep="")
		mod = pcg_build_mod_from_feats(mod, resp = cg_trace$lk27_1k, 
								formu = formu, feats= feats, 
								f_train = f_train, f_test = f_test, 
								f_train_all = f_train_all, f_test_all = f_test_all, tag = "base")
	}

	if(rebuild_atac) {
		feats_atac = cbind(feats, gw$feats_iqdn, atac = log2(10+cg_trace$atac))
		formu_atac = paste("resp ~ ", paste(colnames(feats_atac),collapse="+"))
		formu_atac = paste(formu_atac, " + GC*CG",sep="")
		mod = pcg_build_mod_from_feats(mod, resp = cg_trace$lk27_1k, 
								formu = formu_atac, feats= feats_atac,
								f_train = f_train, f_test = f_test, 
								f_train_all = f_train_all, f_test_all = f_test_all, tag = "base_atac")
	}
	return(mod)
}

pcg_build_mod_from_feats = function(mod, resp, formu, feats, f_train, f_test, f_train_all, f_test_all, tag)
{
	set.seed(42)
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

# train linear and extract QC

	feats_r = as.data.frame(feats)
	feats_r$resp = resp
	lmod = lm(formu, feats_r[f_train,])
	
	lm_pred = predict(lmod, feats_r)
	lm_pred = pmax(lm_pred, 1)
	lm_stat = c(lm_r2_test = cor(resp[f_test], lm_pred[f_test])**2,
					lm_r2_train = cor(resp[f_train], lm_pred[f_train])**2,
					lm_r2_test_all = cor(resp[f_test_all], lm_pred[f_test_all])**2,
					lm_r2_train_all = cor(resp[f_train_all], lm_pred[f_train_all])**2,
					lm_auc_65_train = comp_auc(lm_pred, ifelse(resp > 6.5,1,0), f_train_all)$auc,
					lm_auc_65_test = comp_auc(lm_pred, ifelse(resp > 6.5,1,0), f_test_all)$auc,
					lm_auc_6_train = comp_auc(lm_pred, ifelse(resp > 6,1,0), f_train_all)$auc,
					lm_auc_6_test = comp_auc(lm_pred, ifelse(resp > 6,1,0), f_test_all)$auc,
					lm_auc_5_train = comp_auc(lm_pred, ifelse(resp > 5,1,0), f_train_all)$auc,
					lm_auc_5_test = comp_auc(lm_pred, ifelse(resp > 5,1,0), f_test_all)$auc)

	dtrain = xgb.DMatrix(as.matrix(feats)[f_train,], label = resp[f_train])
	dtest = xgb.DMatrix(as.matrix(feats)[f_test,], label = resp[f_test])

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
					xg_auc_65_train = comp_auc(xg_pred, ifelse(resp > 7.5,1,0), f_train_all)$auc,
					xg_auc_65_test = comp_auc(xg_pred, ifelse(resp > 7.5,1,0), f_test_all)$auc,
					xg_auc_6_train = comp_auc(xg_pred, ifelse(resp > 7,1,0), f_train_all)$auc,
					xg_auc_6_test = comp_auc(xg_pred, ifelse(resp > 7,1,0), f_test_all)$auc,
					xg_auc_5_train = comp_auc(xg_pred, ifelse(resp > 6.5,1,0), f_train_all)$auc,
					xg_auc_5_test = comp_auc(xg_pred, ifelse(resp > 6.5,1,0), f_test_all)$auc,
					xg_n_under4 = sum((xg_pred < 5 & resp > 6.5)[f_test_all | f_train_all]),
					xg_n_overr4 = sum((xg_pred > 6.5 & resp < 5)[f_test_all | f_train_all]))

	mod$gw$xg_mods[[tag]] = xgb_mod
	mod$gw$lm_mods[[tag]] = lmod
	mod$gw$xg_pred[[tag]] = xg_pred
	mod$gw$lm_pred[[tag]] = lm_pred
	mod$gw$mods_stat[[tag]] = c(formu=formu, lm_stat, xg_stat)
	return(mod)

}

peaks_eb = function(){
peaks_epi = readRDS('./data/peaks_atac_CRJK_0373_atac_wt_to_wt_epi_eb_d5.rds')
peaks_epi = peaks_epi$peaks
dim(peaks_epi)

peaks_meso = readRDS('./data/peaks_atac_CRJK_0377_atac_wt_to_wt_meso_eb_d5.rds')
peaks_meso = peaks_meso$peaks
dim(peaks_meso)

peaks_endo = readRDS('./data/peaks_atac_CRJK_0376_atac_wt_to_wt_endo_eb_d5.rds')
peaks_endo = peaks_endo$peaks
dim(peaks_endo)

peaks_eb = rbind(peaks_epi,peaks_meso) %>% rbind(.,peaks_endo)

peaks_eb_can = gintervals.canonic(peaks_eb)

peaks_eb_cen = peaks_eb_can

peaks_eb_cen$start = peaks_eb_can$start + (peaks_eb_can$end-peaks_eb_can$start)/2

peaks_eb_cen$end = peaks_eb_cen$start + 150

peaks_eb_cen$start = peaks_eb_cen$start -150
    
return(peaks_eb_cen)
    }
	
gen_atac_mat = function(tracks_atac,atac_peaks){


    for (i in 1:length(tracks_atac)) {
            tr = as.character(tracks_atac[i])
            nm = names(tracks_atac[i])
            gvtrack.create(nm, tr, 
                "sum")
            gvtrack.iterator(nm, sshift = -140, 
                eshift = 140)
        
    }

atac = atac_peaks
#atac = peaks_eb_cen
atac = atac[ !atac$chr %in% c('chrY','chrM','chrX'),] %>% arrange(chrom,start)

ge_a = gextract(names(tracks_atac),intervals = atac,iterator = 20)

ge_a[is.na(ge_a)] = 0

atac_mat = cbind(atac,tgs_matrix_tapply(t(ge_a[,names(tracks_atac)]),ge_a$intervalID,max))
return(atac_mat)
}

fig5_ribo_cov_atac = function(mod,tracks_atac){
ribo = mod$tss[ grepl('Rps|Rpl',mod$tss$geneSymbol),]

ribo$start = ribo$start-500
ribo$end = ribo$end+500


ribo_c = gintervals.canonic(ribo)



ribo_ge = gextract(names(tracks_atac) ,intervals = ribo_c,iterator = 20)

ribo_ge[is.na(ribo_ge)] = 0

ribo_cov = colSums(ribo_ge[,names(tracks_atac)])
return(ribo_cov)
}
fig5_track_cov = function(){
if (file.exists('./data/tracks_cov.rds')){
    tracks_cov = readRDS('./data/tracks_cov.rds')
    cov_vect = as.numeric(tracks_cov$cov_vect)
    names(cov_vect) = rownames(tracks_cov)
}else{
    cov_vect = c()
    for ( i in mod$epi_tracks$short_name){
        cov_vect = c(cov_vect,gsummary(i)[5])
    }



    names(cov_vect) = mod$epi_tracks$short_name
}
return(cov_vect)
}
fig5_gext_all_tracks = function(mod=mod){
cgd = mod$cgdom_ann
gext = gextract(names_k27, mod$cgdom_ann, 
        iterator = 20)
all_seed_k27 = as.data.frame((tgs_matrix_tapply(t(gext[,names_k27]),gext$intervalID,mean)))
all_seed_k27 = cbind(cgd[,c('chrom','start','end')],all_seed_k27)
all_seed_k27$pred = cgd$pred
all_seed_k27$pred_seed_k27 = cgd$pred_seed_k27
all_seed_k27$pred_seed_k4 = cgd$pred_seed_k4
return(all_seed_k27)
}

fig5_create_plt_mat = function(all_seed_k27,nm1,nm2,quant){

more = 5.5#5.5
delta =1.5
#    len = all_seed_k27$end - all_seed_k27$start
#    len_norm = 300/len
track1_r = as.numeric(all_seed_k27[,nm1])
track2_r = as.numeric(all_seed_k27[,nm2])
#q995_nm1 = as.numeric(quant[nm1,5])
#q995_nm2 = as.numeric(quant[nm2,5])
#ratio = as.numeric(cov_vect[nm1])/as.numeric(cov_vect[nm2])
#ratio = q995_nm1 - q995_nm2
#if ( ratio > 0 ){
#    norm_1 = 0
#    norm_2 = ratio
#}else {
#    norm_1 = abs(ratio)
#    norm_2 = 0    
#}
#message(norm_1)
#message(norm_2)
track1  = track1_r #+ norm_1
track2 = track2_r #+ norm_2
q98_nm1 = as.numeric(quant[nm1,3])
q98_nm2 = as.numeric(quant[nm2,3])

plt = data.frame(track1=track1,track2 = track2,pred=all_seed_k27$pred,
pred_seed_k27=all_seed_k27$pred_seed_k27,pred_seed_k4=all_seed_k27$pred_seed_k4)

cons = plt[,1] > q98_nm1 & plt[,2] > q98_nm2
track2_loss = plt[,1] > q98_nm1 & (plt[,1] - plt[,2])>delta

track1_loss = plt[,2] > q98_nm2 &  (plt[,2] - plt[,1])>delta

plt$col = ifelse(track2_loss, paste0(nm2,"_loss"),'else')
#plt$col = ifelse(cons, 'cons',plt$col)
plt$col = ifelse(track1_loss, paste0(nm1,"_loss"),plt$col)
df_tss = all_seed_k27
df_tss$col = plt$col
df_tss_neigh = gintervals.neighbors1(df_tss,mod$tss)
print(paste0(nm2,"_loss"))
df_tss_neigh[ df_tss_neigh$col==paste0(nm2,"_loss") & df_tss_neigh$dist==0,]

print(paste0(nm1,"_loss"))
df_tss_neigh[ df_tss_neigh$col==paste0(nm1,"_loss") & df_tss_neigh$dist==0,]

loss = df_tss_neigh[ df_tss_neigh$col==paste0(nm2,"_loss") & df_tss_neigh$dist==0,]

cons2 = plt[,1] > q98_nm1 & plt[,2] > q98_nm2
track2_loss = plt[,1] > q98_nm1 & (plt[,1] - plt[,2])>delta

#track1_loss = plt[,2] > 6 & plt[,1] < 5
plt$col2 = ifelse(cons2, 'cons','else')
plt$col2 = ifelse(track2_loss, paste0(nm2,"_loss"),plt$col2)

plt$col2 = ifelse(track1_loss, paste0(nm1,"_loss"),plt$col2)

plt$col2 = ifelse(plt$track1<5 & plt$track2<5, 'k4',plt$col2)
return(plt)
}

fig5_norm_atac = function(atac_mat,nm1_a,nm2_a,ribo_cov=ribo_cov){
track1_a = as.numeric(atac_mat[,nm1_a])
track2_a = as.numeric(atac_mat[,nm2_a])
#ratio_a = (as.numeric(cov_vect_atac[ nm1_a])/as.numeric(cov_vect_atac[ nm2_a]))
ratio_a = (as.numeric(ribo_cov[ nm1_a])/as.numeric(ribo_cov[ nm2_a]))
if ( ratio_a > 1 ){
    norm_1_a = 1
    norm_2_a = ratio_a
}else {
    norm_1_a = ratio_a
    norm_2_a = 1    
}
track1_a_n  = log2(4+(track1_a * norm_1_a))
track2_a_n = log2(4+(track2_a * norm_2_a))


plt_a = data.frame(track1_a=track1_a_n,track2_a = track2_a_n)

atac_normed = cbind(atac_mat[,c('chrom','start','end')],plt_a)
return(atac_normed)
    }