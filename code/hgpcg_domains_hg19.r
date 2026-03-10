library("misha")
library("dplyr")
library("zoo")
library("tglkmeans")
library("tgstat")
library("misha.ext")
library("data.table")

#source("scripts/hgpcg_cgdd_classify.r")
plt_cg_histone_trends_hg19 = function(mod)
{
    genome = gintervals.all()
#foc_chroms = genome[ genome$chrom %in% c('chr1','chr2','chr18','chr10','chr11'),]

foc_chroms = genome[ genome$chrom %in% c('chr1'),]
    tracks = c(
    'jk.henikoff_n19.h1_h3k27me3_1',
    #'jk.epipcg.pcg.CRJK_0403_k27me3_wt_to_wt_eb_d3_a70ls',
    'jk.henikoff_n19.h1_h3k4me3_1')
    short_nms = c('k27','k4')
	for(i in tracks) {
		gvtrack.create(gsub('jk.henikoff_n19.','',i), i, "sum") 
		gvtrack.iterator(gsub('jk.henikoff_n19.','',i), sshift = -140, eshift= 140)
	}

	gvtrack.create("CGh", "seq.CG", "sum")
	gvtrack.create("GCh", "seq.G_or_C", "sum")

	dst = gextract(gsub('jk.henikoff_n19.','',tracks), iterator=200, intervals=foc_chroms,colnames = c('k27','k4'))
	dst = dst[dst$start > 3e+6,] #remove prefix of NNNN
	lk4 = log2(2+dst$k4)
	lk27 = log2(2+dst$k27)
	stats = matrix(nrow=0,ncol=10)


	for(h in seq(0,2500,50)) { 
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
plot_cor_to_cg_hg19 = function(vars = c('c_g_k27','c_g_k4'),stats){

stats_gg = stats %>% as.data.frame()%>%rownames_to_column(var = 'scale')

stats_gg$scale = factor(levels = stats_gg$scale,stats_gg$scale)



options( repr.plot.width=15,repr.plot.height=10)
gg = stats_gg[,c('scale',vars)] %>% pivot_longer(-1) %>% ggplot(aes(x=scale, y = (value), group = name,col=name))+
coord_cartesian(ylim = c(0.25,.62))+scale_color_manual(values = c('darkblue','darkred'))+
geom_line(size=.7)+theme_bw()+theme(axis.text.x = element_text(angle = 90, hjust = 1), 
                                         panel.grid.major = element_blank(),
  panel.grid.minor = element_blank())
    print(gg)
    return(gg)
    }
set.seed(42)
hg19_finalize_hcgdds = function(mod){


mod$cgdom_ann$ID = 1:nrow(mod$cgdom_ann)

cgdd_iq = fread(here("output/iq-hg19-model_preds-diff-preds.tsv"), sep = "\t")%>%
arrange(chrom,start)
mod$cgdom_ann$ID = 1:nrow(mod$cgdom_ann)


mod$cgdom_ann$pred = as.numeric(cgdd_iq$pred)
mod$cgdom_ann$pred5mc = as.numeric(cgdd_iq$pred_5mc)

mod$seqmod_loc_pred = mod$cgdom_ann$pred

mod$seqmod_loc_pred5mc = mod$cgdom_ann$pred5mc
mod$cgdom_ann_orig = mod$cgdom_ann
cgd = mod$cgdom_ann



	hpcg_gen_k4_vt(mod)
	hpcg_gen_k27_vt(mod)
	
	hpcg_gen_atac_vt()
  



cgd = mod$cgdom_ann
cgd_data <- gextract(c("h1_1", 'h1_1_k4'), intervals = cgd, iterator = cgd) %>% arrange(intervalID) %>% 
    mutate(
        h1_cnt_k4 = h1_1_k4/(end - start) * 300, 
        h1_cnt = h1_1/(end - start) * 300
    )

cgd <- cgd %>% 
    mutate(
        h1_k4_sum = cgd_data$h1_cnt_k4,
        h1_k27_sum = cgd_data$h1_cnt,
        l10_h1_k4_sum = log2(10 + h1_k4_sum),
        l10_h1_k27_sum = log2(10 + h1_k27_sum)
    )

head(cgd)
mod$cgdom_ann = cgd
saveRDS(mod,'./data//files_hg19/chache_mod_hg19.rds')
return(mod)
    }
hpcg_init = function()
{
	options(gmax.data.size=1e+9)

	cache_file <- "./data/files_hg19/cache_mod_hg19.rds"

	if (file.exists(cache_file)) {
		mod <- readRDS(cache_file)
	} else {
	
	mod = list()

	mod$test_chrom =  c("chr2", "chr8", "chr12", "chr18")

	mod = hpcg_init_tss(mod)

	mod = hpcg_init_epi_track_lib(mod)
	#mod = hpcg_update_track_q_thresh(mod)
	mod = hpcg_gen_cg_h_doms(mod)
	mod = hg19_finalize_hcgdds(mod)
	return(mod)
	}
}

hpcg_init_tss = function(mod)
{
	mod$tss = gintervals.load("intervs.global.tss")
	return(mod)
}
plot_dense_scatter_ylim_legend_hg19=function (data, x, y, lim_d = 3.3, lim_u = 8.5) 
{
    options(repr.plot.width = 10, repr.plot.height = 10)
    my_pal <- colorRampPalette(c("darkgray", "blue3", "red", 
        "yellow"))
    temp = data[, c(x, y)]
    colnames(temp) = c("x", "y")
    kd = MASS::kde2d(temp$x, temp$y, n = 200)
    ix = findInterval(temp$x, kd$x)
    iy = findInterval(temp$y, kd$y)
    ix[ix < 1] = 1
    ix[ix > length(kd$x)] = length(kd$x)
    iy[iy < 1] = 1
    iy[iy > length(kd$y)] = length(kd$y)
    temp$density = kd$z[cbind(ix, iy)]
    gg = ggplot(temp, aes(x = x, y = y, color = density)) + geom_point(size = 0.8, 
        alpha = 1) + scale_color_gradientn(colours = my_pal(200), 
        name = "Density") + labs(x = x, y = y) + theme_bw() + 
        theme(axis.text.x = element_text(angle = -45, hjust = 0)) + 
        xlim(c(0, 8.5)) + ylim(c(0, 8.5))
    print(gg)
}
hpcg_init_epi_track_lib = function(mod)
{
	mod$epi_tracks = as.data.frame(fread("data/files_hg19/index_tracks.txt"))
	mod$epi_tracks$short_name_k4 = paste(mod$epi_tracks$short_name,"k4", sep="_")
	mod$epi_tracks_all = mod$epi_tracks
	mod$epi_tracks = mod$epi_tracks[mod$epi_tracks$use_comb==1,]
	mod$k4_track_thresh = read.table("data/files_hg19/track_thresh_k4.txt", 
        sep = "\t", stringsAsFactors = F, header = T)
	mod$k27_track_thresh = read.table("data/files_hg19/track_thresh.txt", 
        sep = "\t", stringsAsFactors = F, header = T)

	return(mod)
}
hpcg_update_track_q_thresh = function(mod)
{
	hpcg_gen_k27_vt(mod)
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
	hpcg_gen_k4_vt(mod)
	tnm_k4 = mod$epi_tracks$
	tnm_k4 =  mod$epi_tracks$short_name_k4
	tnm_k4 = tnm_k4[!is.na(tnm_k4)]
	for(tnm in tnm_k4) {
		thresh_k4[[tnm]] = gquantiles(tnm, c(0.5, 0.8, 0.9, 0.98, 0.985,0.99,0.995))
	}
	th_mat_k4 = do.call('rbind',thresh_k4)
	rownames(th_mat_k4) = tnm_k4
	write.table(th_mat_k4, file="data/track_thresh_k4.txt", quote=F, sep="\t")
	mod$k4_track_thresh = th_mat_k4
	return(mod)
}

hpcg_gen_k27_vt = function(mod, win_len=20) 
{
	wl2 = win_len/2
	ndx = mod$epi_tracks
	for(i in 1:nrow(ndx)) {
		if(!is.na(ndx$track_k27[i])) {
			gvtrack.create(ndx$short_name[i], ndx$track_k27[i], "sum") 
			gvtrack.iterator(ndx$short_name[i], sshift = -150+wl2, eshift= 150-wl2)
		}
	}
}

hpcg_gen_k4_vt = function(mod, win_len=20) 
{
	wl2 = win_len/2
	ndx = mod$epi_tracks
	for(i in 1:nrow(ndx)) {
		if(!is.na(ndx$track_k4[i])) {
			gvtrack.create(ndx$short_name_k4[i], ndx$track_k4[i], "sum") 
			gvtrack.iterator(ndx$short_name_k4[i], sshift = -150+wl2, eshift= 150-wl2)
		}
	}
}
hpcg_gen_atac_vt = function (tname = "marginal", w_ext = 140) 
{
    message("set vt for ", tname)
    nm = "meers_mol_cell_2019.H1_ATACseq_dense"
    gvtrack.create("atac_ext", nm, "sum")
    gvtrack.iterator("atac_ext", sshift = -w_ext, eshift = w_ext)
}


hpcg_gen_cg_h_doms = function(mod, force_update=F)
{
	if(!force_update & file.exists("data/files_hg19/cgdom_ann.RDS")) {
		mod$cgdom_ann = readRDS("data/files_hg19/cgdom_ann.RDS")
		mod$cgd_rpt = readRDS("data/files_hg19/cgd_rpt.RDS")
		return(mod)
	}

	hpcg_gen_k4_vt(mod)
	hpcg_gen_k27_vt(mod)
	
	hpcg_gen_atac_vt()
	cgd = gscreen("seq.CG_500_mean>0.04")
	cgd = cgd[!cgd$chrom %in% c("chrM", "chrY"),]
	cgd$l = cgd$end-cgd$start
	cgd = cgd[cgd$l>300,]

	rownames(cgd) = 1:nrow(cgd)
	
#annotate: max CG, max k27, k4, GC
	prof = gextract(c("h1_1", "jk.henikoff_n19.h1_h3k27me3_2reps", "h1_1_k4", "atac_ext","seq.CG_500_mean", 'seq.GC_500_mean'), iterator=20, intervals=cgd)
	prof[is.na(prof)]= 0
	cgd$es_k27_max = tapply(prof[,"h1_1"], prof$intervalID, max)
	cgd$es_k4_max = tapply(prof[,"h1_1_k4"], prof$intervalID, max)

	T_k27_995 = mod$k27_track_thresh["h1_1",7]
	T_k4_995 = mod$k4_track_thresh["h1_1_k4",7]
	T_k27_99 = mod$k27_track_thresh["h1_1",6]
	T_k4_99 = mod$k4_track_thresh["h1_1_k4",6]
	up_k27 = ifelse(prof[,"h1_1"]>T_k27_995,1,0)
	up_k4 = ifelse(prof[,"h1_1_k4"]>T_k4_995,1,0)
	cgd$es_k27_cov = tapply(up_k27 * (1-up_k4), prof$intervalID, mean)
	cgd$es_k4_cov = tapply((1-up_k27) * up_k4, prof$intervalID, mean)
	cgd$es_biv_cov = tapply(up_k27 * up_k4, prof$intervalID, mean)

	up_k27_1 = ifelse(prof[,"h1_1"]>T_k27_99,1,0)
	up_k4_1 = ifelse(prof[,"h1_1_k4"]>T_k4_99,1,0)
	cgd$es_k27_cov_1 = tapply(up_k27_1 * (1-up_k4_1), prof$intervalID, mean)
	cgd$es_k4_cov_1 = tapply((1-up_k27_1) * up_k4_1, prof$intervalID, mean)
	cgd$es_biv_cov_1 = tapply(up_k27_1 * up_k4_1, prof$intervalID, mean)

	cgd$atac_max = tapply(prof[,"atac_ext"], prof$intervalID, max)
	cgd$cg_max = tapply(prof[,"seq.CG_500_mean"], prof$intervalID, max)
	cgd$gc_max = tapply(prof[,"seq.GC_500_mean"], prof$intervalID, max)

	fcov = (prof$h1_1_k4+prof$h1_1)>20
	prof[!fcov,"h1_1_k4"] = NA
	prof[!fcov,"h1_1"] = NA
	cgd$es_k4_rat = tapply(log2((10+prof[,"h1_1_k4"])/(10+prof[,"h1_1"])), prof$intervalID, max, na.rm=T)
	cgd$es_k27_rat = tapply(log2((10+prof[,"h1_1"])/(10+prof[,"h1_1_k4"])), prof$intervalID, max, na.rm=T)

	cgd_tss = gintervals.neighbors(cgd, mod$tss)
	cgd$tss_dist = cgd_tss$dist
	cgd$tss_strand = cgd_tss$strand
	cgd$tss_start = cgd_tss$start
	cgd$tss_gene = cgd_tss$geneSymbol
	
	gvtrack.create("exon_d", "intervs.global.exon", "distance")
	gvtrack.create("line_d", "intervs.global.rmsk_line", "distance")
	gvtrack.create("ltr_d", "intervs.global.rmsk_ltr", "distance")
	gvtrack.create("simp_d", "intervs.global.rmsk_simple_repeat", "distance")
	gvtrack.create("sine_d","intervs.global.rmsk_sine", "distance")
	gvtrack.create("lowcomplex_d", "intervs.global.rmsk_low_complexity", "distance")
	gvtrack.create("blacklist_d", 'ENCODE.blacklist', "distance")
	gvtrack.create("bgc_d", 'intervs.global.bgc', "distance")

	cgd_rptdat = gextract(c("ifelse(exon_d==0,1,0)",
									"ifelse(line_d==0,1,0)", 
									"ifelse(ltr_d==0,1,0)", 
									"ifelse(sine_d==0,1,0)", 
									"ifelse(lowcomplex_d==0,1,0)", 
									"ifelse(!is.na(blacklist_d) & blacklist_d==0,1,0)", 
									"ifelse(!is.na(bgc_d) & bgc_d==0,1,0)", 
									"ifelse(simp_d==0,1,0)"), 
									iterator=20, intervals=cgd, 
									colnames=c("exon","line","ltr","sine","low_complex", "simple", "blacklist", "bgc"))
	cgd_rpts = data.frame(
						exon = tapply(cgd_rptdat$exon, cgd_rptdat$intervalID, mean),
						line = tapply(cgd_rptdat$line, cgd_rptdat$intervalID, mean),
						sine = tapply(cgd_rptdat$sine, cgd_rptdat$intervalID, mean),
						ltr = tapply(cgd_rptdat$ltr, cgd_rptdat$intervalID, mean),
						simp = tapply(cgd_rptdat$simp, cgd_rptdat$intervalID, mean),
						bgc = tapply(cgd_rptdat$bgc, cgd_rptdat$intervalID, mean),
						blacklist = tapply(cgd_rptdat$blacklist, cgd_rptdat$intervalID, mean),
						low_complex = tapply(cgd_rptdat$low_complex, cgd_rptdat$intervalID, mean))
	cgd_rpts$tot_rpt = rowSums(cgd_rpts[,c("line","sine","ltr","simp","low_complex")])
	
	cgd_mapab = gextract(c('mapab.length_50',
                         'mapab.wgEncodeDukeMapabilityUniqueness20bp',
                         'mapab.wgEncodeDukeMapabilityUniqueness35bp'
       ), iterator = cgd, intervals = cgd,
        colnames = c('mapab50', 'map_duke_20','map_duke_35'))

	cgd = cbind(cgd, cgd_mapab[,c("map_duke_20", "map_duke_35", "mapab50")])

	cgd = cbind(cgd, cgd_rpts)
	f_bad = cgd$line > 0.5 | cgd$ltr > 0.5 
	f_bad1 =  !f_bad & cgd$tot_rpt > 0.5 
	f_lowmapability = !f_bad & !f_bad1 & (cgd$es_k27_max+cgd$es_k4_max) == 0
	f_mask = f_bad | f_bad1 | cgd_mapab$map_duke_20 < 0.5
	#cgd$l = cgd$l+500

	mod$cgdom_ann = cgd[!f_mask,]
	mod$cgd_all = cgd
	
	saveRDS(mod$cgdom_ann, "data/files_hg19/cgdom_ann.RDS")
	saveRDS(mod$cgd_all, "data/files_hg19/cgd_rpt.RDS")

	return(mod)
}

hpcg_add_region_cg_data = function(mod)
{
	cgd = mod$cgdom_ann

	gvtrack.create("d_cgd",cgd[,1:3], "distance")
	gvtrack.create("d_cgdnotss",cgd[abs(cgd$tss_dist)>1000,1:3], "distance")

	all_c = gintervals.all()
	all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]
	cg_trace = gextract(c("ifelse(!is.na(seq.CG_500_mean) & seq.CG_500_mean > 0.04, 1, 0)"),
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
#regional CG elements
	cgd_cent = cgd[,c("chrom","start","end")]
	cgd_cent$start = (cgd_cent$start+cgd_cent$end)/2
	cgd_cent$end = cgd_cent$start + 1
	cgd_reg = gintervals.neighbors(cgd_cent, cgd_regstat)
	mod$cgdom_cg_regstat = cgd_reg
	return(mod)
}

old_scale50_code = function(mod)
{
	cgd = mod$cgdom_ann
	####seed clust
	all_c = gintervals.all()
	all_c = all_c[!all_c$chrom %in% c("chrM", "chrY"),]
	pall = gextract(c("h1_1"), iterator=500, intervals = all_c, colnames=c("k27"))
	T_k27 = mod$k27_track_thresh["h1_1","0.99"]
	pall$dom = ifelse(pall$k27 > T_k27, 1, 0)
	pall$dom50 = rollmean(pall$dom,100,f='e')	
	all_to50 = gintervals.neighbors(cgd, pall)	
	cgd$scale50k = all_to50$dom50

}

old_classes_code = function(mod)
{
	
	###divide into clusters
	d = mod$cgdom_ann
	d$type = ifelse(log2(1+d$es_k27_max) >= 6 & log2(1+d$es_k4_max)<6 ,'cl1','cl5')

	d$type = ifelse(log2(1+d$es_k27_max) >= 6 & log2(1+d$es_k4_max)>6 ,'cl2',d$type)

	d$type = ifelse(log2(1+d$es_k27_max) < 6 & log2(1+d$es_k27_max) > 4 & log2(1+d$es_k4_max)>5.5,'cl3',d$type)

	d$type = ifelse( log2(1+d$es_k27_max) < 4 & log2(1+d$es_k4_max)>5.5,'cl4',d$type)
	d$type_50k = ifelse(d$scale50k<.06,'seed','clust')
	mod$cgdom_ann = d	
	
	return(mod)
}
