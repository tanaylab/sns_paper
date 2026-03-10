create_sns_conf_mm = function(){
conf_mm = list()
sns = list()
conf_mm$random_seed = 19
conf_mm$root = "mm10"
conf_mm$k27_track_nm = "jk.epipcg.pcg.CRJK_0364_k27me3_eb_j1_d4_a"
conf_mm$k4_track_nm = "jk.epipcg.pcg.CRJK_0411_k4me3_wt_to_wt_eb_d3"
conf_mm$atac_track_nm = "jk.epipcg.atac.CRJK_0389_atac_wt_to_wt_eb_d3"
conf_mm$cg_track_nm = "seq.CG_500_mean_new"
conf_mm$gc_track_nm = "seq.GC20_bin20/20"

conf_mm$tss_interv_nm = "intervs.global.tss"

conf_mm$line_interv_nm = "intervs.global.rmsk_line"
conf_mm$ltr_interv_nm = "intervs.global.rmsk_ltr"
conf_mm$sine_interv_nm = "intervs.global.rmsk_sine"
conf_mm$simple_repeat_interv_nm = "intervs.global.rmsk_simple_repeat"
conf_mm$low_complex_interv_nm = "intervs.global.rmsk_low_complexity"

conf_mm$test_chrom = c("chr4", "chr10", "chr14", "chr15")
conf_mm$remove_chrom = c("chrX")
sns$conf = conf_mm
return(sns)
    }
sns_gen_cnt_norm_track = function(out_tn, cnt_tn, atac_tn, 
													cnt_k_reg = NA, atac_k_reg = NA,
													win = 300,
													k_add = 5,
													line_interv_nm = "intervs.global.rmsk_line", 
													ltr_interv_nm = "intervs.global.rmsk_ltr")
{
	normalizer = sns_gen_cnt_normalizer(cnt_tn = cnt_tn, atac_tn = atac_tn, 
							cnt_k_reg = cnt_k_reg, 
							atac_k_reg = atac_k_reg,
							win = win,
							line_interv_nm = line_interv_nm, 
							ltr_interv_nm = ltr_interv_nm)

	gtrack.create(out_tn, 
		description=sprintf("%s normalized by atac %s, win = %s", cnt_tn, atac_tn, win),
		expr = "k_add + sns_norm_exp_simple(cnt, atac, normalizer, zero_cap = F)", iterator=20)

}

sns_gen_cnt_normalizer = function(cnt_tn, atac_tn, 
													cnt_k_reg = NA, atac_k_reg = NA,
													win = 300,
													line_interv_nm = "intervs.global.rmsk_line", 
													ltr_interv_nm = "intervs.global.rmsk_ltr")
{
	win_ext = (win-20)/2
	gvtrack.create("cnt", cnt_tn, "sum")
	gvtrack.iterator("cnt", sshift = -win_ext, eshif = win_ext)

	gvtrack.create("atac", atac_tn, "sum")
	gvtrack.iterator("atac", sshift = -win_ext, eshif = win_ext)

	gvtrack.create("line_d", line_interv_nm, "distance")
	gvtrack.create("ltr_d", ltr_interv_nm, "distance")

	prof = gextract(c("cnt","atac","line_d", "ltr_d"), iterator=20, 
							intervals=gintervals.all()[c(2,6,12),])

	prof$atac[is.na(prof$atac)] = 0
	prof$cnt[is.na(prof$cnt)] = 0
	prof = prof[prof$line_d != 0 & prof$ltr_d != 0,]

	if(is.na(cnt_k_reg)) {
		n1 = sum(prof$cnt == 1)
		n01 = sum(prof$cnt <= 1)
#cnt regularization is roughly the probabiltiy of 1 given count is smaller or equal to 1 - which represent some upper bound on the probabiltiy of 1 given 0. If there are no 0 count this is going to be just 0
		cnt_k_reg = n1/(1e-5+n01)
	}
	if(is.na(atac_k_reg)) {	
#we force atac regularization to be at least 1, no interest in smaller effects
		atac_k_reg = 1
	}

	l_cnt = log2(cnt_k_reg + prof$cnt)

	platac = log2(atac_k_reg + prof$atac)
	max_atac_bin = floor(4*quantile(platac, 1-1000/nrow(prof)))
	platac_bin = pmin(floor(4*platac),max_atac_bin)

	trend = tapply(l_cnt, platac_bin, mean, na.rm=T)

	normalizer = list()
	normalizer$atac_tn = atac_tn
	normalizer$cnt_tn = cnt_tn
	normalizer$trend = trend
	normalizer$max_atac_bin = max_atac_bin
	normalizer$min_atac_bin = min(as.numeric(names(trend)))
	normalizer$cnt_k_reg = cnt_k_reg
	normalizer$atac_k_reg = atac_k_reg
	normalizer$win = win
	return(normalizer)
	
}

sns_extract_norm_cnt = function(normalizer, intervs, geomean=F)
{
	cnt_tn = normalizer$cnt_tn
	atac_tn = normalizer$atac_tn

	win_ext = (normalizer$win-20)/2
	gvtrack.create("cnt", cnt_tn, "sum")
	gvtrack.iterator("cnt", sshift = -win_ext, eshif = win_ext)

	gvtrack.create("atac", atac_tn, "sum")
	gvtrack.iterator("atac", sshift = -win_ext, eshif = win_ext)

	prof = gextract(c("cnt", "atac"), iteravals = intervs, iterator=20)
		
	prof$atac[is.na(prof$atac)] = 0
	prof$cnt[is.na(prof$cnt)] = 0

	platac = log2(normalizer$atac_k_reg + prof$atac)
	max_atac_bin = normalizer$max_atac_bin
	min_atac_bin = normalizer$min_atac_bin
	platac_bin = as.character(pmax(pmin(floor(4*platac),max_atac_bin), min_atac_bin))
	
	prof_expected = normalizer$trend[platac_bin]

	if(geomean) {
		vals = log2(normalizer$cnt_k_reg + prof$cnt)
		norm_cnt = tapply(vals - prof_expected, prof$intervalID, mean)
	} else {
		tot_cnt = tapply(prof$cnt, prof$intervalID, mean)
		tot_pred = tapply((2**prof_expected)-normalizer$cnt_k_reg, prof$intervalID, mean)
		
		norm_cnt = log2(tot_cnt + normalizer$cnt_k_reg) - log2(tot_pred + normalizer$cnt_k_reg)
	}

	return(norm_cnt)	
}

sns_norm_exp_simple = function(cnt, atac, normalizer, zero_cap = T)
{
	platac = log2(normalizer$atac_k_reg + ifelse(is.na(atac), 0, atac))
	max_atac_bin = normalizer$max_atac_bin
	min_atac_bin = normalizer$min_atac_bin
	platac_bin = as.character(pmax(pmin(floor(4*platac),max_atac_bin), min_atac_bin))
	
	prof_log_expected = normalizer$trend[platac_bin]

	#norm_cnt = 2**(log2(cnt+ normalizer$cnt_k_reg) - prof_expected)
	norm_cnt = log2(cnt + normalizer$cnt_k_reg) - prof_log_expected
	if(zero_cap) {
		norm_cnt = pmax(norm_cnt,0)
	}

	return(norm_cnt)
}
