comp_capped_regional_atac = function(track_nm, q=0.95, win_local = 300, win_regional = 20000)
{
	gvtrack.create("cov_local", track_nm, "sum")
	local_ext = (win_local-20)/2
	local_n = win_local/20
	gvtrack.iterator("cov_local", sshift = -local_ext, eshift = local_ext)

	T_q = gquantiles("cov_local",q)

	cap_q = function(x) { x[is.na(x)] = 0; return(ifelse(x > T_q, T_q, x)/local_n) }
	gtrack.smooth(sprintf("%s_regcap", paste0(track_nm)), 
						desc= "capped regional avg", 
						expr = "cap_q(cov_local)",
						winsize = win_regional,
						iterator=20)

}


comp_peaks_over_punc_region = function(track_nm, norm_track_nm, peak_ext, eps=3, T_q = 0.98)
{
	peak_l = peak_ext*2+20
	peak_n = peak_l/20
	gvtrack.create("cov_loc", track_nm, "sum")
	gvtrack.iterator("cov_loc", sshift=-peak_ext, eshift=peak_ext)	
	norm_expr = sprintf("(%s+cov_loc)/(%s+peak_n*%s)", eps, eps, norm_track_nm)
	T_rat = gquantiles(norm_expr, T_q, iterator=20)
	message("ratio thresh ", T_rat)
	peaks = gscreen(sprintf("%s > %s", norm_expr, T_rat), iterator=20)

	return(peaks)
}

center_intervs_on_max = function(interv, track, width)
{
	a = gextract(c(track), iterator=20, intervals=interv, colnames=c("v"))
	x_cent = (interv$start+interv$end)/2
	d_cent = abs((a$start+a$end)/2-x_cent[a$intervalID])
	b = a[order(-a$v+0.01*d_cent),]
	b = b[!duplicated(b$intervalID),]
	b = b[order(b$intervalID),]
	cent_int = data.frame(chrom = b$chrom, start = b$start - width, end=b$end + width, v=b$v)
	
	return(cent_int)
}

center_and_union_intervs = function(interv, track, width, min_space)
{
	interv1 = center_intervs_on_max(interv, track, width)

	interv2 = interv1
	interv2$start = interv2$start - min_space/2	
	interv2$end = interv2$end + min_space/2
	interv3 = gintervals.canonic(interv2)
	interv3$start = interv3$start + min_space/2
	interv3$end = interv3$end - min_space/2

#center again (only affect unified intervals)
	interv4 = center_intervs_on_max(interv3, track, width)

	return(interv4)
}
add_peak_cggc = function(npeaks,tss) 
{	
	npeaks = npeaks
	gvtrack.create("GC500", "seq.G_or_C", "sum")
	all_cg = gextract(c("seq.CG_500_mean_new", "GC500"), iterator=npeaks,
								intervals=npeaks, colnames=c("CG", "GC"))
	all_cg$GC = all_cg$GC/(all_cg$end - all_cg$start)
	npeaks$CG = all_cg$CG
	npeaks$GC = all_cg$GC

	a = gintervals.neighbors(npeaks, tss[,c("chrom", "start", "end", "strand", "geneSymbol")])
	npeaks$d_tss = a$dist
	npeaks$gene_tss  = a$geneSymbol
	npeaks$gene_strand  = a$strand
	#mod$npeaks = npeaks
	return(npeaks)
#find the max energy per interval - realign
}
