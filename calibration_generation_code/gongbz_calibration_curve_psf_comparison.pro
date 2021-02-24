function gongbz_calibration_curve_psf_comparison, tau01, zconst, use_taucore, zlvl, tablestrs, titlestr, fluxnorm=fluxnorm, reverse=reverse, nterms=nterms

	bzerr = 10.0

	if(tau01 eq 0) then taustr = '1'
	if(tau01 eq 1) then taustr = '0.1'
	if(zconst eq 1) then zconststr = 'Y'
	if(zconst eq 0) then zconststr = 'N'
	if(use_taucore eq 0) then taucorestr = 'N'
	if(use_taucore eq 1) then taucorestr = 'Y'

	zlvlstr = strtrim(string(round(zlvl)),2)


	bmodlabel = 'Model B!Dz!N at '
	if(zconst) then bmodlabel = bmodlabel+'median '
	bmodlabel = bmodlabel+'tau='+taustr+'+'+zlvlstr+'km'

	plotfilename='plots/cal_curves/GONG_calibration_curve_fit_'+titlestr+'_'
	if(zconst eq 1) then plotfilename = plotfilename+'med'
	plotfilename = plotfilename+'z'+zlvlstr+'_tau'+strjoin(strsplit(taustr,'.',/extract))
	if(use_taucore eq 0) then plotfilename = plotfilename+'_cont'
	plotfilename = plotfilename+'_PSF_comparison'

	set_plot,'z'
	plotfile_suffix = '.png'
	device,set_resolution=[2304,768]
	bgcolor = 0
	fgcolor = 255
	
	print,'tau='+taustr+' zconst:'+zconststr+' taucore:'+taucorestr+' zlvl='+zlvlstr+' fluxnorm:'+string(keyword_set(fluxnorm))+' reverse:'+string(keyword_set(reverse))+' nterms:'+strtrim(string(nterms),2)
	
	print,'With PSF:'	
	get_calibration_curves, tau01, 0, zconst, use_taucore, zlvl, bzobs0, bzmod0, bzcal0, qsflags0, arflags0, tile_levels0, ondisk0, qs_modfit, qs_obsfit, ar_modfit, ar_obsfit, reverse=reverse, fluxnorm=fluxnorm, arcoeffs=arcoeffs, qscoeffs=qscoeffs, nterms=nterms	
	
	print,'Without PSF:'
	get_calibration_curves, tau01, 1, zconst, use_taucore, zlvl, bzobs_nopsf0, bzmod_nopsf0, bzcal_nopsf0, qsflags_nopsf0, arflags_nopsf0, tile_levels0, ondisk0, qs_modfit_nopsf, qs_obsfit_nopsf, ar_modfit_nopsf, ar_obsfit_nopsf, reverse=reverse, fluxnorm=fluxnorm, arcoeffs=arcoeffs_nopsf, qscoeffs=qscoeffs_nopsf, nterms=nterms

	chi2eqstr = '$\chi^2 = \frac{1}{N}\sum_{i=1}{N} \frac{(B_{z,\textrm{NO PSF}}^i-B_{z,\textrm{Cal}}^i)^2}{\sigma^2}$'
	
	tablestrs = ['','\begin{frame}{PSF/Scale Reduced $\chi^2$ Comparison}',chi2eqstr,titlestr]
	tablestrs = [tablestrs,'{\footnotesize\begin{tabular}{ l | c | c | c | c | c | c }']
	tablestrs = [tablestrs,'Scale & Uncal AR & Uncal QS & AR Cal no PSF & QS Cal no PSF & AR cal w/PSF & QS Cal w/PSF \\']
	tablestrs = [tablestrs,'\hline']

	nxa=[1024,512,256]
	nya=[900,450,225]
	pmina=[1000,500,250]
	pmaxa=[3000,1500,750]
	scalestrs = ['Full','Half','Quarter']
	nscales = n_elements(scalestrs)
	
	for i=0,nscales-1 do begin
	
		nx=nxa[i]
		ny=nya[i]
		pmin=pmina[i]
		pmax=pmaxa[i]
		bzobs = rebin(bzobs0,nx,ny)
		bzmod = rebin(bzmod0,nx,ny)
		bzcal = rebin(bzcal0,nx,ny)
		bzobs_nopsf = rebin(bzobs_nopsf0,nx,ny)
		bzmod_nopsf = rebin(bzmod_nopsf0,nx,ny)
		bzcal_nopsf = rebin(bzcal_nopsf0,nx,ny)
		qsflags_nopsf = round(rebin(qsflags_nopsf0,nx,ny))
		arflags_nopsf = round(rebin(arflags_nopsf0,nx,ny))
		tile_levels = round(rebin(tile_levels0,nx,ny))
		ondisk = round(rebin(ondisk0,nx,ny))

		qspx_gongres = where(qsflags_nopsf*ondisk*(tile_levels eq 1))
		arpx_gongres = where(arflags_nopsf*ondisk*(tile_levels eq 1))
		
		!p.multi=[0,3,0]
		plot,[-pmin,pmax],[-pmin,pmax],xrange=[-pmin,pmax],yrange=[-pmin,pmax], xtitle = bmodlabel+' (NO PSF)', ytitle='Simulated GONG B!Dz!N',xstyle=1,ystyle=1,linestyle=1, title='Magnetogram Calibration w/o MHD Model PSF', charsize=2,charthick=2
		oplot,bzmod_nopsf[arpx_gongres],bzcal_nopsf[arpx_gongres],psym=3,color=128
		oplot,bzmod_nopsf[qspx_gongres],bzcal_nopsf[qspx_gongres],psym=3
		oplot,[-pmin,pmax],0.5*[-pmin,pmax],linestyle=3
		ssw_legend,['y=x line', 'y=x/2 line'], linestyle=[1,3], thick=[1,1], charsize=2,charthick=2

		plot,[-pmin,pmax],[-pmin,pmax],xrange=[-pmin,pmax],yrange=[-pmin,pmax], xtitle = bmodlabel+' (NO PSF)', ytitle='Simulated GONG B!Dz!N',xstyle=1,ystyle=1,linestyle=1, title='Magnetogram values before Calibration', charsize=2,charthick=2
		oplot,bzmod_nopsf[arpx_gongres],bzobs[arpx_gongres],psym=3,color=128
		oplot,bzmod_nopsf[qspx_gongres],bzobs[qspx_gongres],psym=3
		oplot,[-pmin,pmax],0.5*[-pmin,pmax],linestyle=3
		oplot,qs_modfit, qs_obsfit, thick=7, color = 255
		oplot,qs_modfit, qs_obsfit, thick=3, color = 0
		oplot,ar_modfit, ar_obsfit, thick=7, color = 255
		oplot,ar_modfit, ar_obsfit, thick=3, color = 0, linestyle=2
		ssw_legend,['y=x line', 'y=x/2 line','QS fit', 'AR fit'], linestyle=[1,3,0,0], thick=[1,1,7,7], charsize=2,charthick=2
		ssw_legend,['y=x line', 'y=x/2 line','QS fit', 'AR fit'], linestyle=[1,3,0,2], thick=[1,1,3,3], charsize=2,charthick=2

		plot,[-pmin,pmax],[-pmin,pmax],xrange=[-pmin,pmax],yrange=[-pmin,pmax], xtitle = bmodlabel+' (NO PSF)', ytitle='Simulated GONG B!Dz!N',xstyle=1,ystyle=1,linestyle=1, title='Magnetogram Calibration with MHD model PSF', charsize=2,charthick=2
		oplot,bzmod_nopsf[arpx_gongres],bzcal[arpx_gongres],psym=3,color=128
		oplot,bzmod_nopsf[qspx_gongres],bzcal[qspx_gongres],psym=3
		oplot,[-pmin,pmax],0.5*[-pmin,pmax],linestyle=3
		ssw_legend,['y=x line', 'y=x/2 line'], linestyle=[1,3], thick=[1,1], charsize=2,charthick=2
		
		chi2_ar_uncal = strtrim(string(total((bzmod_nopsf[arpx_gongres]-bzobs_nopsf[arpx_gongres])^2/bzerr^2,/double)/n_elements(arpx_gongres)),2)
		chi2_qs_uncal = strtrim(string(total((bzmod_nopsf[qspx_gongres]-bzobs_nopsf[qspx_gongres])^2/bzerr^2,/double)/n_elements(qspx_gongres)),2)
		chi2_ar_cal = strtrim(string(total((bzmod_nopsf[arpx_gongres]-bzcal[arpx_gongres])^2/bzerr^2,/double)/n_elements(arpx_gongres)),2)
		chi2_qs_cal = strtrim(string(total((bzmod_nopsf[qspx_gongres]-bzcal[qspx_gongres])^2/bzerr^2,/double)/n_elements(qspx_gongres)),2)
		chi2_ar_calnoPSF = strtrim(string(total((bzmod_nopsf[arpx_gongres]-bzcal_nopsf[arpx_gongres])^2/bzerr^2,/double)/n_elements(arpx_gongres)),2)
		chi2_qs_calnoPSF = strtrim(string(total((bzmod_nopsf[qspx_gongres]-bzcal_nopsf[qspx_gongres])^2/bzerr^2,/double)/n_elements(qspx_gongres)),2)
		
		tablestrs = [tablestrs,scalestrs[i]+' & '+chi2_ar_uncal+' & '+chi2_qs_uncal+' & '+chi2_ar_calnoPSF+' & '+chi2_qs_calnoPSF+' & '+chi2_ar_cal+' & '+chi2_qs_cal+' \\']
		write_png,plotfilename+'_'+scalestrs[i]+plotfile_suffix,tvrd()
		
	endfor
	
	tablestrs = [tablestrs,'\end{tabular}}\end{frame}']
	
	set_plot,'x'
	
	return, {bzmod:bzmod0, bzmod_nopsf:bzmod_nopsf0, bzobs:bzobs0, qs_modfit:qs_modfit, qs_obsfit:qs_obsfit, ar_modfit:ar_modfit, $
			ar_obsfit:ar_obsfit, qs_modfit_nopsf:qs_modfit_nopsf, qs_obsfit_nopsf:qs_obsfit_nopsf, ar_modfit_nopsf:ar_modfit_nopsf, $
			ar_obsfit_nopsf:ar_obsfit_nopsf, qsflags:qsflags0, arflags:arflags0, qsflags_nopsf:qsflags_nopsf0, $ arflags_nopsf:arflags_nopsf0, tile_levels:tile_levels0, ondisk:ondisk0, arcoeffs:arcoeffs, arcoeffs_nopsf:arcoeffs_nopsf, $
			qscoeffs:qscoeffs, qscoeffs_nopsf:qscoeffs_nopsf, ondisk:ondisk0, tile_levels:tile_levels0}

end