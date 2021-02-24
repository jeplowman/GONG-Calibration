pro spectrum_restore, specfile, zcore, zcont, atmos, atmos_geometry, geometry, ray

	sobj=obj_new('IDL_Savefile',specfile)
	snames = sobj->Names()
	if(total(snames eq 'RAY')) then sobj->Restore,'ray'
	sobj->Restore,'zcore'
	sobj->Restore,'zcont'
	sobj->Restore,'atmos'
	if(total(snames eq 'ATMOS_GEOMETRY')) then sobj->Restore,'atmos_geometry'
	sobj->Restore,'geometry'
	sobj->cleanup

	;restore, specfile

end

pro get_calibration_curves, tau01, nopsf, zconst, use_taucore, zlvl, bzobs, bzmod, bzcal, qsflags, arflags, tile_levels, ondisk, qs_modfit, qs_obsfit, ar_modfit, ar_obsfit, nterms=nterms, ar_nterms=ar_nterms, xchang=xchang, fluxnorm=fluxnorm, arcoeffs=arcoeffs, qscoeffs=qscoeffs, binfit=binfit, modelfile=modelfile, simfile=simfile, darkfile=darkfile, flatfile=flatfile, gongitot=gongitot, deconvolve=deconvolve, theta_ray = theta_ray, blos_cube_inc=blos_cube_inc, blos_cube_vert=blos_cube_vert, z_fid=z_fid, zarr=zarr, ray=ray, use_contrib=use_contrib, qs_chi2=qs_chi2, ar_chi2=ar_chi2, bmax_qs=bmax_qs, bmax_ar=bmax_ar, phi_ray=phi_ray, geometry=geometry, atmos_geometry=atmos_geometry, gong_geo=gong_geo, bz_tauone=bz_tauone, dopp=dopp, post_deconvolve=post_deconvolve, tlvl0=tlvl0, dc_niter=dc_niter, nodeconv_bzobs=nodeconv_bzobs, zr0_inv=zr0_inv, sample=sample, atmo_w=atmo_w, use_pc_thold=use_pc_thold, wpsf=wpsf, spectrum_directory=spectrum_directory

	if(tau01 eq 0) then taustr = '1'
	if(tau01 eq 1) then taustr = '0.1'
	if(nopsf eq 0) then psfstr = 'Y'
	if(nopsf eq 1) then psfstr = 'N'
	if(zconst eq 1) then zconststr = 'Y'
	if(zconst eq 0) then zconststr = 'N'
	if(use_taucore eq 0) then taucorestr = 'N'
	if(use_taucore eq 1) then taucorestr = 'Y'
	
	if(n_elements(tlvl0) eq 0) then tlvl0 = 1
	if(n_elements(bmax_qs) eq 0) then begin
		bmax_qs = 400*sqrt(tlvl0)
		if(keyword_set(theta_ray)) then bmax_qs *= cos(theta_ray)
	endif
	if(n_elements(bmax_ar) eq 0) then bmax_ar = 4000

	if(n_elements(modelfile) eq 0) then modelfile = 'save/gongrh_ARQScube_111317.sav'
	if(n_elements(simfile) eq 0) then simfile = '/data/jplowman/GONG_sim_data_rh_ARQS_111417/terif170402/TE170402154716.fits.gz'

	if(n_elements(darkfile) eq 0) then darkfile = '/data/jplowman/GONG_data/tecqb170402/tecqb170402t1452/avgdrk.fits.gz'
	if(n_elements(flatfile) eq 0) then flatfile = '/data/jplowman/GONG_data/tecqb170402/tecqb170402t1452/flat.fits.gz'

	if(n_elements(dc_niter) eq 0) then dc_niter = 20
	
	restore,modelfile
;	if(savstr.specfile eq '/data/jplowman/rh_runs/run_rhsc3d_ARQS_60_combined/run_rhsc3d_ARQS_spectrum.sav') then savstr.specfile = '/data/jplowman/rh_runs/run_rhsc3d_ARQS_60_combined/run_rhsc3d_ARQS_spectrum.sav'
	specfile = savstr.specfile
	
	; Fix absolute file path reference in the saved spectrum file path. Actual runs were done before this change,
	; it is untested.
	if(n_elements(spectrum_directory) eq 1) then specfile = spectrum_directory+file_basename(specfile)
	spectrum_restore, specfile, zcore, zcont, atmos, atmos_geometry, geometry, ray
	;restore,specfile

	if(n_elements(atmos_geometry) eq 0) then atmos_geometry=geometry

	if(tag_exist(savstr,'npbin')) then begin
		npbin=savstr.npbin
	endif else begin
		npbin = 1200
	endelse
	
	if(tag_exist(savstr,'atmo_w')) then atmo_w = savstr.atmo_w
	if(tag_exist(savstr,'wpsf')) then wpsf = savstr.wpsf
	
	bzobs = gong_dirty_raw_process(simfile, flatfile, darkfile, rawhdr=rawhdr, npbin=npbin, dopp=dopp, inten=gongitot, ndarkframes=540.0, deconvolve=deconvolve, dc_niter = dc_niter)

	if(keyword_set(post_deconvolve)) then begin
		nodeconv_bzobs = bzobs
		bzobs = gong_deconvolve_magnetogram(bzobs, gongitot, ondisk=ondisk, inten=gongitot, dc_niter=dc_niter)
	endif
		
	if(keyword_set(theta_ray) gt 0) then begin
		if(n_elements(phi_ray) eq 0) then phi_ray = 90.0*!pi/180.0
		; Additional rotation of coordinate system  about x axis (for inclined rays):
		bangle_fac = (cos(phi_ray)*cos(atmos.chi_b) + sin(phi_ray)*sin(atmos.chi_b))*sin(theta_ray)*sin(atmos.gamma_b) + cos(theta_ray)*cos(atmos.gamma_b)
;		bangle_fac = (cos(phi_ray)*cos(atmos.chi_b) + sin(phi_ray)*sin(atmos.chi_b))*sin(theta_ray)*sin(atmos.gamma_b) + cos(theta_ray)*cos(atmos.gamma_b)
		thetay=theta_ray
	endif else begin
		bangle_fac = cos(atmos.gamma_b)
		;nw = n_elements(savstr.wavelengths)
		;ray = get_vertical_ray([1,round(0.5*nw)-1,nw-2])
	endelse

	if(n_elements(ray) eq 0) then ray = 'None'

	; Change this here to see if we're on an inclined ray and compute the z level from inclined tau if so. Might also consider
	; updating to use new equivalent code for the vertical tau. Could also update to use contribution functions.
	if(is_struct(ray)) then begin
		if(keyword_set(tau01)) then fcontrib = 0.1
		if(keyword_set(use_contrib)) then begin
			print,'Using Contribution function depth'
			if(keyword_set(tau01) eq 0) then fcontrib = 0.5
			depth = get_contrib_function_depth_reduced(ray, fcontrib, zray=z_fid)
		endif else begin
			print,'Using tau depth'
			if(keyword_set(tau01) eq 0) then fcontrib = 1.0			
			depth = get_contrib_function_depth_reduced(ray, fcontrib, zray=z_fid, /taudepth)
		endelse
		if(use_taucore) then z_fid = reform(z_fid[*,*,1])
		if(use_taucore eq 0) then z_fid = reform(mean(z_fid[*,*,[0,0]],dimension=3))
		z_fid = interpol(geometry.z,lindgen(geometry.nz),z_fid)
	endif else begin
		if(use_taucore eq 0) then z_fid=zcont*1.e3
		if(use_taucore eq 1) then z_fid=zcore*1.e3
		if(tau01 eq 1 and use_taucore eq 1) then z_fid = zcore01*1.e3
		if(tau01 eq 1 and use_taucore eq 0) then z_fid = zcont01*1.e3
	endelse
	if(zconst eq 1) then z_fid = median(z_fid)+0.0*z_fid
	if(n_elements(zlvl) eq 1) then z_fid=z_fid+zlvl*1.0e3

	gong_geo = gong_sim_image_geometry(rawhdr, savstr, geometry, arflags=arflags)
	zarr = atmos_geometry.z
	
	print,'BLOS factor range=',min(bangle_fac),max(bangle_fac)
	blos_cube_vert=atmos.b*bangle_fac
	
	bzmod = bz_tautile_samp(gong_geo, geometry, blos_cube_vert, z_fid, nopsf=nopsf, thetay=thetay, atmos_geometry=atmos_geometry, bz_cube=blos_cube_inc, ray=ray, bz_tauone=bz_tauone, zr0_inv=zr0_inv, sample=sample, atmo_w=atmo_w, wpsf=wpsf);, fcontrib=0.5)

;	bzmod = bz_tautile_samp(gong_geo, geometry, blos_cube_vert, z_fid, nopsf=nopsf, thetay=thetay, atmos_geometry=atmos_geometry, bz_cube=blos_cube_inc, ray=ray, bz_tauone=bz_tauone, zr0_inv=zr0_inv, fcontrib=0.5, atmo_w=atmo_w, wpsf=wpsf)
	
	ondisk = gong_geo.ondisk
	tile_levels = gong_geo.tile_levels

	if(n_elements(dlarge) eq 0) then dlarge = 30*tlvl0
	if(n_elements(dsmall) eq 0) then dsmall = 3*tlvl0
	;izi_ms_large = median(gongitot,dlarge)
	;izi_ms_small = median(gongitot,dsmall)
	izi_ms_large = median2(gongitot,dlarge)
	izi_ms_small = median2(gongitot,dsmall)
	;arflags = (gongitot lt 0.98*izi_ms_large) or (abs(bzobs) gt bmax_qs)
	if(keyword_set(use_pc_thold)) then arflags = get_ar_thold(gongitot,bzobs,bmax_qs=bmax_qs)
	qsflags = (arflags eq 0) and (abs(bzobs) lt bmax_qs) and ondisk
	arflags = (gongitot lt 0.98*izi_ms_large) and (qsflags eq 0) and ondisk

	bz_thold = 0.0
	it_thold = 0.92
	s_dilate = replicate(1,7,7)
	s_erode = replicate(1,3,3)
	gongi_qs = median(gongitot[where(ondisk and tile_levels eq tlvl0)])
	gongi_thold = dilate(erode(gongitot le it_thold*gongi_qs,s_erode),s_dilate)
	;qsflags = (gongi_thold eq 0)*(abs(bzobs) lt bmax_qs)
	qsflags *= ondisk
	;arflags = (gongi_thold eq 1)*(abs(bzobs) lt bmax_ar)
	arflags *= ondisk
	qspx = where(qsflags)
	arpx = where(arflags,arcount)
	gongres = ondisk*(tile_levels eq tlvl0)
	qspx_gongres = where(gongres*qsflags*(abs(bzobs) gt bz_thold),qscount)
	arpx_gongres = where(gongres*arflags,arcount)
	errs = sqrt(25^2+bzobs^2)

	help,bzobs
	help,qspx_gongres
	help,arpx_gongres

	if(keyword_set(binfit)) then begin
		if(n_elements(nterms) eq 0) then nterms = 50
		if(n_elements(ar_nterms) eq 0) then ar_nterms=nterms
		if(qscount gt 1) then begin
			qs_modfit = gong_curve_fluxcons_binfit(bzmod[qspx_gongres], bzobs[qspx_gongres], errs[qspx_gongres], qs_obsfit, unsigned=1, xchang=xchang, coeffs=qscoeffs, domean=1, npts=nterms)
		endif else begin
			qs_modfit = -2e4+1e4*findgen(nterms)/(nterms-1.0)
			qs_obsfit = qs_modfit
			qscoeffs = qs_modfit
		endelse
		if(arcount gt 1) then begin
			ar_modfit = gong_curve_fluxcons_binfit(bzmod[arpx_gongres], bzobs[arpx_gongres], errs[arpx_gongres], ar_obsfit, unsigned=1, xchang=xchang, coeffs=arcoeffs, domean=1, npts=ar_nterms)
		endif else begin
			ar_modfit = -2e4+1e4*findgen(nterms)/(nterms-1.0)
			ar_obsfit = ar_modfit
			arcoeffs = ar_modfit
		endelse
	endif else begin
		if(n_elements(nterms) eq 0) then nterms = 5	
		if(n_elements(ar_nterms) eq 0) then ar_nterms=nterms
		if(qscount gt 1) then begin
			qs_modfit = gong_curve_linfit(bzmod[qspx_gongres], bzobs[qspx_gongres], errs[qspx_gongres], qs_obsfit, unsigned=1, nterms=nterms, coeffs=qscoeffs, xchang=xchang, fluxnorm=fluxnorm)
		endif else begin
			qs_modfit = -2e4+1e4*findgen(nterms)/(nterms-1.0)
			qs_obsfit = qs_modfit
			qscoeffs = qs_modfit
		endelse
		if(arcount gt 1) then begin
			ar_modfit = gong_curve_linfit(bzmod[arpx_gongres], bzobs[arpx_gongres], errs[arpx_gongres], ar_obsfit, unsigned=1, nterms=ar_nterms, coeffs=arcoeffs, xchang=xchang, fluxnorm=fluxnorm)
		endif else begin
			ar_modfit = -2e4+1e4*findgen(nterms)/(nterms-1.0)
			ar_obsfit = ar_modfit
			arcoeffs = ar_modfit
		endelse
	endelse
	
	nx = n_elements(bzobs[*,0])
	ny = n_elements(bzobs[0,*])
	bzcal = fltarr(nx,ny)
	bzmodel_qs2 = spl_init(qs_obsfit,qs_modfit,/double)
	if(n_elements(qspx) gt 1) then bzcal[qspx] = spl_interp(qs_obsfit,qs_modfit,bzmodel_qs2,bzobs[qspx],/double)

	bzmodel_ar2 = spl_init(ar_obsfit,ar_modfit,/double)
	if(n_elements(arpx) gt 1) then bzcal[arpx] = spl_interp(ar_obsfit,ar_modfit,bzmodel_ar2,bzobs[arpx],/double)
	
	residuals = ((bzmod-bzcal)/errs)^2.0
	
	print,'Quiet sun Coefficients: ', qscoeffs
	print,'Active region Coefficients: ', arcoeffs

	qspx_gongres = where(gongres*qsflags)

	
end
