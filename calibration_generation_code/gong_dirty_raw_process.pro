function gong_dirty_raw_process, rawfile, flatfile, darkfile, rawdat=rawdat, flatdat=flatdat, darkdat=darkdat, rawhdr=rawhdr, flathdr=flathdr, darkhdr=darkhdr, ipermi=ipermi, npbin=npbin, dithermean=dithermean, dopp=dopp, gongphi1=gongphi1, gongphi2=gongphi2, inten=inten, ndarkframes=ndarkframes, flatindex=flatindex, deconvolve=deconvolve, ondisk=ondisk, psf=psf, dc_niter=dc_niter, atmo_gaussian=atmo_gaussian, rebin_fac=rebin_fac, rebin_range=rebin_range

	if(n_elements(darkdat) eq 0) then darkdat = transpose(double(readfits(darkfile,darkheader)))
	if(n_elements(flatdat) eq 0) then flatdat = transpose(double(readfits(flatfile,flatheader)),[1,0,2])
	if(n_elements(rawdat) eq 0) then rawdat = transpose(double(readfits(rawfile,rawhdr)),[1,0,2])
	
;	rawdat += 108600L*600L/8L - 108600L/8L
	
	nx = n_elements(rawdat[*,0,0])
	ny = n_elements(rawdat[0,*,0])
	
	if(n_elements(flatindex) eq 0) then flatindex = 1
	if(n_elements(ndarkframes) eq 0) then begin
		ndarkframes = 1
		if(n_elements(darkheader) gt 0) then ndarkframes = sxpar(darkheader,'ACCUM')
	endif
	if(n_elements(npbin) eq 0) then npbin = 540; sxpar(rawhdr,'ACCUM')
	if(n_elements(npbin) eq 1) then npbin = npbin+fltarr(6)
	;print,'npbin=',strtrim(string(npbin),2)
	;print,'dark ACCUM='+strtrim(string(sxpar(darkheader,'ACCUM')),2)+', raw ACCUM='+strtrim(string(sxpar(rawhdr,'ACCUM')),2)

	if(n_elements(ipermi) eq 0) then ipermi = [1,0,2]
	
	; The value of the dither does not appear to be reflected in the GONG data - darks or otherwise,
	; so we assume it's already been removed as part of initial raw file generation. Otherwise, it
	; should be 182.0/8.0 per 1/60 sec frame...
	if(n_elements(dithermean) eq 0) then dithermean = 0.0; 182.0/8.0

	for i=0,5 do rawdat[*,*,i] -= darkdat*npbin[i]/ndarkframes+dithermean*npbin[i]
	flatdat = flatdat[*,*,flatindex]/median(flatdat[*,*,flatindex])
;	for i=0,5 do flatdat[*,*,i] -= darkdat
	for i=0,5 do rawdat[*,*,i] /= flatdat

	if(nx gt 100 and ny gt 100) then begin
		diskcen_rawval = median(rawdat[(round(0.5*nx)-50):(round(0.5*nx)+50),(round(0.5*ny)-50):(round(0.5*ny)+50),*])
		ondisk = label_region(max(rawdat gt 0.05*diskcen_rawval,dimension=3))
		ondisk = ondisk eq ondisk[round(0.5*nx),round(0.5*ny)]
	endif

	if(n_elements(dc_niter) eq 0) then dc_niter = 20
	if(keyword_set(deconvolve)) then rawdat = gong_sim_deconvolve_psf(rawdat, ondisk, psf=psf, niter=dc_niter, atmo_gaussian=atmo_gaussian)
	
	if(n_elements(rebin_fac) ne 0) then begin
		rawdat = rebin(rawdat[rebin_range[0]:rebin_range[1],rebin_range[2]:rebin_range[3],*],(rebin_range[1]-rebin_range[0]+1)/rebin_fac,(rebin_range[3]-rebin_range[2]+1)/rebin_fac,6)
	endif
	
	gongphi1 = atan(sqrt(3.0d0)*(rawdat[*,*,ipermi[1]]-rawdat[*,*,ipermi[2]]),(rawdat[*,*,ipermi[1]]+rawdat[*,*,ipermi[2]]-2.0d0*rawdat[*,*,ipermi[0]]))
	gongphi2 = atan(sqrt(3.0d0)*(rawdat[*,*,ipermi[1]+3]-rawdat[*,*,ipermi[2]+3]),(rawdat[*,*,ipermi[1]+3]+rawdat[*,*,ipermi[2]+3]-2.0d0*rawdat[*,*,ipermi[0]+3]))

	dopp1 = 2151.86d0*gongphi1
	dopp2 = 2151.86d0*gongphi2

	dopp = 0.5d0*(dopp2+dopp1)
	inten = total(rawdat,3)
	
	return, 0.5d0*(dopp2-dopp1)*0.704d0

end


function gong_deconvolve_magnetogram, bz_in, inten_in, dopp_in=dopp_in, ipermi=ipermi, modamp=modamp, ondisk=ondisk, dopp_out=dopp_out, inten=inten, dc_niter=dc_niter, phase0 = phase0

	if(n_elements(modamp) eq 0) then modamp = 0.85
	if(n_elements(dopp_in) eq 0) then dopp_in = 0.0
	if(n_elements(ipermi) eq 0) then ipermi = [1,0,2]
	if(n_elements(phase0) eq 0) then phase0 = -!pi/3.0
	
	bz=bz_in
	inten=inten_in
;	bz = readfits(bzfile,bzheader)
;	inten = readfits(intenfile,bzheader)
	
	nx = n_elements(inten[*,0])
	ny = n_elements(inten[0,*])
	
	if(n_elements(ondisk) eq 0) then ondisk = inten gt 0.1*median(inten[round(nx/2-100):round(nx/2+100),round(ny/2-100):round(ny/2+100)])
	
	bzphase = bz/0.704/2151.86
	phase_offset = dopp_in/2151.86
	
	ipvphase = phase0+phase_offset-bzphase
	imvphase = phase0+phase_offset+bzphase
	
	rawdat = dblarr(nx,ny,6)
	rawdat[*,*,0] = inten*(1.0+modamp*cos(ipvphase))
	rawdat[*,*,3] = inten*(1.0+modamp*cos(imvphase))
	rawdat[*,*,1] = inten*(1.0+modamp*cos(2.*!pi/3.0-ipvphase))
	rawdat[*,*,4] = inten*(1.0+modamp*cos(2.*!pi/3.0-imvphase))
	rawdat[*,*,2] = inten*(1.0+modamp*cos(4.*!pi/3.0-ipvphase))
	rawdat[*,*,5] = inten*(1.0+modamp*cos(4.*!pi/3.0-imvphase))

	if(n_elements(dc_niter) eq 0) then dc_niter = 20.0
	rawdat = gong_sim_deconvolve_psf(rawdat, ondisk, psf=psf, niter=dc_niter, atmo_gaussian=atmo_gaussian)
	
	gongphi1 = atan(sqrt(3.0)*(rawdat[*,*,ipermi[1]]-rawdat[*,*,ipermi[2]]),(rawdat[*,*,ipermi[1]]+rawdat[*,*,ipermi[2]]-2.0*rawdat[*,*,ipermi[0]]))
	gongphi2 = atan(sqrt(3.0)*(rawdat[*,*,ipermi[1]+3]-rawdat[*,*,ipermi[2]+3]),(rawdat[*,*,ipermi[1]+3]+rawdat[*,*,ipermi[2]+3]-2.0*rawdat[*,*,ipermi[0]+3]))

	dopp1 = 2151.86*gongphi1
	dopp2 = 2151.86*gongphi2

	dopp_out = 0.5*(dopp2+dopp1)
	inten = total(rawdat,3)
	
	return, 0.5*(dopp2-dopp1)*0.704

end