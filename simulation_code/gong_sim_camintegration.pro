;+
; function gong_sim_camintegration, specfluxipv_in, specfluximv_in, dark, flat, waves, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe,
;		planckfac=planckfac, seed=seed, llo=llo, do_plot=do_plot
; 
; Given an input pair of datacubes (nx*ny*nw) for the I+V and I-V polarization states of GONG, produce 
; the six planes of the GONG raw data. There are three of these for each of the two polarization states,
; one for each of the 3 sets of binned similar halfwave plate angles (ori3={0,90,180,270}; ori3={30,120,210,300};
; and ori3={60,150,240,330} in gong_wave_resp). Takes into account wavelength response, dark & flat fields,
; detector response curve (currently linear response is assumed), and the variation in noise characteristics
; as a function of wavelength (and therefore photon energy; this effect is very subtle due to the narrow bandpass,
; and its use here is most likely overkill).
;
; [NOTES: (1) The wavelength response functions currently just use the instantaneous angle, but the wavelength
; response functions will need to be modified so that they are summed over the range of angles, since the
; cameras integrate continuously. (2) We will need to determine whether the nominal 0 in the GONG integration
; refers to the angle at start or the average angle during the integration. (3) Is there an overall phase
; difference between the ori3 definition used here and that for GONG? Does it matter? (4) In the longer term
; we may need to take into account the fact that the integration time intervals vary somewhat over the image
; as a result of the software 'shutter' employed.]
;
; Inputs:
;	specfluxipv_in: Array (nx*ny*nw) of spectral fluxes (e.g., in erg/s/cm) falling on each pixel for the I+V modulator state.
;	specfluxipv_in: Array (nx*ny*nw) of spectral fluxes (e.g., in erg/s/cm) falling on each pixel for the I-V modulator state.
;	dark: Array (nx*ny) giving the dark field for the image.
;	flat: Array (nx*ny) giving the flat field for the image.
;	waves: Wavelengths corresponding to the spectral fluxes.
;	exptime: The exposure time (default 1/60 second)
;   gaintabx: Input axis for the gain table (i.e., energy input in detector).
;   gaintaby: Output axis for the gain table (DNs for a given energy input).
;	rdn: Read noise of the detector.
;	qe: Quantum efficiency (default = 0.16)
;   planckfac: Constant factor which defines the conversion from photons to energy at a given wavelength.
;         That is, energy = n_phot*planckfac. Equal to h*c. All units must be consistent with this
;         default is 1.987821e-16, which assumes that all input units are centimeters, grams, and seconds.
;   seed: Seed for random number generation - uses IDL's randomu by way of poidev.pro
;	llo: An additional tunable parameter for the wavelength response (see gong_wave_resp.pro)
;	do_plot: If set, makes plots of the image as it's being integrated
;	plotfile: Name of file to use for output plot (otherwise will use whatever device settings are current)
;	nw_hi: Number of high-resolution wavelengths - datacubes are resampled to this resolution when multiplying against the wavelength
;		response functions. Set to 80 by default, should be a multiple of nw_lo.
;	nw_lo: Number of low-resolution wavelengths - datacubes are rebinned to this resolution after multiplying against the wavelength
;		response functions. This is done because the wavelength response functions require high spectral resolution, but the succeeding
;		computations do not. Set to 8 by default.
;	wmin: Minimum of wavelength range. Default: 676.6785e-7 cm
;	wmax: Maximum of wavelength range. Default: 676.8785e-7 cm
;	gongarea: GONG's effective area. Default: 0.35 cm^2 (2.8 cm aperture, ~3.5% throughput) 
;				(Changed 102617 - was 6%) Changed back for testing 102717
;-
function gong_sim_camintegration, specfluxipv_in, specfluximv_in, dark, flat, waves, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, planckfac=planckfac, seed=seed, llo=llo, do_plot=do_plot, plotfile=plotfile, ori3_offset=ori3_offset, rebin_wav=rebin_wav, gongarea=gongarea, nonoise=nonoise, npbin=npbin, dithermeans=dithermeans, nodisc=nodisc
; Energy in each pixel is sum(lambda_i*energyfac)


	; Find array dimensions:
	nx = n_elements(specfluxipv_in[*,0,0])
	ny = n_elements(specfluxipv_in[0,*,0])
	nw = n_elements(specfluxipv_in[0,0,*])

	
	if(n_elements(wmin) eq 0) then wmin = 676.25d-7
	if(n_elements(wmax) eq 0) then wmax = 677.25d-7
;	if(n_elements(wmin) eq 0) then wmin = double(676.6785d-7)
;	if(n_elements(wmax) eq 0) then wmax = double(676.8785d-7)
	if(n_elements(nw_lo) eq 0) then nw_lo = 2L
	if(n_elements(nw_hi) eq 0) then nw_hi = 1024L
	if(n_elements(gongarea) eq 0) then gongarea = 0.6d0
	waves_hi = wmin+(wmax-wmin)*dindgen(nw_hi)/(nw_hi-1.0d0)
;	waves_hi = waves
;	nw_hi = n_elements(waves)
	waves_lo = rebin(waves_hi,nw_lo)
	wcont = where((waves_hi gt max(waves)) or (waves_hi lt min(waves)))

	if(n_elements(ori1) eq 0) then ori1 = 0.0d0
	if(n_elements(ori2) eq 0) then ori2 = 0.0d0
	if(n_elements(ori3min) eq 0) then ori3min = 0.0d0
;	if(n_elements(ori3max) eq 0) then ori3max = 360*60*5.0d0
	if(n_elements(nori3) eq 0) then nori3 = 36000L ; Assumes 10 minute integration...
	if(n_elements(exptime) eq 0) then exptime = 1.0d0/60.0d0
	if(n_elements(ori3max) eq 0) then ori3max = (nori3-1.0d0)*exptime*360.0d0*5.0d0
	if(n_elements(ori3_offset) eq 0) then ori3_offset = 0.0d0
	
	; Note that this is a periodic function, so we don't go all the way to the max (which should be equivalent to the min).
	ori3s = ori3min+(ori3max-ori3min)*dindgen(nori3)/(nori3-1.0d0) - ori3_offset
	uori3s = ori3s[uniq(ori3s mod 360,sort(ori3s mod 360))] mod 360
	nuori3 = n_elements(uori3s)
	specfluxipv = dblarr(nx,ny,nuori3,nw_lo)
	specfluximv = dblarr(nx,ny,nuori3,nw_lo)
	for k=0,nuori3-1 do begin
		print,'Precomputing spectral fluxes for unique tuning number',k
		wresp = double(gong_wave_resp2(ori1,ori2,uori3s[k],waves_hi,llo=llo))
		for i=0,nx-1 do begin
			for j=0,ny-1 do begin
				ic_ipv = mean(specfluxipv_in[i,j,[0,1,nw-1,nw-2]])
				ic_imv = mean(specfluximv_in[i,j,[0,1,nw-1,nw-2]])
				sp_ipv = interpol(reform(double(specfluxipv_in[i,j,*])),double(waves-wmin),double(waves_hi-wmin),/nan)
				sp_imv = interpol(reform(double(specfluximv_in[i,j,*])),double(waves-wmin),double(waves_hi-wmin),/nan)
				sp_ipv[wcont] = ic_ipv
				sp_imv[wcont] = ic_imv
				specfluxipv[i,j,k,*] = gongarea*rebin(sp_ipv*wresp,nw_lo)
				specfluximv[i,j,k,*] = gongarea*rebin(sp_imv*wresp,nw_lo)
;				specfluxipv[i,j,k,*] = gongarea*rebin(reform(double(specfluxipv_in[i,j,*]))*wresp,nw_lo)
;				specfluximv[i,j,k,*] = gongarea*rebin(reform(double(specfluximv_in[i,j,*]))*wresp,nw_lo)
			endfor
		endfor
	endfor
				
;	message,'Breakpoint!'
				
	lcm=0
	wresps = fltarr(nori3,nw)
	images = dblarr(nx,ny,6)
;	specfluxipv = specfluxipv_in
;	specfluximv = specfluximv_in
	npbin = lonarr(6)
	lastimg_ipv = fltarr(nx,ny)
	lastimg_imv = fltarr(nx,ny)
	dither_offset = 32.0
	dithermeans = fltarr(6)
	for i=0,nori3-1 do begin
;;		if(((i mod 3600) ne 0) and ((i mod 12) eq 0) and ((i mod 3600) lt 1800)) then dither_offset = dither_offset + 2.0d0
;;		if(((i mod 3600) ne 0) and ((i mod 12) eq 0) and ((i mod 3600) gt 1800)) then dither_offset = dither_offset - 2.0d0

		if(((i mod 12) eq 0) and ((i mod 3600) lt 1800)) then dither_offset = dither_offset + 2.0d0
		if(((i mod 12) eq 0) and ((i mod 3600) ge 1800)) then dither_offset = dither_offset - 2.0d0

;		wresp = gong_wave_resp2(ori1,ori2,ori3s[i],waves,llo=llo)
;		for j=0,nw-1 do begin
;			specfluxipv[*,*,j] = specfluxipv_in[*,*,j]*wresp[j]
;			specfluximv[*,*,j] = specfluximv_in[*,*,j]*wresp[j]
;		endfor
		uori_index = where((ori3s[i] mod 360) eq uori3s)
		bin = floor(ori3s[i]/30) mod 3
		secs = floor(i*exptime)
		; 102617: Moved dither_offset to be before bit truncation, where it should be. Most straightforward to implement by tacking onto dark...
		if((secs mod 2) eq 0) then begin
;			img_ipv = gong_sim_detectorplane(reform(specfluxipv[*,*,uori_index,*]), dark+dither_offset/8.0d0, reform(flat[*,*,bin]), waves_lo, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, planckfac=planckfac, seed=seed, nonoise=nonoise)
			img_ipv = gong_sim_detectorplane(reform(specfluxipv[*,*,uori_index,*]), dark, reform(flat[*,*,bin]), waves_lo, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, planckfac=planckfac, seed=seed, nonoise=nonoise, nodisc=nodisc)+dither_offset
			images[*,*,bin] += img_ipv
			npbin[bin]++
			dithermeans[bin] += dither_offset
;			lastimg_ipv = img_ipv			
		endif else begin
;			img_imv = gong_sim_detectorplane(reform(specfluximv[*,*,uori_index,*]), dark+dither_offset/8.0d0, reform(flat[*,*,bin+3]), waves_lo, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, planckfac=planckfac, seed=seed, nonoise=nonoise)
			img_imv = gong_sim_detectorplane(reform(specfluximv[*,*,uori_index,*]), dark, reform(flat[*,*,bin+3]), waves_lo, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, planckfac=planckfac, seed=seed, nonoise=nonoise, nodisc=nodisc)+dither_offset
			images[*,*,bin+3] += img_imv
			npbin[bin+3]++
			dithermeans[bin+3] += dither_offset
;			lastimg_imv = img_imv
		endelse
;		img_ipv = gong_sim_detectorplane(specfluxipv, dark, reform(flat[*,*,bin]), waves_hi, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, planckfac=planckfac, seed=seed);+dither_offset
;		img_imv = gong_sim_detectorplane(specfluximv, dark, reform(flat[*,*,bin+3]), waves_hi, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, planckfac=planckfac, seed=seed);+dither_offset
		if(keyword_set(do_plot) and min(npbin) gt 0 and min(npbin) eq max(npbin)) then begin
			if(n_elements(plotfile) gt 0) then begin
				set_plot,'ps'
				device,/encapsulated,bits=8,/inches,xsize=8,ysize=4,filename=plotfile[0]
			endif
			!p.multi=[0,3,0]
			imgcor = images
			for j=0,5 do imgcor[*,*,j]-=dark*npbin[j]
			imgcor = imgcor/flat
			imgcor[where(flat lt 0.85)]=1.0
			imean = total(imgcor,3,/double)
			vmean = total(imgcor[*,*,3:5]-imgcor[*,*,0:2],3,/double)
			stokesi = imgcor[*,*,0:2]+imgcor[*,*,3:5]
			itanphi = reform(sqrt(3.0)*(stokesi[*,*,1]-stokesi[*,*,2])/(stokesi[*,*,1]+stokesi[*,*,2]-2.0*stokesi[*,*,0]))
			vdiff = imgcor[*,*,3]-imgcor[*,*,0] - (imgcor[*,*,5]-imgcor[*,*,2])
;			plot_image,images[*,*,0]
;			plot_image,images[*,*,1]
;			plot_image,images[*,*,2]
			plot_image,imean;,min=0.9*median(imean),max=1.2*median(imean)
			plot_image,itanphi,min=-1,max=1
			plot_image,vdiff;/imean,min=-0.025,max=0.025
			if(n_elements(plotfile) gt 0) then device,/close
			print,i,ori3s[i],bin,max(images[*,*,0]),min(images[*,*,0]),median(images[*,*,0]), dither_offset
			print,npbin
		endif
	endfor
	
	dithermeans = dithermeans/npbin
	
	print,'Done with integration; dither means = '+string(dithermeans)
	print,'npbin = '+string(npbin)
	
	if(n_elements(discfac) eq 0) then discfac = 1.0D0
	if(keyword_set(nodisc)) then begin
		return,images/discfac/8.0
	endif else begin
		return,ishft(ulong(images/discfac),-3)
	endelse
	
end