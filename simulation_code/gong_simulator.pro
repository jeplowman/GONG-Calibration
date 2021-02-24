;+
; pro gong_simulator, datacube, wavelengths, x0, y0, dx0, dy0, dx, dy, nx, ny, tunings, psf=psf, wpsf=wpsf
; 
; Preliminary 'strawman' version of program to simulate the GONG magnetograph instruments. Given an input datacube,
; specifying Stokes vector components as a function of sky position and wavelength (units of power per unit
; square angle and wavelength), forms a datacube corresponding to the GONG observations. These datacubes contain
; both of the polarizations that GONG observes in a single array, with wavelengths of individual subimages being integrated
; over a range of angles of the rotating GONG half-wave plate. The individual exposures (taken at 60 Hz) impose a finer-grained
; integration over angles (the rotation frequency of the HWP is 5 Hz, so each is integrated over 1/12, or 30 degrees, of the 
; Half-wave plate). The polarization state is switched every second, so that every two seconds a complete set of wavelengths
; and polarization states is collected. These are then binned again, over 1 minute intervals, into 3 wavelength bins and the
; two polarization states, so that the rawest data recorded contains 3 frames - the first 3 are the wavelengths at one polarization
; state, the next 4 at the other. This is meant to replicate the data collection of the actual GONG instrument. The image is assumed
; to be static across the integration period.
;
; Inputs:
;    datacube: Array (dimensions [nx,ny,nwave,nstokes]; nstokes is always 4) giving the data at each Stokes vector and wavelength.
;    wavelengths: Vector (dimension nwave) giving the wavelengths in datacube.
;    x0: x center of datacube relative to origin of output image (should be negative to avoid clipping)
;    y0: y center of datacube relative to origin of output image (should be negative to avoid clipping)
;    dx0: x spacing of points in datacube image slices (nominally arcseconds)
;    dy0: y spacing of points in datacube image slices (nominally arcseconds)
;    dx: x spacing of points in output image (nominally arcseconds). Default: 2.5
;    dy: y spacing of points in output image (nominally arcseconds). Default: 2.5
;    nx: x size of output image. Default: 1024
;    ny: y size of output image. Default: 1024
;
; Optional Keyword Input/Output:
;    psf: Point spread function. If not assigned, will use a Gaussian with a width equal to wpsf.
; Optional Keyword Inputs:
;    wpsf: Width of psf. Ignored if psf is passed in.
;    mueller: Array (dimensions [2x4x4]) containing two mueller matrices for the two polarization states
function gong_simulator, datacube, wavelengths, x0, y0, dx0, dy0, dx, dy, nx, ny, psf=psf, wpsf=wpsf, mueller=mueller, dark=dark, flat=flat, seed=seed, llo=llo

	nw = n_elements(wavelengths)

	print,'Running GONG simulator: wmin='+string(min(wavelengths))+', wmax='+string(max(wavelengths))+', dw='+string(min(deriv(wavelengths)))+', nw='+string(nw)
	
	if(n_elements(mueller) eq 0) then begin
		mueller = dblarr(2,4,4)
		mueller[0,*,*] = diag_matrix([1,0,0,1]) ; I+V assumed to pass I & V unchanged, zero everything else.
		mueller[1,*,*] = diag_matrix([1,0,0,-1]) ; I-V assumed to pass I unchanged, flip sign of V, zero everything else.
	endif

	; First apply Mueller matrices, since it's a quick simplifying operation:
	specint_ipv = gong_sim_polarimetry(datacube,wavelengths,reform(mueller[0,*,*]))
	specint_imv = gong_sim_polarimetry(datacube,wavelengths,reform(mueller[1,*,*]))
	
	specflux_ipv = fltarr(nx,ny,nw)
	specflux_imv = fltarr(nx,ny,nw)
	; Now Apply PSFs and pixelization:
	for i=0,nw-1 do begin
		specflux_ipv[*,*,i] = gong_sim_pixel_psf(reform(specint_ipv[*,*,i]),x0,y0,dx0,dy0,dx,dy,nx,ny,psf=psf,wpsf=wpsf)
		specflux_imv[*,*,i] = gong_sim_pixel_psf(reform(specint_imv[*,*,i]),x0,y0,dx0,dy0,dx,dy,nx,ny,psf=psf,wpsf=wpsf)
	endfor

	; Detector plane / Cam integration happens last:
	cube = gong_sim_camintegration(specflux_ipv, specflux_imv, dark, flat, wavelengths, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, seed=seed, llo=llo)

	return, cube
	
end