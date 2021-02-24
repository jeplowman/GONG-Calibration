;+
; gong_sim_detectorplane, photons, dark, flat, gaintabx, gaintaby, rdn=rdn, qe=qe
; 
; Converts spectral flux incident on detector plane (nominally GONG's, but form
; is generic) to DN counts recorded by instrument, accounting for gain table,
; flat field and dark field, (wavelength dependent) shot noise. read noise, etc.
;
; Inputs:
;
;   specflux (nx*ny*nw): Spectral fluxes (e.g,. in erg/s/cm) onto each pixel as a function of wavelength.
;   dark: Dark field for the detector. Can be a scalar or an (nx*ny) image.
;   flat: Flat field for the detector. Can be a scalar or an (nx*ny) image.
;	waves (nw): The wavelengths in specflux, gaintabx, gaintaby, and qe (e.g., in cm).
;   exptime: The exposure time (e.g., in seconds).
;   gaintabx: Input axis for the gain table (i.e., energy input in detector).
;   gaintaby: Output axis for the gain table (DNs for a given energy input).
;
; Optional inputs:
;   rdn: The read noise of the detector, in DN units. Can be scalar or image (nx*ny). Default: 1.0
;   qe (nw): The quantum efficiency of the detector, as a function of wavelength. Default: 0.75
;   planckfac: Constant factor which defines the conversion from photons to energy at a given wavelength.
;         That is, energy = n_phot*planckfac. Equal to h*c. All units must be consistent with this
;         default is 1.987821e-16, which assumes that all input units are centimeters, grams, and seconds.
;   seed: Seed for random number generation - uses IDL's randomu by way of poidev.pro
;
; Returns: Simulated image (nx*ny) produced by the detector given the inputs.
;
; Author: Joseph Plowman
; 03/17/17
;-
function gong_sim_detectorplane, specflux, dark, flat, waves, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, planckfac=planckfac, seed=seed, adu_fac=adu_fac, nonoise=nonoise, nodisc=nodisc


	; Find array dimensions:
	nx = n_elements(specflux[*,0,0])
	ny = n_elements(specflux[0,*,0])
	nw = n_elements(specflux[0,0,*])
	if(n_elements(adu_fac) eq 0) then adu_fac = 33.0d0 ; 1 ADU is 33 photons/electrons, per Jack...

	
	; Default value for Planck factor, hc (assumes all inputs units are literally CGS, including wavelength):
	if(n_elements(planckfac) eq 0) then planckfac=1.987821d-16

	if(n_elements(gaintabx) eq 0 or n_elements(gaintaby) eq 0) then begin
		dnmax = 2.0^23.0-1 ; Assumes 16-bit ADC...
		gaintaby = dnmax*dindgen(10)/9.0d0
		gaintabx = gaintaby ; Assumes gain is linear...
		emax = dnmax
	endif else begin
		dnmax = max(gaintaby)
		emax = min(gaintabx[where(gaintaby eq dnmax)])
	endelse
		
	
	; Default values for read noise and quantum efficiency:
	if(n_elements(rdn) eq 0) then rdn = 0.6*8*adu_fac ; Based on standard deviation of a GONG average dark, plus the 3-bit shift
	if(n_elements(qe) ne nw) then qe = 0.16+fltarr(nw) ; Per email from Jack...

;	print,'Computing Photon Fluxes:'
;	tstart=systime(1)
	; Compute average photon fluxes:
;	dw = deriv(waves) ; The wavelength bins to convert from spectral fluxes to fluxes at each wavelength.
	dw = waves[1]-waves[0] + dblarr(nw) ; Assumes uniform bins in wavelength...

	photflux = specflux
	for i=0,nw-1 do photflux[*,*,i] *= dw[i]*qe[i]*waves[i]/planckfac
;	print,min(photflux)*exptime,max(total(photflux,3))*exptime,' '+strtrim(string(systime(1)-tstart),2)+' Seconds elapsed'
;	tstart=systime(1)
	
;	print,'Generating Poisson counts:'
	meancount = photflux*exptime > 0.0d0
	meancount = meancount < (emax*adu_fac+10.0*sqrt(emax*adu_fac)) ; sanitize flues according to the ADU limit
	for i=0,nw-1 do meancount[*,*,i] *=	flat
	if(keyword_set(nonoise)) then begin
		dn = meancount
		dn = (total(dn,3,/double) < emax) + dark*8.0d0*adu_fac
	endif else begin		
		; Generate Poisson photon counts and convert them to energies:
		dn = poidev(meancount > 0.0d0,seed=seed)
		; Apply flat and add read noise:
;		dn = (total(dn,3,/double) < emax) + randomn(seed,nx,ny)*rdn
		dn = (total(dn,3,/double) < emax) + randomn(seed,nx,ny)*rdn + poidev(dark*adu_fac*8.0d0 > 0.0d0, seed=seed)
	endelse
		
	; Apply `gain table' (i.e,. converts from incident energy to DN), bit shift ('compression'), add dark frame:
;	return, ishft(ulong(interpol(gaintaby,gaintabx,(dn+dark*8.0d0)>0.0d0)),-3)
;	return, interpol(gaintaby,gaintabx,(dn+dark*8.0d0)>0.0d0)/8.0D0
	if(keyword_set(nodisc)) then begin
		return, dn/adu_fac
	endif else begin
		return, round(dn/adu_fac)
	endelse
	
	; shouldn't darks be randomly Poisson distributed?
	
end

