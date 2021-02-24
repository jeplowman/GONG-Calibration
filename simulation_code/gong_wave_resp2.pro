;+
; function gong_wave_resp2,ori1,ori2,ori3, waves
;
; Computes wavelength response for specified tuning of Lyot filters and
; Michelson. Inputs:
;
;   ori1: Tuning of second Lyot element in degrees (GONG fixes this to 0 or 90)
;   ori2: Tuning of third Lyot element in degrees (GONG fixes this to 0 or 90)
;   ori3: Tuning of Michelson in degrees
;
; Input/output:
;
;   waves: Wavelegth abscissa of response function, in cm (default is 676 to 677 *10^-7).
;
; Returns: wavelength response at specified tuning at each wavelength.
;
; Authors: Valentin M Pillet, Joseph Plowman
;
; 03/20/17
;-
function gong_wave_resp2,ori1,ori2,ori3, waves, llo=llo, dori3=dori3

; ori1 is  tuning of second lyot element in degrees 
; ori2 is  tuning of third lyot element in degrees 
; ori3 is  tuning of Michelson 
; 
; run like this for regular spectral line tunning
; IDL>for ind=0,180,10 do gong_fil_sim,0.,0.,float(ind)
; run like this for continuum line tunning
; IDL> for ind=0,180,10 do gong_fil_sim,90.,90.,float(ind)

	if(n_elements(tune_angle) ne 1) then tune_angle = 8.2d0

	if(n_elements(waves) ne 0) then begin
		lambda = double(waves)*1.d8 ; cm to Angstrom conversion.
	endif else begin	
		lambda=6762.5d0+(dindgen(8000)+1)/800.0d0
		waves = lambda*1.d-8 ; Angstrom to cm conversion.
	endelse
	lambdamu=lambda/1.d4

	if(n_elements(llo) eq 0) then llo=6767.781d0
	if(n_elements(dori3) eq 0) then dori3 = 30.0d0

	; prefilter
	fwhm=5.0d0 ; ansgtroms interf. filter
	;filter=1.0d0/(1.0d0+((lambda-llo)/(fwhm/2.0d0))^4.0d0) ; double cavity filter

	; Alternative prefilter profile from GONG instrument memo 86-4:
	c2 = 1.6569d0
	filter = 1.0d0/(1.0d0+c2*((lambda-llo)/(fwhm))^2.0d0)^2.0d0

	; Lyot:
	delth=608.d0 ; from GONG Report_7 (microns)

	elem1=(cos((!pi*delth/lambdamu)-67.7d0*!pi/180.0d0))^2
	elem2=(cos((!pi*2.d0*delth/lambdamu)+45.3d0*!pi/180.0d0-2.0d0*ori1*!pi/180.0d0))^2
	elem3=(cos((!pi*4.d0*delth/lambdamu)-93.0d0*!pi/180.0d0-2.0d0*ori2*!pi/180.0d0))^2

	lyot=elem1*elem2*elem3

	; Michelson (mysterious factor of 2 on ori3 is needed to match wiki description):
	fsr=0.305d0 ; Angstroms
	betathic=llo^2/fsr/1.d4
;	mic=(cos((!pi*betathic/lambdamu)-2.0d0*ori3*!pi/180.d0))^2.0d0
	arg = (!pi*betathic/lambdamu)+tune_angle*!pi/180.0d0-2.0d0*ori3*!pi/180.0d0
	twodelt = 2.0d0*dori3*!pi/180.0d0
	mic = (2.0d0*twodelt-(sin(2.0d0*(arg-twodelt))-sin(2.0d0*arg)))/(4.0d0*twodelt)
;	mic = (4.0d0*delt-(sin((arg-delt))-sin(arg)))/(8.0d0*delt)
	
	return,filter*lyot*mic
	
end

function gong_wr_wings,ori1,ori2,ori3, waves_in, llo=llo, dori3=dori3, wcen=wcen, wavrange=wavrange

		if(n_elements(wcen) eq 0) then wcen = 676.7781d-7
		if(n_elements(bandw) eq 0) then wavrange = 2.0d-7
		wmin = wcen-0.5d0*wavrange
		wmax = wcen+0.5d0*wavrange

		wmin_in = min(waves_in)
		wmax_in = max(waves_in)

		nw = 4000

		wing_resp = 0.0d0
		if(wmin_in gt wmin) then begin
			waves = wmin+(wmin_in-wmin)*dindgen(nw)/(nw-1.0d0)
			wr = gong_wave_resp2(ori1,ori2,ori3,waves,llo=llo,dori3=dori3)
			wing_resp += int_tabulated(waves,wr,/double)
		endif
		if(wmax_in lt wmax) then begin
			waves = wmax_in+(wmax-wmax_in)*dindgen(nw)/(nw-1.0d0)
			wr = gong_wave_resp2(ori1,ori2,ori3,waves,llo=llo,dori3=dori3)
			wing_resp += int_tabulated(waves,wr,/double)
		endif

		return,wing_resp

end

function gong_wave_resp2_0, ori1, ori2, ori3, waves, llo=llo, dori3=dori3

	wr = gong_wave_resp(ori1,ori2,ori3,waves,llo=llo)
	if(n_elements(dori3) eq 0) then dori3 = 30.0

	ndori3 = 100
	nw = n_elements(waves)

	dori3a = dori3*findgen(ndori3)/(ndori3-1.0)	
	
	for i=1,ndori3-1 do wr += gong_wave_resp(ori1,ori2,ori3+dori3a[i],waves,llo=llo)
	
	return,wr/ndori3
	
end
