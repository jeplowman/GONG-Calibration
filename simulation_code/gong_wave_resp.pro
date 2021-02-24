;+
; function gong_wave_resp,ori1,ori2,ori3, waves
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
function gong_wave_resp,ori1,ori2,ori3, waves, llo=llo

; ori1 is  tuning of second lyot element in degrees 
; ori2 is  tuning of third lyot element in degrees 
; ori3 is  tuning of Michelson 
; 
; run like this for regular spectral line tunning
; IDL>for ind=0,180,10 do gong_fil_sim,0.,0.,float(ind)
; run like this for continuum line tunning
; IDL> for ind=0,180,10 do gong_fil_sim,90.,90.,float(ind)

	if(n_elements(waves) ne 0) then begin
		lambda = waves*1.e8 ; cm to Angstrom conversion.
	endif else begin	
		lambda=6762.5d0+(findgen(8000)+1)/800.
		waves = lambda*1.e-8 ; Angstrom to cm conversion.
	endelse
	lambdamu=lambda/1.d4

	if(n_elements(llo) eq 0) then llo=6767.781d0

	; prefilter
	fwhm=5.0d0 ; ansgtroms interf. filter
	filter=1./(1.+((lambda-llo)/(fwhm/2.))^4) ; double cavity filter

	; Lyot:
	delth=608.d0 ; from GONG Report_7 (microns)

	elem1=(cos((!pi*delth/lambdamu)-67.7*!pi/180.))^2
	elem2=(cos((!pi*2.d0*delth/lambdamu)+45.3*!pi/180.-2.*ori1*!pi/180.))^2
	elem3=(cos((!pi*4.d0*delth/lambdamu)-93.*!pi/180.-2.*ori2*!pi/180.))^2

	lyot=elem1*elem2*elem3

	; Michelson (mysterious factor of 2 on ori3 is need to match wiki description):
	fsr=0.305
	betathic=llo^2/fsr/1.e4
	mic=(cos((!pi*betathic/lambdamu)+8.2*!pi/180.-2.0*ori3*!pi/180.))^2
	
	return,filter*lyot*mic
	
end
