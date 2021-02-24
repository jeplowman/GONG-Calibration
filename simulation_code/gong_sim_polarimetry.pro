;+
; gong_sim_polarimetry, specint_stokes, waves, mueller, fullout=fullout
;
; Applies a Mueller matrix to an input datacube. Inputs:
;
;   specint_stokes (nx*ny*nw*nstokes): The spectral intensities (e.g., in erg/s/cm/asec^2) as a
;   function of position and wavelength for each Stokes component.
;
;   waves: The wavelengths corresponding to specint_stokes.
;
;   mueller: The mueller matrix/matrices. Can be a function of position & wavelength
;   (dimensions nx*ny*nw*nstokes*nstokes), position (dimensions nx*ny*nstokes*nstokes), or
;   constant across the image (dimensions nstokes*nstokes).
;
; Optional inputs:
;   
;   fullout: Keyword; if set, Returns an output Stokes Vector for each pixel and wavelength, 
;   rather than summing over Stokes components (as is the case when counts are recorded on
;   a sensor).
;
; Returns: spectral intensity at each pixel and wavelength, summed over Stokes components
;    unless fullout keyword is set.
;
; Author: Joseph Plowman
; 03/22/17
;-
function gong_sim_polarimetry, specint_stokes, waves, mueller, fullout=fullout

	nx = n_elements(specint_stokes[*,0,0,0])
	ny = n_elements(specint_stokes[0,*,0,0])
	nw = n_elements(waves)
	nstokes = 4
;	specint_out = fltarr(nx,ny,nw,nstokes)
	specint_out = fltarr(nx,ny,nw)
	msize = size(mueller,/n_dimensions)
	
	print,'Running gong_sim_polarimetry; wavelength info:'+string([min(waves),max(waves),min(deriv(waves)),nw])
	
	; Multiply the Mueller matrices against the input spectral intensity Stokes vectors:
	if(msize eq 5) then begin
		; (Mueller matrix is fully specified for every wavelength & pixel):
		specint_out = specint_stokes*mueller
		for i=0,nstokes-1 do begin
			for j=0,nstokes-1 do begin
;				specint_out[*,*,*,i]+=specint_stokes[*,*,*,j]*mueller[*,*,*,i,j]
				specint_out+=specint_stokes[*,*,*,j]*mueller[*,*,*,i,j]
			endfor
		endfor
	endif else if (msize eq 4) then begin
		; (Mueller matrix is specified at every pixel):
		for i=0,nstokes-1 do begin
			for j=0,nstokes-1 do begin
;				for k=0,nw-1 do specint_out[*,*,k,i]+=specint_stokes[*,*,k,j]*mueller[*,*,i,j]
				for k=0,nw-1 do specint_out[*,*,k]+=specint_stokes[*,*,k,j]*mueller[*,*,i,j]
			endfor
		endfor
	endif else begin
		; (Mueller is constant across image)
		for i=0,nstokes-1 do begin
			for j=0,nstokes-1 do begin
;				specint_out[*,*,*,i] = specint_out[*,*,*,i]+specint_stokes[*,*,*,j]*mueller[i,j]
				specint_out += reform(specint_stokes[*,*,*,j]*mueller[i,j])
			endfor
		endfor
	endelse

;	; Sum over Stokes vectors:
;	if(keyword_set(fullout) eq 0) then specint_out = total(specint_out,4,/double)
	
	return,specint_out
	
end