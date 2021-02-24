;+
; This function applies median smoothing to an image. Unlike the built-in IDL median function
; it uses a circular, not square, window, which avoids some of the artifacts that a
; retangular window can produce (because it is anisotropic). Edge treatment is equivalent
; to IDL's edge_truncate for convolution/smoothing. Inputs:
;
;	datain: The image to be median smoothed.
;	widthin: The width of the median smoothing window. This is forced to be an odd integer, which ensures
;			that the window is centered. 
;
; Returns a median smoothed image.
;
;-
function median2, datain, widthin
		
	nx = n_elements(datain[*,0])
	ny = n_elements(datain[0,*])
	
	pad = floor(widthin*0.5)
	width = pad*2+1
	
	data = dblarr(nx+2*pad,ny+2*pad)
	data[pad:-pad-1,pad:-pad-1]=datain

	for i=0,pad-1 do begin
		data[i,*] = data[pad,*]
		data[nx+2*pad-i-1,*] = data[-pad-1,*]
	endfor
	for i=0,pad-1 do begin &$
		data[*,i] = data[*,pad]
		data[*,ny+2*pad-i-1] = data[*,-pad-1]
	endfor
	
	xa = (lindgen(width))#transpose(1+lonarr(width))-pad
	ya = (1+lonarr(width))#transpose(lindgen(width))-pad
	rada = sqrt(xa*xa+ya*ya)
	xa_in = xa[where(floor(rada) le pad)]+pad
	ya_in = ya[where(floor(rada) le pad)]+pad
	dataout = dblarr(nx,ny)
	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			dataout[i,j] = median(data[xa_in+i,ya_in+j])
		endfor
	endfor

	return,dataout
	
end
