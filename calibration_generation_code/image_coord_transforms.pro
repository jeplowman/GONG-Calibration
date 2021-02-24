function rotmat_2d, theta

	return, [[cos(theta),-sin(theta)],[sin(theta),cos(theta)]]

end

; Convert from image pixel coordinates (x_pix_in, y_pix_in) to arcseconds (x_arcsecs, y_arcsecs), 
; given the following input parameters:
; 
;	center: The image center in arcseconds
;	scale: The pixel scale of the image
;	rollcen: The roll center of the image (in arcseconds)
;	rota: The rotation angle of the image about rollcen, in degrees.
;	nx: The number of pixels along the x direction (i.e., scale[0])
;	ny: The number of pixels along the y direction (i.e., scale[1])
pro pix_to_arcsecs, center, scale, rollcen, rota, nx, ny, x_pix_in, y_pix_in, x_arcsecs, y_arcsecs

	tilt = rota*2.0d0*acos(0.0d0)/180.0d0

	; Translate from corner coords:
	x_pix = x_pix_in - 0.5d0*(nx-1.0d0)
	y_pix = y_pix_in - 0.5d0*(ny-1.0d0)
	
	; Convert to arcseconds and translate from image center to roll center:
	x_arcsecs0 = x_pix*scale[0] + center[0] - rollcen[0]
	y_arcsecs0 = y_pix*scale[1] + center[1] - rollcen[1]
	
	; Apply roll rotation:
	x_arcsecs = x_arcsecs0*cos(tilt) - y_arcsecs0*sin(tilt)
	y_arcsecs = x_arcsecs0*sin(tilt) + y_arcsecs0*cos(tilt)

	; Translate back from roll center:
	x_arcsecs = x_arcsecs+rollcen[0]
	y_arcsecs = y_arcsecs+rollcen[1]

end

; Convert from arcseconds (x_arcsecs_in, y_arcsecs_in) to image pixel coordinates (x_pix, y_pix), 
; given the following input parameters:
; 
;	center: The image center in arseconds
;	scale: The pixel scale of the image
;	rollcen: The roll center of the image (in arcseconds)
;	rota: The rotation angle of the image about rollcen, in degrees.
;	nx: The number of pixels along the x direction (i.e., scale[0])
;	ny: The number of pixels along the y direction (i.e., scale[1])
pro arcsecs_to_pix, center, scale, rollcen, rota, nx, ny, x_arcsecs_in, y_arcsecs_in, x_pix, y_pix

	tilt = -rota*2.0d0*acos(0.0d0)/180.0d0

	; Translate to roll center:
	x_arcsecs = x_arcsecs_in-rollcen[0]
	y_arcsecs = y_arcsecs_in-rollcen[1]
	
	; Apply roll rotation:
	x_arcsecs0 = x_arcsecs*cos(tilt) - y_arcsecs*sin(tilt)
	y_arcsecs0 = x_arcsecs*sin(tilt) + y_arcsecs*cos(tilt)

	; Translate from roll center to image center and convert to pixels:
	x_pix = (x_arcsecs0-center[0]+rollcen[0])/scale[0]
	y_pix = (y_arcsecs0-center[1]+rollcen[1])/scale[1]
	
	; Translate to corner coords:
	x_pix += 0.5d0*(nx-1.0d0)
	y_pix += 0.5d0*(ny-1.0d0)
		
end

; px_refval must be defined if rollcen is not defined. If px_rollcen is undefined, 
; will assume it's the same as cen_px. If refval is undefined, assume it's at the origin of the (non-px)
; coordinate system.
function find_image_center, cen_px, scale, px_refval, refval=refval, rollcen=rollcen, px_rollcen=px_rollcen, rota=rota

	if(n_elements(rota) eq 0) then rota = 0.0d0
	if(n_elements(rollcen) eq 0 and n_elements(px_rollcen) eq 0 and rota mod 360.0d0 eq 0) then rollcen = [0.0d0,0.0d0]
	theta = rota*!pi/180.0d0
	rmat = rotmat_2d(theta)
	rmat_inv = rotmat_2d(-theta)

	; If one or the other of the rollcens are zero, then set them using the refvals:
	if(n_elements(px_rollcen) eq 0 or n_elements(rollcen) eq 0) then begin
		; Assume reference position is origin if it and at least one roll center is unassigned:
		if(n_elements(refval) eq 0) then refval = [0.0d0,0.0d0]
		; If BOTH roll center coordinates are undefined, have to assume rollcen is at image center:
		if(n_elements(px_rollcen) eq 0 and (n_elements(rollcen) eq 0 or n_elements(px_refval) eq 0)) then px_rollcen=cen_px
		if(n_elements(rollcen) eq 0) then rollcen = refval-scale*(rmat_inv#(px_refval-px_rollcen))
		if(n_elements(px_rollcen) eq 0) then px_rollcen = px_refval - (1.0d0/scale)*(rmat#(refval-rollcen))
	endif  
	
	; Both rollcen and px_rollcen must be initialized if px_refval is not:
	if(n_elements(px_refval) eq 0) then begin
		if(n_elements(refval) eq 0) then refval = [0.0d0,0.0d0]
		px_refval = (1.0d0/scale)*(rmat#(refval-rollcen))+px_rollcen
	endif
	if(n_elements(refval) eq 0) then begin
		refval = scale*(rmat_inv#(px_refval-px_rollcen))+rollcen
	endif

	return, scale*(rmat_inv#(cen_px-px_rollcen))+rollcen
	
end
	
;+
; Resample an image whose orientation and scale are specified by the structure index. The orientation and scale
; of the new image are defined in terms of the following parameters:
;	
;	center: The image center for the new image ([x,y])
;	scale: The pixel scale for the new image ([dx,dy])
;	nx: The number of pixels in the x direction
;	ny: The number of pixels in the y direction
;	rollcen: The roll center for the image rotation ([x,y]) (optional keyword)
;	rota: The roll angle for the image rotation, in degrees (optional keyword)
;
; The index structure must also specify these parameters, via elements xc,yc (for center), dx, dy (for scale)
; roll_center (for rollcen), and roll_angle (for rota). Any resample parameters not specified on input
; are set to be the same as in the index structure. The resampling is done using IDL's 'interpolate'
; function with cubic=-0.5
; 
; An output index containing the new parameters defined above, along the exposure time (but nothing else)
; will be returned in the optional keyword indexout.
;-
function resample_image, image, index, center, scale, nx, ny, rollcen=rollcen, rota=rota, smooth_radius = smooth_radius, indexout=indexout, missing=missing, sample=sample 

	; Orientation, scale, and dimensions of input image:
	nx0 = n_elements(image[*,0])
	ny0 = n_elements(image[0,*])
	if(where(tag_names(index) eq 'XC') eq -1) then begin
		scale0 = [index.cdelt1,index.cdelt2]
		center0 = find_image_center(0.5d0*[nx0,ny0],scale0,[index.crpix1,index.crpix2],refval=[index.crval1,index.crval2])
		pix_to_arcsecs, center0, scale0, [0.0d0,0.0d0], 0.0d0, nx0, ny0, index.crpix1, index.crpix2, rollcenx, rollceny
		rollcen0 = [rollcenx,rollceny]
		rota0 = index.crota2
	endif else begin		
		center0 = [index.xc,index.yc]
		scale0 = [index.dx, index.dy]
		rollcen0 = index.roll_center
		rota0 = index.roll_angle
	endelse
	
	; Any resample parameters not set on input are assigned to be the same as the input image:
	if(n_elements(center) eq 0) then center = center0
	if(n_elements(scale) eq 0) then scale = scale0
	if(n_elements(nx) eq 0) then nx = nx0
	if(n_elements(ny) eq 0) then ny = ny0
	if(n_elements(rollcen) eq 0) then rollcen = rollcen0
	if(n_elements(rota) eq 0) then rota = rota0

	; Arrays containing pixel coordinates of the new image:
	xpix = dblarr(nx,ny)
	ypix = dblarr(nx,ny)

	for i=0,ny-1 do xpix[*,i] = dindgen(nx)
	for i=0,nx-1 do ypix[i,*] = dindgen(ny)

	; Compute arcsecond coordinates of each pixel in new image:
	pix_to_arcsecs, center, scale, rollcen, rota, nx, ny, xpix, ypix, xa, ya
	
	; Convert these to the pixel coordinates WRT the original image
	arcsecs_to_pix, center0, scale0, rollcen0, rota0, nx0, ny0, xa, ya, x_pix0, y_pix0
	
	; Create the output index structure if it has been specified:
	if(n_elements(indexout) gt 0) then begin
		exptimes = interpolate(index.exptimes,x_pix0,y_pix0,cubic=-0.5d0)
		indexout = {nx:nx, ny:ny, xc:center[0], yc:center[1], dx:scale[0], dy:scale[1], roll_center:rollcen, roll_angle:rota, exptimes:exptimes, exptime0:index.exptime0}
	endif
	
	if(keyword_set(sample)) then begin
		image_out = image[floor(x_pix0)<(nx0-1),floor(y_pix0)<(ny0-1)]
	endif else begin
		image_out = interpolate(image,x_pix0,y_pix0,missing=missing,/double)
	endelse

	; Create and return the resampled image:
	return, image_out

end
