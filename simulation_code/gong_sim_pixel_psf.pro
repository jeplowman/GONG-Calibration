function gong_sim_pixelize, image_in, origin0, dx0, dy0, origin, dx, dy, nx, ny, oversamp_fac=oversamp_fac, rota0=rota0, rota2=rota2, xc=xc, yc=yc, sample=sample, average=average, missing=missing

	nx0 = n_elements(image_in[*,0])
	ny0 = n_elements(image_in[0,*])
	if(n_elements(rota0) eq 0) then rota0 = 0.0d0
;	if(n_elements(rota2) eq 0) then rota2 = 0.0
	
	if(n_elements(oversamp_fac) eq 0) then oversamp_fac = 5.0d0
	nx2 = nx*oversamp_fac
	ny2 = ny*oversamp_fac
	dx2 = dx/oversamp_fac
	dy2 = dy/oversamp_fac
	
	x0 = origin0[0]
	y0 = origin0[1]
	
	index = {crval1:x0,crval2:y0,crpix1:0,crpix2:0,cdelt1:dx0,cdelt2:dy0,crota2:rota0}
;	index = {xc:x0,yc:y0,dx:dx0,dy:dy0,roll_center:[0.0,0.0],roll_angle:0.0}

	center = find_image_center(0.5d0*[nx2,ny2], [dx2,dy2], [0.0d0,0.0d0], refval=origin, rollcen=rollcen, rota=rota2)
	
	image_out2 = resample_image(image_in, index, center, [dx2,dy2], nx2, ny2, rota=rota2, missing=missing, sample=sample)
	
	if(keyword_set(average) eq 0) then image_out2 *= dx*dy
	
	return, rebin(double(image_out2), nx, ny, sample=sample)
	
end


function gong_gaussian_psf, dx0, dy0, wpsf=wpsf, nxhalf=nxhalf, nyhalf=nyhalf

	if(n_elements(wpsf) eq 0) then wpsf = 2.0D ; A Gaussian width of 2.0 arcseconds is assumed unless otherwise specified
	
	if(n_elements(nxhalf) eq 0) then nxhalf = 7+2*ceil(wpsf/dx0)
	if(n_elements(nyhalf) eq 0) then nyhalf = 7+2*ceil(wpsf/dy0)
	nx = 2*nxhalf+1
	ny = 2*nyhalf+1
	print,nxhalf,nyhalf,nx,ny
	
	xa = dx0*(dindgen(nx)-nxhalf)
	ya = dy0*(dindgen(ny)-nyhalf)
	xa2 = xa#(1.0d0+dblarr(ny))
	ya2 = (1.0d0+dblarr(nx))#ya
	
	ra = sqrt(xa2*xa2+ya2*ya2)
	
	gausspsf = dx0*dy0*exp(-0.5D0*(ra/(0.5D0*wpsf))^2.0D0)/sqrt(!pi*2.0D0*(0.5D0*wpsf)^2.0D0)^2.0D0
	
	return, gausspsf/total(gausspsf,/double)
	
end


function gong_airy_psf, dx0, dy0, lambda0=lambda0, ap=ap, atmo_w=atmo_w, wpsf=wpsf, nxhalf=nxhalf, nyhalf=nyhalf

	arcsec_conv = 3600.0d0*180.0d0/!pi
	rayleigh_fac = 1.22d0
	if(n_elements(wpsf) eq 0 and n_elements(ap) eq 0) then wpsf = 6.1d0
	if(n_elements(lambda0) eq 0) then lambda0 = 676.7785d-7
	; If the aperture is not specified, it's assumed consistent with 5.0 arcsecond resolution per Rayleigh criterion:
	if(n_elements(ap) eq 0) then ap = arcsec_conv*rayleigh_fac*lambda0/(wpsf)
	
	wpsf = arcsec_conv*rayleigh_fac*lambda0/ap
	if(keyword_set(atmo_w)) then begin
		wpsf = sqrt(wpsf*wpsf+atmo_w*atmo_w)
		ap = arcsec_conv*rayleigh_fac*lambda0/wpsf
	endif
	
	if(n_elements(nxhalf) eq 0) then nxhalf = 5+2*ceil(wpsf/dx0)
	if(n_elements(nyhalf) eq 0) then nyhalf = 5+2*ceil(wpsf/dy0)
	nx = 2*nxhalf+1
	ny = 2*nyhalf+1
	
	xa = dx0*(dindgen(nx)-nxhalf)
	ya = dy0*(dindgen(ny)-nyhalf)
	xa2 = xa#(1.0d0+dblarr(ny))
	ya2 = (1.0d0+dblarr(nx))#ya	
	ra = sqrt(xa2*xa2+ya2*ya2)*!pi/180.0d0/3600.0d0
	
	kafunc = (!pi*ap/lambda0)*sin(ra)
	airyfunc = (2.0d0*beselj(kafunc,1)/kafunc)^2.0d0
	airyfunc[where(ra eq 0)] = 1.0d0
		
	return, airyfunc/total(airyfunc)
	
end

function gong_sim_apply_psf, image_in, dx0, dy0, psf, wpsf=wpsf, atmo_w=atmo_w, oversamp_fac=oversamp_fac, nopsf=nopsf, sample=sample, nxhalf=nxhalf, nyhalf=nyhalf

	nx0 = n_elements(image_in[*,0])
	ny0 = n_elements(image_in[0,*])
	if(n_elements(oversamp_fac) eq 0) then oversamp_fac = 1.0d0
	dx = dx0/oversamp_fac
	dy = dy0/oversamp_fac
	nx = nx0*oversamp_fac
	ny = ny0*oversamp_fac
	image = rebin(image_in,nx,ny,/sample)

	;print,n_elements(psf),keyword_set(atmo_w),keyword_set(nopsf)

	help,image

	if(keyword_set(nopsf) eq 0) then begin
		if(n_elements(psf) eq 0) then psf = gong_airy_psf(dx,dy, wpsf=wpsf, nxhalf=nxhalf, nyhalf=nyhalf);, atmo_w=atmo_w)
		help,psf
		if(keyword_set(atmo_w)) then begin
			print,'Applying '+strtrim(string(wpsf),2)+' arcsecond PSF with '+strtrim(string(atmo_w),2)+' arcsecond atmospheric seeing:'
			atmo_psf = gong_gaussian_psf(dx,dy,wpsf=atmo_w, nxhalf=nxhalf, nyhalf=nyhalf)
			image = convol(image,atmo_psf,/edge_truncate)
		endif else begin
			print,'Applying '+strtrim(string(wpsf),2)+' arcsecond PSF without atmospheric seeing:'
		endelse
		image = convol(image,psf,/edge_truncate)
	endif

	return, rebin(image, nx0, ny0, sample=sample)

end

function gong_sim_pixel_psf, image_in, origin0, dx0, dy0, origin, dx, dy, nx, ny, oversamp_fac=oversamp_fac, psf=psf, wpsf=wpsf, nopsf=nopsf, average=average, missing=missing, rota0=rota0, atmo_w=atmo_w, sample=sample, nxhalf=nxhalf, nyhalf=nyhalf

	nx0 = n_elements(image_in[*,0])
	ny0 = n_elements(image_in[0,*])

	if(keyword_set(nopsf) eq 0) then begin
		image = gong_sim_apply_psf(image_in, dx0, dy0, psf, wpsf=wpsf, atmo_w=atmo_w, sample=sample, nxhalf=nxhalf, nyhalf=nyhalf)
	endif else begin
		image = image_in
	endelse
	
	return, gong_sim_pixelize(image, origin0, dx0, dy0, origin, dx, dy, nx, ny, oversamp_fac=oversamp_fac, sample=sample, average=average, missing=missing, rota0=rota0)
	
end
