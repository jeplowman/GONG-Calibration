function bz_tautile_samp, gong_geometry, rh_geometry, bz_in, zarr, nopsf=nopsf, scales=scales, psf=psf, thetay=thetay, atmo_w=atmo_w, atmos_geometry=atmos_geometry, bz_cube=bz_cube, ray=ray, fcontrib=fcontrib, bz_tauone=bz_tauone, zr0_inv=zr0_inv, sample=sample, wpsf=wpsf

	au_m_toasec = (3600.0*180/!pi)/1.496e11
	dx_asec = rh_geometry.dx*au_m_toasec
	dy_asec = rh_geometry.dy*au_m_toasec		
	
	if(n_elements(scales) eq 0 and tag_exist(gong_geometry,'scales')) then scales = gong_geometry.scales
	if(n_elements(ray) eq 0) then ray = 'None'
	
	ondisk = gong_geometry.ondisk
	nx = gong_geometry.nx
	ny = gong_geometry.ny
	nz = n_elements(bz_in[0,0,*])
	dx = gong_geometry.dx
	dy = gong_geometry.dy

	bz_cube = bz_in
	if(keyword_set(thetay) eq 1 and (n_elements(fcontrib) eq 0 and is_struct(ray) eq 0)) then begin
		; Would be better to interpolate rather than integer shift, but the MHD simulation has much higher resolution
		; than GONG, so this should suffice:
		zar = reverse(atmos_geometry.z);-0.5*(atmos_geometry.z[0]-atmos_geometry.z[1])
		for i=0,nz-1 do bz_cube[*,*,i] = shift(bz_in[*,*,i],0,round(zar[i]/atmos_geometry.dy/cos(thetay)))
	endif

	gongorigin = gong_geometry.gongorigin
	rhorigin = gong_geometry.rhorigin
	xb_fac = gong_geometry.xb_fac
	yb_fac = gong_geometry.yb_fac
	
	nz2 = 51
	dz = 3.0e3
	za2 = dz*(findgen(nz2)-0.5*(nz2-1))
	
	if(is_struct(ray) eq 0) then begin
		print,'Using simple vertical ray path for 3D MHD height sampling:'
		print,'Averaging over the following heights:',median(zarr)+za2
		bz_tauone = fltarr(rh_geometry.nx,rh_geometry.ny)
		for i=0,rh_geometry.nx-1 do for j=0,rh_geometry.ny-1 do bz_tauone[i,j] = mean(interpol(reform(bz_cube[i,j,*]),rh_geometry.z,zarr[i,j]+za2))
	endif else begin
		print,'Using complex path for 3D MHD height sampling:'
		xr = fltarr(rh_geometry.nx,rh_geometry.ny,nz2)
		yr = fltarr(rh_geometry.nx,rh_geometry.ny,nz2)
		zr = fltarr(rh_geometry.nx,rh_geometry.ny,nz2)
		xr0 = reform(ray.xray[*,*,*,0])
		yr0 = reform(ray.yray[*,*,*,0])
		zr0 = reform(ray.zray[*,*,*,0])
		zr0_inv = reform(ray.zray[*,*,*,0])
		if(tag_exist(ray,'yray_good') eq 0) then begin
			yr0 = xr0
			for i=0,rh_geometry.nx-1 do xr0[i,*,*]=i
		endif
		; Change this here so it can use height input (zarr) if no fcontrib:
		if(n_elements(fcontrib) ne 0) then begin
			print,'Using contrib function depth'
			depth = get_contrib_function_depth(ray, fcontrib, contray=contray, spectind = [0,floor(n_elements(ray.nspect)/2)])
;			depth = mean(depth[*,*,0:2],dimension=3)
;			depth = mean(depth,dimension=3)
;			for i=0,rh_geometry.nx-1 do for j=0,rh_geometry.ny-1 do begin
;				xr[i,j] = interpol(reform(xr0[i,j,*]),ray.raylength,depth[i,j])
;				yr[i,j] = interpol(reform(yr0[i,j,*]),ray.raylength,depth[i,j])
;				zr[i,j] = interpol(reform(zr0[i,j,*]),ray.raylength,depth[i,j])
;			endfor
;			contray = (total(contray,4)+reform(contray[*,*,*,1]))/4.0
;			contray = shift(reform(contray[*,*,*,1]),0,0,0)
;			bz_cube = interpolate(bz_cube,xr0,yr0,zr0)
;			bz_tauone = total(bz_cube*contray,3)/total(contray,3)
			help,contray
			help,ray
			bz_tauone = total(interpolate(bz_cube,xr0,yr0,zr0)*reform(contray[*,*,*,1]),3)/total(reform(contray[*,*,*,1]),3)*5.0/5.0
			;bz_tauone += total(interpolate(bz_cube,xr0,yr0,zr0)*reform(contray[*,*,*,2]),3)/total(reform(contray[*,*,*,2]),3)*0.0/5.0
		endif else begin
			print,'Using input height:'
			print,'Averaging over the following heights:',median(zarr)+za2
			help,xr,yr,zr,xr0,yr0,zr0
			for i=0,rh_geometry.nx-1 do for j=0,rh_geometry.ny-1 do begin
				zr0_inv[i,j,*] = interpol(rh_geometry.z,findgen(rh_geometry.nz),reform(zr0[i,j,*]))
				xr[i,j,*] = interpol(reform(xr0[i,j,*]),zr0_inv[i,j,*],zarr[i,j]+za2)
				yr[i,j,*] = interpol(reform(yr0[i,j,*]),zr0_inv[i,j,*],zarr[i,j]+za2)
				zr[i,j,*] = interpol(reform(zr0[i,j,*]),zr0_inv[i,j,*],zarr[i,j]+za2)
			endfor
			print,reform(zr0[0,0,*])
			print,reform(zr0_inv[0,0,*])
			bz_tauone = interpolate(bz_cube,xr,yr,zr)
			;message,'Breakpoint!'
			if(nz2 gt 1) then bz_tauone = total(bz_tauone,3)/float(nz2)
		endelse
	endelse
		
	bz_tauone_tile = multi_tile_2(bz_tauone,dx_asec,dy_asec,xb_fac=xb_fac,yb_fac=yb_fac,dx_out=dx0,dy_out=dy0,scales=scales, sample=sample)
	bz_tauone_tile_lo = 1e4*gong_sim_pixel_psf(bz_tauone_tile,rhorigin,dx0,dy0,gongorigin,dx,dy,nx,ny,nopsf=nopsf,psf=psf, sample=sample, atmo_w=atmo_w, wpsf=wpsf)*ondisk/dx/dy

	return, bz_tauone_tile_lo
	
end
