function ray_interpolate_3D_rh_func, g0, mu, f0, x_out, y_out, z_out, raylength=raylength, xdist=xdist, zdist=zdist, xaxis=xaxis, yaxis=yaxis

	if(keyword_set(xaxis) eq 0 and keyword_set(yaxis) eq 0) then message, 'Please specify xaxis or yaxis in ray_interpolate_3D_rh_func'

	if(keyword_set(yaxis)) then begin
		g2 = {nrays:g0.nrays, nx:g0.ny, ny:g0.nx, nz:g0.nz, angleset:g0.angleset, xmu:mu, ymu:g0.xmu, dx:g0.dy, dy: g0.dx, $
				z:g0.z, vx:transpose(g0.vy,[1,0,2]), vy:transpose(g0.vx,[1,0,2]), vz:transpose(g0.vz,[1,0,2]), x:g0.dy*dindgen(g0.ny)}
		f2 = transpose(f0,[1,0,2])
	endif
	
	if(keyword_set(xaxis)) then begin
		g2 = create_struct(g0,'x',g0.dx*dindgen(g0.nx))
		g2.xmu = mu
		f2 = f0
	endif
	
	nx = g2.nx
	ny = g2.ny
	nz = g2.nz
	nf = n_elements(f0[0,0,0,*])
	
	xa = fltarr(nx,ny,nz)
	ya = fltarr(nx,ny,nz)
	za = fltarr(nx,ny,nz)
	for i=0,nx-1 do xa[i,*,*] = i
	for i=0,ny-1 do ya[*,i,*] = i
	for i=0,nz-1 do za[*,*,i] = i
	
	f_out = fltarr(nx,ny,nz)
		
	for i=0,ny-1 do begin
		for j=0,nx-1 do begin
			rt = raytrace(g2,0,j,xdist=xdist,zdist=zdist)
			raylength = sqrt(xdist*xdist+zdist*zdist)
			if(i eq 0 and j eq 0) then begin
				nlos = n_elements(raylength)
				f_out = fltarr(nx,ny,nlos,nf)
				x_out = fltarr(nx,ny,nlos)
				y_out = fltarr(nx,ny,nlos)
				z_out = fltarr(nx,ny,nlos)
				raylength0 = raylength
				xdist0 = xdist
				zdist0 = zdist
			endif
			for k=0,nf-1 do begin
				fji_out = rayinterpolate(reform(f2(*,i,*,k)),rt)
				; raytrace doesn't always give the same ray interpolation points, which is a problem when we need
				; to return a 3D array of uniformly gridded values. So we reinterpolate using arc length. This should
				; be OK as long as each ray has consistent start and end points, which appears to be the case.
				f_out(j,i,*,k) = interpol(fji_out,raylength,raylength0)
			endfor
			xji_out = rayinterpolate(reform(xa(*,i,*)),rt)
			yji_out = rayinterpolate(reform(ya(*,i,*)),rt)
			zji_out = rayinterpolate(reform(za(*,i,*)),rt)
			x_out(j,i,*) = interpol(xji_out,raylength,raylength0)
			y_out(j,i,*) = interpol(yji_out,raylength,raylength0)
			z_out(j,i,*) = interpol(zji_out,raylength,raylength0)
			message,'Breakpoint!'
		endfor
	endfor
	
;	message,'Breakpoint!'

	if(keyword_set(yaxis)) then begin
		if(nf gt 1) then begin
			f_out = transpose(f_out,[1,0,2,3])
			x_out = transpose(x_out,[1,0,2,3])
			y_out = transpose(y_out,[1,0,2,3])
			z_out = transpose(z_out,[1,0,2,3])
		endif else begin
			f_out = transpose(f_out,[1,0,2])
			x_out = transpose(x_out,[1,0,2])
			y_out = transpose(y_out,[1,0,2])
			z_out = transpose(z_out,[1,0,2])
		endelse
	endif
	
	xdist=xdist0
	zdist=zdist0
	raylength=raylength0
	
	return, f_out
	
end