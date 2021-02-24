function refactor_inclined_ray, g0, ray, xaxis=xaxis, yaxis=yaxis

	s0 = ray.s
	chi0 = ray.chi

	if(keyword_set(xaxis) eq 0 and keyword_set(yaxis) eq 0) then message, 'Please specify xaxis or yaxis in refactor_inclined_ray'

	if(keyword_set(yaxis)) then begin
		mu = ray.muy
		g2 = {nrays:g0.nrays, nx:g0.ny, ny:g0.nx, nz:g0.nz, angleset:g0.angleset, xmu:mu, ymu:g0.xmu, dx:g0.dy, dy: g0.dx, $
				z:g0.z, vx:transpose(g0.vy,[1,0,2]), vy:transpose(g0.vx,[1,0,2]), vz:transpose(g0.vz,[1,0,2]), x:g0.dy*dindgen(g0.ny)}
		s2 = transpose(s0,[1,0,2,3])
		chi2 = transpose(chi0,[1,0,2,3])
	endif
	
	if(keyword_set(xaxis)) then begin
		mu = ray.mux
		g2 = create_struct(g0,'x',g0.dx*dindgen(g0.nx))
		g2.xmu = mu
		s2 = s0
		chi2 = chi0
	endif
	
	nx = g2.nx
	ny = g2.ny
	nz = g2.nz
	nf = n_elements(ray.nspect)
	
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
				s_out = fltarr(nx,ny,nlos,nf)
				chi_out = fltarr(nx,ny,nlos,nf)
				x_out = fltarr(nx,ny,nlos)
				y_out = fltarr(nx,ny,nlos)
				z_out = fltarr(nx,ny,nlos)
				raylength0 = raylength
				xdist0 = xdist
				zdist0 = zdist
			endif
			for k=0,nf-1 do begin
				sji_out = rayinterpolate(reform(s2(*,i,*,k)),rt)
				chiji_out = rayinterpolate(reform(chi2(*,i,*,k)),rt)
				; raytrace doesn't always give the same ray interpolation points, which is a problem when we need
				; to return a 3D array of uniformly gridded values. So we reinterpolate using arc length. This should
				; be OK as long as each ray has consistent start and end points, which appears to be the case.
				s_out(j,i,*,k) = interpol(sji_out,raylength,raylength0)
				chi_out(j,i,*,k) = interpol(chiji_out,raylength,raylength0)
			endfor
			xji_out = rayinterpolate(reform(xa(*,i,*)),rt)
			yji_out = rayinterpolate(reform(ya(*,i,*)),rt)
			zji_out = rayinterpolate(reform(za(*,i,*)),rt)
			x_out(j,i,*) = interpol(xji_out,raylength,raylength0)
			y_out(j,i,*) = interpol(yji_out,raylength,raylength0)
			z_out(j,i,*) = interpol(zji_out,raylength,raylength0)
		endfor
	endfor
	
;	message,'Breakpoint!'

	if(keyword_set(yaxis)) then begin
		if(nf gt 1) then begin
			s_out = transpose(s_out,[1,0,2,3])
			chi_out = transpose(chi_out,[1,0,2,3])
			x_out = transpose(x_out,[1,0,2])
			y_out = transpose(y_out,[1,0,2])
			z_out = transpose(z_out,[1,0,2])
		endif else begin
			s_out = transpose(s_out,[1,0,2])
			chi_out = transpose(chi_out,[1,0,2])
			x_out = transpose(x_out,[1,0,2])
			y_out = transpose(y_out,[1,0,2])
			z_out = transpose(z_out,[1,0,2])
		endelse
	endif
	
	xdist=xdist0
	zdist=zdist0
	raylength=raylength0
	
	return, {mux:ray.mux, muy:ray.muy, I:ray.I, nspect:ray.nspect, chi:ray.chi, s:ray.s, chi2:chi_out, s2:s_out, stokes_q:ray.stokes_q, stokes_u:ray.stokes_u, stokes_v:ray.stokes_v, xray:x_out, yray:y_out, zray:z_out, raylength:raylength, nlos:nlos}
	
end