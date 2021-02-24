function ray_resample_atmos, atmos, ray, geometry, geometry_resamp=geometry_resamp

	nx = geometry.nx
	ny = geometry.ny
	nz = geometry.nz
	nlos = ray.nlos
	nhydr = atmos.nhydr
	
	xr0 = reform(ray.xray[*,*,*,0])
	yr0 = reform(ray.yray[*,*,*,0])
	zr0 = reform(ray.zray[*,*,*,0])
	if(tag_exist(ray,'yray_good') eq 0) then begin
		yr0 = xr0
		for i=0,geometry.nx-1 do xr0[i,*,*]=i
	endif

	geometry_resamp = {nrays:geometry.nrays, nx:geometry.nx, ny:geometry.ny, nz:nlos, angleset:geometry.angleset, xmu:ray.mux, $
			ymu:ray.muy, wmu:geometry.wmu, dx:geometry.dx, dy:geometry.dy, z:reverse(ray.raylength), $
			vx:interpolate(geometry.vx,xr0,yr0,zr0), vy:interpolate(geometry.vy,xr0,yr0,zr0), vz:interpolate(geometry.vz,xr0,yr0,zr0)}
		
	atmos_resamp = {nhydr:atmos.nhydr, nelem:atmos.nelem, moving:atmos.moving, id:atmos.id, elements:atmos.elements, $
			stokes:atmos.stokes, backgrflags:atmos.backgrflags, backgrrecno:atmos.backgrrecno, nh:fltarr(nx,ny,nlos,nhydr), $
			t:interpolate(atmos.t,xr0,yr0,zr0), b:interpolate(atmos.b,xr0,yr0,zr0), gamma_b:interpolate(atmos.gamma_b,xr0,yr0,zr0), $
			chi_b:interpolate(atmos.chi_b,xr0,yr0,zr0), n_elec:interpolate(atmos.n_elec,xr0,yr0,zr0)}

	for i=0,nhydr-1 do atmos_resamp.nh[*,*,*,i] = interpolate(reform(atmos.nh[*,*,*,i]),xr0,yr0,zr0)

	return, atmos_resamp
	
end