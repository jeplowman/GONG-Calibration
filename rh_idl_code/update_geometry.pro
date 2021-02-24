pro update_geometry,i,j,k,geometry,geometry_out, info_struc, xbin=xbin, ybin=ybin, zbin=zbin, lbin=lbin

	if(n_elements(xbin) eq 0) then xbin=1
	if(n_elements(ybin) eq 0) then ybin=1
	if(n_elements(zbin) eq 0) then zbin=1
	if(n_elements(lbin) eq 0) then lbin=1

	nx0 = info_struc.n0
	ny0 = info_struc.n2
	nz0 = info_struc.n1
	nl0 = info_struc.nlambdas
	
	nx = nx0/round(xbin)
	ny = ny0/round(ybin)
	nz = nz0/round(zbin)
	nl = nl0/round(lbin)
	
	xcrop = info_struc.xcrop/xbin
	ycrop = info_struc.ycrop/ybin
	nperlambdabin = info_struc.nperlambdabin/lbin

	nlambdas = info_struc.nlambdas/lbin
	
	if(n_elements(geometry_out) eq 0) then begin
		geometry_out = {nrays:geometry.nrays, nx:nx, ny:ny, nz:nz, angleset:geometry.angleset, xmu:geometry.xmu, ymu:geometry.ymu, $
				wmu:geometry.wmu, dx:geometry.dx*xbin, dy:geometry.dy*ybin, z:rebin(geometry.z,nz), vx:fltarr(nx,ny,nz), vy:fltarr(nx,ny,nz), vz:fltarr(nx,ny,nz)}
	endif
	
	ilo = i*xcrop
	ihi = (i+1)*xcrop-1
	jlo = j*ycrop
	jhi = (j+1)*ycrop-1

	geometry_out.vx[ilo:ihi,jlo:jhi,*] = rebin(geometry.vx,xcrop,ycrop,nz)
	geometry_out.vy[ilo:ihi,jlo:jhi,*] = rebin(geometry.vy,xcrop,ycrop,nz)
	geometry_out.vz[ilo:ihi,jlo:jhi,*] = rebin(geometry.vz,xcrop,ycrop,nz)
	
end