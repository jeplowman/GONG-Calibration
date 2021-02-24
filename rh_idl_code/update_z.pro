pro update_z,i,j,k,z,z_out,info_struc, xbin=xbin, ybin=ybin, zbin=zbin, lbin=lbin

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
	
	if(n_elements(z_out) eq 0) then z_out = fltarr(nx,ny)
	
	ilo = i*xcrop
	ihi = (i+1)*xcrop-1
	jlo = j*ycrop
	jhi = (j+1)*ycrop-1

	z_out[ilo:ihi,jlo:jhi] = rebin(z,xcrop,ycrop)
	
end