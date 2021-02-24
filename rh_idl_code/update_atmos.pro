pro update_atmos,i,j,k,atmos,atmos_out,info_struc, xbin=xbin, ybin=ybin, zbin=zbin, lbin=lbin

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
	
	
	if(n_elements(atmos_out) eq 0) then begin
		atmos_out = {nhydr:atmos.nhydr, nelem:atmos.nelem, moving:atmos.moving, T:fltarr(nx,ny,nz), nh:fltarr(nx,ny,nz,atmos.nhydr), id:atmos.id, elements:atmos.elements, stokes:atmos.stokes, B:fltarr(nx,ny,nz), gamma_b:fltarr(nx,ny,nz), chi_b:fltarr(nx,ny,nz), n_elec:fltarr(nx,ny,nz), backgrflags:replicate(atmos.backgrflags[0],nl), backgrrecno:atmos.backgrrecno}
		for ii=0,nl-1 do begin
			atmos_out.backgrflags[ii].hasline = 1
			atmos_out.backgrflags[ii].ispolarized = 1
		endfor
	endif
	
	ilo = i*xcrop
	ihi = (i+1)*xcrop-1
	jlo = j*ycrop
	jhi = (j+1)*ycrop-1
	klo = k*nperlambdabin
	khi = (k+1)*nperlambdabin-1
	
	atmos_out.t[ilo:ihi,jlo:jhi,*] = rebin(atmos.t,xcrop,ycrop,nz)
	atmos_out.n_elec[ilo:ihi,jlo:jhi,*] = rebin(atmos.n_elec,xcrop,ycrop,nz)
;	atmos_out.vturb[ilo:ihi,jlo:jhi,*] = rebin(atmos.vturb,xcrop,xcrop,nz)
	atmos_out.nh[ilo:ihi,jlo:jhi,*,*] = rebin(atmos.nh,xcrop,ycrop,nz,atmos.nhydr)
	atmos_out.b[ilo:ihi,jlo:jhi,*] = rebin(atmos.b,xcrop,ycrop,nz)
	atmos_out.gamma_b[ilo:ihi,jlo:jhi,*] = rebin(atmos.gamma_b,xcrop,ycrop,nz)
	atmos_out.chi_b[ilo:ihi,jlo:jhi,*] = rebin(atmos.chi_b,xcrop,ycrop,nz)

end
