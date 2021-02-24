pro update_spectrum,i,j,k,spectrum,spectrum_out,info_struc, xbin=xbin, ybin=ybin, zbin=zbin, lbin=lbin

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
	lambdas = rebin(info_struc.lambdas,nperlambdabin)

	if(n_elements(spectrum_out) eq 0) then begin
		spectrum_out = {nspect:nlambdas, lambda:lambdas, i:fltarr(nx,ny,nl), vacuum_to_air:spectrum.vacuum_to_air, air_limit:spectrum.air_limit, $
				stokes_q:fltarr(nx,ny,nl), stokes_u:fltarr(nx,ny,nl), stokes_v:fltarr(ny,ny,nl), as_rn:lonarr(nl)}
	endif
	
	ilo = i*xcrop
	ihi = (i+1)*xcrop-1
	jlo = j*ycrop
	jhi = (j+1)*ycrop-1
	klo = k*nperlambdabin
	khi = (k+1)*nperlambdabin-1

	spectrum_out.I[ilo:ihi,jlo:jhi,klo:khi] = rebin(spectrum.I,xcrop,ycrop,nperlambdabin)
	spectrum_out.stokes_Q[ilo:ihi,jlo:jhi,klo:khi] = rebin(spectrum.stokes_Q,xcrop,ycrop,nperlambdabin)
	spectrum_out.stokes_U[ilo:ihi,jlo:jhi,klo:khi] = rebin(spectrum.stokes_U,xcrop,ycrop,nperlambdabin)
	spectrum_out.stokes_V[ilo:ihi,jlo:jhi,klo:khi] = rebin(spectrum.stokes_V,xcrop,ycrop,nperlambdabin)
	spectrum_out.as_rn[klo:khi] = rebin(spectrum.as_rn,nperlambdabin)
	
end

