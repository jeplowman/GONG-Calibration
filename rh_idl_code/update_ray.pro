pro update_ray, i, j, k, ray, ray_out, info_struc, xbin=xbin, ybin=ybin, zbin=zbin, lbin=lbin

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
	nlos = ray.nlos
	nspect = n_elements(ray.nspect)

	xcrop = info_struc.xcrop/xbin
	ycrop = info_struc.ycrop/ybin
	nperlambdabin = info_struc.nperlambdabin/lbin

	if(n_elements(ray_out) eq 0) then begin
		ray_out = {mux:ray.mux, muy:ray.muy, I:fltarr(nx,ny,nl), nspect:ray.nspect, chi:fltarr(nx,ny,nz,nspect), s:fltarr(nx,ny,nz,nspect), stokes_q:fltarr(nx,ny,nl), stokes_u:fltarr(nx,ny,nl), stokes_v:fltarr(nx,ny,nl), s2:fltarr(nx,ny,nlos,nspect), chi2:fltarr(nx,ny,nlos,nspect), raylength:ray.raylength, xray:fltarr(nx,ny,nlos,nspect), yray:fltarr(nx,ny,nlos,nspect), zray:fltarr(nx,ny,nlos,nspect), nlos:nlos}
	endif
	
	ilo = i*xcrop
	ihi = (i+1)*xcrop-1
	jlo = j*ycrop
	jhi = (j+1)*ycrop-1
	klo = k*nperlambdabin
	khi = (k+1)*nperlambdabin-1

	ray_out.I[ilo:ihi,jlo:jhi,klo:khi] = rebin(ray.I,xcrop,ycrop,nperlambdabin)
	ray_out.stokes_Q[ilo:ihi,jlo:jhi,klo:khi] = rebin(ray.stokes_Q,xcrop,ycrop,nperlambdabin)
	ray_out.stokes_U[ilo:ihi,jlo:jhi,klo:khi] = rebin(ray.stokes_U,xcrop,ycrop,nperlambdabin)
	ray_out.stokes_V[ilo:ihi,jlo:jhi,klo:khi] = rebin(ray.stokes_V,xcrop,ycrop,nperlambdabin)
	
	ray_out.chi[ilo:ihi,jlo:jhi,*,*] = rebin(ray.chi,xcrop,ycrop,nz,nspect)
	ray_out.s[ilo:ihi,jlo:jhi,*,*] = rebin(ray.s,xcrop,ycrop,nz,nspect)

	ray_out.chi2[ilo:ihi,jlo:jhi,*,*] = rebin(ray.chi2,xcrop,ycrop,nlos,nspect)
	ray_out.s2[ilo:ihi,jlo:jhi,*,*] = rebin(ray.s2,xcrop,ycrop,nlos,nspect)
	ray_out.xray[ilo:ihi,jlo:jhi,*] = rebin(ray.xray,xcrop,ycrop,nlos)/xbin
	ray_out.yray[ilo:ihi,jlo:jhi,*] = rebin(ray.yray,xcrop,ycrop,nlos)/ybin+jlo
	ray_out.zray[ilo:ihi,jlo:jhi,*] = rebin(ray.zray,xcrop,ycrop,nlos)

end	
