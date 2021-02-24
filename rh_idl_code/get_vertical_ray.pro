function get_vertical_ray, nspect

	@files.common
	@geometry.common
	@spectrum.common
	@opacity.common

	result = openJ('J.dat')
	backgroundFile = 'background.dat'
	result = openOpacity('opacity.out')

	KM_TO_M = 1.0E3
	MIN_FRACTION = 1.0E-6
	
	z = geometry.z
	dz = deriv(z)
;	dz = z - shift(z, 1)
;	dz[0] = dz[1]
	
	nx = geometry.nx
	ny = geometry.ny
	nz = n_elements(z)
	
	nlamb = n_elements(nspect)
	
	nlos = n_elements(z)
	s_out = fltarr(nx,ny,nlos,nlamb)
	tau_out = s_out
	chi_out = s_out
	x_out = fltarr(nx,ny,nlos)
	y_out = fltarr(nx,ny,nlos)
	z_out = fltarr(nx,ny,nlos)
	for i=0,nx-1 do x_out[i,*,*] = i
	for i=0,ny-1 do y_out[*,i,*] = i
	for i=0,nz-1 do z_out[*,*,i] = i
	
	for la = 0, nlamb-1 do begin
		readJ, nspect[la]
		readOpacity, nspect[la]
		for i=0,nx-1 do begin
			for k=0,ny-1 do begin
				tau = -total(dz*(chi_as[i,k,*]+chi_c[i,k,*]),/cumulative)
				dtau = deriv(tau)
				s = (eta_c[i,k,*]+J[i,k,*]*scatt[i,k,*])/(chi_c[i,k,*]+chi_as[i,k,*])
				s_out[i,k,*,la] = s
				chi_out[i,k,*,la] = chi_as[i,k,*]+chi_c[i,k,*]
				tau_out[i,k,*,la] = tau
			endfor
		endfor
	endfor

	free_lun, Junit, opacUnit, backgroundUnit
	
	return, {mux:0.0, muy:0.0, I:spectrum.I, nspect:nspect, chi:chi_out, s:s_out, chi2:chi_out, s2:s_out, stokes_q:spectrum.stokes_q, stokes_u:spectrum.stokes_u, stokes_v:spectrum.stokes_v, xray:x_out, yray:y_out, zray:z_out, raylength:reverse(z), nlos:nlos}
	
	
end