function get_lambdacg, spec, lambda_in

	lambda0 = min(lambda_in)
	lambda = lambda_in-lambda0
	nw = n_elements(lambda)

	ic = interpol([spec[0],spec[nw-1]],[lambda[0],lambda[nw-1]],lambda)
	lambdacen = int_tabulated(lambda,lambda*(ic-spec),/double)/int_tabulated(lambda,ic-spec,/double)

	return,lambda0+lambdacen

end


pro get_xz, xray, zray, slice, index, x, z, nz
		
		x = reform(xray[slice,*,index])
		z = reform(zray[slice,*,index])
		ny = n_elements(x)
		keeps = where(abs(deriv(x)) lt ny/3)
		x = x[keeps]
		z = z[keeps]
		ishift = 0
		while((max(abs(deriv(shift(x,-ishift)))) gt ny/3) and (ishift lt 2*ny)) do ishift++
		print,ishift
		x = shift(x,-ishift)
		z = shift(z,-ishift)
		z = nz-1-z

end

pro get_xz0, x, z
		
		ny = n_elements(x)
		keeps = where(abs(deriv(x)) lt ny/3)
		x = x[keeps]
		z = z[keeps]
		ny2 = n_elements(x)
		;for ishift=0,ny2-1 do begin
		;	xshift = shift(x,-ishift)
		;	if(max(abs(xshift[0:ny2-2]-xshift[1:ny2-1])) lt ny/2) then break
		;endfor
		ishift = 0
		while((max(abs(deriv(shift(x,-ishift)))) gt ny/3) and (ishift lt 2*ny)) do ishift++
		print,ishift
		x = shift(x,-ishift)
		z = shift(z,-ishift)

end

pro get_xz_index, atmos, geometry, ray, spectrum, slice, index, xray, zray, x, z

	nx = n_elements(xray[0,*,0])
	nz = n_elements(geometry.z)
	nspect = n_elements(ray.nspect)
	ray_lambdas = spectrum.lambda[ray.nspect]
	dlambda = spectrum.lambda[1]-spectrum.lambda[0]

	x = dblarr(nx)
	z = dblarr(nx)
	for i=0,nx-1 do begin
		lambdacg = double(get_lambdacg(reform(spectrum.i[slice,i,*]),spectrum.lambda))
		x[i] = interpol(xray[slice,i,*],ray_lambdas-lambdacg,index*dlambda)
		z[i] = nz-1-interpol(zray[slice,i,*],ray_lambdas-lambdacg,index*dlambda)
	endfor

	get_xz0,x,z

end

pro make_contrib_function_plot, atmos, geometry, ray, spectrum, tau_index, contrib_indices, xray1, zray1, xray2, zray2, xray3, zray3, slice=slice,min=min,max=max, yr=yr, linestyles=linestyles, titlestr=titlestr, nominal_center_index = nominal_center_index

	nx = geometry.nx
	ny = geometry.ny
	nz = geometry.nz
	nspect = n_elements(ray.nspect)
	if(n_elements(slice) eq 0) then slice = round(ceil(nx/2.0)-1.0)
	if(n_elements(yr) ne 2) then yr = [0,ny-1]
	ny2 = yr[1]-yr[0]+1
	if(n_elements(linestyles) eq 0) then linestyles = [0,1,2]
	if(n_elements(nominal_center_index) eq 0) then nominal_center_index = 1 
	

	b = atmos.b[slice,yr[0]:yr[1],*]
	chi = atmos.chi_b[slice,yr[0]:yr[1],*]
	gamma = atmos.gamma_b[slice,yr[0]:yr[1],*]
	phi_ray = 90.0*!pi/180.0
	theta_ray = asin(ray.muy)
	blos = reform(b*((cos(phi_ray)*cos(chi) + sin(phi_ray)*sin(chi))*sin(theta_ray)*sin(gamma) + cos(theta_ray)*cos(gamma)))
	t = reform(atmos.t[slice,yr[0]:yr[1],*])
	n = reform(atmos.n_elec[slice,yr[0]:yr[1],*])

	dy = geometry.dy/1000000.0
	dz = (geometry.z[0]-geometry.z[1])/1000000.0

	;plot_image,reverse(blos,2),min=min,max=max
	;plot_image,reverse(t,2),min=min,max=max,orgin=[yr[0]*dy,0],scale=[dy,dz]
	greektau_octal = "164B
	greektau = '!9'+string(greektau_octal)+'!X'
	plot_image,reverse(alog10(n),2),min=min,max=max,origin=[yr[0]*dy,0.0],scale=[dy,dz],ytitle='Height (Mm)',title=titlestr
	tau_wavelength = 'Continuum ' + greektau + '=1'; strtrim(string(spectrum.lambda[ray.nspect[tau_index]]),2) + ' nm'
	get_xz, xray1, zray1, slice, tau_index, x, z, nz
	oplot,x*dy,z*dz,linestyle=linestyles[0]
	print,yr[0]*dy

	ncontrib_indices = n_elements(contrib_indices)
	if(ncontrib_indices gt 0) then begin
		contrib_waves = strtrim(string(spectrum.lambda[ray.nspect[contrib_indices]]),2) + ' nm'
		for i=0,ncontrib_indices-1 do begin
			;if(contrib_indices[i] eq 0 or contrib_indices[i] eq nspect-1) then begin
			;	get_xz, xray2, zray2, slice, contrib_indices[i], x2, z2, nz
			;	get_xz, xray3, zray3, slice, contrib_indices[i], x3, z3, nz
			;endif else begin
			;	get_xz_index,atmos,geometry,ray,spectrum,slice,contrib_indices[i]-nominal_center_index,xray2,zray2,x2,z2
			;	get_xz_index,atmos,geometry,ray,spectrum,slice,contrib_indices[i]-nominal_center_index,xray3,zray3,x3,z3
			;endelse
			get_xz, xray2, zray2, slice, contrib_indices[i], x2, z2, nz
			get_xz, xray3, zray3, slice, contrib_indices[i], x3, z3, nz
			oplot,x2*dy,z2*dz,linestyle=linestyles[i+1]
			oplot,x3*dy,z3*dz,linestyle=linestyles[i+1]
		endfor
		ssw_legend,[tau_wavelength,contrib_waves],linestyle=linestyles,/bottom,/left,charsize=1.0/sqrt(product(!p.multi[[1,2]] > 1))
	endif

end