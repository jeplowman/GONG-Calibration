function get_contrib_function_depth, ray, fcontrib, xray=xray, yray=yray, zray=zray, tauray=tauray, contray=contray, taudepth=taudepth, spectinds=spectinds

	if(n_elements(spectinds) eq 0) then spectinds = indgen(n_elements(ray.nspect))
	nx = n_elements(ray.xray[*,0,0,0])
	ny = n_elements(ray.xray[0,*,0,0])
	nlos = n_elements(ray.raylength)
	nw = n_elements(spectinds)
	
	tauray = fltarr(nx,ny,nlos,nw)
	contray = fltarr(nx,ny,nlos,nw)
	depth = fltarr(nx,ny,nw)
	xray = fltarr(nx,ny,nw)
	yray = fltarr(nx,ny,nw)
	zray = fltarr(nx,ny,nw)

	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			for k0=0,nw-1 do begin
				k = spectinds[k0]
				dl = deriv(ray.raylength)
				tau = total(reform(ray.chi2[i,j,*,k])*dl,/cumulative)
				cont = ray.s2[i,j,*,k]*exp(-(tau < 50.0))*ray.chi2[i,j,*,k]
				ilambda = total(cont*dl,/cumulative)
				if(keyword_set(taudepth) eq 0) then depth[i,j,k0] = interpol(ray.raylength,ilambda,fcontrib*max(ilambda))
				if(keyword_set(taudepth)) then depth[i,j,k] = interpol(ray.raylength,tau,fcontrib)
				xray[i,j,k0] = interpol(ray.xray[i,j,*],ray.raylength,depth[i,j,k0])
				yray[i,j,k0] = interpol(ray.yray[i,j,*],ray.raylength,depth[i,j,k0])
				zray[i,j,k0] = interpol(ray.zray[i,j,*],ray.raylength,depth[i,j,k0])
				tauray[i,j,*,k0] = tau
				contray[i,j,*,k0] = cont
			endfor
		endfor
	endfor

	return,depth
	
end


function get_contrib_function_depth_reduced, ray, fcontrib, zray=zray, taudepth=taudepth

	nx = n_elements(ray.s2[*,0,0,0])
	ny = n_elements(ray.s2[0,*,0,0])
	nlos = n_elements(ray.raylength)
	nw = n_elements(ray.nspect)
	
	depth = fltarr(nx,ny,nw)
	zray = fltarr(nx,ny,nw)

	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			for k=0,nw-1 do begin
				dl = deriv(ray.raylength)
				tau = total(reform(ray.chi2[i,j,*,k])*dl,/cumulative)
				cont = ray.s2[i,j,*,k]*exp(-(tau < 50.0))*ray.chi2[i,j,*,k]
				ilambda = total(cont*dl,/cumulative)
				if(keyword_set(taudepth) eq 0) then depth[i,j,k] = interpol(ray.raylength,ilambda,fcontrib*max(ilambda))
				if(keyword_set(taudepth)) then depth[i,j,k] = interpol(ray.raylength,tau,fcontrib)
				zray[i,j,k] = interpol(ray.zray[i,j,*],ray.raylength,depth[i,j,k])
			endfor
		endfor
	endfor

	return,depth
	
end
