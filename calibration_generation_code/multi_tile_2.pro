function multi_tile_2, image_in, dx_in, dy_in, xb_fac=xb_fac, yb_fac=yb_fac, dx_out=dx_out, dy_out=dy_out, rsun=rsun, total=total, scales=scales, xsmin=xsmin, ysmin=ysmin, tile_levels = tile_levels, sample=sample

	nx_in = n_elements(image_in[*,0])
	ny_in = n_elements(image_in[0,*])

	if(n_elements(scales) eq 0) then scales = 2^([[0,1],[2,3]])
	if(n_elements(rsun) eq 0) then rsun = 981 ; Upper limit of solar size in arseconds
	if(n_elements(xsmin) eq 0) then xsmin = 2*rsun
	if(n_elements(ysmin) eq 0) then ysmin = 2*rsun
	
	if(n_elements(xb_fac) eq 0) then xb_fac = 1
	if(n_elements(yb_fac) eq 0) then yb_fac = 1
	
	dx_out = dx_in*xb_fac
	dy_out = dy_in*yb_fac
	
	nx_tile = round(nx_in/xb_fac)
	ny_tile = round(ny_in/yb_fac)
		
	ntile_x = max(scales)*ceil(0.5*xsmin/(nx_in*dx_in)/max(scales))
	ntile_y = max(scales)*ceil(0.5*ysmin/(ny_in*dy_in)/max(scales))
		
	nx_out = nx_tile*ntile_x*2
	ny_out = ny_tile*ntile_y*2
	
	
	image_out = fltarr(nx_out,ny_out)
	tile_levels = fltarr(nx_out,ny_out)
	
	for i=0,1 do begin
		for j=0,1 do begin
			ix = i*nx_tile*ntile_x
			jy = j*ny_tile*ntile_y
			nt_x = ntile_x/scales[i,j]
			nt_y = ntile_y/scales[i,j]
			nx_t = nx_tile*scales[i,j]
			ny_t = ny_tile*scales[i,j]
			for k=0,nt_x-1 do begin
				for l=0,nt_y-1 do begin
					kx = ix+k*nx_t
					ly = jy+l*ny_t
					image_out[kx:(kx+nx_t-1),ly:(ly+ny_t-1)] = rebin(image_in,nx_t,ny_t,sample=sample)
					tile_levels[kx:(kx+nx_t-1),ly:(ly+ny_t-1)] = scales[i,j]
					if(keyword_set(total)) then image_out[kx:(kx+nx_t-1),ly:(ly+ny_t-1)] *= 1.0*nx_in*ny_in/(1.0*nx_t*ny_t)
				endfor
			endfor
		endfor
	endfor
	
	return,image_out
	
end