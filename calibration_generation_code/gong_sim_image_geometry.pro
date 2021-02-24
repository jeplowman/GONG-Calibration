function gong_sim_image_geometry, gongheader, savstr, geometry, r_ondisk_frac=r_ondisk_frac, arflags=arflags

	au_cm = 1.49e13
	rsun_cm = 6.957e10
	au_m_toasec = (3600.0*180/!pi)/1.496e11

	nx = sxpar(gongheader,'NAXIS2') ; These are flipped because the gong raws have north-south on the x axis,
	ny = sxpar(gongheader,'NAXIS1') ; but we transpose them when we read the data.

	xc = 0.5*nx
	yc = 0.5*ny
	dx = 2.5
	dy = 2.5

	xa_px = (findgen(nx)#(1.0+fltarr(ny)))
	ya_px = ((1.0+fltarr(nx))#findgen(ny))

	xa = dx*(findgen(nx)#(1.0+fltarr(ny))-xc)
	ya = dy*((1.0+fltarr(nx))#findgen(ny)-yc)
	xa_cm = (xa*au_cm)*(!pi/180/3600)
	ya_cm = (ya*au_cm)*(!pi/180/3600)

	ra_cm = sqrt(xa_cm*xa_cm+ya_cm*ya_cm)

	if(n_elements(r_ondisk_frac) eq 0) then r_ondisk_frac = 0.95
	ondisk = ra_cm le rsun_cm*r_ondisk_frac

	phi = asin(xa_cm/rsun_cm)*180/!pi
	theta = acos((ya*au_cm/rsun_cm)*(!pi/180/3600))*180/!pi-90

	dr = diff_rot(1.0/24/3600,theta)*!pi/180

	vel = xa_cm*dr

	xb_fac = savstr.xb_fac
	yb_fac = savstr.yb_fac

	dx_asec = geometry.dx*au_m_toasec
	dy_asec = geometry.dy*au_m_toasec

	xa0_px = findgen(geometry.nx)#(1+fltarr(geometry.ny))
	ya0_px = (1+fltarr(geometry.nx))#findgen(geometry.ny)
	ra0_px = sqrt((xa0_px-0.5*(geometry.nx-1))^2+(ya0_px-0.5*(geometry.ny-1))^2)

	if(tag_exist(savstr,'scales')) then begin
		scales = savstr.scales
	endif else begin
		scales = 2^([[0,1],[2,3]])
	endelse
	
	xasim = (findgen(geometry.nx)-0.5*geometry.nx)#(1+fltarr(geometry.ny))
	yasim = (1+fltarr(geometry.nx))#(findgen(geometry.ny)-0.5*geometry.ny)
	rasim = sqrt(xasim*xasim+yasim*yasim)

	arflags = multi_tile_2(rasim lt 266,dx_asec,dy_asec,xb_fac=xb_fac,yb_fac=yb_fac,dx_out=dx0,dy_out=dy0,tile_levels=tile_levels,scales=scales)

	nx0 = n_elements(arflags[*,0])
	ny0 = n_elements(arflags[0,*])

	xref = 0.0
	yref = 0.0
	ixref = 0.5*nx0
	iyref = 0.5*ny0

	rhscale = [dx0,dy0]
	rhorigin = find_image_center([0.0,0.0],rhscale,[ixref,iyref],refval=[xref,yref])

	gongny = sxpar(gongheader,'NAXIS1')
	gongnx = sxpar(gongheader,'NAXIS2')
	yc_pix = sxpar(gongheader,'X_CENTER')
	xc_pix = sxpar(gongheader,'Y_CENTER')
	scale = [2.5,2.5]
	rollcen = [0.0,0.0]
	rota = 0.0
	px_refval = 0.5*[gongnx,gongny]; [xc_pix,yc_pix] ; These center values don't make sense as sun center pixels - they're way off!?
	gongorigin = find_image_center([0,0], scale, px_refval, refval=refval, rollcen=rollcen, px_rollcen=px_rollcen, rota=rota)

	tile_levels_lo = gong_sim_pixelize(tile_levels,rhorigin,dx0,dy0,gongorigin,dx,dy,nx,ny)/dx/dy
	arflags = gong_sim_pixelize(arflags,rhorigin,dx0,dy0,gongorigin,dx,dy,nx,ny)/dx/dy

	gong_geometry = {gongorigin:gongorigin, rhorigin:rhorigin, dx:dx, dy:dy, nx:nx, ny:ny, ondisk:ondisk, xb_fac:xb_fac, yb_fac:yb_fac, tile_levels:tile_levels_lo, scales:scales}

	rdilate=8
	nxdilate = 2*rdilate+1
	nydilate = 2*rdilate+1
	xadilate = (findgen(nxdilate)-rdilate)#(1+fltarr(nydilate))
	yadilate = (1+fltarr(nxdilate))#(findgen(nydilate-rdilate))
	yadilate = (1+fltarr(nxdilate))#(findgen(nydilate)-rdilate)
	s_dilate = sqrt(xadilate^2+yadilate^2) lt rdilate
	arflags = dilate(arflags,s_dilate)


	return, gong_geometry

end