function do_mc_cal_inner, img, imgflags, cal_meas, cal_gt, calflags

	imgvals = img[imgflags]
	n_imgvals = n_elements(imgvals)
	imgvals += 0.1*randomn(seed,n_imgvals)
	cal_measvals = [-cal_meas[calflags],cal_meas[calflags]]
	cal_gtvals = [-cal_gt[calflags],cal_gt[calflags]]
	n_calvals = n_elements(cal_measvals)
	msort = sort(cal_measvals)
	indices = floor(interpol(lindgen(n_calvals),cal_measvals[msort],imgvals)+randomu(seed,n_imgvals))
	meas_cal = cal_gtvals[msort[indices]]

	print,n_elements(cal_meas),n_elements(cal_gt),n_elements(imgvals),min(meas_cal),max(meas_cal),min(cal_gtvals),max(cal_gtvals)

	nbins = 40
	binbs = floor((n_calvals-1)*findgen(nbins+1)/nbins)
	binlocs = 0.5*(binbs[0:nbins-1]+binbs[1:nbins])
	deltas = fltarr(nbins)
	for i=0,nbins-1 do deltas[i] = 0.1+ 1.0*(stdev(cal_gtvals[msort[binbs[i]:binbs[i+1]]])^2+stdev(cal_measvals[msort[binbs[i]:binbs[i+1]]])^2)/sqrt(n_calvals)
	
	return, meas_cal+(interpol(deltas,binlocs,imgvals) < 1.0)*randomn(seed,n_imgvals)

end

function do_mc_cal, image, qsflags, arflags, angles, cal_measimages, cal_gtimages, cal_qsflags, cal_arflags, cal_angles, bzobs_qs, bzmodel_qs, bzobs_ar, bzmodel_ar, angles_mc

	n_angles = n_elements(cal_angles)
	nx = n_elements(image[*,0])
	ny = n_elements(image[0,*])

	img_cal = fltarr(nx,ny)
	ang_dithimg = randomu(seed,nx,ny)
	angle_indices0 = interpol(indgen(n_angles),cal_angles,angles_mc) < n_angles
	angle_indices = floor(ang_dithimg + angle_indices0) < ceil(angle_indices0)

	qs_px = where(qsflags)
	ar_px = where(arflags)

	bzi_cals = fltarr(nx,ny,n_angles)
	for i=0,n_angles-1 do begin
		bzi_cal = fltarr(nx,ny)	
		bzi_cal[qs_px] = interpol(bzmodel_qs[i,*],bzobs_qs[i,*],image[where(qsflags)])		
		bzi_cal[ar_px] = interpol(bzmodel_ar[i,*],bzobs_ar[i,*],image[where(arflags)])
		bzi_cals[*,*,i] = bzi_cal
	endfor

	for i=0,nx-1 do begin
		for j=0,ny-1 do begin
			img_cal[i,j] = interpol(bzi_cals[i,j,*],cal_angles,angles[i,j])
		endfor
	endfor
	
	for i=0,n_angles-1 do begin
		calimg_mc = fltarr(nx,ny)
		resid_img = fltarr(nx,ny)
		cal_measimg = cal_measimages[*,*,i]
		cal_gtimg = cal_gtimages[*,*,i]
		cal_wqs = where(cal_qsflags[*,*,i])
		cal_war = where(cal_arflags[*,*,i])
		resid_img[cal_wqs] = cal_gtimg[cal_wqs]-interpol(bzmodel_qs[i,*],bzobs_qs[i,*],cal_measimg[cal_wqs])
		resid_img[cal_war] = cal_gtimg[cal_war]-interpol(bzmodel_ar[i,*],bzobs_ar[i,*],cal_measimg[cal_war])
		calimg_mc[where(qsflags)] = do_mc_cal_inner(image, where(qsflags), cal_measimg, resid_img, cal_wqs)
		calimg_mc[where(arflags)] = do_mc_cal_inner(image, where(arflags), cal_measimg, resid_img, cal_war)
		img_cal[where(angle_indices eq i)]+=calimg_mc[where(angle_indices eq i)]
		print,n_elements(where(qsflags)),n_elements(where(arflags)),n_elements(where(angle_indices eq i)),min(calimg_mc),max(calimg_mc)
	endfor

	return, img_cal

end

function apply_calibration_curve2, bzobs_qs, bzmodel_qs, bzobs_ar, bzmodel_ar, bzifile, izifile, dlarge=dlarge, dsmall=dsmall, out_dir=out_dir, plot_dir=plot_dir, prange=prange, transpose=transpose, post_deconvolve=post_deconvolve, qspts_bzmod=qspts_bzmod, arpts_bzmod=arpts_bzmod, qspts_bzobs=qspts_bzobs, arpts_bzobs=arpts_bzobs, angles=angles, cal_measimages=cal_measimages, cal_gtimages=cal_gtimages, cal_qsflags=cal_qsflags, cal_arflags=cal_arflags, bzi_uncal=bzi_uncal, do_mc=do_mc, mc_max_sunspot_angle=mc_max_sunspot_angle

	;rsun_cm = 6.957e10
	;au_cm_toasec = (3600.0*180/!pi)/1.496e13
	;rsun_asec = rsun_cm/au_cm_toasec
	
	if(n_elements(dlarge) eq 0) then dlarge = 30
	if(n_elements(dsmall) eq 0) then dsmall = 3

	if(n_elements(bmax_qs) eq 0) then bmax_qs = max(min(abs(bzobs_qs),dimension=1))
	if(n_elements(bmax_ar) eq 0) then bmax_ar = 4000

	izi = readfits(izifile,izihdr)
	bzi = readfits(bzifile,bzihdr)
	datestr = strtrim(sxpar(bzihdr,'DATE-OBS'),2)
	timestr = strtrim(sxpar(bzihdr,'TIME-OBS'),2)
	timestr = strmid(timestr,0,strpos(timestr,'.')+2)
	bzi_uncal = bzi
	if(keyword_set(post_deconvolve)) then begin
		bzi0 = bzi
		bzi = gong_deconvolve_magnetogram(bzi, izi, inten=izi)
	endif else begin
		bzi0 = bzi
	endelse

	zerovals = where(bzi0 eq 0)
	bzi[zerovals] = 1.0e-10
	
	;izi_ms_large = median(izi,dlarge)
	;izi_ms_small = median(izi,dsmall)
	izi_ms_large = median2(izi,dlarge)
	izi_ms_small = median2(izi,dsmall)

	rsun_asec = sxpar(izihdr,'RADIUS')*3600.0D0*180.0D0/!PI
	
	nx = n_elements(izi[*,0])
	ny = n_elements(izi[0,*])

	xa = lindgen(nx)#(1+lonarr(ny))
	ya = (1+lonarr(nx))#lindgen(ny)

	regions = label_region(izi lt 0.5*max(izi_ms_large),/all_neighbors,/ulong)
	offdisk = regions eq regions[2,2] or (xa eq 0) or (xa eq nx-1) or (ya eq 0) or (ya eq ny-1)
	ondisk = offdisk eq 0

	xc = mean(xa[where(ondisk)])
	yc = mean(ya[where(ondisk)])
;	xc = sxpar(izihdr,'CRPIX1') ; Axes are flipped since we transposed image to make north up
;	yc = sxpar(izihdr,'CRPIX2')
	dx = 2.5*1.021;sxpar(izihdr,'CDELT1')
	dy = 2.5*1.021;sxpar(izihdr,'CDELT2')

	ra = sqrt((xa-xc)^2+(ya-yc)^2)
	rsun = rsun_asec/sqrt(dx*dy); max(ra[where(ondisk)])
	
	ondisk = ra lt rsun*0.9875
	offdisk = ra gt rsun*0.9875
	
	rmedimg = fltarr(nx,ny)
	rlvls = fltarr(nx,ny)
	nrbins = 50
	binno = findgen(nrbins)
	rbins = 0.98*rsun*sqrt(findgen(nrbins)/(nrbins-1.0))
	binno2 = spl_init(rbins,binno,/double)
	rlvls[where(ondisk)] = round(spl_interp(rbins,binno,binno2,ra[where(ondisk)]) < (nrbins -1))
	rlvls[where(offdisk or ra ge 0.98*rsun)]=-1
	
	for i=0,max(rlvls) do rmedimg[where(rlvls eq i)] = median(izi[where(rlvls eq i)])

	arflg = (izi lt 0.98*izi_ms_large)*(ra lt 0.95*rsun) or (ra lt 0.95*rsun)*(abs(bzi0) gt bmax_qs); 0*bzi0;
	qsflg = (arflg eq 0)*(ra lt 0.995*rsun)
	qs_px = where(qsflg)
	ar_px = where(arflg)
	n_angles = n_elements(angles)

	bzi_cals = fltarr(nx,ny,n_angles)
	;bzi[qs_px] = bzi[qs_px] < max(bzobs_qs)
	;bzi[qs_px] = bzi[qs_px] > min(bzobs_qs)
	;bzi[ar_px] = bzi[ar_px] < max(bzobs_ar)
	;bzi[ar_px] = bzi[ar_px] > min(bzobs_ar)
	img_angles = (abs(asin(ra/rsun))) < max(angles)
	img_angles_mc = img_angles
	if(n_elements(mc_max_sunspot_angle) eq 1) then begin
		img_angles_mc[ar_px] = img_angles_mc[ar_px] < mc_max_sunspot_angle
	endif
	if(keyword_set(do_mc)) then begin
		bzi_cal = do_mc_cal(bzi, qsflg, arflg, img_angles, cal_measimages, cal_gtimages, cal_qsflags, cal_arflags, angles, bzobs_qs, bzmodel_qs, bzobs_ar, bzmodel_ar, img_angles_mc)
	endif else begin
		for i=0,n_angles-1 do begin
			bzi_cal = fltarr(nx,ny)	
			bzi_cal[qs_px] = interpol(bzmodel_qs[i,*],bzobs_qs[i,*],bzi[qs_px])		
			bzi_cal[ar_px] = interpol(bzmodel_ar[i,*],bzobs_ar[i,*],bzi[ar_px])
			bzi_cals[*,*,i] = bzi_cal
		endfor
		
		bzi_cal = fltarr(nx,ny)
		testimg = fltarr(nx,ny)
		for i=0,nx-1 do begin
			for j=0,ny-1 do begin
				bzi_cal[i,j] = interpol(bzi_cals[i,j,*],angles,img_angles[i,j])
				testimg[i,j] = interpol(indgen(n_angles),angles,img_angles[i,j])
			endfor
		endfor
	endelse
	
	bzi_cal[where(offdisk)] = bzi[where(offdisk)]
		
	writefits, out_dir+file_basename(bzifile), bzi_cal, bzihdr, /compress
	
	print,'Calibrated net flux = ',total(bzi_cal[where(ondisk)]), ' Uncalibrated net flux = ',total(bzi0[where(ondisk)])
	print,'Calibrated AR net flux = ',total(bzi_cal[ar_px]), ' Uncalibrated AR net flux = ',total(bzi0[ar_px])
	print,'Calibrated QS net flux = ',total(bzi_cal[qs_px]), ' Uncalibrated QS net flux = ',total(bzi0[qs_px])

	if(n_elements(plot_dir) eq 1) then begin
		if(n_elements(prange) eq 0) then prange = max(abs([bzi0*ondisk,bzi_cal*ondisk]))
		set_plot,'z'
		device,set_resolution=round([3.9*nx,3.9*ny])
		;device,/encapsulated,bits=8,xsize=18,ysize=6,filename=plot_dir+file_basename(bzifile,'.fits.gz')+'.eps'
		!p.multi=[0,2,2]
		plot_image,bzi0*ondisk,min=-prange,max=prange,title='Original GONG image',charsize=2.5
		plot_image,bzi_cal*ondisk,min=-prange,max=prange,title='Curve-corrected GONG image',charsize=2.5
		cal_ratio = abs(bzi_cal)/abs(bzi)
		plot_image,cal_ratio,title='Ratio of corrected to uncorrected',min=0,max=2,charsize=2.5
		plot_image,qsflg,title='QS flags',charsize=2.5
		write_png,plot_dir+file_basename(bzifile,'.fits.gz')+'.png',tvrd()
		;device,/close
		xl = min(xa[where(ondisk)])
		xh = max(xa[where(ondisk)])
		yl = min(ya[where(ondisk)])
		yh = max(ya[where(ondisk)])
		xl += floor(0.5*(xh-xl))
		nx2 = (xh-xl+1)
		ny2 = (yh-yl+1)
		device,set_resolution=round([2*nx2,ny2])
		img1out = bzi_cal*ondisk+prange*(ondisk eq 0)
		img2out = bzi0*ondisk+prange*(ondisk eq 0)
		imgout = ([img2out[xl:xh,yl:yh],img1out[xl:xh,yl:yh]] < prange > (-prange))
		imgout[(nx2-1):nx2,*] = -prange
		;imgout = ([[bzi_cal[*,yl:yh]*ondisk[*,yl:yh]+prange*(ondisk[*,yl:yh] eq 0)],[bzi0[*,yl:yh]*ondisk[*,yl:yh]+prange*(ondisk[*,yl:yh] eq 0)]] < prange) > (-prange)
		;imgout[*,(yh-yl):(yh-yl+1)] = -prange
		tvscl,imgout
		xyouts,nx2-5,ny2-30,'Uncalibrated',charsize=2,charthick=2,/device,color=0,alignment=1.0
		xyouts,nx2-5,ny2-65,'GONG',charsize=2,charthick=2,/device,color=0,alignment=1.0
		xyouts,2*nx2-5,ny2-30,'Calibrated',charsize=2,charthick=2,/device,color=0,alignment=1.0
		xyouts,2*nx2-5,ny2-65,'GONG',charsize=2,charthick=2,/device,color=0,alignment=1.0
		xyouts,nx2-5,5,datestr,charsize=2,charthick=2,/device,color=0,alignment=1.0
		xyouts,2*nx2-5,5,timestr,charsize=2,charthick=2,/device,color=0,alignment=1.0
		;xyouts,5,yh-yl+1-35,'Calibrated GONG',charsize=3,charthick=3,/device,color=0
		;xyouts,5,2*(yh-yl+1)-35,'Uncalibrated GONG',charsize=3,charthick=3,/device,color=0
		;xyouts,5,5,file_basename(bzifile,'.fits.gz'),charsize=3,charthick=3,/device,color=0
		write_png,plot_dir+file_basename(bzifile,'.fits.gz')+'_2panel.png',tvrd()
		set_plot,'x'
		
		
	endif
	
	bzi_cal[zerovals] = 0.0

	return,bzi_cal
	
end
