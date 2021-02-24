function cal_curve_fit_function, indeps, parms
	
	nparms = n_elements(parms)
	nindeps = n_elements(indeps)
		
	output = fltarr(nindeps,nparms+1)
		
	for i=0,nparms-1 do begin
		output[*,i+1] = indeps^i
		output[*,0] += parms[i]*output[*,i+1]
	endfor
	
	return, output
	
end

function gong_cal_curve_lmfit, bz_model, bz_obs, errs, xvals, nterms=nterms, guess=guess

	if(n_elements(nterms) eq 0) then nterms=8
	if(n_elements(guess) eq 0) then begin
		guess = fltarr(nterms)
		guess[1] = 1.0
	endif
	
	coefs = lmfit(bz_model,bz_obs,guess,measure_errors=errs,function_name = 'cal_curve_fit_function',/double)
	
	ymin = min(bz_model)
	ymax = max(bz_model)
	if(n_elements(nout) eq 0) then nout = 101.0
	yvals = ymin+(ymax-ymin)*findgen(nout)/(nout-1.0)
	
	fitvals = cal_curve_fit_function(yvals,coefs)
	
	xvals = fitvals[*,0]
	
	return,yvals
	
end

function gong_curve_linfit, bz_model_in, bz_obs_in, errs, xvals, nterms=nterms, nout=nout, coeffs=coeffs, zeroterm=zeroterm, xchang=xchang, fluxnorm=fluxnorm, unsigned=unsigned
	
	if(keyword_set(xchang)) then begin
		bz_model=bz_obs_in
		bz_obs=bz_model_in
	endif else begin
		bz_model=bz_model_in
		bz_obs=bz_obs_in
	endelse

	if(keyword_set(unsigned)) then begin
;		bz_model = abs(bz_model)
;		bz_obs = abs(bz_obs)
		bz_model = [-bz_model,bz_model]
		bz_obs = [-bz_obs,bz_obs]
		errs = [errs,errs]
	endif
	
	bomin = min(bz_obs)
	bomax = max(bz_obs)
	bmmin = -max(abs(bz_model))
	bmmax = max(abs(bz_model))
	bz_model_norm = 1; max(abs([bmmin,bmmax]))
	bz_model = bz_model/bz_model_norm
	
	if(n_elements(nterms) eq 0) then nterms = 5

	npts = n_elements(bz_model)
	
	mmat = dblarr(nterms,npts)
	bvec = dblarr(nterms)
	
	for i=0,nterms-1 do begin
		if(keyword_set(zeroterm)) then begin
			fac = (2*i-1) > 0
		endif else begin
			fac = 2*i+1
		endelse
		mmat[i,*] = bz_model^fac/errs
		bvec[i] = total(mmat[i,*]*bz_obs/errs)
	endfor
	
	amat = mmat#transpose(mmat)
	
	coeffs = invert(amat,/double)#bvec
	
	if(n_elements(nout) eq 0) then nout = 101
	yvals = (bmmin+(bmmax-bmmin)*findgen(nout)/(nout-1.0))/bz_model_norm
	xvals = fltarr(nout)
	for i=0,nterms-1 do begin
		if(keyword_set(zeroterm)) then begin
			fac = (2*i-1) > 0
		endif else begin
			fac = 2*i+1
		endelse
		xvals+=coeffs[i]*yvals^fac
	endfor
	
	if(keyword_set(fluxnorm)) then begin
		coeffs = total(bz_obs)/total(bz_model)
		xvals = coeffs*yvals
	endif
	
	if(keyword_set(xchang)) then begin
		yvals_temp = yvals
		yvals = xvals
		xvals = yvals_temp
	endif
	
	return,yvals*bz_model_norm
	
end
	
function gong_curve_binfit, bz_model_in, bz_obs_in, errs, bincen, npts=npts, unsigned=unsigned, xchang=xchang, coeffs=coeffs, mincount=mincount, domean=domean, bomin=bomin, bomax=bomax, reweight_prior=reweight_prior

	if(n_elements(mincount) ne 1) then mincount = 100

	if(keyword_set(xchang)) then begin
		bz_model=bz_obs_in
		bz_obs=bz_model_in
	endif else begin
		bz_model=bz_model_in
		bz_obs=bz_obs_in
	endelse
	
	if(keyword_set(unsigned)) then begin
;		bz_model = abs(bz_model)
;		bz_obs = abs(bz_obs)
		bz_model = [-bz_model,bz_model]
		bz_obs = [-bz_obs,bz_obs]
	endif

	ntot = n_elements(bz_obs)
	if(n_elements(npts) eq 0) then npts = round(0.5*n_elements(bz_obs)/mincount)
	npts = npts < floor(ntot/mincount)

	if(n_elements(bomin) eq 0) then bomin = min(bz_obs)
	if(n_elements(bomax) eq 0) then bomax = max(bz_obs)
	bmmin = min(bz_model)
	bmmax = max(bz_model)
	if(keyword_set(reweight_prior)) then begin
		bzm_srta = sort(bz_model)
		bzm_chist_x = bz_model[bzm_srta]
		bzm_chist_y = findgen(ntot)/(ntot-1.0)
		ibzma = bzm_srta[ntot*(lindgen(npts))/(npts-1)]
		bzmlin = bz_model[ibzma]
		bzm_chist = interpol(bzm_chist_y,bzm_chist_x,bzmlin)
		bzm_marg0 = fltarr(npts)
		bzm_marg0[1:npts-2] = (bzm_chist[2:npts-1]-bzm_chist[0:npts-3])/(bzmlin[2:npts-1]-bzmlin[0:npts-3])
		bzm_marg0[0] = (bzm_chist[1]-bzm_chist[0])/(bzmlin[1]-bzmlin[0])
		bzm_marg0[npts-1] = (bzm_chist[npts-1]-bzm_chist[npts-2])/(bzmlin[npts-1]-bzmlin[npts-2])
		bzm_marg = interpol(bzm_marg0,bzmlin,bz_model)
		bzm_prior = 1.0+0.0*bz_model; exp(-0.5*((bz_model)/stdev(bz_model))^2)/sqrt(2.0*!pi*stdev(bz_model)^2)
		weights = bzm_prior/bzm_marg
	endif else begin
		weights = 1.0+fltarr(ntot)
	endelse

;	message,'Breakpoint!'
	
	bzo_srta = sort(bz_obs)
	ibzoa = bzo_srta[ntot*lindgen(npts)/npts]
;	binlo = bz_obs[ibzoa]
;	binhi = bz_obs[[ibzoa[1:npts-1],bzo_srta[ntot-1]]]
	
	bomin = bz_obs[bzo_srta[1599]]
	bomax = bz_obs[bzo_srta[ntot-1600]]
	bomin -= (bomax-bomin)/npts
	bomax += (bomax-bomin)/npts

	;binlo = bomin+(bomax-bomin)*findgen(npts)/npts	
	;binhi = bomin+(bomax-bomin)*(findgen(npts)+1)/npts
	;bincen = 0.5*(binlo+binhi)

	binrange = min(abs([bomin,bomax]))
	bincen = 2.0*binrange*dindgen(npts)/(npts-1.0) - binrange
	bincen[where(abs(bincen) lt 1.0d-7)] = 0.0
	binlo = bincen-(bincen[1]-bincen[0])*0.5
	binhi = bincen+(bincen[1]-bincen[0])*0.5

	binvals = fltarr(npts)
	missing = intarr(npts)
		
	for i=0,npts-1 do begin
		wvals = where(bz_obs ge binlo[i] and bz_obs lt binhi[i], count)
		if(i eq 0) then wvals = where(bz_obs lt binhi[i],count)
		if(i eq npts-1) then wvals = where(bz_obs ge binlo[i],count)
;		if(count gt mincount) then begin
			
			if(keyword_set(domean) eq 0) then binvals[i] = median(bz_model[wvals])
;			if(keyword_set(domean)) then binvals[i] = mean(bz_model[wvals])
;			if(keyword_set(domean)) then binvals[i] = total(weights[wvals]*bz_model[wvals])/total(weights[wvals])
			if(keyword_set(domean)) then binvals[i] = bincen[i]*total(weights[wvals]*bz_model[wvals])/total(weights[wvals]*bz_obs[wvals])
			;if(keyword_set(domean) eq 0) then bincen[i] = median(bz_obs[wvals])
			;if(keyword_set(domean)) then bincen[i] = mean(bz_obs[wvals])			
;		endif else begin
;			missing[i] = 1
;		endelse
	endfor
	
	wgood = where(missing eq 0)
	binvals = binvals[wgood]
	bincen = bincen[wgood]
	
;	binvals[wmiss] = interpol(binvals[wgood],bincen[wgood],bincen[wmiss])
	
;	if(keyword_set(unsigned)) then begin
;		binvals = [-reverse(binvals),binvals]
;		bincen = [-reverse(bincen),bincen]
;	endif
	
	if(keyword_set(xchang)) then begin
		retval = bincen
		bincen = binvals
	endif else begin
		retval=binvals
	endelse
	
	coeffs=retval
	
	return, retval
	
end

; This doesn't work - will require binning the data to a uniform grid...
function gong_cal_curve_savgol, bz_model, bz_obs, xvals, nright=nright, nleft=nleft, order=order, degree=degree, double=double, delta=delta, nout=nout


	bomin = min(bz_obs)
	bomax = max(bz_obs)
	bmmin = min(bz_model)
	bmmax = max(bz_model)

	npts = n_elements(bz_obs)
	if(n_elements(order) eq 0) then order = 0
	if(n_elements(degree) eq 0) then degree = 3

	if(n_elements(nright) eq 0 or n_elements(nleft) eq 0) then begin
		if(n_elements(delta) eq 0) then delta = 100.0
		nright = npts*delta/(bomax-bomin)
		nleft = npts*delta/(bomax-bomin)
	endif
	
	filter = savgol(nleft,nright,order,degree,double=double)
	bzo_smooth = CONVOL(bz_obs, filter, /EDGE_TRUNCATE)
	
	return,bzo_smooth
	
end
