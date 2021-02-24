function fluxcons_binfit, xa, ya, bincens, nbins=nbins, binsz_thold=binsz_thold, bins_in=bins_in, unsigned=unsigned

	if(n_elements(binsz_thold) eq 0) then binsz_thold = 400

	ntot = n_elements(xa)
	xmin = min(xa)-1.0e-5
	xmax = max(xa)+1.0e-5

	xa_srta = sort(xa)
	if(n_elements(bins_in) gt 0) then begin
		nbins = n_elements(bins_in)
		bincens = bins_in
		binbounds = [xmin,0.5*(bincens[0:nbins-2]+bincens[1:nbins-1]),xmax]
		npts = nbins+1
	endif else begin
		npts = nbins+1
		bomin = xa[xa_srta[binsz_thold-1]]
		bomax = xa[xa_srta[ntot-binsz_thold]]
		if(keyword_set(unsigned)) then begin
			binrange = min(abs([bomin,bomax]))
			binbounds = [xmin,2.0*binrange*dindgen(npts-2)/(npts-3.0) - binrange,xmax]
		endif else begin
			binbounds = [xmin,bomin+(bomax-bomin)*dindgen(npts-2)/(npts-3.0),xmax]
		endelse
		bincens = 0.5*(binbounds[0:npts-2]+binbounds[1:npts-1])
	endelse

	binlo_m = binbounds[0:nbins-1]
	binhi_m = bincens
	binlo_p = bincens
	binhi_p = binbounds[1:nbins]

	solmat = dblarr(nbins,nbins)
	solvec = dblarr(nbins)

	for i=0,nbins-1 do begin
		wvals_m = where(xa le binhi_m[i] and xa gt binlo_m[i],countm)
		wvals_p = where(xa le binhi_p[i] and xa gt binlo_p[i],countp)
		if(countm gt 0) then begin
			solvec[i] += total(ya[wvals_m],/double)
			if(i eq 0) then	begin
				solmat[i,i] += total((bincens[i+1]-xa[wvals_m])/(bincens[i+1]-bincens[i]),/double)
				solmat[i,i+1] += total((xa[wvals_m]-bincens[i])/(bincens[i+1]-bincens[i]),/double)
			endif
			if(i eq nbins-1) then begin
				solmat[i,i] += total((xa[wvals_m]-bincens[i-1])/(bincens[i]-bincens[i-1]),/double)
				solmat[i,i-1] += total((bincens[i]-xa[wvals_m])/(bincens[i]-bincens[i-1]),/double)
			endif
			if(i gt 0 and i lt nbins-1) then begin
				solmat[i,i] += total((xa[wvals_m]-bincens[i-1])/(bincens[i]-bincens[i-1]),/double)
				solmat[i,i-1] += total((bincens[i]-xa[wvals_m])/(bincens[i]-bincens[i-1]),/double)
			endif
		endif
		if(countp gt 0) then begin
			solvec[i] += total(ya[wvals_p],/double)
			if(i eq 0) then	begin
				solmat[i,i] += total((bincens[i+1]-xa[wvals_p])/(bincens[i+1]-bincens[i]),/double)
				solmat[i,i+1] += total((xa[wvals_p]-bincens[i])/(bincens[i+1]-bincens[i]),/double)
			endif
			if(i eq nbins-1) then begin
				solmat[i,i] += total((xa[wvals_p]-bincens[i-1])/(bincens[i]-bincens[i-1]),/double)
				solmat[i,i-1] += total((bincens[i]-xa[wvals_p])/(bincens[i]-bincens[i-1]),/double)
			endif
			if(i gt 0 and i lt nbins-1) then begin
				solmat[i,i] += total((bincens[i+1]-xa[wvals_p])/(bincens[i+1]-bincens[i]),/double)
				solmat[i,i+1] += total((xa[wvals_p]-bincens[i])/(bincens[i+1]-bincens[i]),/double)
			endif
		endif
	endfor

	LUDC, solmat, INDEX, /column, /double
	solution = lusol(solmat,index,solvec, /column, /double)

	ya_inv = interpol(solution,bincens,xa)
	
	yainv_binned = dblarr(nbins)

	for i=0,nbins-1 do begin
		wvals_m = where(xa le binhi_m[i] and xa gt binlo_m[i],countm)
		wvals_p = where(xa le binhi_p[i] and xa gt binlo_p[i],countp)
		if(countm gt 0) then yainv_binned[i] += total(ya_inv[wvals_m],/double)
		if(countp gt 0) then yainv_binned[i] += total(ya_inv[wvals_p],/double)
	endfor

	print,binbounds

	return, solution

end

function fluxcons_binfit_v0, xa, ya, bincens, nbins=nbins, binsz_thold=binsz_thold

	npts = nbins+1
	if(n_elements(binsz_thold) eq 0) then binsz_thold = 1600

	ntot = n_elements(xa)
	xmin = min(xa)-1.0e-5
	xmax = max(xa)+1.0e-5

	xa_srta = sort(xa)
	ixsrta = xa_srta[round(ntot*lindgen(npts)/(npts-1.0))]
	bomin = xa[xa_srta[binsz_thold-1]]
	bomax = xa[xa_srta[ntot-binsz_thold]]
	binrange = min(abs([bomin,bomax]))
	binbounds = [xmin,2.0*binrange*dindgen(npts-2)/(npts-3.0) - binrange,xmax]
	bincens = 0.5*(binbounds[0:npts-2]+binbounds[1:npts-1])
	bincens[where(abs(bincens) lt 1.0d-5)] = 0.0
	
;	binlo_p = [xmin,bincens[1:nbins-1]]
;	binhi_p = binbounds[1:nbins]
;	binlo_m = binhi_p
;	binhi_m = [bincens[1:nbins-1],xmax]

	binlo_p = [xmin,binbounds[1:nbins-1]]
	binhi_p = bincens
	binlo_m = binhi_p
	binhi_m = [binbounds[1:nbins-1],xmax]

	solmat = dblarr(nbins,nbins)
	solvec = dblarr(nbins)

	for i=0,nbins-1 do begin
		if(i gt 0) then begin
			wvals_m = where(xa lt binhi_m[i-1] and xa ge binlo_m[i-1])
			solmat[i,i-1] = total((bincens[i]-xa[wvals_m])/(bincens[i]-bincens[i-1]))
			solmat[i,i] += total((xa[wvals_m]-bincens[i-1])/(bincens[i]-bincens[i-1]))
			solvec[i] += total(ya[wvals_m])
		endif
		if(i lt nbins-1) then begin
			wvals_p = where(xa lt binhi_p[i] and xa ge binlo_p[i])
			solmat[i,i+1] = total((xa[wvals_p]-bincens[i])/(bincens[i+1]-bincens[i]))
			solmat[i,i] += total((bincens[i+1]-xa[wvals_p])/(bincens[i+1]-bincens[i]))
			solvec[i] += total(ya[wvals_p])
		endif
	endfor

	;LUDC, solmat, INDEX, /column, /double
	;solution = lusol(solmat,index,solvec, /column, /double)

	SVDC, solmat, W, U, V, /column, /double
	; Compute the solution and print the result:
	wlo = where(abs(w/max(w)) lt 1.0e-5, count)
	if(count gt 0) then w[wlo] = 0.0
	solution = SVSOL(U, W, V, solvec, /column, /double)

	print,xmin,xmax,min(ya),max(ya)
	print,'bincens: ',bincens
	print,'binlo_p: ',binlo_p
	print,'binhi_p: ',binhi_p
	print,'binlo_m: ',binlo_m
	print,'binhi_m: ',binhi_m

	solution[where(bincens eq 0.0)] = 0.0

	return, solution

end


function fluxcons_binfit_nplusone, xa, ya, bincens, npts=npts, binsz_thold=binsz_thold

	nbins = npts-1
	if(n_elements(binsz_thold) eq 0) then binsz_thold = 1600

	ntot = n_elements(xa)
	xmin = min(xa)
	xmax = max(xa)

	xa_srta = sort(xa)
	ixsrta = xa_srta[round(ntot*lindgen(npts)/(npts-1.0))]
	bomin = xa[xa_srta[binsz_thold-1]]
	bomax = xa[xa_srta[ntot-binsz_thold]]
	binrange = min(abs([bomin,bomax]))
	bincens = [xmin,2.0*binrange*dindgen(npts-2)/(npts-3.0) - binrange,xmax]
	bincens[where(abs(bincens) lt 1.0d-5)] = 0.0
	
	binlo_p = bincens[0:npts-2]
	binhi_p = 0.5*(bincens[0:npts-2]+bincens[1:npts-1])
	binlo_m = binhi_p
	binhi_m = bincens[1:npts-1]

	solmat = dblarr(npts,npts)
	solvec = dblarr(npts)

	for i=0,nbins do begin
		if(i gt 0) then begin
			wvals_m = where(xa lt binhi_m[i-1] and xa ge binlo_m[i-1])
			solmat[i,i-1] = total((bincens[i]-xa[wvals_m])/(bincens[i]-bincens[i-1]))
			solmat[i,i] += total((xa[wvals_m]-bincens[i-1])/(bincens[i]-bincens[i-1]))
			solvec[i] += total(ya[wvals_m])
		endif
		if(i lt nbins) then begin
			wvals_p = where(xa lt binhi_p[i] and xa ge binlo_p[i])
			solmat[i,i+1] = total((xa[wvals_p]-bincens[i])/(bincens[i+1]-bincens[i]))
			solmat[i,i] += total((bincens[i+1]-xa[wvals_p])/(bincens[i+1]-bincens[i]))
			solvec[i] += total(ya[wvals_p])
		endif
	endfor

	LUDC, solmat, INDEX, /column, /double
	solution = lusol(solmat,index,solvec, /column, /double)

	solution[where(bincens eq 0.0)] = 0.0

	return, solution

end

function gong_curve_fluxcons_binfit, bz_model_in, bz_obs_in, errs, bincen, npts=npts, unsigned=unsigned, xchang=xchang, coeffs=coeffs, mincount=mincount, domean=domean, bomin=bomin, bomax=bomax, reweight_prior=reweight_prior, binsz_thold=binsz_thold, bins_in=bins_in

	if(n_elements(mincount) ne 1) then mincount = 100
	if(n_elements(binsz_thold) eq 0) then binsz_thold = 400

	if(keyword_set(xchang)) then begin
		bz_model=bz_obs_in
		bz_obs=bz_model_in
	endif else begin
		bz_model=bz_model_in
		bz_obs=bz_obs_in
	endelse
	
	if(keyword_set(unsigned)) then begin
		bz_model = [-bz_model,bz_model]
		bz_obs = [-bz_obs,bz_obs]
	endif

	binvals = fluxcons_binfit(bz_obs, bz_model, bincen, nbins=npts, binsz_thold=binsz_thold, unsigned=unsigned, bins_in=bins_in)
		
	if(keyword_set(xchang)) then begin
		retval = bincen
		bincen = binvals
	endif else begin
		retval=binvals
	endelse

	if(keyword_set(unsigned)) then begin
		zerobins = where(abs(bincen) lt 1.0d-5,count)
		if(count gt 0) then begin
			bincen[zerobins] = 0.0
			retval[zerobins] = 0.0
		endif
	endif

	coeffs=retval
	
	return, retval
	
end
