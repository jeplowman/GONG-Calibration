pro reldiff_plot, x1, y1_in, x2, y2_in, npts=npts, xrange=xrange, yrange=yrange, title=title, xtitle=xtitle, ytitle=ytitle, overplot=overplot, linestyle=linestyle, color=color, thick=thick

	xa = [x1,x2]
	xa = xa[UNIQ(xa, SORT(xa))]
	
	y1 = interpol(y1_in,x1,xa)
	y2 = interpol(y2_in,x2,xa)
	
	reldiff = 1-y2/y1
	w0 = where(y1 eq 0 or ~finite(y1) or ~finite(y2), count)
	if(count gt 0) then reldiff(w0) = 0.0
	
	if(keyword_set(overplot)) then begin
		plot, xa, reldiff, linestyle=linestyle, color=color, thick=thick
	endif else begin
		plot, xa, reldiff, xrange=xrange, yrange=yrange, title=title, xtitle=xtitle, ytitle=ytitle, linestyle=linestyle, color=color, thick=thick
	endelse
	
end