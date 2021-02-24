function get_ar_thold, gongitot, bzobs, dlarge=dlarge, dsmall=dsmall, bmax_qs=bmax_qs, rerode=rerode, rdilate=rdilate, it_thold=it_thold

	if(n_elements(dlarge) eq 0) then dlarge=50
	if(n_elements(dsmall) eq 0) then dsmall=3
	if(n_elements(it_thold) eq 0) then it_thold = 0.97
	if(n_elements(bmax_qs) eq 0) then bmax_qs = 300

	izi_ms_large = median(gongitot,dlarge)
	izi_ms_small = median(gongitot,dsmall)


	rerode=4
	nxerode = 2*rerode+1
	nyerode = 2*rerode+1
	xaerode = (findgen(nxerode)-rerode)#(1+fltarr(nyerode))
	yaerode = (1+fltarr(nxerode))#(findgen(nyerode-rerode))
	yaerode = (1+fltarr(nxerode))#(findgen(nyerode)-rerode)
	s_erode = sqrt(xaerode^2+yaerode^2) lt rerode

	rdilate=11
	nxdilate = 2*rdilate+1
	nydilate = 2*rdilate+1
	xadilate = (findgen(nxdilate)-rdilate)#(1+fltarr(nydilate))
	yadilate = (1+fltarr(nxdilate))#(findgen(nydilate-rdilate))
	yadilate = (1+fltarr(nxdilate))#(findgen(nydilate)-rdilate)
	s_dilate = sqrt(xadilate^2+yadilate^2) lt rdilate
	gongi_thold = dilate(erode(gongitot le it_thold*izi_ms_large or abs(bzobs) gt bmax_qs,s_erode),s_dilate)

	return,gongi_thold

end
