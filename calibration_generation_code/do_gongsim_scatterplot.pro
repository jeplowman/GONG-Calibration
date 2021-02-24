pro do_gongsim_scatterplot, bmin, bmax, bz_obstitles, bz_modtitles, bz_obs, bz_model, obsfit, modfit, label, position=position, noerase=noerase

	plot,[bmin,bmax],[bmin,bmax],xrange=[bmin,bmax],yrange=[bmin,bmax], xtitle=bz_obstitles, ytitle=bz_modtitles, xstyle=1,ystyle=1,thick=2,title=label, charsize=1.0,position=position,noerase=noerase
	oplot,[0,0],[bmin,bmax],linestyle=2
	oplot,[bmin,bmax],[0,0],linestyle=2
	oplot,bz_obs,bz_model, psym=3, color=128
	oplot,-bz_obs,-bz_model, psym=3, color=128
	oplot,obsfit,modfit,linestyle=0,thick=8
	oplot,obsfit,modfit,linestyle=2,thick=4,color=255
	if(n_elements(position) eq 2) then begin
		ssw_legend,['Calibration Curve','y=x line'],linestyle=[0,0],thick=[8,2],charsize=1.0, charthick=2, position=[position[1],position[3]]-0.05,/noerase, clear=255
		ssw_legend,['Calibration Curve','y=x line'],linestyle=[2,0],thick=[4,2],color=[255,0], charthick=2, charsize=1.0, position=[position[1],position[3]]-0.05,/noerase
	endif else begin
		ssw_legend,['Calibration Curve','y=x line'],linestyle=[0,0],thick=[8,2],charsize=1.0,charthick=2, clear=255
		ssw_legend,['Calibration Curve','y=x line'],linestyle=[2,0],thick=[4,2],color=[255,0],charsize=1.0,charthick=2
	endelse

end
