pro make_rh_Bfile, b, inc, azi, filename

	GAUSS_TO_TESLA = 1.0d-4
	openw,unit,/xdr,filename,/get_lun
	writeu,unit,b*GAUSS_TO_TESLA,inc,azi
	free_lun,unit
	
end