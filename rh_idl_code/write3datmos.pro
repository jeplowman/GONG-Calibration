pro write3datmos, atmos, outputfile, bfile=bfile

	tags = tag_names(atmos)
	openw, lun, outputfile, /GET_LUN, /XDR

;	nx = long(atmos.nx)
;	ny = long(atmos.ny)
;	nz = long(atmos.nz)
;	nhydr = long(atmos.nhydr)
;	print,nx,ny,nz,nhydr
	writeu, lun, atmos.nx, atmos.ny, atmos.nz, atmos.nhydr
	writeu, lun, long(atmos.boundary)
	writeu, lun, double(atmos.dx), double(atmos.dy), double(atmos.z)

	writeu, lun, double(atmos.t)
	writeu, lun, double(atmos.n_elec)     ;; n_elec
	writeu, lun, double(atmos.vturb)     ;; vturb

	;; Remember: second dimension in the original cube is the vertical

	writeu, lun, double(atmos.vx)
	writeu, lun, double(atmos.vy)
	writeu, lun, double(atmos.vz)

	writeu, lun, double(atmos.nh)

	free_lun, lun
	
	
	
	if(n_elements(bfile) eq 1 and total(tags eq 'B') gt 0) then begin
		openw,lun,bfile,/get_lun,/xdr
		bdata = dblarr(atmos.nx,atmos.ny,atmos.nz,3)
		bdata[*,*,*,0] = atmos.B
		bdata[*,*,*,1] = atmos.gamma
		bdata[*,*,*,2] = atmos.chi
		writeu,lun,bdata
		free_lun,lun
	endif

end
