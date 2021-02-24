pro rebin_atmos, infile, outfile, nx=nx, ny=ny, nz=nz, bfile=bfile, boutfile=boutfile, atmos=atmos

	atmos = read3datmos(infile,bfile=bfile)
	
	nx0 = atmos.nx
	ny0 = atmos.ny
	nz0 = atmos.nz
	
	if(n_elements(nx) eq 0) then nx = nx0
	if(n_elements(ny) eq 0) then ny = ny0
	if(n_elements(nz) eq 0) then nz = nz0
	
	nt = n_tags(atmos)
	names = tag_names(atmos)
	
	atmos2 = {nx:long(nx),ny:long(ny),nz:long(nz)}
	
	for i=0,nt-1 do begin
		if(names(i) ne 'NX' and names(i) ne 'NY' and names(i) ne 'NZ') then begin
			elem = atmos.(i)
			sze = [size(elem),-1] ; Augument the size array with a -1 so the if statement below doesn't crash for scalars...
			if(sze[0] eq 3 and sze[1] eq nx0 and sze[2] eq ny0 and sze[3] eq nz0) then begin
				if(names[i] eq 'T' or names[i] eq 'N_ELEC' or names[i] eq 'NH') then begin
					sanitize_array,elem,/nonzero
				endif else begin
					sanitize_array,elem
				endelse
				elem = rebin(elem,nx,ny,nz)
			endif else if(sze[3] eq nz0) then begin
				elem = rebin(elem,nz)
			endif
			atmos2 = create_struct(atmos2,names[i],elem)
		endif
	endfor
	
	atmos2.dx = atmos.dx*round(nx0/nx)
	atmos2.dy = atmos.dx*round(ny0/ny)

	write3datmos,atmos2,outfile,bfile=boutfile

;	write3datmos0, long(atmos2.nx), long(atmos2.ny), long(atmos2.nz), double(atmos2.dx),  double(atmos2.dy), double(atmos2.z[1]-atmos2.z[0]), atmos2.T, atmos2.vx, atmos2.vy, atmos2.vz, atmos2.nh, outfile, bfile=boutfile
	
end
	
	
	