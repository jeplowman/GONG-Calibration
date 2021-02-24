pro write3datmos0, N0, N2, NZNEW, DX, DY, DZ, T, vx, vy, vz, nHtot, outputfile, bfile=bfile


	openw, lun, outputFile, /GET_LUN, /XDR

	lBoundVal = [1L,  2L]

	NHYDR = 1L
	NX    = N0
	NY    = N2
	Nz    = NZNEW

	z  = reverse(dindgen(Nz) * DZ)

	writeu, lun, NX, NY, NZ, NHYDR
	writeu, lun, lBoundVal
	writeu, lun, DX, DY, z

	writeu, lun, double(T)
	writeu, lun, dblarr(Nx, Ny, Nz)     ;; n_elec
	writeu, lun, dblarr(Nx, Ny, Nz)     ;; vturb

	;; Remember: second dimension in the original cube is the vertical

	writeu, lun, double(vx)
	writeu, lun, double(vy)
	writeu, lun, double(vz)

	writeu, lun, double(nHtot)

	free_lun, lun

end