FUNCTION readvarcube, N0, N1, N2, filename

  dummy = fltarr(N0, N1, N2)

  openr, lun, /GET_LUN, filename
  readu, lun, dummy
  free_lun, lun

  return, dummy
END

pro readrempelcube, vfiles, bfiles, tfile, dfile, dx, dy, dz, n0, n1, n2, outfilehead, outputFile=outputFile, B_outputFile=B_outputFile, irradiated=irradiated, zero=zero, planck=planck, kmin=kmin, kmax=kmax, outrange_lo=outrange_lo, outrange_hi=outrange_hi

;; Convert Mathias' input files to RH input atmosphere

  if(n_elements(outrange_lo) eq 0) then outrange_lo = long([0,0,0])
  if(n_elements(outrange_hi) eq 0) then outrange_hi = long([n0,n2,n1])-1L

  if(n_elements(outfilehead) ne 0) then begin
	  outputFile   = outfilehead + '_' + strjoin(strtrim(string(outrange_lo),2),'_') + '_' + strjoin(strtrim(string(outrange_hi),2),'_')+'.atmos'
	  B_outputFile = outfilehead + '_' + strjoin(strtrim(string(outrange_lo),2),'_') + '_' + strjoin(strtrim(string(outrange_hi),2),'_')+'.B'
  endif

  if(n_elements(irradiated) eq 0) then IRRADIATED = 0L
  if(n_elements(zero) eq 0) then ZERO = 1L
  if(n_elements(PLANCK) eq 0) then PLANCK = 2L
  
  ;; Physical constants

  WGHT_PER_H   = 1.4271D0
  AVG_MOL_WGHT = 1.2981D0
  AMU          = 1.6605402D-27
  KBOLTZMANN   = 1.380658D-23

  ;; Unit conversions

  CM_TO_M   = 1.0D-2
  KM_TO_M   = 1.0D+3
  G_TO_KG   = 1.0D-3
  DYNE_TO_N = 1.0D-5
  GAUSS_TO_TESLA = 1.0D-4

  imin = outrange_lo[0]
  imax = outrange_hi[0]
  jmin = outrange_lo[1]
  jmax = outrange_hi[1]
  kmin = outrange_lo[2]
  kmax = outrange_hi[2]

  nxnew = imax - imin + 1
  nynew = jmax - jmin + 1
  NZNEW = KMAX - KMIN + 1

  v_files   = vfiles
  B_files   = bfiles
  T_file    = tfile
  dens_file = dfile

  Bx = (reverse(transpose(readvarcube(N0, N1, N2, B_files[0]), [0, 2, 1]), 3, $
                /OVERWRITE))[imin:imax, jmin:jmax, KMIN:KMAX] * GAUSS_TO_TESLA
  By = (reverse(transpose(readvarcube(N0, N1, N2, B_files[2]), [0, 2, 1]), 3, $
                /OVERWRITE))[imin:imax, jmin:jmax, KMIN:KMAX] * GAUSS_TO_TESLA
  Bz = (reverse(transpose(readvarcube(N0, N1, N2, B_files[1]), [0, 2, 1]), 3, $
                /OVERWRITE))[imin:imax, jmin:jmax, KMIN:KMAX] * GAUSS_TO_TESLA

  B     = sqrt(Bx^2 + By^2 + Bz^2)
  gamma = acos(Bz / B)
  chi   = atan(By, Bx)

  openw, lun, /GET_LUN, /XDR, B_outputFile
  writeu, lun, B, gamma, chi
  free_lun, lun
  
  Bx=0
  By=0
  Bz=0
  B=0
  gamma=0
  chi=0

  vx = (reverse(transpose(readvarcube(N0, N1, N2, v_files[0]), [0, 2, 1]), 3, $
                /OVERWRITE))[imin:imax, jmin:jmax, KMIN:KMAX] * (CM_TO_M / KM_TO_M)
  vy = (reverse(transpose(readvarcube(N0, N1, N2, v_files[2]), [0, 2, 1]), 3, $
                /OVERWRITE))[imin:imax, jmin:jmax, KMIN:KMAX] * (CM_TO_M / KM_TO_M)
  vz = (reverse(transpose(readvarcube(N0, N1, N2, v_files[1]), [0, 2, 1]), 3, $
                /OVERWRITE))[imin:imax, jmin:jmax, KMIN:KMAX] * (CM_TO_M / KM_TO_M)
 
  T = (reverse(transpose(readvarcube(N0, N1, N2, T_file[0]), [0, 2, 1]), 3, $
              /OVERWRITE))[imin:imax, jmin:jmax, KMIN:KMAX]
  rho = (reverse(transpose(readvarcube(N0, N1, N2, dens_file), $
                           [0, 2, 1]), 3, /OVERWRITE))[imin:imax, jmin:jmax, KMIN:KMAX]
  nHtot = (rho * G_TO_KG / CM_TO_M^3) / (WGHT_PER_H * AMU)


  openw, lun, outputFile, /GET_LUN, /XDR

  lBoundVal = [ZERO,  PLANCK]

  NHYDR = 1L
  NX    = nxnew
  NY    = Nynew
  Nz    = NZNEW

  z  = reverse(dindgen(Nz) * double(DZ))

  writeu, lun, long(NX), long(NY), long(NZ), long(NHYDR)
  writeu, lun, long(lBoundVal)
  writeu, lun, double(DX), double(DY), double(z)

  writeu, lun, double(T)
  writeu, lun, dblarr(Nx, Ny, Nz)     ;; n_elec
  writeu, lun, dblarr(Nx, Ny, Nz)     ;; vturb

  ;; Remember: second dimension in the original cube is the vertical

  writeu, lun, double(vx)
  writeu, lun, double(vy)
  writeu, lun, double(vz)

  writeu, lun, double(nHtot)

  free_lun, lun

END
