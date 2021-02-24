KM_TO_M = 1.0E3

@files.common
@geometry.common
@atmos.common
@opacity.common
@spectrum.common

  metalFile      = "metals.out"
  moleculeFile   = "molecules.out"
  opacFile       = 'opacity.out'
  backgroundFile = 'background.dat'

  result = readInput('input.out')
  result = readgeometry('geometry.out')
  result = readatmos('atmos.out')
  result = readspectrum('spectrum.out')
  result = openOpacity(opacFile)

  minpos = intarr(geometry.Nx, geometry.Ny)
  lambda = findgen(spectrum.Nspect)
  FOR m=0, geometry.Ny-1 DO BEGIN
    FOR n=0, geometry.Nx-1 DO BEGIN
      minpos[n, m] = round(parmin(lambda, reform(spectrum.I[n, m, *])))
    ENDFOR
  ENDFOR
  themins = unique_elements(minpos)
;  themins = minpos[uniq(minpos,sort(minpos))]
  Nmin = n_elements(themins)
  print, "  Found ", Nmin, " possible minimum wavelength positions"

  zeff   = fltarr(geometry.Nx, geometry.Ny, 2)
  zcont  = fltarr(geometry.Nx, geometry.Ny)
  zcore  = fltarr(geometry.Nx, geometry.Ny)

  if(n_elements(tauvalue) ne 0) then print,tauvalue
  
  readOpacity, 0, 0
  FOR m=0, geometry.Ny-1 DO $
   zeff[*, m, 0] = tauone(geometry.z, reform(chi_c[*, m, *]),TAUVALUE=tauvalue)

  FOR n=0, Nmin-1 DO BEGIN
    loc = where(minpos EQ themins[n], Nloc)
    readOpacity, themins[n], 0
   
    FOR l=0, Nloc-1 DO BEGIN
      ix = loc[l] MOD geometry.Nx
      iy = loc[l] / geometry.Nx
      zeff[ix, iy, 1] = tauone(geometry.z, reform(chi_c[ix, iy, *]),TAUVALUE=tauvalue)
    ENDFOR
  ENDFOR

  FOR m=0, geometry.Ny-1 DO BEGIN
    FOR n=0, geometry.Nx-1 DO BEGIN
      zcont[n, m] = $
       interpolate(geometry.z, zeff[n, m, 0], CUBIC=-0.5) / KM_TO_M
      zcore[n, m] = $
       interpolate(geometry.z, zeff[n, m, 1], CUBIC=-0.5) / KM_TO_M
    ENDFOR
  ENDFOR

END
