cd,'/rhdata/jplowman/rh_v2/rhf1d/idl_codes_jp/'

.compile /rhdata/jplowman/rh_v2/rhf1d/idl_codes_jp/unique_elements.pro
.compile /rhdata/jplowman/rh_v2/rhf1d/idl_codes_jp/gettauonesurf.pro
.compile /rhdata/jplowman/rh_v2/rhf1d/idl_codes_jp/parmin.pro

cd,'/rhdata/jplowman/rh_v2/rhsc3d/run_rhsc3d_plage'

@initrh
.r readall
.r /rhdata/jplowman/rh_v2/rhf1d/idl_codes_jp/gettauonesurf.pro

save, atmos, spectrum, geometry, zcont, zcore, filename='rh_3dspectrum.sav'