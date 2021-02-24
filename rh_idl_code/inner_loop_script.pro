.compile /data/jplowman/rh_idl_code_jp/refactor_inclined_ray.pro
restore,'inner_loop_script_inputs.sav'

print, rhdir
cd, rhdir
@initrh_noprompt
.r readall
tauvalue = 0.1
.r /data/jplowman/rh_idl_code_jp/gettauonesurf.pro
zcont01 = zcont
zcore01 = zcore
tauvalue = 0.3
.r /data/jplowman/rh_idl_code_jp/gettauonesurf.pro
zcont03 = zcont
zcore03 = zcore
tauvalue = 0.7
.r /data/jplowman/rh_idl_code_jp/gettauonesurf.pro
zcont07 = zcont
zcore07 = zcore
tauvalue = 2.0
.r /data/jplowman/rh_idl_code_jp/gettauonesurf.pro
zcont20 = zcont
zcore20 = zcore
tauvalue = 3.0
.r /data/jplowman/rh_idl_code_jp/gettauonesurf.pro
zcont30 = zcont
zcore30 = zcore
tauvalue = 1.0
.r /data/jplowman/rh_idl_code_jp/gettauonesurf.pro

if(n_elements(rayfilename) gt 0) then if(rayfilename ne '') then begin &$
	ray = readray(rayfilename) &$
	spectrum.i = ray.i &$
	spectrum.stokes_q = ray.stokes_q &$
	spectrum.stokes_u = ray.stokes_u &$
	spectrum.stokes_v = ray.stokes_v &$
;	s2 = ray_interpolate_3D_rh_func(geometry,ray.muy,ray.s,/yaxis, xray, yray, zray, raylength=raylength, xdist=xdist, zdist=zdist) &$
;	chi2 = ray_interpolate_3D_rh_func(geometry,ray.muy,ray.chi,/yaxis, xray, yray, zray, raylength=raylength, xdist=xdist, zdist=zdist) &$
	ray = refactor_inclined_ray(geometry, ray, /yaxis) &$
endif

save, atmos, spectrum, geometry, zcont, zcore, zcont01, zcore01, zcont03, zcore03, zcont07, zcore07, zcont20, zcore20, zcont30, zcore30, ray, filename='spectrum.sav'