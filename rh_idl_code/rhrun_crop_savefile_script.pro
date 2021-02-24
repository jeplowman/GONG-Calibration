resize_fac = 300
nx = 64
ny = 128

run_label = 'run_rhsc3d_ARQS'
infilename = '/data/jplowman/rh_runs/run_rhsc3d_ARQS_60_slice/'+run_label+'_spectrum.sav'
outfilename = '/data/jplowman/rh_runs/run_rhsc3d_ARQS_60_slice/'+run_label+'crop_spectrum.sav'

restore,infilename

atmos_crop = {nhydr:atmos.nhydr, nelem:atmos.nelem, moving:atmos.moving, t:atmos.t[0:nx-1,0:ny-1,*], $
	n_elec:atmos.n_elec[0:nx-1,0:ny-1,*], vturb:atmos.vturb[0:nx-1,0:ny-1,*], nh:atmos.nh[0:nx-1,0:ny-1,*,*], $
	id:atmos.id, elements:atmos.elements, B:atmos.b[0:nx-1,0:ny-1,*], gamma_b:atmos.gamma_b[0:nx-1,0:ny-1,*], $
	chi_b:atmos.chi_b[0:nx-1,0:ny-1,*], backgrflags:atmos.backgrflags, backgrrecno:atmos.backgrrecno}

atmos_geometry_crop = {nrays:atmos_geometry.nrays, nx:nx, ny:ny, nz:atmos_geometry.nz, angleset:atmos_geometry.angleset, $
	xmu:atmos_geometry.xmu, ymu:atmos_geometry.ymu, wmu:atmos_geometry.wmu, dx:atmos_geometry.dx*resize_fac, $
	dy:atmos_geometry.dy*resize_fac, z:atmos_geometry.z, vx:atmos_geometry.vx[0:nx-1,0:ny-1,*], $
	vy:atmos_geometry.vy[0:nx-1,0:ny-1,*], vz:atmos_geometry.vz[0:nx-1,0:ny-1,*]}

geometry_crop = {nrays:geometry.nrays, nx:nx, ny:ny, nz:geometry.nz, angleset:geometry.angleset, $
	xmu:geometry.xmu, ymu:geometry.ymu, wmu:geometry.wmu, dx:geometry.dx*resize_fac, dy:geometry.dy*resize_fac, $
	z:geometry.z, vx:geometry.vx[0:nx-1,0:ny-1,*], vy:geometry.vy[0:nx-1,0:ny-1,*], $
	vz:geometry.vz[0:nx-1,0:ny-1,*]}

ray_crop = {mux:ray.mux, muy:ray.muy, I:ray.I[0:nx-1,0:ny-1,*], nspect:ray.nspect, chi:ray.chi[0:nx-1,0:ny-1,*,*], $
	s:ray.s[0:nx-1,0:ny-1,*,*], chi2:ray.chi2[0:nx-1,0:ny-1,*,*], s2:ray.s2[0:nx-1,0:ny-1,*,*], $
	stokes_q:ray.stokes_q[0:nx-1,0:ny-1,*], stokes_u:ray.stokes_u[0:nx-1,0:ny-1,*], $
	stokes_v:ray.stokes_v[0:nx-1,0:ny-1,*], xray:ray.xray[0:nx-1,0:ny-1,*], yray:ray.yray[0:nx-1,0:ny-1,*], $
	zray:ray.zray[0:nx-1,0:ny-1,*], raylength:ray.raylength, nlos:ray.nlos}

spectrum_crop = {nspect:spectrum.nspect, lambda:spectrum.lambda, I:spectrum.I[0:nx-1,0:ny-1,*], $
	vacuum_to_air:spectrum.vacuum_to_air, air_limit:spectrum.air_limit, stokes_q:spectrum.stokes_q[0:nx-1,0:ny-1,*], $
	stokes_u:spectrum.stokes_u[0:nx-1,0:ny-1,*], stokes_v:spectrum.stokes_v[0:nx-1,0:ny-1,*], as_rn:spectrum.as_rn}

spectrum=spectrum_crop
ray=ray_crop
geometry=geometry_crop
atmos_geometry=atmos_geometry_crop
atmos=atmos_crop

zcont = zcont[0:nx-1,0:ny-1]
zcont01 = zcont01[0:nx-1,0:ny-1]
zcont03 = zcont03[0:nx-1,0:ny-1]
zcont07 = zcont07[0:nx-1,0:ny-1]
zcont20 = zcont20[0:nx-1,0:ny-1]
zcont30 = zcont30[0:nx-1,0:ny-1]

zcore = zcore[0:nx-1,0:ny-1]
zcore01 = zcore01[0:nx-1,0:ny-1]
zcore03 = zcore03[0:nx-1,0:ny-1]
zcore07 = zcore07[0:nx-1,0:ny-1]
zcore20 = zcore20[0:nx-1,0:ny-1]
zcore30 = zcore30[0:nx-1,0:ny-1]

save, atmos, spectrum, geometry, atmos_geometry, zcont, zcore, zcont01, zcore01, zcont03, zcore03, zcont07, zcore07, zcont20, zcore20, zcont30, zcore30, ray, filename=outfilename
