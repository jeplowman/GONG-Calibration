;cd,'/mnt/usb-WD_easystore_25FB_3753474C37303543-0:0-part1/jplowman/rh_runs/run_rhsc3d_ARQS_15_0_0'
cd,'/mnt/usb-WD_easystore_25FB_3753474C37303543-0:0-part1/jplowman/rh_runs/runs_rhsc3d_ARQS_6768_2/run_rhsc3d_ARQS_15_0_0'

@initrh
.r readall

ray = readray('spectrum_0.00_0.94')

;s2 = ray_interpolate_3D_rh_func(geometry,ray.muy,ray.s[*,*,*,0],/yaxis, x_out, y_out, z_out, raylength=raylength, xdist=xdist, zdist=zdist)
;chi2 = ray_interpolate_3D_rh_func(geometry,ray.muy,ray.chi[*,*,*,0],/yaxis, x_out, y_out, z_out, raylength=raylength, xdist=xdist, zdist=zdist)


.compile /data/jplowman/rh_v2/idl/raytrace.pro
.compile /data/jplowman/rh_idl_code_jp/ray_interpolate_3D_rh_func.pro
.compile /data/jplowman/rh_idl_code_jp/refactor_inclined_ray.pro
.compile /data/jplowman/rh_idl_code_jp/get_vertical_ray.pro


ray2 = refactor_inclined_ray(geometry,ray,/yaxis)
vray = get_vertical_ray(ray.nspect)
vray2 = get_vertical_ray([19,26])

;KM_TO_M = 1.0E3

;@files.common
;@geometry.common
;@atmos.common
;@opacity.common
;@spectrum.common

;  metalFile      = "metals.out"
;  moleculeFile   = "molecules.out"
;  opacFile       = 'opacity.out'
;  backgroundFile = 'background.dat'

;  result = readInput('input.out')
;  result = readgeometry('geometry.out')
;  result = readatmos('atmos.out')
;  result = readspectrum('spectrum.out')
;  result = openOpacity(opacFile)

;  readOpacity, 0, 0


cd,'/data/jplowman/rh_idl_code_jp'
