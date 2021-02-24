.compile get_contrib_function_depth.pro
.compile make_contrib_function_plot.pro

idlcode_dir = '/data/jplowman/rh_idl_code_jp'
rhrun_dir = '/mnt/usb-WD_easystore_25FB_3753474C37303543-0:0-part1/jplowman/rh_runs/runs_rhsc3d_ARQS_6768_2/run_rhsc3d_ARQS_15_0_0/'
plotfile_vertical = '/data/jplowman/GONG_Simulator/plots/GONG_simulator_paper1_plots/contrib_plot_0degree.eps'

anglestr0 = strtrim(string(70),2)
rayfile_inclined = rhrun_dir+'spectrum_0.00_0.94'
plotfile_inclined = '/data/jplowman/GONG_Simulator/plots/GONG_simulator_paper1_plots/contrib_plot_'+anglestr0+'degree.eps'
anglestr = anglestr0 + ' degree'

cd,rhrun_dir
@initrh
.r readall
cd,idlcode_dir

.compile /data/jplowman/rh_v2/idl/raytrace.pro
.compile /data/jplowman/rh_idl_code_jp/ray_interpolate_3D_rh_func.pro
.compile /data/jplowman/rh_idl_code_jp/refactor_inclined_ray.pro
.compile /data/jplowman/rh_idl_code_jp/get_vertical_ray.pro

if(is_struct(ray) eq 0) then ray = readray(rayfile_inclined)
if(is_struct(ray2) eq 0) then ray2 = refactor_inclined_ray(geometry,ray,/yaxis)
if(is_struct(vray) eq 0) then vray = get_vertical_ray(ray.nspect)

fcontrib = 0.5
if(keyword_set(raycoord_complete) eq 0) then begin &$
	depth05 = get_contrib_function_depth(ray2, fcontrib, xray=xray05, yray=yray05, zray=zray05, tauray=tauray, contray=contray) &$
	depth01 = get_contrib_function_depth(ray2, 0.1, xray=xray01, zray=zray01) &$
	depth09 = get_contrib_function_depth(ray2, 0.9, xray=xray09, zray=zray09) &$
	depthtau1 = get_contrib_function_depth(ray2, 1.0, xray=xraytau1, zray=zraytau1,/taudepth) &$
	vdepth05 = get_contrib_function_depth(vray, fcontrib, yray=vxray05, xray=vyray05, zray=vzray05, tauray=vtauray, contray=vcontray, /taudepth) &$
	vdepth01 = get_contrib_function_depth(vray, 0.1, yray=vxray01, xray=vyray01, zray=vzray01) &$
	vdepth09 = get_contrib_function_depth(vray, 0.9, yray=vxray09, xray=vyray09, zray=vzray09) &$
	vdepthtau1 = get_contrib_function_depth(vray, 1.0, yray=vxraytau1, zray=vzraytau1,/taudepth) &$
	raycoord_complete = 1 &$
endif

set_plot,'ps'
!p.font=0
device,/times,bits=8,/inches,xsize=6,ysize=4
!p.multi=[0,1,3]
yr = [100,400]

device,/times,bits=8,filename=plotfile_inclined,/encapsulated
make_contrib_function_plot, atmos, geometry, ray2, spectrum, [0], [1,2], xraytau1, zraytau1, xray01, zray01, xray09, zray09, slice=slice, min=min, max=max, linestyles=[0,1,2], yr=yr, titlestr = anglestr+' contribution, near line core'
make_contrib_function_plot, atmos, geometry, ray2, spectrum, [0], [4,5], xraytau1, zraytau1, xray01, zray01, xray09, zray09, slice=slice, min=min, max=max, linestyles=[0,1,2], yr=yr, titlestr = anglestr+' contribution, near line wings'
make_contrib_function_plot, atmos, geometry, ray2, spectrum, [0], [0], xraytau1, zraytau1, xray01, zray01, xray09, zray09, slice=slice, min=min, max=max, linestyles=[0,2], yr=yr, titlestr = anglestr+' contribution, continuum'
device,/close

device,/times,bits=8,filename=plotfile_vertical,/encapsulated
make_contrib_function_plot, atmos, geometry, vray, spectrum, [0], [1,2], vxraytau1, vzraytau1, vxray01, vzray01, vxray09, vzray09, slice=slice, min=min, max=max, linestyles=[0,1,2], yr=yr, titlestr = 'Vertical contribution, near line core'
make_contrib_function_plot, atmos, geometry, vray, spectrum, [0], [4,5], vxraytau1, vzraytau1, vxray01, vzray01, vxray09, vzray09, slice=slice, min=min, max=max, linestyles=[0,1,2], yr=yr, titlestr = 'Vertical contribution, near line wings'
make_contrib_function_plot, atmos, geometry, vray, spectrum, [0], [0], vxraytau1, vzraytau1, vxray01, vzray01, vxray09, vzray09, slice=slice, min=min, max=max, linestyles=[0,2], yr=yr, titlestr = 'Vertical contribution, continuum'
device,/close
