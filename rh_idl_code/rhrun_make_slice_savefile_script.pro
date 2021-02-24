.compile update_atmos.pro
.compile update_spectrum.pro
.compile update_ray.pro
.compile update_geometry.pro
.compile update_z.pro

infofilename = '/mnt/usb-WD_easystore_25FB_3753474C37303543-0:0-part1/jplowman/rh_runs/run_rhsc3d_ARQS_info.sav'
run_label = 'run_rhsc3d_ARQS'
outfilename = '/data/jplowman/rh_runs/run_rhsc3d_ARQS_60_slice/'+run_label+'_spectrum.sav'

rayfilename = 'spectrum_0.00_0.87'
ray_fac_dx = 1.0
ray_fac_dy = 0.5

restore,infofilename

master_path = file_dirname(infofilename)
rh_rundirectoryhead = master_path+'/'+file_basename(rh_rundirectoryhead)+'/'

firstiterate=1
recompute_save=1
i=0
j=0
k=0

rhdir = master_path + '/' + run_label + '_' + strtrim(string(i),2) + '_' + strtrim(string(j),2)+'_'+strtrim(string(k),2)+'/'
savefile = rhdir+'spectrum.sav'
$export SHELL=/bin/csh
save,rhdir,master_path,rayfilename,filename='inner_loop_script_inputs.sav'
if(recompute_save or (file_test(savefile) eq 0)) then spawn,'idl -e @inner_loop_script.pro'
$export SHELL=/bin/bash
restore, savefile

if(n_elements(ray_out) eq 0) then ray_out = 'None'

atmos_geometry = geometry
geometry.dx = geometry.dx*ray_fac_dx
geometry.dy = geometry.dy*ray_fac_dy

save, atmos, spectrum, geometry, atmos_geometry, zcont, zcore, zcont01, zcore01, zcont03, zcore03, zcont07, zcore07, zcont20, zcore20, zcont30, zcore30, ray, filename=outfilename
