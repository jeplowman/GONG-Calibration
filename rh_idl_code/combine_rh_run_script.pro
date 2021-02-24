.compile update_atmos.pro
.compile update_spectrum.pro
.compile update_ray.pro
.compile update_geometry.pro
.compile update_z.pro

;infofilename = '/data/jplowman/rh_runs/run_rhsc3d_ARQS_info.sav'
;run_label = 'run_rhsc3d_ARQS'
;outfilename = master_path+'/'+run_label+'_spectrum.sav'

infofilename = '/mnt/usb-WD_easystore_25FB_3753474C37303543-0:0-part1/jplowman/rh_runs/run_rhsc3d_ARQS_info.sav'
run_label = 'run_rhsc3d_ARQS'
outfilename = '/data/jplowman/rh_runs/run_rhsc3d_ARQS_75_combined/'+run_label+'_spectrum.sav'

rayfilename = 'spectrum_0.00_0.97'
ray_fac_dx = 1.0
ray_fac_dy = 0.258819

restore,infofilename

master_path = file_dirname(infofilename)
rh_rundirectoryhead = master_path+'/'+file_basename(rh_rundirectoryhead)+'/'

info_struc = {lambdamin:lambdamin, lambdamax:lambdamax, nlambdas:nlambdas, nlambdabins:nlambdabins, nperlambdabin:nperlambdabin, lambdas:lambdas, $
		dx:dx, dy:dy, dz:dz, n0:n0, n1:n1, n2:n2, xcrop:xcrop, ycrop:ycrop, nxcrop:nxcrop, nycrop:nycrop, data_directory:data_directory, $
		rh_rundirectoryhead:rh_rundirectoryhead}

xbin=2
ybin=2

firstiterate=1
recompute_save=1
for i=0,nxcrop-1 do begin &$
	for j=0,nycrop-1 do begin &$
		for k=0,nlambdabins-1 do begin &$
			rhdir = master_path + '/' + run_label + '_' + strtrim(string(i),2) + '_' + strtrim(string(j),2)+'_'+strtrim(string(k),2)+'/' &$
			savefile = rhdir+'spectrum.sav' &$
			$export SHELL=/bin/csh &$
			save,rhdir,master_path,rayfilename,filename='inner_loop_script_inputs.sav' &$
			if(recompute_save or (file_test(savefile) eq 0)) then spawn,'idl -e @inner_loop_script.pro' &$
			$export SHELL=/bin/bash &$
			restore, savefile &$
			update_atmos,i,j,k,atmos,atmos_out,info_struc,xbin=xbin,ybin=ybin &$
			update_spectrum,i,j,k,spectrum,spectrum_out,info_struc,xbin=xbin,ybin=ybin &$
			if(n_elements(rayfilename) eq 1) then update_ray,i,j,k,ray,ray_out,info_struc,xbin=xbin,ybin=ybin &$
			update_geometry,i,j,k,geometry,geometry_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcont,zcont_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcore,zcore_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcont01,zcont01_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcore01,zcore01_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcont03,zcont03_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcore03,zcore03_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcont07,zcont07_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcore07,zcore07_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcont20,zcont20_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcore20,zcore20_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcont30,zcont30_out,info_struc,xbin=xbin,ybin=ybin &$
			update_z,i,j,k,zcore30,zcore30_out,info_struc,xbin=xbin,ybin=ybin &$
			print,'Done with ',i,j,k &$
		endfor &$
	endfor &$
endfor

if(n_elements(ray_out) eq 0) then ray_out = 'None'

atmos = atmos_out
delvar,atmos_out
spectrum = spectrum_out
delvar,spectrum_out
ray = ray_out
delvar,ray_out
geometry = geometry_out
delvar,geometry_out
zcont = zcont_out
delvar,zcont_out
zcore = zcore_out
delvar,zcore_out
zcore01=zcore01_out
delvar,zcore01_out
zcont01=zcont01_out
delvar,zcont01_out
zcore03=zcore03_out
delvar,zcore03_out
zcont03=zcont03_out
delvar,zcont03_out
zcore07=zcore07_out
delvar,zcore07_out
zcont07=zcont07_out
delvar,zcont07_out
zcore20=zcore20_out
delvar,zcore20_out
zcont20=zcont20_out
delvar,zcont20_out
zcore30=zcore30_out
delvar,zcore30_out
zcont30=zcont30_out
delvar,zcont30_out

atmos_geometry = geometry
geometry.dx = geometry.dx*ray_fac_dx
geometry.dy = geometry.dy*ray_fac_dy

save, atmos, spectrum, geometry, atmos_geometry, zcont, zcore, zcont01, zcore01, zcont03, zcore03, zcont07, zcore07, zcont20, zcore20, zcont30, zcore30, ray, filename=outfilename
