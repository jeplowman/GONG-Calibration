.compile ../rh_idl_code/get_contrib_function_depth.pro
.compile ../rh_idl_code/ray_resample_atmos.pro
.compile ../rh_idl_code/get_vertical_ray.pro
.compile median2.pro
.compile get_ar_thold.pro
.compile image_coord_transforms.pro
.compile gong_sim_pixel_psf.pro
.compile multi_tile_2.pro
.compile gong_dirty_raw_process.pro
.compile gong_sim_image_geometry.pro
.compile bz_tautile_samp.pro
.compile gong_cal_curve_lmfit.pro
.compile gong_curve_fluxcons_binfit.pro
.compile get_calibration_curves.pro
.compile gongbz_calibration_curve_psf_comparison.pro
.compile apply_calibration_curve2.pro
.compile do_gongsim_scatterplot.pro

tau01=0
nopsf=0
zconst=0
use_taucore=0
zlvl=175.0
binfit=1
nterms_in_qs=[35,21,15,9,7,7]
nterms_in_ar=[21,9,7,7,7,7]
nterms = 1001
xchang=0
outres1 = 1024
brange = 1000
brange_diff = 100
deconvolve=0
post_deconvolve=0
use_contrib=0
atmo_w = 3.6 ; Smoothing due to atmospheric seeing.
active_tile_level = 1
do_mc=0
tlstr = ''
if(active_tile_level ne 1) then tlstr = '_'+strtrim(string(active_tile_level),2)+'x'
max_sunspot_angle_index = 1

angles_deg = [0.0,25.0,45.0,60.0,70.0,75.0]
angles = angles_deg*!pi/180.0
anglestrs = strtrim(string(round(angles_deg)),2)
n_angles = n_elements(angles)
zlvls = zlvl+dblarr(n_angles)

; Original runs were made with these as absolute paths; change to relative path name is untested.

gong_sim_run_directory = '../gong_sim_runs/'
gong_data_dir = '../GONG_data/'

simfile_dir = '../GONG_sim_data_rh/'
simfile_angledirs = ['ARQScube2/','ARQScube2_25degrees/','ARQScube2_45degrees/', $
		'ARQScube2_60degrees/','ARQScube2_70degrees/','ARQScube_75degrees/']
simfile_name = 'terif170402/TE170402154716.fits.gz'
simfiles = simfile_dir+simfile_angledirs+simfile_name
modelfiles = ['gongrh_ARQScube_2_012820.sav', 'gongrh_ARQScube2_25degrees_012820.sav', $
		'gongrh_ARQScube2_45degrees_012820.sav','gongrh_ARQScube_2_60degrees_012820.sav', $ 
		'gongrh_ARQScube2_70degrees_012820.sav','gongrh_ARQScube_75degrees_012820.sav']
modelfiles = gong_sim_run_directory+modelfiles
dailyfile_directory = gong_data_dir+'daily_files_100608_42/'

spectrum_directory = '../rh_runs/' ; Used to fix path to spectrum. This is untested.
out_dir = dailyfile_directory+'exp_cal_allangles_binfit/'
plot_dir = dailyfile_directory+'cal_images_allangles_binfit/'
paper_plots_dir = '../plots/contribfn_testplots_040320/'
paper_plot_names = paper_plots_dir+'gongcal_plot_'+anglestrs+tlstr+'.eps'
paper_arscatter_plot_names = paper_plots_dir+'gongcal_plot_sunspot_'+anglestrs+tlstr+'.eps'

bmin = -200
bmax = 200
set_plot,'ps'
!p.multi=0; [0,3,0]
device,/encapsulated,/inches,xsize=4.5,ysize=3,bits=8
brangestr = ' (+- '+strtrim(string(round(bmax)),2)+' Gauss)'
bz_obstitles = anglestrs+'!Uo!N magnetogram'
bz_modtitles = anglestrs+'!Uo!N ground truth'

net_model_fluxes = fltarr(n_angles)
net_obs_fluxes = fltarr(n_angles)
net_cal_fluxes = fltarr(n_angles)

smooth_range = 2.0*brange*dindgen(1000)-brange
smooth_width = 2.5/(smooth_range[0]-smooth_range[1])

qs_slopes = dblarr(n_angles)
ar_slopes = dblarr(n_angles)

qs_modfits = fltarr(n_angles,nterms)
qs_obsfits = fltarr(n_angles,nterms)
ar_modfits = fltarr(n_angles,nterms)
ar_obsfits = fltarr(n_angles,nterms)
for i=0,n_angles-1 do begin &$
	print,'Computing curves for '+anglestrs[i]+' degrees:' &$
	get_calibration_curves, tau01, nopsf, zconst, use_taucore, zlvls[i], bz_obs, bz_model, bz_cal, qsflag, arflag, tile_level, ondisk, qs_modfit0, qs_obsfit0, ar_modfit0, ar_obsfit0, binfit=binfit, nterms=nterms_in_qs[i], ar_nterms=nterms_in_ar[i], xchang=xchang, simfile=simfiles[i], modelfile=modelfiles[i], deconvolve=deconvolve, theta_ray=angles[i], use_contrib=use_contrib, post_deconvolve=post_deconvolve, atmo_w=atmo_w, tlvl0=active_tile_level, gongitot=gongitot, spectrum_directory=spectrum_directory &$
	if(i eq 0) then begin &$
		nx = n_elements(bz_obs[*,0]) &$
		ny = n_elements(bz_obs[0,*]) &$
		bz_gt = fltarr(nx,ny,n_angles) &$
		bz_meas = fltarr(nx,ny,n_angles) &$
		qsflags = fltarr(nx,ny,n_angles) &$
		arflags = fltarr(nx,ny,n_angles) &$
		tile_levels = fltarr(nx,ny,n_angles) &$
		ondisks = fltarr(nx,ny,n_angles) &$
		gongitots = fltarr(nx,ny,n_angles) &$
		xa = lindgen(nx)#(1+lonarr(ny)) &$
		ya = (1+lonarr(nx))#lindgen(ny) &$
	endif &$
	bz_gt[*,*,i]=bz_model &$
	bz_meas[*,*,i]=bz_obs &$
	qsflags[*,*,i]=qsflag &$
	arflags[*,*,i]=arflag &$
	tile_levels[*,*,i]=tile_level &$
	ondisks[*,*,i] = ondisk &$
	gongitots[*,*,i] = gongitot &$
	qspx_gongres = where(qsflag and (tile_level eq active_tile_level) and ondisk) &$
	arpx_gongres = where(arflag and (tile_level eq active_tile_level) and ondisk) &$
	xl = min(xa[qspx_gongres]) &$
	xh = max(xa[qspx_gongres]) &$
	yl = min(ya[qspx_gongres]) &$
	yh = max(ya[qspx_gongres]) &$
	qs_interp_range = max(abs(bz_obs[qspx_gongres])) &$
	ar_interp_range = max(abs(bz_obs[arpx_gongres])) &$
	qs_obsfit = 2*qs_interp_range*dindgen(nterms)/(nterms-1.0)-qs_interp_range &$
	ar_obsfit = 2*ar_interp_range*dindgen(nterms)/(nterms-1.0)-ar_interp_range &$
	qs_modfit = interpol(qs_modfit0,qs_obsfit0,qs_obsfit) &$
	ar_modfit = interpol(ar_modfit0,ar_obsfit0,ar_obsfit) &$
	qs_slope = deriv(qs_obsfit,qs_modfit) &$
	foo = min(abs(qs_obsfit),qsfit_ix0) &$
	qs_slope = qs_slope[qsfit_ix0] &$
	qslabel='Non-sunspot calibration, slope='+strtrim(string(qs_slope),2) &$
	ar_slope = deriv(ar_obsfit,ar_modfit) &$
	foo = min(abs(ar_obsfit),arfit_ix0) &$
	ar_slope = ar_slope[arfit_ix0] &$
	arlabel='Sunspot calibration, slope='+strtrim(string(ar_slope),2) &$
	net_model_fluxes[i] = total(bz_model[qspx_gongres]) &$
	net_obs_fluxes[i] = total(bz_obs[qspx_gongres]) &$
	net_cal_fluxes[i] = total(bz_cal[qspx_gongres]) &$
	device,filename = paper_plot_names[i], /inches, xsize=9,ysize=6 &$
	plot_image,bz_model[xl:xh,yl:yh]*ondisk[xl:xh,yl:yh],ytitle=bz_modtitles[i]+brangestr,min=bmin, max=bmax, charsize=1.0, position=[0.075, 0.54, 0.375, 0.99], thick=2, /noadjust &$
	plot_image,bz_obs[xl:xh,yl:yh]*ondisk[xl:xh,yl:yh],ytitle=bz_obstitles[i]+brangestr,min=bmin, max=bmax, charsize=1.0, position=[0.075, 0.04, 0.375, 0.49], thick=2, /noerase, /noadjust &$
	!p.multi=0 &$
	do_gongsim_scatterplot,bmin,bmax,bz_obstitles[i],bz_modtitles[i],bz_obs[qspx_gongres],bz_model[qspx_gongres], qs_obsfit, qs_modfit, qslabel,position=[0.475,0.08,0.975,0.95],/noerase &$
	device,/close &$
	device,filename = paper_arscatter_plot_names[i], /inches, xsize=4.5,ysize=4.5 &$
	!p.multi=0 &$
	do_gongsim_scatterplot, bmin*20, bmax*20, bz_obstitles[i], bz_modtitles[i], bz_obs[arpx_gongres], bz_model[arpx_gongres], ar_obsfit, ar_modfit, arlabel &$
	device,/close &$
	qs_slopes[i] = qs_slope &$
	ar_slopes[i] = ar_slope &$
	qs_modfits[i,*] = qs_modfit &$
	qs_obsfits[i,*] = qs_obsfit &$
	ar_modfits[i,*] = ar_modfit &$
	ar_obsfits[i,*] = ar_obsfit &$
endfor

set_plot,'ps'
!p.multi=0
cmin = 64
cmax = 191
colors = cmin+(cmax-cmin)*findgen(n_angles)/(n_angles-1)
device,/inches,xsize=5,ysize=5,filename = paper_plots_dir+'calibration_curve_plots_combined'+tlstr+'.eps'
plot,[bmin,bmax],[bmin,bmax],xrange=[bmin,bmax],yrange=[bmin,bmax], xtitle='GONG simulator magnetogram', ytitle='Ground truth magnetogram', xstyle=1,ystyle=1,thick=1,title='Calibration curves, all angles', charsize=1,linestyle=2
oplot,[0,0],[bmin,bmax],linestyle=1 &$
oplot,[bmin,bmax],[0,0],linestyle=1 &$
for i=0,n_angles-1 do oplot,qs_obsfits[i,*],qs_modfits[i,*],linestyle=0,thick=3,color=colors[i]
ssw_legend,['y=x line',anglestrs+' degrees'],linestyle=[0,0+intarr(n_angles)],color=[0,colors],thick=[1,3+fltarr(n_angles)],/top,/left
device,/close

for i=0,n_angles-1 do begin &$
	if(angles_deg[i] gt angles_deg[max_sunspot_angle_index]+1) then begin &$
		ar_modfits[i,*] = ar_modfits[max_sunspot_angle_index,*] &$
		ar_obsfits[i,*] = ar_obsfits[max_sunspot_angle_index,*] &$
	endif &$
endfor

bzifiles = findfile(dailyfile_directory+'bbzqa*.fits.gz')
izifiles = findfile(dailyfile_directory+'bbiqa*.fits.gz')
ndays = n_elements(bzifiles)

cal_qsflags = qsflags*(tile_levels eq 1)
cal_arflags = arflags*(tile_levels eq 1)

if(active_tile_level eq 1) then for i=0,ndays-1 do bzical = apply_calibration_curve2(qs_obsfits, qs_modfits, ar_obsfits, ar_modfits, bzifiles[i], izifiles[i], out_dir=out_dir, plot_dir=plot_dir, prange=100, post_deconvolve=post_deconvolve, angles=angles, bzi_uncal=bzi_uncal, cal_measimages=bz_meas, cal_gtimages=bz_gt, cal_qsflags=cal_qsflags, cal_arflags=cal_arflags,do_mc=do_mc,mc_max_sunspot_angle=0.0*!pi/180.0)

cal_measimages = bz_meas
cal_gtimages = bz_gt

save,angles,angles_deg,anglestrs,qs_obsfits,qs_modfits,ar_obsfits,ar_modfits,qs_slopes,ar_slopes,filename='save/calibration_curves_tl'+strtrim(string(active_tile_level),2)+'.sav'

if(active_tile_level eq 1) then  save, cal_gtimages, cal_measimages, qsflags, arflags, tile_levels, ondisks,gongitots,angles,angles_deg,anglestrs,qs_obsfits,qs_modfits,ar_obsfits,ar_modfits,filename='save/calibration_data_040320.sav'
