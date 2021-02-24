.compile get_contrib_function_depth.pro

fcontrib = 0.5
depth = get_contrib_function_depth(ray2, fcontrib, xray=xray, yray=yray, zray=zray, tauray=tauray, contray=contray)
depth01 = get_contrib_function_depth(ray2, 0.1, xray=x01, zray=z01)
depth09 = get_contrib_function_depth(ray2, 0.9, xray=x09, zray=z09)

vdepth = get_contrib_function_depth(vray, fcontrib, xray=vy, yray=vx, zray=vz, tauray=vtauray, contray=vcontray, /taudepth)
vdepth01 = get_contrib_function_depth(vray, 0.1, yray=vx01, zray=vz01)
vdepth09 = get_contrib_function_depth(vray, 0.9, yray=vx09, zray=vz09)

tauray_cont = mean(tauray[*,*,*,[0,2]],dimension=4)
vtauray_cont = mean(vtauray[*,*,*,[0,2]],dimension=4)

tauray_core = reform(tauray[*,*,*,1])
vtauray_core = reform(vtauray[*,*,*,1])


vz_cont = mean(vz[*,*,[0,2]],dimension=3)
vz09_cont = mean(vz09[*,*,[0,2]],dimension=3)
vz01_cont = mean(vz01[*,*,[0,2]],dimension=3)
vx_cont = mean(vx[*,*,[0,2]],dimension=3)
x05_cont = mean(xray[*,*,[0,2]],dimension=3)
z05_cont = mean(zray[*,*,[0,2]],dimension=3)
vx05_cont = mean(vx[*,*,[0,2]],dimension=3)
vz05_cont = mean(vz[*,*,[0,2]],dimension=3)
depth01_cont = mean(depth01[*,*,[0,2]],dimension=3)
depth09_cont = mean(depth09[*,*,[0,2]],dimension=3)
x01_cont = mean(x01[*,*,[0,2]],dimension=3)
x09_cont = mean(x09[*,*,[0,2]],dimension=3)
z01_cont = mean(z01[*,*,[0,2]],dimension=3)
z09_cont = mean(z09[*,*,[0,2]],dimension=3)

depth01_core = reform(depth01[*,*,1])
depth09_core = reform(depth09[*,*,1])
vz09_core = reform(vz09[*,*,1])
vz01_core = reform(vz01[*,*,1])
x01_core = reform(x01[*,*,1])
x09_core = reform(x09[*,*,1])
z01_core = reform(z01[*,*,1])
z09_core = reform(z09[*,*,1])

x05_core = reform(xray[*,*,1])
z05_core = reform(zray[*,*,1])

vx05_core = reform(vx[*,*,1])
vz05_core = reform(vz[*,*,1])


islice = 44
;jslice1=1013
;jslice2=1028
jslice1=1018
jslice2=1018
jslice = 0.5*(jslice1+jslice2)

nx = n_elements(depth[*,0,0])
ny = n_elements(depth[0,*,0])
nz = n_elements(geometry.z)

dz = geometry.z[0]-geometry.z[1]

xa = findgen(nx)#(1.0+fltarr(ny))
ya = (1.0+fltarr(nx))#findgen(ny)

distance = sqrt((x05_cont-jslice)^2+(z05_cont-vz_cont[islice,jslice])^2)
dmin = min(distance,imin)

offset = round(jslice-ya[imin])

;offset = -round(4.5*median(ya[31,*]-xray[31,*]))

if(n_elements(ray0) eq 0) then ray0 = ray
if(n_elements(ray0) ne 0) then ray = ray0

ray.i = shift(ray.i,0,offset,0)
ray.stokes_q = shift(ray.stokes_q,0,offset,0)
ray.stokes_u = shift(ray.stokes_u,0,offset,0)
ray.stokes_v = shift(ray.stokes_v,0,offset,0)

nw = n_elements(spectrum.i[0,0,*])
icont = 0.5*(spectrum.i[*,*,0]+spectrum.i[*,*,nw-1])
icont_ray = 0.5*(ray.i[*,*,0]+ray.i[*,*,nw-1])
icore = min(spectrum.i,dimension=3);/icont
qtot = mean(spectrum.stokes_q,dimension=3);/icore
utot = mean(spectrum.stokes_u,dimension=3);/icore
vabstot = mean(abs(spectrum.stokes_v),dimension=3);/icore

qprange = 0.1*max(abs(qtot))
uprange = 0.1*max(abs(utot))
icoreprange = 1.2*median(icore)
icontprange = 1.2*median(icont)


specnorm = fltarr(nx,ny,nw)
for i=0,nw-1 do specnorm[*,*,i] = icont
raynorm = fltarr(nx,ny,nw)
for i=0,nw-1 do raynorm[*,*,i] = icont_ray


set_plot,'z'
device,set_resolution=[2560,1280]
outdir = '/data/jplowman/GONG_Simulator/plots/'

!p.multi=[0,1,5]
plot_image, transpose(icont),title = 'Vertical I!DContinuum!N', charsize=4,/nosquare, min=0, max=icontprange
oplot,[0,ny],[islice,islice],linestyle=2
plot_image, transpose(icore),title = 'Vertical I!DCore!N', charsize=4,/nosquare, min=0, max=icoreprange
oplot,[0,ny],[islice,islice],linestyle=2
plot_image, transpose(qtot),title = 'Vertical Q!Dtot!N', charsize=4, min=-qprange, max=qprange,/nosquare
oplot,[0,ny],[islice,islice],linestyle=2
plot_image, transpose(utot),title = 'Vertical U!Dtot!N', charsize=4, min=-uprange, max=uprange, /nosquare
oplot,[0,ny],[islice,islice],linestyle=2
plot_image, transpose(vabstot),title = 'Vertical |V!Dtot!N|', charsize=4,/nosquare
oplot,[0,ny],[islice,islice],linestyle=2


write_tiff,outdir+'spectrum_vertical_strip_averages.tiff',reverse(tvrd(),2)

icont_60 = 0.5*(ray.i[*,*,0]+ray.i[*,*,nw-1])
icore_60 = min(ray.i,dimension=3);/icont
qtot_60 = mean(ray.stokes_q,dimension=3);/icore
utot_60 = mean(ray.stokes_u,dimension=3);/icore
vabstot_60 = mean(abs(ray.stokes_v),dimension=3);/icore

!p.multi=[0,1,5]
plot_image, transpose(icont_60),title = '60 degree I!DContinuum!N', charsize=4, min=0, max=icontprange, scale=[0.5,1],/nosquare
oplot,[0,ny],[islice,islice]
plot_image, transpose(icore_60),title = '60 degree I!DCore!N', charsize=4, min=0, max=icoreprange, scale=[0.5,1],/nosquare
oplot,[0,ny],[islice,islice]
plot_image, transpose(qtot_60),title = '60 degree Q!Dtot!N', charsize=4, min=-qprange, max=qprange,scale=[0.5,1],/nosquare
oplot,[0,ny],[islice,islice]
plot_image, transpose(utot_60),title = '60 degree U!Dtot!N', charsize=4, min=-uprange, max=uprange,scale=[0.5,1],/nosquare
oplot,[0,ny],[islice,islice]
plot_image, transpose(vabstot_60),title = '60 degree |V!Dtot!N|', charsize=4, min=0, max=max(vabstot),scale=[0.5,1],/nosquare
oplot,[0,ny],[islice,islice]

write_tiff,outdir+'spectrum_60degree_strip_averages.tiff',reverse(tvrd(),2)

iprange = 1.2*median(abs(spectrum.i[islice,*,*])); max(abs([spectrum.i[islice,*,*],ray.i[islice,*,*]]))
qprange = max(abs([spectrum.stokes_q[islice,*,*]/specnorm[islice,*,*],ray.stokes_q[islice,*,*]/raynorm[islice,*,*]]))
uprange = max(abs([spectrum.stokes_u[islice,*,*]/specnorm[islice,*,*],ray.stokes_u[islice,*,*]/raynorm[islice,*,*]]))
vprange = max(abs([spectrum.stokes_v[islice,*,*]/specnorm[islice,*,*],ray.stokes_v[islice,*,*]/raynorm[islice,*,*]]))

!p.multi=[0,1,4]
plot_image, reform(spectrum.i[islice,*,*]), title = 'Vertical Stokes I/I!Dcont!N', charsize=4,/nosquare, min=0, max=iprange
oplot,[jslice1,jslice1],[0,nw],linestyle=2
oplot,[jslice2,jslice2],[0,nw],linestyle=2
plot_image, reform(spectrum.stokes_q[islice,*,*]/specnorm[islice,*,*]),min=-qprange,max=qprange, title = 'Vertical Q/I!Dcont!N', charsize=4,/nosquare
oplot,[jslice1,jslice1],[0,nw],linestyle=2
oplot,[jslice2,jslice2],[0,nw],linestyle=2
plot_image, reform(spectrum.stokes_u[islice,*,*]/specnorm[islice,*,*]),min=-uprange,max=uprange, title = 'Vertical U/I!Dcont!N', charsize=4,/nosquare
oplot,[jslice1,jslice1],[0,nw],linestyle=2
oplot,[jslice2,jslice2],[0,nw],linestyle=2
plot_image, reform(spectrum.stokes_v[islice,*,*]/specnorm[islice,*,*]),min=-vprange,max=vprange, title = 'Vertical V/I!Dcont!N', charsize=4,/nosquare
oplot,[jslice1,jslice1],[0,nw],linestyle=2
oplot,[jslice2,jslice2],[0,nw],linestyle=2

write_tiff,outdir+'spectrum_vertical_strip_slice.tiff',reverse(tvrd(),2)

!p.multi=[0,1,4]
plot_image, reform(ray.i[islice,*,*]), title = '60 degree Stokes I', charsize=4, /nosquare, scale=[0.5,1], min=0, max=iprange
oplot,[jslice1,jslice1]*0.5,[0,nw],linestyle=2
oplot,[jslice2,jslice2]*0.5,[0,nw],linestyle=2
plot_image, reform(ray.stokes_q[islice,*,*]/raynorm[islice,*,*]),min=-qprange,max=qprange, title = '60 degree Q/I!Dcont!N', charsize=4,/nosquare, scale=[0.5,1]
oplot,[jslice1,jslice1]*0.5,[0,nw],linestyle=2
oplot,[jslice2,jslice2]*0.5,[0,nw],linestyle=2
plot_image, reform(ray.stokes_u[islice,*,*]/raynorm[islice,*,*]),min=-uprange,max=uprange, title = '60 degree U/I!Dcont!N', charsize=4,/nosquare, scale=[0.5,1]
oplot,[jslice1,jslice1]*0.5,[0,nw],linestyle=2
oplot,[jslice2,jslice2]*0.5,[0,nw],linestyle=2
plot_image, reform(ray.stokes_v[islice,*,*]/raynorm[islice,*,*]),min=-vprange,max=vprange, title = '60 degree V/I!Dcont!N', charsize=4,/nosquare, scale=[0.5,1]
oplot,[jslice1,jslice1]*0.5,[0,nw],linestyle=2
oplot,[jslice2,jslice2]*0.5,[0,nw],linestyle=2

write_tiff,outdir+'spectrum_60degree_strip_slice.tiff',reverse(tvrd(),2)

speciplot = mean(spectrum.i[islice,jslice1:jslice2,*],dimension=2)
rayiplot = mean(ray.i[islice,jslice1:jslice2,*],dimension=2)
ipltmax = max([speciplot,rayiplot])

specqplot = mean(spectrum.stokes_q[islice,jslice1:jslice2,*]/specnorm[islice,jslice1:jslice2,*],dimension=2)
rayqplot = mean(ray.stokes_q[islice,jslice1:jslice2,*]/raynorm[islice,jslice1:jslice2,*],dimension=2)
qpltrng = [min([specqplot,rayqplot]),max([specqplot,rayqplot])]

specuplot = mean(spectrum.stokes_u[islice,jslice1:jslice2,*]/specnorm[islice,jslice1:jslice2,*],dimension=2)
rayuplot = mean(ray.stokes_u[islice,jslice1:jslice2,*]/raynorm[islice,jslice1:jslice2,*],dimension=2)
upltrng = [min([specuplot,rayuplot]),max([specuplot,rayuplot])]

specvplot = mean(spectrum.stokes_v[islice,jslice1:jslice2,*]/specnorm[islice,jslice1:jslice2,*],dimension=2)
rayvplot = mean(ray.stokes_v[islice,jslice1:jslice2,*]/raynorm[islice,jslice1:jslice2,*],dimension=2)
vpltrng = [min([specvplot,rayvplot]),max([specvplot,rayvplot])]


!p.multi=[0,2,2]
plot,spectrum.lambda,speciplot,title='Stokes I', charsize=2, yrange=[0,1.2*ipltmax], ystyle=1
oplot,spectrum.lambda,rayiplot,linestyle=2
ssw_legend,['Vertical','60 degree'],linestyle=[0,2]

plot,spectrum.lambda,specqplot,title='Stokes Q/I!Dcont!N', charsize=2, yrange=qpltrng
oplot,spectrum.lambda,rayqplot,linestyle=2
ssw_legend,['Vertical','60 degree'],linestyle=[0,2]

plot,spectrum.lambda,specuplot,title='Stokes U/I!Dcont!N', charsize=2, yrange=upltrng
oplot,spectrum.lambda,rayuplot,linestyle=2
ssw_legend,['Vertical','60 degree'],linestyle=[0,2]

plot,spectrum.lambda,specvplot,title='Stokes V/I!Dcont!N', charsize=2, yrange=vpltrng
oplot,spectrum.lambda,rayvplot,linestyle=2
ssw_legend,['Vertical','60 degree'],linestyle=[0,2]

write_tiff,outdir+'spectrum_vertical_vs_60degree_compare.tiff',reverse(tvrd(),2)

!p.multi=[0,2,2]
tray = interpolate(atmos.t,reform(ray2.yray[islice,jslice-offset,*]),reform(ray2.xray[islice,jslice-offset,*]),reform(ray2.zray[islice,jslice-offset,*]))
tvert = reform(atmos.t[islice,jslice,*])
tauray_plt = reform(mean(tauray[islice,jslice-offset,*,[0,2]],dimension=4))
contray_plt = reform(mean(contray[islice,jslice-offset,*,[0,2]],dimension=4))
vtauray_plt = reform(mean(vtauray[islice,jslice,*,[0,2]],dimension=4))
vcontray_plt = reform(mean(vcontray[islice,jslice,*,[0,2]],dimension=4))

plot,tauray_plt, dz*reform(ray2.zray[islice,jslice-offset,*]), title='Vertical depth vs. tau', charsize=2, xrange=[1.e-4,10], xtitle='tau', ytitle='Vertical depth', /xlog
oplot, vtauray_plt, reverse(geometry.z), linestyle=2
ssw_legend, ['60 degree ray','Vertical ray'], linestyle=[0,2]

plot, tauray_plt, ray2.raylength, title='Path length vs. tau', charsize=2, xrange=[1.e-4,10], xtitle='tau', ytitle='length', /xlog
oplot, vtauray_plt, reverse(geometry.z), linestyle=2
ssw_legend, ['60 degree ray','Vertical ray'], linestyle=[0,2]

plot, tauray_plt, contray_plt, title='Contribution function vs. tau', charsize=2, xrange=[1.e-4,10], xtitle='tau', ytitle='Contribution function', yrange=[0,max([contray_plt,vcontray_plt])], /xlog
oplot, vtauray_plt, vcontray_plt, linestyle=2
ssw_legend, ['60 degree ray','Vertical ray'], linestyle=[0,2]

plot, tauray_plt, tray, title='Temperature vs. tau', charsize=2, xrange=[1.e-4,10], xtitle='tau', ytitle='Temperature', yrange=[0,8000], ystyle=1, /xlog
oplot, vtauray_plt, tvert, linestyle=2
ssw_legend, ['60 degree ray','Vertical ray'], linestyle=[0,2]

write_tiff,outdir+'distance_temp_vs_tau.tiff',reverse(tvrd(),2)

ray = ray0

!p.multi=0

plot_image,reverse(reform(atmos.t[islice,*,*]),2),/nosquare, charsize=2
oplot,x09_cont[islice,*],(nz-1)-z09_cont[islice,*],linestyle=2
oplot,x01_cont[islice,*],(nz-1)-z01_cont[islice,*],linestyle=2
oplot,(nz-1)-vz01_cont[islice,*]
oplot,(nz-1)-vz09_cont[islice,*]
oplot,[jslice1,jslice1],[0,100]
oplot,[jslice2,jslice2],[0,100]
oplot,ray2.xray[islice,jslice1-offset,*],(nz-1)-ray2.zray[islice,jslice1-offset,*]
oplot,ray2.xray[islice,jslice2-offset,*],(nz-1)-ray2.zray[islice,jslice2-offset,*]

write_tiff,outdir+'temperature_vertical_field_compare_cont.tiff',reverse(tvrd(),2)

plot_image,reverse(reform(atmos.t[islice,*,*]),2),/nosquare, charsize=2
oplot,x09_core[islice,*],(nz-1)-z09_core[islice,*],linestyle=2
oplot,x01_core[islice,*],(nz-1)-z01_core[islice,*],linestyle=2
oplot,(nz-1)-vz01_core[islice,*]
oplot,(nz-1)-vz09_core[islice,*]
oplot,[jslice1,jslice1],[0,100]
oplot,[jslice2,jslice2],[0,100]
oplot,ray2.xray[islice,jslice1-offset,*],(nz-1)-ray2.zray[islice,jslice1-offset,*]
oplot,ray2.xray[islice,jslice2-offset,*],(nz-1)-ray2.zray[islice,jslice2-offset,*]

write_tiff,outdir+'temperature_vertical_field_compare_core.tiff',reverse(tvrd(),2)

set_plot,'x'
jb_low = 0
jb_high = 255
gammab = reverse(rebin(reform(atmos.gamma_b[islice,jb_low:jb_high,*]),jb_high+1,48),2)
chib = reverse(rebin(reform(atmos.chi_b[islice,jb_low:jb_high,*]),jb_high+1,48),2)
bz = cos(gammab)
by = sin(gammab)*sin(chib)
y2 = geometry.dy*findgen(jb_high+1)
z2 = rebin(reverse(geometry.z),48)-(geometry.z[0]-geometry.z[1])
plot_image,reverse(reform(abs(atmos.gamma_b[islice,jb_low:jb_high,*])),2),scale=[geometry.dx,(geometry.z[0]-geometry.z[1])], charsize=2
velovect,by,bz,y2,z2,/overplot,color=0


lambda0 = spectrum.lambda[ray2.nspect[0]]
t_tau1ray = interpol(tray,tauray_plt,1)
t_tau1vert = interpol(tvert,vtauray_plt,1)
planck_tau1vert = planck(t_tau1vert,lambda0,/hz)
planck_tau1ray = planck(t_tau1ray,lambda0,/hz)
intens_tau1vert = speciplot[ray2.nspect[0]]
intens_tau1ray = rayiplot[ray2.nspect[0]]

print,'Continuum Planck function at vertical tau=1: '+string(planck_tau1vert)
print,'Continuum Planck function at inclined tau=1: '+string(planck_tau1ray)

print,'Continuum intensity at vertical tau=1: '+string(intens_tau1vert)
print,'Continuum intensity at inclined tau=1: '+string(intens_tau1ray)

print,'Vertical to inclined Planck ratio: '+string(planck_tau1vert/planck_tau1ray)
print,'Vertical to inclined intens ratio: '+string(intens_tau1vert/intens_tau1ray)