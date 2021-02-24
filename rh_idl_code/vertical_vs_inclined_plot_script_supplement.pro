depth_tau1 = get_contrib_function_depth(ray2, 1.0, /taudepth, xray=yrtau, yray=yxtau, zray=zrtau)
vdepth_tau1 = get_contrib_function_depth(vray, 1.0, /taudepth, xray=yxrtau, yray=vxrtau, zray=vzrtau)

xrtau_core = reform(xrtau[*,*,1])
zrtau_core = reform(zrtau[*,*,1])
xrtau_cont = mean(xrtau[*,*,[0,2]],dimension=3)
zrtau_cont = mean(zrtau[*,*,[0,2]],dimension=3)
vxrtau_core = reform(vxrtau[*,*,1])
vzrtau_core = reform(vzrtau[*,*,1])
vxrtau_cont = mean(vxrtau[*,*,[0,2]],dimension=3)
vzrtau_cont = mean(vzrtau[*,*,[0,2]],dimension=3)

set_plot,'z'
!p.multi=[0,1,4]

plot_image,reverse(reform(atmos.t[islice,*,*]),2),/nosquare, charsize=2, scale=[geometry.dy,dz]/1000.0, title='Vertical Continuum contribution function and optical depth heights'
oplot,vx05_cont[islice,*]*geometry.dy/1000.0,((nz-1)-vz05_cont[islice,*])*dz/1000.0
oplot,vxrtau_cont[islice,*]*geometry.dy/1000.0,((nz-1)-vzrtau_cont[islice,*])*dz/1000.0, linestyle=2
ssw_legend,['Contribution function 50%','Optical Depth Unity'], linestyle=[0,2],/right

plot_image,reverse(reform(atmos.t[islice,*,*]),2),/nosquare, charsize=2, scale=[geometry.dy,dz]/1000.0, title='Vertical core contribution function and optical depth heights'
oplot,vx05_core[islice,*]*geometry.dy/1000.0,((nz-1)-vz05_core[islice,*])*dz/1000.0
oplot,vxrtau_core[islice,*]*geometry.dy/1000.0,((nz-1)-vzrtau_core[islice,*])*dz/1000.0, linestyle=2
ssw_legend,['Contribution function 50%','Optical Depth Unity'], linestyle=[0,2],/right

plot_image,reverse(reform(atmos.t[islice,*,*]),2),/nosquare, charsize=2, scale=[geometry.dy,dz]/1000.0, title='Inclined Continuum contribution function and optical depth heights'
oplot,x05_cont[islice,*]*geometry.dy/1000.0,((nz-1)-z05_cont[islice,*])*dz/1000.0
oplot,xrtau_cont[islice,*]*geometry.dy/1000.0,((nz-1)-zrtau_cont[islice,*])*dz/1000.0, linestyle=2
ssw_legend,['Contribution function 50%','Optical Depth Unity'], linestyle=[0,2],/right

plot_image,reverse(reform(atmos.t[islice,*,*]),2),/nosquare, charsize=2, scale=[geometry.dy,dz]/1000.0, title='Inclined core contribution function and optical depth heights'
oplot,x05_core[islice,*]*geometry.dy/1000.0,((nz-1)-z05_core[islice,*])*dz/1000.0
oplot,xrtau_core[islice,*]*geometry.dy/1000.0,((nz-1)-zrtau_core[islice,*])*dz/1000.0, linestyle=2
ssw_legend,['Contribution function 50%','Optical Depth Unity'], linestyle=[0,2],/right

write_tiff,outdir+'height_comparison_plots.tiff',reverse(tvrd(),2)


print,'Inclined line core median contribution function height = '+strtrim(string(interpol(geometry.z,findgen(geometry.nz),median(z05_core[*,0:500]))/1000.0),2)+' km'
print,'Inclined line core median tau=1 height = '+strtrim(string(interpol(geometry.z,findgen(geometry.nz),median(zrtau_core[*,0:500]))/1000.0),2)+' km'

print,'Inclined continuum median contribution function height = '+strtrim(string(interpol(geometry.z,findgen(geometry.nz),median(z05_cont[*,0:500]))/1000.0),2)+' km'
print,'Inclined continuum median tau=1 height = '+strtrim(string(interpol(geometry.z,findgen(geometry.nz),median(zrtau_cont[*,0:500]))/1000.0),2)+' km'

print,'Vertical line core median contribution function height = '+strtrim(string(interpol(geometry.z,findgen(geometry.nz),median(vz05_core[*,0:500]))/1000.0),2)+' km'
print,'Vertical line core median tau=1 height = '+strtrim(string(interpol(geometry.z,findgen(geometry.nz),median(vzrtau_core[*,0:500]))/1000.0),2)+' km'

print,'Vertical continuum median contribution function height = '+strtrim(string(interpol(geometry.z,findgen(geometry.nz),median(vz05_cont[*,0:500]))/1000.0),2)+' km'
print,'Vertical continuum median tau=1 height = '+strtrim(string(interpol(geometry.z,findgen(geometry.nz),median(vzrtau_cont[*,0:500]))/1000.0),2)+' km'



