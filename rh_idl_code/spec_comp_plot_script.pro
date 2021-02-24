;savfile1 = '../rh_example_files_FALC_noB/ivspec.sav'
;savfile2 = '../rh_example_files_FALC_677nonLTE_noB/ivspec.sav'
;plotfilehead = 'stokesIVray_2file_noB_'

;savfile1 = '../rh_example_files_FALC/ivspec.sav'
;savfile2 = '../rh_example_files_FALC_677nonLTE/ivspec.sav'
;plotfilehead = 'stokesIVray_2file_'

savfile1 = '../rh_example_files_FALP/ivspec.sav'
savfile2 = '../rh_example_files_FALP_677nonLTE/ivspec.sav'
plotfilehead = 'stokesIVray_2file_'

restore,savfile2
spectrum2=spectrum
ray2=ray
atom_files2=atom_files
atmos2=atmos

restore,savfile1

lambdamin=676.7
lambdamax=676.85

ftsread,data,lambdamin*10,lambdamax*10,xlam=xlam,sdir='.'
spec_atlas = interpol(data,xlam*0.1,spectrum.lambda)
normfac = median(spec_atlas/ray.I)

active_atoms = atom_files
nactive = n_elements(active_atoms)

if(active_atoms ne '') then begin &$
	for i=0,nactive-1 do begin &$
		atomname = strsplit(active_atoms[i],'.',/extract) &$
		active_atoms[i] = atomname[1] &$
	endfor &$
end

plotname = plotfilehead+atmos.id+strjoin(strtrim(active_atoms,2),'_')+'_'+strjoin(strtrim(string([lambdamin,lambdamax]),2),'_')+'.eps'

set_plot,'ps'
device,/inches,/encapsulated,filename=plotname,xsize=8,ysize=8,bits=8

!p.multi=[0,2,2]
plot,spectrum.lambda,normfac*ray.I,xrange=[lambdamin,lambdamax], title='Stokes I', xtitle='wavelength (nm)'
oplot,spectrum.lambda,normfac*ray.I, color=128, thick=11
oplot,spectrum2.lambda,normfac*ray2.I, thick=1
oplot,0.1*xlam,data,linestyle=0,thick=3
oplot,0.1*xlam,data,linestyle=0,thick=1,color=255
ssw_legend,['LTE','Non LTE','FTS Atlas'],/bottom,linestyle=[0,0,0],color=[128,0,0],thick=[11,1,3]
ssw_legend,['LTE','Non LTE','FTS Atlas'],/bottom,linestyle=[0,0,0],color=[128,0,255],thick=[11,1,1]

plot,spectrum.lambda,normfac*ray.stokes_q,xrange=[lambdamin,lambdamax], title='Stokes Q', xtitle='wavelength (nm)'
oplot,spectrum.lambda,normfac*ray.stokes_q, color=128, thick=11
oplot,spectrum2.lambda,normfac*ray2.stokes_q, thick=1
ssw_legend,['LTE','Non LTE'],/bottom,linestyle=[0,0],color=[128,0],thick=[11,1]

plot,spectrum.lambda,normfac*ray.stokes_u,xrange=[lambdamin,lambdamax], title='Stokes U', xtitle='wavelength (nm)'
oplot,spectrum.lambda,normfac*ray.stokes_u, color=128, thick=11
oplot,spectrum2.lambda,normfac*ray2.stokes_u, thick=1
ssw_legend,['LTE','Non LTE'],/bottom,linestyle=[0,0],color=[128,0],thick=[11,1]

plot,spectrum.lambda,normfac*ray.stokes_v,xrange=[lambdamin,lambdamax], title='Stokes V', xtitle='wavelength (nm)'
oplot,spectrum.lambda,normfac*ray.stokes_v, color=128, thick=11
oplot,spectrum2.lambda,normfac*ray2.stokes_v, thick=1
ssw_legend,['LTE','Non LTE'],/bottom,linestyle=[0,0],color=[128,0],thick=[11,1]

device,/close
