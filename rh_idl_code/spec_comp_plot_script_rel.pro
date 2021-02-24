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

active_atoms = atom_files
nactive = n_elements(active_atoms)

if(active_atoms ne '') then begin &$
	for i=0,nactive-1 do begin &$
		atomname = strsplit(active_atoms[i],'.',/extract) &$
		active_atoms[i] = atomname[1] &$
	endfor &$
end

plotname = plotfilehead+atmos.id+strjoin(strtrim(active_atoms,2),'_')+'_'+strjoin(strtrim(string([lambdamin,lambdamax]),2),'_')+'_rel.eps'

set_plot,'ps'
device,/inches,/encapsulated,filename=plotname,xsize=8,ysize=8,bits=8


reldiff_plot,spectrum.lambda,ray.i,spectrum2.lambda,ray2.i,xrange=[lambdamin,lambdamax], title='1-I!DnLTE!N/I!DLTE!N', xtitle='wavelength (nm)'

reldiff_plot,spectrum.lambda,ray.stokes_q,spectrum2.lambda,ray2.stokes_q,xrange=[lambdamin,lambdamax], title='1-Q!DnLTE!N/Q!DLTE!N', xtitle='wavelength (nm)',yrange=[-0.25,0.25]

reldiff_plot,spectrum.lambda,ray.stokes_u,spectrum2.lambda,ray2.stokes_u,xrange=[lambdamin,lambdamax], title='1-U!DnLTE!N/U!DLTE!N', xtitle='wavelength (nm)',yrange=[-0.25,0.25]

reldiff_plot,spectrum.lambda,ray.stokes_v,spectrum2.lambda,ray2.stokes_v,xrange=[lambdamin,lambdamax], title='1-V!DnLTE!N/V!DLTE!N', xtitle='wavelength (nm)',yrange=[-0.25,0.25]

device,/close
