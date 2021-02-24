cd, '/data/jplowman/rh_runs/run_rhsc3d_ARQS_0_0_0'
@initrh_noprompt
.r readall
spec00 = total(spectrum.i,3)
cd, '/data/jplowman/rh_runs/run_rhsc3d_ARQS_0_1_0'
@initrh_noprompt
.r readall
spec01 = total(spectrum.i,3)
cd, '/data/jplowman/rh_runs/run_rhsc3d_ARQS_0_2_0'
@initrh_noprompt
.r readall
spec02 = total(spectrum.i,3)
cd, '/data/jplowman/rh_runs/run_rhsc3d_ARQS_0_3_0'
@initrh_noprompt
.r readall
spec03 = total(spectrum.i,3)
cd, '/data/jplowman/rh_runs/run_rhsc3d_ARQS_1_0_0'
@initrh_noprompt
.r readall
spec10 = total(spectrum.i,3)
cd, '/data/jplowman/rh_runs/run_rhsc3d_ARQS_1_1_0'
@initrh_noprompt
.r readall
spec11 = total(spectrum.i,3)
cd, '/data/jplowman/rh_runs/run_rhsc3d_ARQS_1_2_0'
@initrh_noprompt
.r readall
spec12 = total(spectrum.i,3)
cd, '/data/jplowman/rh_runs/run_rhsc3d_ARQS_1_3_0'
@initrh_noprompt
.r readall
spec13 = total(spectrum.i,3)

set_plot,'x'
loadct,0

!p.multi=[0,2,4]
plot_image,spec03,max=max(spec00),min=0
plot_image,spec13,max=max(spec00),min=0
plot_image,spec02,max=max(spec00),min=0
plot_image,spec12,max=max(spec00),min=0
plot_image,spec01,max=max(spec00),min=0
plot_image,spec11,max=max(spec00),min=0
plot_image,spec00,max=max(spec00),min=0
plot_image,spec10,max=max(spec00),min=0
