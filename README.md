# GONG-Calibration-
Data and software used in the GONG magnetogram calibration effort described in Plowman &amp; Berger, Solar Physics, (2020 a,b,c)

These contain the code to run RH on the MURaM cube, create synthetic GONG data from the outputs, compute calibration curves by comparing the synthetic GONG data with the MURaM simulations, and apply the curves to the GONG data.

Contents are as follows:

rh_idl_code: rh_multicube_run_script.pro will split a MURaM cube into slices and run rh on it. rhrun_make_slice_savefile_script.pro will turn those slices into IDL save files. combine_rh_run_script.pro will downsample and combine the slices. Note that running RH is a very computationally intensive process, and that these scripts are not turn-key -- they will require tinkering and rewriting to adapt them to their computing environment. Note also that some of them may be final version of these codes. This directory also contains some miscellaneous pieces of code used elsewhere.

simulation_code: This contains the code to run the actual GONG simulator using the recombined RH output. gong_sim_script_rh3.pro runs the simulator, while gong_synth_fits_script.pro turns the simulator output into GONG-like fits files.

calibration_generation_code: This contains code to compare the simulator output with the MURaM data cube and produce calibration curves. See calibration_curve_script2.pro

Data and example files may be added at a later date, depending on github storage limits (they're large), or can be provided on request.
