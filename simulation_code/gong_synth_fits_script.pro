; This code turns a simulator run into a GONG file. Not much more than
; a copy/paste of the data into the file, though. It's not strictly
; necessary to do this for all the files in the day directory, really,
; but it is what it is.

gong_data_dir = '../GONG_data/'
simfile_directory = '../GONG_sim_data_rh/'
gong_sim_run_directory = '../gong_sim_runs/'


infilenames = findfile(gong_data_dir+'terif170402/'+'*.fits.gz')

;out_dir = simfile_directory+'ARQScube2/terif170402/'
;restore, gong_sim_run_directory+'gongrh_ARQScube_2_012820.sav'

;out_dir = simfile_directory+'ARQScube2_25degrees/terif170402/'
;restore, gong_sim_run_directory+'gongrh_ARQScube2_25degrees_012820.sav'

;out_dir = simfile_directory+'ARQScube2_45degrees/terif170402/'
;restore, gong_sim_run_directory+'gongrh_ARQScube2_45degrees_012820.sav'

;out_dir = simfile_directory+'ARQScube2_60degrees/terif170402/'
;restore,  gong_sim_run_directory+'gongrh_ARQScube_2_60degrees_012820.sav'

;out_dir = simfile_directory+'ARQScube2_70degrees/terif170402/'
;restore, gong_sim_run_directory+'gongrh_ARQScube2_70degrees_012820.sav'

out_dir = simfile_directory+'ARQScube_75degrees/terif170402/'
restore, gong_sim_run_directory+'gongrh_ARQScube_75degrees_012820.sav'

rsun_cm = 6.957e10
au_cm = 1.496e13
rsun_asec = 3600*(180/!pi)*rsun_cm/au_cm

nx = n_elements(savstr.cube[*,0,0])
ny = n_elements(savstr.cube[0,*,0])
nt = n_elements(savstr.cube[0,0,*])

outfilenames = out_dir+file_basename(infilenames)

nfiles = n_elements(infilenames)

for i=0,nfiles-1 do begin &$
	gongdata = float(readfits(infilenames[i],gongheader)) &$
	writefits, outfilenames[i], transpose(savstr.cube,[1,0,2]), gongheader, /compress &$
endfor
