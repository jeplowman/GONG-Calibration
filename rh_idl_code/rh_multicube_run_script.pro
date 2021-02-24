.compile readrempelcube.pro
;.compile /rhdata/jplowman/rh_v2/rhf1d/idl_codes_jp/read3datmos.pro
.compile make_wavetable_file.pro

idlcode_directory = '/rhdata/jplowman/rh_v2/rhf1d/idl_codes_jp/'
data_directory = '/rhdata/jplowman/rempel_QS_sunspot/'

outfilehead = data_directory+'rempel_largescale_rhcube'
outputFile = data_directory+'rempel_largescale_rhcube.atmos'
B_outputFile = data_directory+'rempel_largescale_rhcube.B'

rh_templatedirectory = '/rhdata/jplowman/rh_v2/rhsc3d/template_run_rhsc3d_ARQS'
rh_rundirectoryhead = '/rhdata/jplowman/rh_v2/rhsc3d/run_rhsc3d_ARQS'
wavetable_filename = 'wavetable.wave'
infofilename = '/rhdata/jplowman/rh_v2/rhf1d/idl_codes_jp/run_rhsc3d_ARQS_info.sav'

vfiles = data_directory+['velx_535000.float','vely_535000.float','velz_535000.float']
bfiles = data_directory+['magx_535000.float','magy_535000.float','magz_535000.float']
tfile = data_directory+'temp_535000.float'
dfile = data_directory+'dens_535000.float'


lambdamin = 676.730
lambdamax = 676.830
nlambdas = 48
nlambdabins = 1
nperlambdabin = nlambdas/nlambdabins

lambdas = lambdamin+(lambdamax-lambdamin)*findgen(nlambdas)/(nlambdas-1.0)

dx=48
dy=48
dz=24

n0=2048
n1=96
n2=2048

xcrop = 1024
ycrop = 512
nxcrop = n0/xcrop
nycrop = n0/ycrop

for i=0,nxcrop-1 do begin &$
	for j=0,nycrop-1 do begin &$
		for k=0,nlambdabins-1 do begin &$
			outrange_lo=[i*xcrop,j*ycrop,0] &$
			outrange_hi=[(i+1)*xcrop,(j+1)*ycrop,n1]-1 &$
			readrempelcube, vfiles, bfiles, tfile, dfile, dx, dy, dz, n0, n1, n2, outputFile=outputFile, B_outputFile=B_outputFile, $
					outrange_lo=outrange_lo, outrange_hi=outrange_hi &$
			rhdir = rh_rundirectoryhead + '_' + strtrim(string(i),2) + '_' + strtrim(string(j),2)+'_'+strtrim(string(k),2)+'/' &$
			print, rhdir &$
			spawn,'cp -r '+rh_templatedirectory+' '+rhdir &$
			cd, rhdir &$
			make_wavetable_file, lambdamin, lambdamax, nlambda, filename=wavetable_filename, lambdas=lambdas[k*nperlambdabin:((k+1)*nperlambdabin-1)] &$
			spawn,'../rhsc3d > rhsc3d_log.txt' &$
			spawn,'../solveray > solveray_log.txt' &$
			cd, idlcode_directory &$
	endfor &$
endfor

save,lambdamin,lambdamax,nlambdas,nlambdabins,nperlambdabin,lambdas,dx,dy,dz,n0,n1,n2,xcrop,ycrop,nxcrop,nycrop,data_directory,rh_rundirectoryhead,filename = infofilename