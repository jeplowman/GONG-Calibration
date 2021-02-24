.compile multi_tile_2.pro
.compile image_coord_transforms.pro
.compile sanitize_array.pro
.compile gong_sim_polarimetry.pro
.compile gong_sim_pixel_psf.pro
.compile gong_wave_resp.pro
.compile gong_wave_resp2.pro
.compile gong_sim_detectorplane.pro
.compile gong_sim_camintegration.pro
.compile gong_simulator.pro

au_m_toasec = (3600.0d0*180.0d0/!pi)/1.496d11
asec_cm = 1.496d13/(3600.0d0*180.0d0/!pi)
;dithermean = 182.0d0
dithermean = 183.0d0
npbin = 1200
rsun_cm = 6.957d10
au_cm = 1.496d13
rsun_asec = 3600.0d0*(180.0d0/!pi)*rsun_cm/au_cm
c=3d10
;atmo_w = 5.0 ; Smoothing due to atmospheric seeing (old value).
atmo_w = 3.6 ; Smoothing due to atmospheric seeing.

cm2tom2 = double(1.0d-4)
ergtoj = double(1.0d-7)
bsun = 1366.0*cm2tom2/ergtoj

; This factor is to sanitize nonphysically large values in the input datacube.
; limit is set to equivalent of roughly 10 times the typical solar value
; for the wavelength bin and pixel size:
fluxlimit_fac = double(10.0d0)

;scales=[[1,1],[1,1]] ; Uncomment this to run the whole disk at GONG native scale

;-----------Most recent standard runs-----------------------;

; GONG simulator runs used hardwired directory structure; these edits to use rh_run_directory
; and gong_sim_run_directory instead are untested.
rh_run_directory = '../rh_runs/'
gong_sim_run_directory = '../gong_sim_runs/'
gong_data_directory = '../GONG_data/'

; Uncomment the corresponding block below to run simulator for given angle:

; 0 degrees.
;specfile = rh_run_directory+'run_rhsc3d_ARQS_mu00_spectrum.sav'
;xb_fac = 4
;yb_fac = 4
;outfilename = gong_sim_run_directory+'gongrh_ARQScube_2_012820.sav'

; 25 degrees.
;specfile = rh_run_directory+'run_rhsc3d_ARQSmu42_spectrum.sav'
;xb_fac = 4
;yb_fac = 4
;outfilename = gong_sim_run_directory+'gongrh_ARQScube2_25degrees_012820.sav'

; 45 degrees.
;specfile = rh_run_directory+'runs_rhsc3d_ARQS_6768_2/run_rhsc3d_ARQSmu71_spectrum.sav'
;xb_fac = 4
;yb_fac = 8
;outfilename = gong_sim_run_directory+'gongrh_ARQScube2_45degrees_012820.sav'

; 60 degrees.
;specfile = rh_run_directory+'run_rhsc3d_ARQS_mu87_spectrum.sav'
;xb_fac = 4
;yb_fac = 8
;outfilename = 'save/gongrh_ARQScube_2_60degrees_012820.sav'

; 70 degrees.
;specfile = rh_run_directory+'run_rhsc3d_ARQS_mu94_spectrum.sav'
;xb_fac = 4
;yb_fac = 16
;outfilename = gong_sim_run_directory+'gongrh_ARQScube2_70degrees_012820.sav'

; 75 degrees. Note that an older rh run is used. IIRC, though, that probably doesn't matter...
specfile = rh_run_directory+'run_rhsc3d_ARQS75degree_spectrum.sav'
xb_fac = 4
yb_fac = 16
outfilename = gong_sim_run_directory+'gongrh_ARQScube_75degrees_012820.sav'

nopsf=0
nonoise=0
nodisc=0

gongfile = gong_data_directory+'terif170402/TE170402154716.fits.gz'
gongdata = float(readfits(gongfile,gongheader))

darkfile = gong_data_directory+'tecqb170402/tecqb170402t1452/avgdrk.fits.gz'
flatfile = gong_data_directory+'tecqb170402/tecqb170402t1452/flat.fits.gz'

 ; Should we be transposing here, or leave in native GONG orientation? Currently transposes
 ; back to native orientation in gong_synth_fits_script.pro...
gongdark = transpose(float(readfits(darkfile,darkheader)))
gongflat = transpose(float(readfits(flatfile,flatheader)),[1,0,2])
gongflat = reform(gongflat[*,*,1])
gongflat = gongflat/median(gongflat)

ndarkframes = float(sxpar(darkheader,'ACCUM'))
; Due to frames being dropped, there are not 600 of each 6 1-minute binned frames, but rather 'ACCUM' from header(?):
; I think dithermean should be subtracted from these darks (like below), but that results in negative darks.
; Maybe there's some other offset in the ADU conversion? Either that, or:
;	They're removed as part of the dark acquisition, before being stored as raws
;	There's no dither applied at all during the raw acquisition
;	ACCUM does not actually reflect the number of dark frames and it is actually smaller (<500)
; Currently I assume they're removed as part of initial raw image acquisition.
;gongdark = gongdark/ndarkframes-dithermean/8.0
gongdark = gongdark/ndarkframes

gongny = sxpar(gongheader,'NAXIS1')
gongnx = sxpar(gongheader,'NAXIS2')
yc_pix = sxpar(gongheader,'X_CENTER')
xc_pix = sxpar(gongheader,'Y_CENTER')
scale = [2.5,2.5]
rollcen = [0.0,0.0]
rota = 0.0
px_refval = 0.5*[gongnx,gongny]; [xc_pix,yc_pix] ; These center values don't make sense as sun center pixels - they're way off!?
gongcenter = find_image_center(0.5*[gongnx,gongny], scale, px_refval, refval=refval, rollcen=rollcen, px_rollcen=px_rollcen, rota=rota)
gongorigin = find_image_center([0,0], scale, px_refval, refval=refval, rollcen=rollcen, px_rollcen=px_rollcen, rota=rota)
seed = 4390682934872

restore,specfile
if(is_struct(spectrum_geometry)) then geometry=spectrum_geometry
delvar,ray
delvar,atmos
delvar,zcont
delvar,zcont01
delvar,zcont03
delvar,zcont07
delvar,zcont20
delvar,zcont30
delvar,zcore
delvar,zcore01
delvar,zcore03
delvar,zcore07
delvar,zcore20
delvar,zcore30

nstokes = 4 ; Always 4 Stokes components.
nw0 = n_elements(spectrum.lambda)
waves = double(1.0d-7*spectrum.lambda)

dx_asec = geometry.dx*au_m_toasec
dy_asec = geometry.dy*au_m_toasec

stokesi=spectrum.i
stokesv=spectrum.stokes_v
delvar,spectrum
im_test = multi_tile_2(reform(stokesi[*,*,0]),dx_asec,dy_asec,xb_fac=xb_fac,yb_fac=yb_fac,dx_out=dx0,dy_out=dy0)
nx0 = n_elements(im_test[*,0])
ny0 = n_elements(im_test[0,*])
datacube = fltarr(nx0,ny0,nw0,nstokes)
;for i=0,nw0-1 do datacube[*,*,i,0] = multi_tile_2(reform(stokesi[*,*,i+1]),dx_asec,dy_asec,xb_fac=xb_fac,yb_fac=yb_fac,dx_out=dx0,dy_out=dy0, scales=scales)
;for i=0,nw0-1 do datacube[*,*,i,3] = multi_tile_2(reform(stokesv[*,*,i+1]),dx_asec,dy_asec,xb_fac=xb_fac,yb_fac=yb_fac,dx_out=dx0,dy_out=dy0, scales=scales)
for i=0,nw0-1 do datacube[*,*,i,0] = multi_tile_2(reform(stokesi[*,*,i]),dx_asec,dy_asec,xb_fac=xb_fac,yb_fac=yb_fac,dx_out=dx0,dy_out=dy0, scales=scales)
for i=0,nw0-1 do datacube[*,*,i,3] = multi_tile_2(reform(stokesv[*,*,i]),dx_asec,dy_asec,xb_fac=xb_fac,yb_fac=yb_fac,dx_out=dx0,dy_out=dy0, scales=scales)
delvar,stokesi
delvar,stokesv

xref = 0.0
yref = 0.0
ixref = 0.5*nx0
iyref = 0.5*ny0
dx0_cm = dx0*asec_cm
dy0_cm = dy0*asec_cm
delvar,geometry

rhscale = [dx0,dy0]
rhorigin = find_image_center([0.0,0.0],rhscale,[ixref,iyref],refval=[xref,yref])

gongw0 = 676.7785d0

nw = nw0

dx = scale[0]; 2*dx0
dy = scale[1]; 2*dy0
nx = gongnx; 512
ny = gongny; 512

dark=gongdark;/npbin

flat = fltarr(nx,ny,6)
for i=0,5 do flat[*,*,i]=gongflat
delvar,gongflat

gongarea = !pi*1.4^2*0.85*0.05 ; GONG effective area

; This factor converts mean spectral radiance in J S^-1 m^-2 Hz^-1 Sr^-1 to Erg s^-2 Asun^-2 cm^-3.
; That is, given a spectral radiance value on the sky, how much power would it produce over a solar radius
; at Earth per unit wavelength and aperture (cm).
fluxconv_tot = (!pi*double(rsun_cm)^2)*double(c)/double(median(waves))^2/double(au_cm)^2*cm2tom2/ergtoj

; This factor converts spectral radiance at the Sun (in J S^-1 m^-2 Hz^-1 Sr^-1) to received on Earth (Erg cm^-3 arcsecond^-2/s).
; That is, power received per unit aperture and wavelength (cm) per square arcsecond patch on the sky.
fluxconv_px = double(asec_cm)*double(asec_cm)*double(c)/double(median(waves))^2/double(au_cm)^2*cm2tom2/ergtoj
print,'Estimaded total flux would be' + string(double(median(datacube[*,*,*,0],/double))*double(fluxconv_tot)*double(median(waves))*1.0d-7*1.0d4) + ' W m^-2'

help,fluxconv_px
help,datacube

datacube *= float(fluxconv_px)

; This limit is to sanitize nonphysically large values in the input datacube.
; limit is set to equivalent of roughly 10 times the typical solar value
; for the wavelength bin and pixel size:
fluxlimit = fluxlimit_fac*fluxconv_px*bsun/fluxconv_tot/median(waves)


wavelengths = double(waves)

print,waves
help,datacube

if(n_elements(mueller) eq 0) then begin &$
	mueller = dblarr(2,4,4) &$
	mueller[0,*,*] = diag_matrix([1,0,0,1]) &$ ; I+V assumed to pass I & V unchanged, zero everything else.
	mueller[1,*,*] = diag_matrix([1,0,0,-1]) &$ ; I-V assumed to pass I unchanged, flip sign of V, zero everything else.
endif

; First apply Mueller matrices, since it's a quick simplifying operation:
specint_ipv = (gong_sim_polarimetry(datacube,wavelengths,reform(mueller[0,*,*])) > 0.0) < fluxlimit
specint_imv = (gong_sim_polarimetry(datacube,wavelengths,reform(mueller[1,*,*])) > 0.0) < fluxlimit
delvar,datacube

xa = dx*(findgen(nx)#(1+fltarr(ny))-0.5*(nx-1))
ya = dy*((1+fltarr(nx))#findgen(ny)-0.5*(ny-1))

ra = sqrt(xa*xa+ya*ya)
mask = ra lt rsun_asec


specflux_ipv = fltarr(nx,ny,nw)
specflux_imv = fltarr(nx,ny,nw)
; Now Apply PSFs and pixelization:
for i=0,nw-1 do begin &$
	print,'Applying pixel PSFs to wavelength '+strtrim(string(i+1),2)+' of '+strtrim(string(nw),2) &$
	specflux_ipv[*,*,i] = mask*gong_sim_pixel_psf(reform(specint_ipv[*,*,i]),rhorigin,dx0,dy0,gongorigin,dx,dy,nx,ny,psf=psf,wpsf=wpsf,nopsf=nopsf, atmo_w=atmo_w) &$
	specflux_imv[*,*,i] = mask*gong_sim_pixel_psf(reform(specint_imv[*,*,i]),rhorigin,dx0,dy0,gongorigin,dx,dy,nx,ny,psf=psf,wpsf=wpsf,nopsf=nopsf, atmo_w=atmo_w) &$
endfor

; Detector plane / Cam integration happens last:
cube = gong_sim_camintegration(specflux_ipv, specflux_imv, dark, flat, wavelengths, exptime, gaintabx, gaintaby, rdn=rdn, qe=qe, seed=seed, llo=llo,/do_plot,plotfile='plots/gong_camintegration_statusplot_rh.eps', nonoise=nonoise, dithermeans=dithermeans, npbin=npbin, nodisc=nodisc)

; Remove the dither offset to be self-consistent:
for i=0,5 do cube[*,*,i] -= npbin[i]*dithermeans[i]/8.0D0

savstr = {cube:cube,wavelengths:wavelengths,dx:dx,dy:dy,dark:dark,flat:flat, xb_fac:xb_fac, yb_fac:yb_fac, specfile:specfile, scales:scales, npbin:npbin, dithermeans:dithermeans, wpsf:wpsf, atmo_w:atmo_w}
save,savstr,filename=outfilename

