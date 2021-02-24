.compile /data/jplowman/rh_idl_code_jp/readrempelcube.pro
.compile /data/jplowman/rh_idl_code_jp/unique_elements.pro
;.compile /data/jplowman/rh_idl_code_jp/gettauonesurf.pro

datadir = '/data/jplowman/rempeldata/'
densfile = datadir+'dens_535000.float'
magxfile = datadir+'magx_535000.float'
magyfile = datadir+'magy_535000.float'
magzfile = datadir+'magz_535000.float'
presfile = datadir+'pres_535000.float'
tempfile = datadir+'temp_535000.float'
velxfile = datadir+'velx_535000.float'
velyfile = datadir+'vely_535000.float'
velzfile = datadir+'velz_535000.float'


nx = 1024
nz = 96
ny = 512

dx=48
dy=48
dz=24

x = dindgen(nx)*double(dx)
y = dindgen(ny)*double(dy)
z = reverse(dindgen(nz) * double(DZ))

;cube = reverse(transpose(readvarcube(nx,nz,ny,densfile),[0,2,1]),3)

if(n_elements(cube) eq 0) then begin &$
	cd,'/data/jplowman/rh_runs/run_rhsc3d_ARQS_1_1_0' &$
	@initrh
	.r readall
	.r /data/jplowman/rh_idl_code_jp/gettauonesurf.pro
endif

cd,'/data/jplowman/rh_idl_code_jp/'

cube = atmos.t
xderiv = (shift(cube,1,0,0)-shift(cube,-1,0,0))*0.5/dx
yderiv = (shift(cube,0,1,0)-shift(cube,0,-1,0))*0.5/dy
zderiv = -(shift(cube,0,0,1)-shift(cube,0,0,-1))*0.5/dz

laplac = sqrt(xderiv*xderiv+yderiv*yderiv+zderiv*zderiv)


cube_tauone = fltarr(nx,ny)
laplac_tauone = fltarr(nx,ny)
for i=0,nx-1 do for j=0,ny-1 do cube_tauone[i,j] = interpol(reform(cube[i,j,*]),z,zcont[i,j])
for i=0,nx-1 do for j=0,ny-1 do laplac_tauone[i,j] = interpol(reform(laplac[i,j,*]),z,zcont[i,j])


zinds0 = fltarr(nx,ny)
zinds = lonarr(nx,ny)
zc2 = fltarr(nx,ny)
for i=0,nx-1 do for j=0,ny-1 do zinds0[i,j] = interpol(findgen(nz),z,zcont[i,j])
zinds = round(zinds0)
cube_tauone2 = fltarr(nx,ny)
for i=0,nx-1 do for j=0,ny-1 do cube_tauone2[i,j] = cube[i,j,zinds[i,j]]

npm = 10
ii = 52
cube_iia = fltarr(2*npm+1,ny)
for i=-npm,npm do for j=0,ny-1 do cube_iia[i+npm,j] = interpol(reform(cube[ii,j,*]),z,zcont[ii,j]-i*dz)
;for i=-npm,npm do for j=0,ny-1 do cube_iia[i+npm,j] = cube[ii,j,zinds[ii,j]+i]

!p.multi=[0,0,2]
plot_image,reverse(alog(reform(cube[ii,*,[25:45]])),2),/nosquare,scale=[dy,dz],origin=[0,z[45]]
oplot,y,zcont[ii,*]
plot_image,reverse(alog(transpose(cube_iia)),2),/nosquare,scale=[dy,dz],origin=[0,-npm*dz]
