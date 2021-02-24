atmos_out = '/rhdata/jplowman/local_dynamo/dyn_ng_8km_053000/Rempel_dyn_053000_rebin.atmos'
B_out = '/rhdata/jplowman/local_dynamo/dyn_ng_8km_053000/Rempel_dyn_053000_rebin.B'

atmos2 = read3datmos(atmos_out, bfile = B_out)

deltan = atmos2.nh[*,*,1:atmos2.nz-1]-atmos2.nh[*,*,0:atmos2.nz-2]
for i=0,atmos2.nz-2 do deltan[*,*,i] = deltan[*,*,i]/(atmos2.z[i+1]-atmos2.z[i])

dnmax = min(deltan,imax,dimension=3,/nan)
dnmax_indices = array_indices(deltan,imax)
ixmax0 = dnmax_indices[0,*]
iymax0 = dnmax_indices[1,*]
izmax0 = dnmax_indices[2,*]

izmax = reform(izmax0,atmos2.nx,atmos2.ny)
zvals = atmos2.z[izmax]

bvals = reform(atmos2.b[ixmax0,iymax0,izmax0],atmos2.nx,atmos2.ny)

set_plot,'z'
device,set_resolution=[atmos2.nx,atmos2.ny]


for i=0,atmos2.nz-1 do begin &$
	tvscl,reform(atmos2.b[*,*,i]*cos(atmos2.gamma[*,*,i])) &$
	write_jpeg,'bvszframes/frame_'+strtrim(string(i),2)+'.jpg',tvrd(),quality=100 &$
endfor