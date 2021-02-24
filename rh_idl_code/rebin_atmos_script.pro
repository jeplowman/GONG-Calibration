.compile write3datmos.pro
.compile write3datmos0.pro
.compile sanitize_array.pro
.compile rebin_atmos.pro

;atmos_in = '/rhdata/serena/local_dynamo/Plage_200G_6x3Mm_8km_ng_018000/Rempel_plage_200G_018000.atmos'
;B_in = '/rhdata/serena/local_dynamo/Plage_200G_6x3Mm_8km_ng_018000/Plage_200G_6x3Mm_8km_ng_018000.B'

;atmos_out = '/rhdata/jplowman/local_dynamo/Plage_200G_6x3Mm_8km_ng_018000/Rempel_plage_200G_018000_rebin4x4x2.atmos'
;B_out = '/rhdata/jplowman/local_dynamo/Plage_200G_6x3Mm_8km_ng_018000/Plage_200G_6x3Mm_8km_ng_018000_rebin4x4x2.B'

atmos_in = '/rhdata/jplowman/local_dynamo/dyn_ng_8km_053000/Rempel_dyn_053000.atmos'
B_in = '/rhdata/jplowman/local_dynamo/dyn_ng_8km_053000/Rempel_dyn_053000.B'

;atmos_out = '/rhdata/jplowman/local_dynamo/dyn_ng_8km_053000/Rempel_dyn_053000_rebin.atmos'
;B_out = '/rhdata/jplowman/local_dynamo/dyn_ng_8km_053000/Rempel_dyn_053000_rebin.B'

atmos_out = '/rhdata/jplowman/local_dynamo/dyn_ng_8km_053000/Rempel_dyn_053000_rebin4x4x2.atmos'
B_out = '/rhdata/jplowman/local_dynamo/dyn_ng_8km_053000/Rempel_dyn_053000_rebin4x4x2.B'

rebin_atmos, atmos_in, atmos_out, nx=192, ny=192, nz=65, bfile=B_in, boutfile = B_out, atmos=atmos

;atmos2 = read3datmos(atmos_out, bfile = B_out)

;atmos3 = read3datmos(atmos_other)