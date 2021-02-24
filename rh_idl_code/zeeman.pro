FUNCTION gLande, S, L, J

  IF (J EQ 0.0) THEN $
   gL = 0.0 $
  ELSE $
   gL = 1.5 + (S*(S + 1.0) - L*(L + 1)) / (2.0*J*(J + 1.0))

  return, gL
END

FUNCTION zeemanstrength, Ju, Mu, Jl, Ml

  q  = long(Ml - Mu);
  dJ = long(Ju - Jl);

  CASE (dJ) OF
    0: BEGIN
      CASE (q) OF
         0: s = 2.0 * Mu^2
        -1: s = (Ju + Mu) * (Ju - Mu + 1.0)
         1: s = (Ju - Mu) * (Ju + Mu + 1.0)
      ENDCASE
    END

    1: BEGIN 
      CASE (q) OF
         0: s = 2.0 * (Ju^2 - Mu^2)
        -1: s = (Ju + Mu) * (Ju + Mu - 1.0)
         1: s = (Ju - Mu) * (Ju - Mu - 1.0)
      ENDCASE
    END 

    -1: BEGIN
      CASE (q) OF
         0: s = 2.0 * (Ju + 1.0)^2 - Mu^2
        -1: s = (Ju - Mu + 1.0) * (Ju - Mu + 2.0)
         1: s = (Ju + Mu + 1.0) * (Ju + Mu + 2.0)
      ENDCASE
    END
  ENDCASE

  return, s
END

FUNCTION zeemanshift, Su, Lu, Ju, Mu, Sl, Ll, Jl, Ml

  gu = gLande(Su, Lu, Ju)
  gl = gLande(Sl, Ll, Jl)

  return, gl*Ml - gu*Mu
END

FUNCTION get_splittings, Su, Lu, Ju, Sl, Ll, Jl

  components = replicate({q: 0L,  shift: 0.0,  strength: 0.0}, $
                         (2*Ju+1) * (2*Jl+1))

  Ncomponent = 0
  FOR Ml=-Jl, Jl, 1 DO BEGIN
    FOR Mu=-Ju, Ju, 1 DO BEGIN
      IF (abs(Mu - Ml) LE 1) THEN BEGIN
        components[Ncomponent].q = Ml - Mu
        components[Ncomponent].shift = $
         zeemanshift(Su, Lu, Ju, Mu, Sl, Ll, Jl, Ml)
        components[Ncomponent].strength = zeemanstrength(Ju, Mu, Jl, Ml)
        Ncomponent = Ncomponent + 1
      ENDIF
    ENDFOR
  ENDFOR

  norm = fltarr(3)
  FOR n=0, Ncomponent-1 DO $
    norm[components[n].q + 1] = $
     norm[components[n].q + 1] + components[n].strength

  FOR n=0, Ncomponent-1 DO BEGIN
    IF (norm[components[n].q + 1] GT 0.0) THEN $
     components[n].strength = $
     components[n].strength / norm[components[n].q + 1]
  ENDFOR

  return, components[0:Ncomponent-1]
END

FUNCTION effectiveLande, g_u, J_u, g_l, J_l

  return,  0.5*(g_u + g_l) + 0.25*(g_u - g_l) * $
   (J_u*(J_u + 1.0) - J_l*(J_l + 1.0))
END

PRO plot_splittings, Su, Lu, Ju, Sl, Ll, Jl

  colors = [16B, 100B, 225B]

  c = get_splittings(Su, Lu, Ju, Sl, Ll, Jl)
  c_pi = c[where(c.q EQ 0)]
  c_sp = c[where(c.q EQ 1)]
  c_sm = c[where(c.q EQ -1)]

  s_max = max(c_sp.strength)
  s_min = max(c_pi.strength)
  bottom = -s_min - (s_min + s_max) / 2.0

  sh_max = max(c.shift)

  plot, /NODATA, c.shift, c.strength, $
   XRANGE=[-sh_max, sh_max]*1.5, XSTYLE=1, $
   YRANGE=[bottom, s_max*1.2], YSTYLE=1, $
   XTITLE='Shift [Larmor units]', YTITLE='Strength'
  oplot, !X.CRANGE, [0.0, 0.0], /LINE
  oplot, [0.0, 0.0], !Y.CRANGE, /LINE

  thick = 3.0

  FOR n=0, n_elements(c_pi)-1 DO $
   oplot, [c_pi[n].shift, c_pi[n].shift], [0.0, c_pi[n].strength] + bottom, $
   COLOR=colors[0], THICK=thick

  FOR n=0, n_elements(c_sp)-1 DO $
   oplot, [c_sp[n].shift, c_sp[n].shift], [0.0, c_sp[n].strength], $
   COLOR=colors[1], THICK=thick
   
  FOR n=0, n_elements(c_sm)-1 DO $
   oplot, [c_sm[n].shift, c_sm[n].shift], [0.0, -c_sm[n].strength], $
   COLOR=colors[2], THICK=thick

  g_eff = effectiveLande(gLande(Su, Lu, Ju), Ju, gLande(Sl, Ll, Jl), Jl)

  arrow, g_eff, -0.3*s_min, g_eff, 0.0, /DATA, /SOLID
  rhannotate, TEXT=string(FORMAT='("gL!Deff!N = ", F6.3)', g_eff), $
   g_eff, -0.32*s_min, CHARSIZE=0.8

END
