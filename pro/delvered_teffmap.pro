;+
;
; DELVERED_TEFFMAP
;
; Return a low-resolution teff map for a set of chips
; (likely belonging to one exposure).
;
; INPUTS:
;  chstr      Input chip structure.  Must have VX/VY arrays.
;
; RETURNS:
;  teffim     Teff map (36x36)
;
; USAGE:
;  IDL>teffim=delvered_teffmap(chstr)
;
; By D. Nidever,  Dec 2023
;-

function delvered_teffmap,chstr

tilenx = 3600
tileny = 3600

;; Need to add them at high resolution
teffim = fltarr(tilenx/10,tileny/10)
for i=0,n_elements(chstr)-1 do begin
  teffim1 = fltarr(tilenx/10,tileny/10)
  ind1d = polyfillv(chstr[i].vx/10,chstr[i].vy/10,tilenx/10,tileny/10)
  (teffim1)(ind1d) = chstr[i].teff
  teffim += teffim1
endfor
;; Rebin to low-res
teffim = rebin(teffim,36,36)

return,teffim

end
