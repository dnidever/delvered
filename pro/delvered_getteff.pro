;+
;
; DELVERED_GETTEFF
;
; Calculate the effective exposure time.
;
; INPUTS:
;  chstr      Input chip structure.
;
; RETURNS:
;  finalstr   Chip structure with tags added.
;
; USAGE:
;  IDL>chstr = delvered_gettau(chstr)
;
; By D. Nidever,  Dec 2023
;-

function delvered_getteff,chstr

nchips = n_elements(chstr)

;; Calculate effective exposure time
;;----------------------------------
;; Neilsen+2015
;; mag_limiting = mag0 + 1.25*log(tau)
;; where tau is t_eff and mag0 is mag_limiting for tau=1
;; this works well, for r, mag0=22.5

;; teff/tau
;; exptime, seeing, zpterm, background

zpterm = chstr.calib_zpterm
gdes = where(chstr.gain lt 2,ngdes)
if ngdes gt 0 then zpterm[gdes] += 1.55

;; band median zpterm  background (ADU/sec)
;; u     +2.03            0.10
;; g     -0.07            1.29
;; r     -0.39            3.77
;; i     -0.33           13.19
;; z     -0.03           23.46
;; Y     +1.08           17.02
ufilters = ['u','g','r','i','z','Y']
medzpterm = [2.03,-0.07,-0.39,-0.33,-0.03,1.08]
medbackground = [0.10,1.29,3.77,13.19,23.46,17.02]
deltazpterm = zpterm * 0
background_fiducial = fltarr(nchips)
for f=0,n_elements(ufilters)-1 do begin
  ind = where(chstr.filter eq ufilters[f],nind)
  if nind gt 0 then begin
    deltazpterm[ind] = zpterm[ind]-medzpterm[f]
    background_fiducial[ind] = medbackground[f]
  endif
endfor
;; positive deltazpterm is "bad"

eta = exp(-deltazpterm)
fwhm_arcsec_fiducial = 0.9
background_fiducial = 1.0
fwhm_arcsec = chstr.fwhm * 0.235
background = chstr.skymode/chstr.exptime
tau = (eta^2)*(fwhm_arcsec/fwhm_arcsec_fiducial)^(-2) * (background/background_fiducial)^(-1)
add_tag,chstr,'eta',0.0,chstr
add_tag,chstr,'background',0.0,chstr
add_tag,chstr,'tau',0.0,chstr
add_tag,chstr,'teff',0.0,chstr
chstr.eta = eta
chstr.background = background
chstr.tau = tau
chstr.teff = tau*chstr.exptime

;; very good correlation between teff and depth

return,chstr

end
