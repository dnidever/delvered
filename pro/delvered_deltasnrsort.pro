;+
;
; DELVERED_DELTASNRSORT
;
; Sort the exposures by how much they
; increase S/N map.
;
; INPUTS:
;  expstr   Exposure structure. Must have teff and teffmap.
;
; RETURNS:
;  expstr   Exposure structure sorted by delta S/N.
;
; USAGE:
;  IDL>expstr = delvered_deltasnrsort(expstr)
;
; By D. Nidever,  Dec 2023
;-

function delvered_deltasnrsort,expstr

nexp = n_elements(expstr)

;; keep track of (S/N)**2 = Sum(teff)
;; Neilsen eq. 1
;; S/N ~ f * sqrt(tau * exptime)
;;     where f is the flux of the star
;; therefore, (S/N)**2 ~ teff
    
;; Keep iterating until all the exposures are
;;  sorted by delta (S/N)
flag = 0
count = 0
totdeltasnr_last = -1
WHILE (flag eq 0) do begin

  ;; Sort the exposures
  if count eq 0 then begin
    ;; Sort them by teff the first time
    si = reverse(sort(expstr.teff))
  endif else begin
    ;; Sort them by deltasnr
    si = reverse(sort(expstr.deltasnr))
  endelse
  expstr = expstr[si]

  ;; Calculate delta S/N for each image
  ;;  delta S/N is the S/N before and after adding the image
  sumteffim = fltarr(36,36)  ;; =(S/N)**2
  nmap = n_elements(sumteffim)
  for j=0,nexp-1 do begin
    snr_old = sqrt(sumteffim)
    snr_new = sqrt(sumteffim+expstr[j].teffmap)
    deltasnr = total(snr_new-snr_old)
    expstr[j].deltasnr = deltasnr
    ;; if this exposure covers NEW area, then the new S/N will add linearly
    ;; coverage in already covered areas will add in quadrature
    ;; adding new area is favored!
    sumteffim += expstr[j].teffmap
    expstr[j].fraccovered = total(sumteffim gt 0)/nmap
  endfor

  ;; sum up the deltasnr values
  totdeltasnr = total(expstr.deltasnr,/integer)
  ;;print,count+1,totdeltasnr

  ;; keep iterating until we have maximized sum(delta S/N)
  if totdeltasnr eq totdeltasnr_last or count ge 10 then flag=1
  totdeltasnr_last = totdeltasnr

  count += 1
ENDWHILE

return,expstr

end
