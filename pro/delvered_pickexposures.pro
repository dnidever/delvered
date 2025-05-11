;+
;
; DELVERED_PICKEXPOSURES
;
; Pick a subset of exposures to process.
;
; INPUTS:
;  chstr          Input chip structure.
;  =maxexposures  The maximum number of exposures to pick.
;                   Default is 50.
;
; RETURNS:
;  finalchstr     Final chip structure to use.
;
; USAGE:
;  IDL>finalchstr = delvered_pickexposures(chstr)
;
; By D. Nidever,  Dec 2023
;-

function delvered_pickexposures,chstr,maxexposures=maxexposures

;; Want best exposure
;; -use depth as that metric
;; -also want good data in each band
;; -maximum number of exposures that we can process
;; -prioritize exposures that completely overlap the brick

nchips = n_elements(chstr)

if n_elements(maxexposure) eq 0 then maxexposures = 50  ;; default

;; dao_depth is instrumental depth
add_tag,chstr,'depth',0.0,chstr
chstr.depth = chstr.dao_depth + 2.5*alog10(chstr.exptime) - chstr.calib_zpterm - chstr.apcor

;; Calculate effective exposure time
;;----------------------------------
chstr = delvered_getteff(chstr)


;; Get fraction of brick overlapped by each chip
nchips = n_elements(chstr)
tilenx = 3600
tileny = 3600
add_tag,chstr,'fracoverlap',0.0,chstr
for i=0,nchips-1 do begin
  vx = chstr[i].vx
  vy = chstr[i].vy
  ;; Calculate the fraction of the brick that this chip overlaps
  if max(vx) ge 0 and min(vx) le tilenx-1 and max(vy) ge 0 and min(vy) le tileny-1 then begin
    ind1d = polyfillv(vx,vy,tilenx,tileny)
    chstr[i].fracoverlap = n_elements(ind1d)/(tilenx*float(tileny))  ;; percent of the brick area overlapped
  endif
endfor

;; diagnostic plot
add_tag,chstr,'mnx',0.0,chstr
add_tag,chstr,'mny',0.0,chstr
chstr.mnx = mean(chstr.vx,dim=1)
chstr.mny = mean(chstr.vy,dim=1)
ind = where(chstr.filter eq 'r')
chstr1 = chstr[ind]
;plot,[0],[0],xr=[-2000,5000],yr=[-2000,5000]
;oplot,[0,3600,3600,0,0],[0,0,3600,3600,0],co=250,thick=5
;colors = scale(chstr1.fracoverlap,[0.0,0.5],[50,250])
;for j=0,n_elements(chstr1)-1 do oplot,[chstr1[j].vx,chstr1[j].vx[0]],[chstr1[j].vy,chstr1[j].vy[0]],co=colors[j]
;for j=0,n_elements(chstr1)-1 do oplot,[chstr1[j].mnx],[chstr1[j].mny],co=colors[j],ps=1

;; Make exposure structure
ue = uniq(chstr.expnum,sort(chstr.expnum))
uexp = chstr[ue].expnum
nexposures = n_elements(uexp)
dum = {expnum:'',filter:'',exptime:0.0,obsdate:'',mjd:0.0d0,nchips:0L,$
       airmass:0.0,zpterm:0.0,nmeasavg:0L,nmeas:0L,eta:0.0,background:0.0,fwhm:0.0,$
       depth:99.99,tau:0.0,teff:0.0,fracoverlap:0.0}
expstr = replicate(dum,nexposures)
;; Loop over the exposures
for i=0,nexposures-1 do begin
  ind = where(chstr.expnum eq uexp[i],nind)
  expstr[i].expnum = uexp[i]
  expstr[i].filter = chstr[ind[0]].filter
  expstr[i].exptime = chstr[ind[0]].exptime
  expstr[i].obsdate = chstr[ind[0]].utdate+'T'+chstr[ind[0]].uttime
  expstr[i].mjd = date2jd(expstr[i].obsdate,/mjd)
  expstr[i].nchips = nind
  expstr[i].airmass = mean([chstr[ind].airmass]) 
  expstr[i].zpterm = mean(chstr[ind].calib_zpterm)
  expstr[i].nmeasavg = mean([chstr[ind].dao_nsources])
  expstr[i].nmeas = total([chstr[ind].dao_nsources],/integer)
  expstr[i].eta = median([chstr[ind].eta])
  expstr[i].background = median([chstr[ind].background])
  expstr[i].fwhm = median([chstr[ind].fwhm])
  expstr[i].depth = median([chstr[ind].depth])
  expstr[i].tau = median([chstr[ind].tau])
  expstr[i].teff = median([chstr[ind].teff])
  expstr[i].fracoverlap = total([chstr[ind].fracoverlap])  
endfor

;; use a low-resolution map
;; for each pixel sum up the "depth metric" for all exposures
;; only add exposures that increase this by a decent amount
;; sort by depth first

;; want to maximize depth but also coverage of the brick

;; pick MAXEXPOSURES per filter
;; then trim down to MAXEXPOSURES total afterwards

;; we also want the BEST exposures ACROSS all the filters!!!
;; maybe save ~5-10 exposures for the absolute best (exclude u-band)


;; Get the teff map for each exposure
print,'Getting teff maps for ',strtrim(nexposures,2),' exposures'
add_tag,expstr,'teffmap',fltarr(36,36),expstr
eindex = create_index(chstr.expnum)
for e=0,nexposures-1 do begin
  chindex = eindex.index[eindex.lo[e]:eindex.hi[e]]
  teffim = delvered_teffmap(chstr[chindex])
  expstr[e].teffmap = teffim
endfor

add_tag,expstr,'deltasnr',0.0,expstr
add_tag,expstr,'fraccovered',0.0,expstr
add_tag,expstr,'picked',0,expstr

;; Find the best exposures overall
gd = where(expstr.filter ne 'u',ngd)
deepexpstr = expstr[gd]
deepexpstr = delvered_deltasnrsort(deepexpstr)

;; Now get the best exposures for each band
uf = uniq(chstr.filter,sort(chstr.filter))
ufilters = chstr[uf].filter
nfilters = n_elements(ufilters)
;; Loop over the filters
for f=0,nfilters-1 do begin
  ind = where(expstr.filter eq ufilters[f],nexp)
  print,strtrim(f+1,2),' ',ufilters[f],' ',strtrim(nexp,2)
  if nind gt 0 then begin
    fexpstr = expstr[ind]
    ;; Sort by delta S/N
    ;;   this adds deltasnr and fraccovered columns
    fexpstr = delvered_deltasnrsort(fexpstr)
    ;; Stuff back in, but in different order
    expstr[ind] = fexpstr
  endif
endfor


;; Put together the final list of exposures
;;-----------------------------------------
;; Add the 10 deepest first
;; then add the 5 best from each filter
;;  if there are left, then continue to add 2 more per filter until we
;;  get to MAXEXPOSURES

;; Add the 10 deepest first
deepexpnum = deepexpstr[0:9<n_elements(deepexpstr)-1].expnum
match,deepexpnum,expstr.expnum,ind1,ind2,/sort
expstr[ind2].picked = 1

flag = 0
count = 0
WHILE (flag eq 0) do begin
  ;; How many to pick
  npick = 1
  
  ;; Filter loop
  for f=0,nfilters-1 do begin
    npicked = total(expstr.picked eq 1,/integer)
    if npicked ge maxexposures then break
    ind = where(expstr.filter eq ufilters[f] and expstr.picked eq 0,nind)
    si = reverse(sort(expstr[ind].deltasnr))  ;; sort by delta S/N
    ind_picked = si[0:npick-1]
    expstr[ind[ind_picked]].picked = 1
  endfor

  npicked = total(expstr.picked eq 1,/integer)
  ;;print,count,npicked
  if npicked ge maxexposures then flag=1

  count += 1
ENDWHILE
ind = where(expstr.picked eq 1,nind)
finalexpstr = expstr[ind]

;; Now find the chipstr elements that belong to these exposures
eindex = create_index(chstr.expnum)
match,finalexpstr.expnum,eindex.value,ind1,ind2,/sort
undefine,chind
for i=0,n_elements(ind1)-1 do begin
  ind = eindex.index[eindex.lo[ind2[i]]:eindex.hi[ind2[i]]]
  push,chind,ind
endfor
finalchstr = chstr[chind]

;; Print information about the final picked exposures
gd = where(expstr.picked eq 1,npicked)
print,strtrim(npicked,2),' exposures picked'
uf = uniq(expstr[gd].filter,sort(expstr[gd].filter))
ufilters = expstr[gd[uf]].filter
nfilters = n_elements(ufilters)
;; Loop over the filters
for f=0,nfilters-1 do begin
  ind = where(expstr[gd].filter eq ufilters[f],nexp)
  print,strtrim(f+1,2),' ',ufilters[f],' ',strtrim(nexp,2)
endfor

;; Delete the extra tags that we added
;;  otherwise this causes problems in delvered_jointbrickcats.pro
todel = ['DEPTH','ETA','BACKGROUND','TAU','TEFF','FRACOVERLAP','MNX','MNY']
finalchstr = remove_tags(finalchstr,todel)

return,finalchstr

end
