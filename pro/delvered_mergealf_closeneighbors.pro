pro delvered_mergealf_closeneighbors,nmg,obj,meas,dcr,newobj,newmeas,logfile=logfile

;; This looks at stars that are close neighbors and merges their results

nobj = n_elements(obj)
nmeas = n_elements(meas)

ndet = obj.ndetu+obj.ndetg+obj.ndetr+obj.ndeti+obj.ndetz+obj.ndety
obj.objid = strtrim(obj.objid,2)
meas.objid = strtrim(meas.objid,2)

if n_elements(logfile) eq 0 then logfile=-1

;; Merge anything within 0.5 X FWHM.  This distance is input as DCR
;; matchall_sph.pro wants dcr in degrees
;matches = MATCHALL_SPH(obj.ra, obj.dec, obj.ra, obj.dec, dcr/3600, nmatches)
matches = MATCHALL_2D(nmg.x,nmg.y,nmg.x,nmg.y,dcr,nmatches)
bd = where(nmatches gt 1,nbd)

;; Loop over the objects that have close matches
mstr = replicate({num:0L,objid:'',primind:0LL,ngroup:0L,objind:lonarr(10)-1,nmeas:lonarr(10)-1,$
                  mag:fltarr(10)+99.99,fluxratio:fltarr(10)+999999.0,snr:fltarr(10)-1.0},nbd)
mstr.num = lindgen(nbd)+1
mstr.objid = obj[bd].objid
objtorem = lonarr(nobj)
prim = lonarr(nobj)
nei = lonarr(nobj)
count = 0LL
for i=0,nbd-1 do begin
  j = bd[i]
  if matches[j+1] ne matches[j] then begin
    ind = matches[matches[j]:matches[j+1]-1]
    ind = ind[sort(nmg[ind].mag)]  ;; sort by magnitude
    nind = n_elements(ind)
    ;; Only merge now if this is the brightest star in the group
    if nmg[j].mag eq min(nmg[ind].mag) then begin
      mstr[i].primind = j
      mstr[i].ngroup = nind
      mstr[i].objind[0:nind-1] = ind
      mstr[i].mag[0:nind-1] = nmg[ind].mag
      mstr[i].nmeas[0:nind-1] = ndet[ind]
      mstr[i].fluxratio[0:nind-1] = 10^((nmg[ind].mag-nmg[ind[0]].mag)/2.5)
      mstr[i].snr[0:nind-1] = 1.087/nmg[ind].err
      objtorem[ind[1:*]] = 1
      prim[ind[0]] = 1
      nei[ind[1:*]] = 1
    endif
  endif
endfor
;; Remove blank ones
torem = where(mstr.ngroup eq 0,ntorem)
if ntorem gt 0 then REMOVE,torem,mstr
nmstr = n_elements(mstr)
printlog,logfile,strtrim(nmstr,2),' objects with close neighbors within ',stringize(dcr,ndec=1),' pixels'

;; what if the neighbor itself has close neighbors?
;; do I need to cluster?
;; could we do it iteratively?

;;; Creating final obj and meas structures
;newobj = obj
;indnei = where(nei eq 1,nindnei)
;printlog,logfile,strtrim(nindnei,2)+' neighbor objects to remove'
;REMOVE,indnei,newobj


;; Merging flux in the measurements
;;  exposure by exposure
objindex = create_index(meas.objid)
MATCH,objindex.value,obj.objid,ind1,ind2,/sort,count=nmatch
si = sort(ind2)  ; sort by obj index
ind1 = ind1[si]
ind2 = ind2[si]
newmeas = meas
meastorem = lonarr(nmeas)
;; Loop over the groups
for i=0,nmstr-1 do begin
  ;print,strtrim(i+1,2),' ',mstr[i].objid,' ',strtrim(mstr[i].ngroup,2)
  ;; Make an array of all the measurments of all objects of this group
  ngrp = mstr[i].ngroup
  ngrpmeas = total(mstr[i].nmeas[0:ngrp-1],/int)
  mind = lonarr(ngrpmeas)
  cnt = 0LL
  for k=0,ngrp-1 do begin
    oind = mstr[i].objind[k]  ;; index into OBJ for this object
    mind1 = objindex.index[objindex.lo[ind1[oind]]:objindex.hi[ind1[oind]]]
    nmind1 = n_elements(mind1)
    mind[cnt:cnt+nmind1-1] = mind1   ;; add to meas index for this group
    cnt += nmind1
  endfor  
  ;; Index the measurements on exposure
  expindex = create_index(meas[mind].exposure)
  nexp = n_elements(expindex.value)
  ;; Loop over exposures
  for k=0,nexp-1 do begin
    eind = mind[expindex.index[expindex.lo[k]:expindex.hi[k]]]
    neind = n_elements(eind)
    ;print,'  exposure ',expindex.value[k],' ',meas[eind[0]].filter 
    ;print,'  ',strjoin(meas[eind].objid,', ')
    ;print,'  ',strjoin(meas[eind].mag,', ')
    if neind gt 1 then begin
      ;; Find the primary
      pind = where(meas[eind].objid eq mstr[i].objid,npind)
      if npind gt 0 then begin
        primind = eind[pind]
        meas1 = meas[primind]
      ;; If the primary is not detected in this exposure, then use the brightest source
      endif else begin
        ;print,'Primary not found'
        si = sort(meas[eind].mag)
        pind = si[0]
        primind = eind[pind]
        meas1 = meas[primind]
        meas1.objid = mstr[i].objid   ; make sure it has the right objid
      endelse
      ;; Now combine the flux,  leave all the other parameters
      ;;  m = -2.5*alog10( 10.^(-m1/2.5) + 10.^(-m2/2.5) )
      imag = meas[eind].imag
      meas1.imag = -2.5*alog10( total( 10^(-imag/2.5) ))
      mag = meas[eind].mag
      meas1.mag = -2.5*alog10( total( 10^(-mag/2.5) ))
      ;; Update in NEWMEAS
      newmeas[primind] = meas1
      ;; Measurments to remove
      meastorem[eind] = 1     ; remove all except for the primary one 
      meastorem[primind] = 0  
    ;; Only ONE measurements, make sure it has the right OBJID
    endif else begin
      ;if newmeas[eind].objid ne mstr[i].objid then print,'Single meas. Primary not found'
      newmeas[eind].objid = mstr[i].objid
    endelse
  ;if i eq 174 then stop,'174'
  endfor
endfor

;stop

;; Removing neighbor measurements
indnei = where(nei eq 1,nindnei)
newobjindex = create_index(newmeas.objid)
MATCH,newobjindex.value,obj[indnei].objid,ind1,ind2,/sort,count=nmatch
;; remove measurements belonging to the removed neighbors
meastorem2 = lonarr(nmeas)
for i=0,nmatch-1 do begin
  ind = newobjindex.index[newobjindex.lo[ind1[i]]:newobjindex.hi[ind1[i]]]
  meastorem2[ind] = 1
endfor
torem2 = where(meastorem2 eq 1,ntorem)
printlog,logfile,strtrim(ntorem,2)+' neighbor measurements to remove'
if ntorem gt 0 then REMOVE,torem,newmeas


;; SOME NEIGHBORS ARE COUNTED TWICE!!!
; 0988m505.1902
; in mstr 140 and 375


;; Average the photometry
DELVERED_SIMPLEAVGMEAS,newmeas,avgobj
;; copy the values to the final object structure
MATCH,newobj.objid,avgobj.objid,ind1,ind2,/sort,count=nmatch
temp = newobj[ind1]
STRUCT_ASSIGN,avgobj[ind2],temp,/nozero
newobj[ind1] = temp
newobj[ind1].uscatter = avgobj[ind2].urms
newobj[ind1].gscatter = avgobj[ind2].grms
newobj[ind1].rscatter = avgobj[ind2].rrms
newobj[ind1].iscatter = avgobj[ind2].irms
newobj[ind1].zscatter = avgobj[ind2].zrms
newobj[ind1].yscatter = avgobj[ind2].yrms
;; copy over other columns from the old object structure
MATCH,newobj.objid,obj.objid,ind1,ind2,/sort,count=nmatch
cols = ['PROB','EBV','MAG_AUTO','MAGERR_AUTO','ASEMI','BSEMI','THETA','ELLIPTICITY','FWHM','BRICKUNIQ']
otags = tag_names(obj)
for i=0,n_elements(cols)-1 do begin
  colind = where(otags eq cols[i],ncolind)
  newobj[ind1].(colind[0]) = obj[ind2].(colind[0])
endfor

stop

end
