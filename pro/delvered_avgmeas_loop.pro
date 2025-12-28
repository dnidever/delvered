;+
;
; DELVERED_AVGMEAS
;
; Average measurements for objects. OBJID must have been set.
;
; INPUTS:
;  expstr  Catalog of exposures and their metadata.
;  meas    Catalog of all measurements with OBJID already filled in.
;
; OUTPUTS:
;  obj     Catalog of unique objects with average photometry,
;            morphology, coordinates, etc.
;
; USAGE:
;  IDL>delvered_avgmeas,meas,obj
;
; By D.Nidever   June 2020
;
;-

pro delvered_avgmeas_loop,expstr,measfiles,measid,obj0,obj

;; Not enough inputs
nmeas = n_elements(measid)
if nmeas eq 0 then begin
  print,'Syntax - delvered_avgmeas,measfiles,measid,obj'
  return
endif
print,strtrim(nmeas,2),' measurements'

;; Create index on OBJID
measid.objid = strtrim(measid.objid,2)
oindex = create_index(measid.objid)
nobj = n_elements(oindex.value)
print,strtrim(nobj,2),' objects'

;; Index from MEAS to OBJ
measobj_index = lon64arr(nmeas)
for i=0,nobj-1 do measobj_index[oindex.index[oindex.lo[i]:oindex.hi[i]]] = i


;; Initialize the final object table, with ALL BANDS
obj_schema = {objid:'',brick:'',depthflag:0,nalfdetiter:0,neimerged:0,ra:0.0d0,dec:0.0d0,rarms:0.0,decrms:0.0,ndet:0,ndetall:0,$
              mlon:0.0d0,mlat:0.0d0,$
              umag:99.99,uerr:9.99,urms:99.99,ndetu:0,$
              gmag:99.99,gerr:9.99,grms:99.99,ndetg:0,$
              rmag:99.99,rerr:9.99,rrms:99.99,ndetr:0,$
              imag:99.99,ierr:9.99,irms:99.99,ndeti:0,$
              zmag:99.99,zerr:9.99,zrms:99.99,ndetz:0,$
              ymag:99.99,yerr:9.99,yrms:99.99,ndety:0,$
              chi:99.99,sharp:99.99,prob:99.99,ebv:99.99,mag_auto:99.99,magerr_auto:9.99,$
              asemi:999999.0,bsemi:999999.0,theta:999999.0,ellipticity:999999.0,fwhm:999999.0,$
              rmsvar:999999.0,madvar:999999.0,iqrvar:999999.0,etavar:999999.0,jvar:999999.0,$
              kvar:999999.0,chivar:999999.0,romsvar:999999.0,variable10sig:-1,nsigvar:999999.0,$
              gaia_match:0,gaia_xdist:999999.0,gaia_sourceid:0LL,gaia_ra:999999.0d0,gaia_ra_error:999999.0,$
              gaia_dec:999999.0d0,gaia_dec_error:999999.0,gaia_parallax:999999.0,gaia_parallax_error:999999.0,gaia_pmra:999999.0,$
              gaia_pmra_error:999999.0,gaia_pmdec:999999.0,gaia_pmdec_error:999999.0,gaia_gmag:999999.0,$
              gaia_gmag_error:999999.0,gaia_bpmag:999999.0,gaia_bpmag_error:999999.0,gaia_rpmag:999999.0,$
              gaia_rpmag_error:999999.0,$
              brickuniq:0B}
tot_schema = {objid:'',ra0:0.0d0,dec0:0.0d0,drawt:0.0d0,ddecwt:0.0d0,coowt:0.0d0,$
              dra:0.0d0,ddec:0.0d0,dra2:0.0d0,ddec2:0.0d0,ndet:0L,$
              mag:dblarr(6),mag2:dblarr(6),relmag:dblarr(6),relmag2:dblarr(6),relmagwt:dblarr(6),$
              magfluxwt:dblarr(6),magwt:dblarr(6),nmag:dblarr(6),$
              chi:0.0d0,sharp:0.0d0,sharpwt:0.0d0,$
              rmsvar:0.0,madvar:0.0,iqrvar:0.0,etavar:0.0,jvar:0.0,$
              kvar:0.0,chivar:0.0,romsvar:0.0,depthflag:0}
totstr = replicate(tot_schema,nobj)
totstr.objid = oindex.value
match,totstr.objid,obj0.objid,ind1,ind2,/sort
totstr[ind1].ra0 = obj0[ind2].ra
totstr[ind1].dec0 = obj0[ind2].dec

; depthflag: 1-allstar, single processing; 2-forced photometry; 3-both
;obj = replicate(obj_schema,nobj)
;obj.objid = oindex.value
;otags = tag_names(obj)

bands = ['u','g','r','i','z','Y']

;; Exposure loop
print,strtrim(n_elements(measfiles),2),' exposures'
for i=0,n_elements(measfiles)-1 do begin

  print,i+1,measfiles[i]
  meas = MRDFITS(measfiles[i],1)
  meas.objid = strtrim(meas.objid,2)

  MATCH,meas.objid,totstr.objid,ind1,ind2,/sort

  band = meas[0].filter
  bandind = where(bands eq band)

  totstr[ind2].ndet += 1

  ;; Filter mean magnitudes
  totstr[ind2].magfluxwt[bandind] += 2.5118864d^meas[ind1].mag * (1.0d0/meas[ind1].err^2)   
  totstr[ind2].magwt[bandind] += 1.0d0/meas[ind1].err^2
  totstr[ind2].nmag[bandind] += 1

  ;; mag rms
  totstr[ind2].mag[bandind] += meas[ind1].mag
  totstr[ind2].mag2[bandind] += meas[ind1].mag^2

  ;; Residual mag relative to the uncertainty
  ;;  set a lower threshold of 0.02 in the uncertainty
  totstr[ind2].relmag[bandind] += meas[ind1].mag/(meas[ind1].err > 0.02)^2
  totstr[ind2].relmag2[bandind] += meas[ind1].mag^2/(meas[ind1].err > 0.02)^2
  totstr[ind2].relmagwt[bandind] += 1.0d0/(meas[ind1].err > 0.02)^2

  ;; Best photometry


  ;; Coordinates, morphology, ebv, chi, sharp, prob
  dra = (meas[ind1].ra-totstr[ind2].ra0)*3600.0d0
  totstr[ind2].dra += dra
  ddec = (meas[ind1].dec-totstr[ind2].dec0)*3600.0d0
  totstr[ind2].ddec += ddec
  ;; Measure rms, RMS
  ;;  sqrt(mean(diff^2)) = sqrt(Sum((m1-avgm)^2)/N) =
  ;;   = sqrt((1/N)*Sum(m1^2-2*m1*avgm+avgm^2))
  ;;  just the Sum
  ;;  Sum(m1^2) - 2*avgm*Sum(m1) + N*avgm^2
  ;;  we just need to keep track of Sum(ra) and Sum(ra^2)
  totstr[ind2].dra2 += dra^2
  totstr[ind2].ddec2 += ddec^2
  ;; S/N weighted coordinates
  coowt = 1.0d0/meas[ind1].err^2
  totstr[ind2].coowt += coowt
  totstr[ind2].drawt += dra*coowt
  totstr[ind2].ddecwt += ddec*coowt

  ;; CHI
  totstr[ind2].chi += meas[ind1].chi

  ;; SHARP
  ;; use "weights" with wt=1 for normal values
  ;; and wt=0.0001 for short, allstar and S/N<5
  ;; Using weights will still use the low S/N
  ;; sharp value if it's the ONLY detection, which
  ;;  is what I want.
  if meas[0].forced eq 0 then begin
    wtsharp = dblarr(n_elements(meas))+1
    snr = 1.087/meas.err
    lowsnrind = where(abs(meas.sharp) lt 1e5 and snr lt 5,nlowsnrind)
    if nlowsnrind gt 0 then wtsharp[lowsnrind]=1e-4
    totstr[ind2].sharp += wtsharp[ind1]*meas[ind1].sharp
    totstr[ind2].sharpwt += wtsharp[ind1]
  endif else begin
    totstr[ind2].sharp += meas[ind1].sharp
    totstr[ind2].sharpwt += 1.0
  endelse

  ;; Depthflag
  ;;  OR combine to depthflag, 1-single-exposure, 2-forced, 3-both
  ;;  within a single exposure we can have some forced and some
  ;;  not, need to do this object by object
  ;; forced=0 -> depthflag=1, forced=1 -> depthflag=2, depthflag=forced+1
  totstr[ind2].depthflag OR= (meas[ind1].forced+1)

  ;;stop
endfor

;; Get mean values

;; weighted coordinates
ra = totstr.drawt/3600.0/totstr.coowt + totstr.ra0
dec = totstr.ddecwt/3600.0/totstr.coowt + totstr.dec0
;ra = totstr.dra/3600.0/totstr.ndet + totstr.ra0
;dec = totstr.ddec/3600.0/totstr.ndet + totstr.dec0
wtra = totstr.drawt/totstr.coowt
wtdec = totstr.ddecwt/totstr.coowt
rarms = sqrt( (totstr.dra2 - 2*wtra*totstr.dra + totstr.ndet*wtra^2)/totstr.ndet )
rarms *= cos(dec/!radeg)
decrms = sqrt( (totstr.ddec2 - 2*wtdec*totstr.ddec + totstr.ndet*wtdec^2)/totstr.ndet )
;rarms = sqrt( (totstr.ra2 - 2*ra*totstr.ra + totstr.ndet*ra^2)/totstr.ndet )
;decrms = sqrt( (totstr.dec2 - 2*dec*totstr.dec + totstr.ndet*dec^2)/totstr.ndet )

;; mean photometry
newmag = dblarr(nobj,6)
newerr = dblarr(nobj,6)
for i=0,5 do begin
  newflux1 = totstr.magfluxwt[i]/totstr.magwt[i]
  newmag1 = 2.50*alog10(newflux1)
  newerr1 = sqrt(1.0/totstr.magwt[i])
  bdmag1 = where(finite(newmag1) eq 0,nbdmag1)
  if nbdmag1 gt 0 then begin
    newmag1[bdmag1] = 99.99
    newerr1[bdmag1] = 9.99
  endif
  newmag[*,i] = newmag1
  newerr[*,i] = newerr1
endfor

;; photometry rms
magrms = dblarr(nobj,6)
;relmagrms = dblarr(nobj,6)
for i=0,5 do begin
  magrms1 = sqrt( (totstr.mag2[i] - 2*newmag[*,i]*totstr.mag[i] + totstr.nmag[i]*newmag[*,i]^2)/totstr.nmag[i] )
  bdrms = where(totstr.nmag[i] eq 0,nbdrms)
  if nbdrms gt 0 then magrms1[bdrms] = 99.99
  magrms[*,i] = magrms1
  ;; relmag
  ;relmagrms1 = sqrt( (totstr.relmag2[i] - 2*newmag[*,i]*totstr.relmag[i] + totstr.nmag[i]*newmag[*,i]^2)/totstr.nmag[i] )
  ;bdrelrms = where(totstr.nmag[i] eq 0,nbdrelrms)
  ;if nbdrelrms gt 0 then relmagrms1[bdrelrms] = 99.99
  ;relmagrms[*,i] = relmagrms1
endfor

;; Make average CHI
avgchi = totstr.chi/totstr.ndet
;gdchi = where(numchi gt 0,ngdchi)
;avgchi = fltarr(nobj)+99.99
;if ngdchi gt 0 then avgchi[gdchi]=totchi[gdchi]/numchi[gdchi]
;obj.chi = avgchi
;; Make average SHARP
;gdsharp = where(totwtsharp gt 0,ngdsharp)
;avgsharp = fltarr(nobj)+99.99
;if ngdsharp gt 0 then avgsharp[gdsharp]=totsharp[gdsharp]/totwtsharp[gdsharp]
;obj.sharp = avgsharp
avgsharp = totstr.sharp/totstr.sharpwt


;; NDET and NDETALL
totstr.ndet = total(totstr.nmag,1)
 
;; depthflag: 1-allstar, single processing; 2-forced photometry; 3-both
obj = replicate(obj_schema,nobj)
struct_assign,totstr,obj,/nozero
otags = tag_names(obj)
obj.ra = ra
obj.dec = dec
obj.rarms = rarms
obj.decrms = decrms
obj.sharp = avgsharp
obj.chi = avgchi
for i=0,5 do begin
  magind = where(otags eq strupcase(bands[i])+'MAG')
  errind = where(otags eq strupcase(bands[i])+'ERR')
  rmsind = where(otags eq strupcase(bands[i])+'RMS')
  numind = where(otags eq 'NDET'+strupcase(bands[i]))
  obj.(magind) = newmag[*,i]
  obj.(errind) = newerr[*,i]
  obj.(rmsind) = magrms[*,i]
  obj.(numind) = totstr.nmag[i]
endfor

;; Variables
newmag0 = newmag
bd = where(newmag0 gt 50)
newmag0[bd] = 0.0
sumresid = dblarr(nobj)
nresid = lonarr(nobj)
sumresidsq = dblarr(nobj)
sumrelresidsq = dblarr(nobj)
for i=0,5 do begin
  nresid += totstr.nmag[i]

  ;; Sum(m1-avgm)
  ;; = Sum(m1)-N*avgm
  sumresid1 = totstr.mag[i]-totstr.nmag[i]*newmag0[*,i]
  sumresid += sumresid1

  ;; Sum( (m1-avgm)^2 )
  ;; = Sum( m1^2 - 2*avgm*m1 + avgm^2 )
  ;; = Sum(m1^2) - 2*avgm*Sum(m1) + N*avgm^2
  sumresidsq1 = totstr.mag2[i] - 2*newmag0[*,i]*totstr.mag[i] + totstr.nmag[i]*newmag0[*,i]^2
  sumresidsq += sumresidsq1

  ;; Sum( (m1-avgm)^2/err1^2 )
  ;; = Sum( (m1/err1)^2 - 2*avgm*m1/err1^2 + (avgm/err1)^2 )
  ;; = Sum((m1/err1)^2) - 2*avgm*Sum(m1/err1^2) + avgm^2 * Sum(1/err1^2)
  sumrelresidsq1 = totstr.relmag2[i] - 2*newmag0[*,i]*totstr.relmag[i] + newmag0[*,i]^2 * totstr.relmagwt[i]
  sumrelresidsq += sumrelresidsq1
  ;; (totstr.mag2[i] - 2*newmag[*,i]*totstr.mag[i] + totstr.nmag[i]*newmag[*,i]^2)
endfor

;; RMS
rmsvar = sqrt((sumresidsq > 0.0)/nresid)
obj.rmsvar = rmsvar

;; Relative variability metrics
chivar = sqrt((sumrelresidsq > 0.0))/nresid
;relresid2 = relresid[gdrelresid]
;pk = relresid2^2-1
;jvar = total( sign(pk)*sqrt(abs(pk)) )/ngdrelresid
;chivar = sqrt(total(relresid2^2))/ngdrelresid
;kdenom = sqrt(total(relresid2^2)/ngdrelresid)
;if kdenom ne 0 then begin
;  kvar = (total(abs(relresid2))/ngdrelresid) / kdenom
;endif else begin
;  kvar = !values.f_nan
;endelse
;; RoMS                                                                                                                                                                                  
;romsvar = total(abs(relresid2))/(ngdrelresid-1)
;obj.jvar = jvar
;obj.kvar = kvar
obj.chivar = chivar
;obj.romsvar = romsvar

stop

;; Fiducial magnitude
tempmag = newmag
tempmag[*,0] = newmag[*,2]  ;; rmag
tempmag[*,1] = newmag[*,1]  ;; gmag
tempmag[*,2] = newmag[*,3]  ;; imag
tempmag[*,3] = newmag[*,4]  ;; zmag
tempmag[*,4] = newmag[*,5]  ;; ymag
tempmag[*,5] = newmag[*,0]  ;; umag
fidmag = fltarr(nobj)
for i=0,nobj-1 do begin
  ind = where(tempmag[i,*] lt 50)
  fidmag[i] = tempmag[i,ind[0]]
endfor

;; Select variables using photometric variability indices.

;; Select Variables
;;  1) Construct fiducial magnitude (done in loop above)
;;  2) Construct median VAR and sigma VAR versus magnitude
;;  3) Find objects that Nsigma above the median VAR line
si = sort(fidmag)   ; NaNs are at end
;;varcol = 'madvar'
varcol = 'rmsvar'
gdvar = where(finite(obj.rmsvar) eq 1 and finite(fidmag) eq 1,ngdvar,comp=bdvar,ncomp=nbdvar)
obj.variable10sig = 0
if ngdvar gt 0 then begin
  binsize = 0.25
  BINDATA,fidmag[gdvar],fidmag[gdvar],fidmagmed_bins,fidmagmed,binsize=binsize,/med
  BINDATA,fidmag[gdvar],fidmag[gdvar],xbin2,numhist,binsize=binsize,/hist
  ;; Fix NaNs in fidmagmed
  bdfidmagmed = where(finite(fidmagmed) eq 0,nbdfidmagmed)
  if nbdfidmagmed gt 0 then fidmagmed[bdfidmagmed] = fidmagmed_bins[bdfidmagmed]
  ;; Median metric
  BINDATA,fidmag[gdvar],obj[gdvar].rmsvar,xbin,varmed,binsize=binsize,/med
  ;; Smooth, it handles NaNs well
  smlen = 5
  smvarmed = gsmooth(varmed,smlen)
  bdsmvarmed = where(finite(smvarmed) eq 0,nbdsmvarmed)
  if nbdsmvarmed gt 0 then smvarmed[bdsmvarmed] = median(smvarmed)
  ;; Interpolate to all the objects
  gv = where(finite(smvarmed) eq 1,ngv,comp=bv,ncomp=nbv)
  if ngv le 1 then begin
    print,'Not enough good RMSVAR values to detect variables'
    obj.variable10sig = 0
    obj.nsigvar = !values.f_nan
    return
  endif
  INTERP,fidmagmed[gv],smvarmed[gv],fidmag[gdvar],objvarmedgd
  objvarmed = dblarr(nobj)
  objvarmed[gdvar] = objvarmedgd
  objvarmed[gdvar] = min(smvarmed[gv]) > objvarmed[gdvar]   ; lower limit
  if nbdvar gt 0 then objvarmed[bdvar]=smvarmed[gv[ngv-1]]  ; objects with bad fidmag, set to last value
  ;; Scatter in metric around median
  ;;  calculate MAD ourselves so that it's around our computed
  ;;  median metric line
  BINDATA,fidmag[gdvar],abs(obj[gdvar].rmsvar-objvarmed[gdvar]),xbin3,varsig,binsize=binsize,/med
  ;;varsig *= 1.4826  ; scale MAD to stddev
  ;; Fix values for bins with few points
  bdhist = where(numhist lt 3,nbdhist,comp=gdhist,ncomp=ngdhist)
  if nbdhist gt 0 then begin
    if ngdhist gt 0 then begin
      varsig[bdhist] = median(varsig[gdhist])
   endif else begin
      varsig[*] = 0.02
   endelse
  endif
  ;; Smooth
  smvarsig = gsmooth(varsig,smlen)
  ;; Interpolate to all the objects
  gv = where(finite(smvarsig) eq 1,ngv,comp=bv,ncomp=nbv)
  if ngv le 1 then begin
    print,'Not enough good RMSVAR values to detect variables'
    obj.variable10sig = 0
    obj.nsigvar = !values.f_nan
    return
  endif
  INTERP,fidmagmed[gv],smvarsig[gv],fidmag[gdvar],objvarsiggd
  objvarsig = dblarr(nobj)
  objvarsig[gdvar] = objvarsiggd
  objvarsig[gdvar] = min(smvarsig[gv]) > objvarsig[gdvar]   ; lower limit
  if nbdvar gt 0 then objvarsig[bdvar]=smvarsig[gv[ngv-1]]  ; objects with bad fidmag, set to last value
  ;; Detect positive outliers
  nsigvarthresh = 10.0
  nsigvar = (obj.rmsvar-objvarmed)/objvarsig
  obj[gdvar].nsigvar = nsigvar[gdvar]
  isvar = where(nsigvar[gdvar] gt nsigvarthresh,nisvar)
  print,strtrim(nisvar,2)+' variables detected'
  if nisvar gt 0 then obj[gdvar[isvar]].variable10sig = 1
endif


;; EBV
glactc,obj.ra,obj.dec,2000.0,glon,glat,1,/deg
obj.ebv = dust_getval(glon,glat,/noloop,/interp)

return



;bd = where(tempmag gt 50)
;tempmag[bd] = !values.f_nan
;magarr = [newobj.rmag,newobj.gmag,newobj.imag,newobj.zmag,newobj.ymag,newobj.umag]
;gfid = where(magarr lt 50,ngfid)
;if ngfid gt 0 then fidmag[i]=magarr[gfid[0]]


stop

;filtindex = create_index(meas1.filter)
;nfilters = n_elements(filtindex.value)
;resid = dblarr(nmeas)+!values.f_nan     ; residual mag
;relresid = dblarr(nmeas)+!values.f_nan  ; residual mag relative to the uncertainty
;for f=0,nfilters-1 do begin
;  filt = strupcase(filtindex.value[f])
;  findx = filtindex.index[filtindex.lo[f]:filtindex.hi[f]]
;  magind = where(otags eq filt+'MAG',nmagind)
;  errind = where(otags eq filt+'ERR',nerrind)
;  gph = where(meas1[findx].mag lt 50,ngph)
;  if ngph gt 1 then begin
;    ;; Residual mag                                                                                                                                                                        
;    resid[findx[gph]] = meas1[findx[gph]].mag-obj[i].(magind)
;    ;; Residual mag relative to the uncertainty                                                                                                                                            
;    ;;  set a lower threshold of 0.02 in the uncertainty                                                                                                                                   
;    relresid[findx[gph]] = sqrt(ngph/(ngph-1.0)) * (meas[findx[gph]].mag-obj[i].(magind))/(meas[findx[gph]].err > 0.02)
;  endif
;endfor


;; Variability metrics
;; Calculate variability indices
gdresid = where(finite(resid) eq 1,ngdresid)
if ngdresid gt 0 then begin
  resid2 = resid[gdresid]
  sumresidsq = total(resid2^2)
  tsi = sort(meas1[gdresid].mjd)
  resid2tsi = resid2[tsi]
  quartiles = cgpercentiles(resid2,percentiles=[0.25,0.50,0.75])
  ;; RMS
  rms = sqrt(sumresidsq/ngdresid)
  ;; MAD
  madvar = 1.4826*median(abs(resid2-quartiles[1]))
  ;; IQR
  iqrvar = 0.741289*(quartiles[2]-quartiles[0])
  ;; 1/eta
  etavar = sumresidsq / total((resid2tsi[1:*]-resid2tsi[0:ngdresid-2])^2)
  obj[i].rmsvar = rms
  obj[i].madvar = madvar
  obj[i].iqrvar = iqrvar
  obj[i].etavar = etavar
endif

;; Calculate variability indices wrt to uncertainties                                                                                                                                      
gdrelresid = where(finite(relresid) eq 1,ngdrelresid)
if ngdrelresid gt 0 then begin
  relresid2 = relresid[gdrelresid]
  pk = relresid2^2-1
  jvar = total( sign(pk)*sqrt(abs(pk)) )/ngdrelresid
  chivar = sqrt(total(relresid2^2))/ngdrelresid
  kdenom = sqrt(total(relresid2^2)/ngdrelresid)
  if kdenom ne 0 then begin
    kvar = (total(abs(relresid2))/ngdrelresid) / kdenom
  endif else begin
    kvar = !values.f_nan
  endelse
  ;; RoMS                                                                                                                                                                                  
  romsvar = total(abs(relresid2))/(ngdrelresid-1)
  obj[i].jvar = jvar
  obj[i].kvar = kvar
  obj[i].chivar = chivar
  obj[i].romsvar = romsvar
endif

;; Fiducial magnitude, used to select variables below                                                                                                                                      
;;  order of priority: r,g,i,z,Y,VR,u                                                                                                                                                      
if obj[i].ndet gt 0 then begin
  magarr = [obj[i].rmag,obj[i].gmag,obj[i].imag,obj[i].zmag,obj[i].ymag,obj[i].umag]
  gfid = where(magarr lt 50,ngfid)
  if ngfid gt 0 then fidmag[i]=magarr[gfid[0]]
endif


stop

;; Exposure names
;;   F23-00912982_15, trim off ccdnum
dum = reform((strsplitter(meas.exposure,'_',/extract))[0,*])
expnum = reform((strsplitter(dum,'-',/extract))[1,*])
eindex = create_index(expnum)
;; match EXPSTR to EINDEX
;;  also some exposures have 0 measurements
MATCH,eindex.value,expstr.expnum,ind1,ind2,/sort
mexpstr = expstr[ind2]   ;; matched expstr, ind1 is in the right order
nexposure = n_elements(mexpstr)
print,strtrim(nexposure,2),' exposures with measurements'
ufilter = mexpstr[uniq(mexpstr.filter,sort(mexpstr.filter))].filter
nufilter = n_elements(ufilter)
print,strtrim(nufilter,2),' unique filters'

;; Loop over unique filters
For f=0,nufilter-1 do begin
  filtind = where(mexpstr.filter eq ufilter[f],nfiltind)
  
  ;; Indices for the magnitude and errors in OBJ
  magind = where(otags eq strupcase(ufilter[f])+'MAGALL')
  errind = where(otags eq strupcase(ufilter[f])+'ERRALL')
  rmsind = where(otags eq strupcase(ufilter[f])+'RMSALL')
  nobsind = where(otags eq 'NDETALL'+strupcase(ufilter[f]))

  ;; Only one exposure for this filter, copy
  if nfiltind eq 1 then begin

    ;; Now copy in the values, MEAS only has "good" detections
    mind = eindex.index[eindex.lo[filtind]:eindex.hi[filtind]]
    obj[measobj_index[mind]].(magind) = meas[mind].mag
    obj[measobj_index[mind]].(errind) = meas[mind].err
    obj[measobj_index[mind]].(nobsind) = 1
    ;; Depthflag
    ;;  OR combine to depthflag, 1-single-exposure, 2-forced, 3-both
    ;;  within a single exposure we can have some forced and some
    ;;  not, need to do this object by object
    ;; forced=0 -> depthflag=1, forced=1 -> depthflag=2, depthflag=forced+1
    obj[measobj_index[mind]].depthflag OR= (meas[mind].forced+1)

  ;; Multiple exposures for this filter to average
  endif else begin

    ;; Loop through all of the exposures and add up the flux, totalwt, etc.
    totalwt = dblarr(nobj)
    totalfluxwt = dblarr(nobj)
    for k=0,nfiltind-1 do begin
      mind = eindex.index[eindex.lo[filtind[k]]:eindex.hi[filtind[k]]]
      totalwt[measobj_index[mind]] += 1.0d0/meas[mind].err^2
      totalfluxwt[measobj_index[mind]] += 2.5118864d^meas[mind].mag * (1.0d0/meas[mind].err^2)
      ;; Depthflag
      obj[measobj_index[mind]].depthflag OR= (meas[mind].forced+1)
    endfor
    newflux = totalfluxwt/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)
    bdmag = where(finite(newmag) eq 0,nbdmag)
    if nbdmag gt 0 then begin
      newmag[bdmag] = 99.99
      newerr[bdmag] = 9.99
    endif

    ;; Measure rms, RMS
    ;;  sqrt(mean(diff^2))
    totaldiff = dblarr(nobj)
    numobs = lonarr(nobj)
    for k=0,nfiltind-1 do begin
      mind = eindex.index[eindex.lo[filtind[k]]:eindex.hi[filtind[k]]]
      totaldiff[measobj_index[mind]] += (newmag[measobj_index[mind]] - meas[mind].mag)^2
      numobs[measobj_index[mind]]++
    endfor
    newrms = sqrt( totaldiff/(numobs>1) )
    if nbdmag gt 0 then newrms[bdmag]=99.99

    ;; Set rms=99.99 for numobs=1
    oneobs = where(numobs eq 1,noneobs)
    if noneobs gt 0 then newrms[oneobs]=99.99

    obj.(magind) = newmag
    obj.(errind) = newerr
    obj.(rmsind) = newrms
    obj.(nobsind) = numobs
  endelse  ; combine multiple exposures for this filter

Endfor  ;; unique filter loop


;; Compute "best" photometry
;;--------------------------------

;; Loop over unique filters
For f=0,nufilter-1 do begin
  filtind = where(mexpstr.filter eq ufilter[f],nfiltind)
  
  ;; Indices for the magnitude and errors in OBJ
  magind = where(otags eq strupcase(ufilter[f])+'MAG')
  errind = where(otags eq strupcase(ufilter[f])+'ERR')
  rmsind = where(otags eq strupcase(ufilter[f])+'RMS')
  nobsind = where(otags eq 'NDET'+strupcase(ufilter[f]))

  ;; Only one exposure for this filter, copy
  if nfiltind eq 1 then begin

    ;; Now copy in the values, MEAS only has "good" detections
    mind = eindex.index[eindex.lo[filtind]:eindex.hi[filtind]]
    obj[measobj_index[mind]].(magind) = meas[mind].mag
    obj[measobj_index[mind]].(errind) = meas[mind].err
    obj[measobj_index[mind]].(nobsind) = 1

  ;; Multiple exposures for this filter to average
  endif else begin

    ;; Loop through all of the exposures and add up the flux, totalwt, etc.
    totalwt = dblarr(nobj)
    totalfluxwt = dblarr(nobj)
    numobs = lonarr(nobj)
    numgdobs = lonarr(nobj)
    for k=0,nfiltind-1 do begin
      mind = eindex.index[eindex.lo[filtind[k]]:eindex.hi[filtind[k]]]
      wt = 1.0d0/meas[mind].err^2
      ;; Severely downweight shallow non-forced photometry with err>0.2 mag
      if mexpstr[filtind[k]].exptime lt 90 then begin
        bd = where(meas[mind].forced eq 0 and meas[mind].err gt 0.2,nbd,comp=gd,ncomp=ngd)
        if nbd gt 0 then wt[bd] *= 1e-4
      endif else gd=lindgen(n_elements(mind))
      numgdobs[measobj_index[mind[gd]]]++
      numobs[measobj_index[mind]]++
      totalwt[measobj_index[mind]] += wt
      totalfluxwt[measobj_index[mind]] += 2.5118864d^meas[mind].mag * wt
    endfor
    newflux = totalfluxwt/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)
    bdmag = where(finite(newmag) eq 0,nbdmag)
    if nbdmag gt 0 then begin
      newmag[bdmag] = 99.99
      newerr[bdmag] = 9.99
    endif
    obj.(nobsind) = numgdobs
    ;; Some objects only have poor shallow exposures, drop them
    ;;   they are not very reliable, and we have the ALL phot for them
    badobs = where(numgdobs eq 0 and newmag lt 50,nbadobs)
    if nbadobs gt 0 then begin
      newmag[badobs] = 99.99
      newerr[badobs] = 9.99
    endif

    ;if ufilter[f] eq 'g' then stop,'gmag 1'

    ;; Measure rms, RMS
    ;;  sqrt(mean(diff^2))
    totaldiff = dblarr(nobj)
    numobs = lonarr(nobj)
    for k=0,nfiltind-1 do begin
      mind = eindex.index[eindex.lo[filtind[k]]:eindex.hi[filtind[k]]]
      if mexpstr[filtind[k]].exptime lt 90 then begin
        bd = where(meas[mind].forced eq 0 and meas[mind].err gt 0.2,nbd,comp=gd,ncomp=ngd)
      endif else begin
        gd = lindgen(n_elements(mind))
        ngd = n_elements(gd)
      endelse
      if ngd gt 0 then begin
        totaldiff[measobj_index[mind[gd]]] += (newmag[measobj_index[mind[gd]]] - meas[mind[gd]].mag)^2
        numobs[measobj_index[mind[gd]]]++
      endif
    endfor
    newrms = sqrt( totaldiff/(numobs>1) )
    if nbdmag gt 0 then newrms[bdmag]=99.99

    ;; Set rms=99.99 for numobs=1
    oneobs = where(numobs eq 1,noneobs)
    if noneobs gt 0 then newrms[oneobs]=99.99

    obj.(magind) = newmag
    obj.(errind) = newerr
    obj.(rmsind) = newrms

    ;if ufilter[f] eq 'g' then stop,'gmag 2'

  endelse  ; combine multiple exposures for this filter

Endfor  ;; unique filter loop


;; NDET and NDETALL
obj.ndet = obj.ndetu+obj.ndetg+obj.ndetr+obj.ndeti+obj.ndetz+obj.ndety
obj.ndetall = obj.ndetallu+obj.ndetallg+obj.ndetallr+obj.ndetalli+obj.ndetallz+obj.ndetally



;; Coordinates, morphology, ebv, etc
;;-----------------------------------

;; Get average chi, sharp, prob
totchi = fltarr(nobj) & numchi = lon64arr(nobj)
totsharp = fltarr(nobj) & totwtsharp = fltarr(nobj)
;totprob = fltarr(nobj) & numprob = lon64arr(nobj)

;; Exposure loop
For i=0,nexposure-1 do begin
  ;; All measurements for this exposure
  mind = eindex.index[eindex.lo[i]:eindex.hi[i]]
  nmind = n_elements(mind)
  ;; CHI
  chi1 = fltarr(nobj)+!values.f_nan
  chi1[measobj_index[mind]] = meas[mind].chi
  gdchi = where(finite(chi1) eq 1 and chi1 lt 1e5,ngdchi)
  if ngdchi gt 0 then begin
    totchi[gdchi] += chi1[gdchi]
    numchi[gdchi]++
  endif
  ;; SHARP
  ;; use "weights" with wt=1 for normal values
  ;; and wt=0.0001 for short, allstar and S/N<5
  ;; Using weights will still use the low S/N
  ;; sharp value if it's the ONLY detection, which
  ;;  is what I want.
  sharp1 = fltarr(nobj)+!values.f_nan
  sharp1[measobj_index[mind]] = meas[mind].sharp
  gdsharp = where(finite(sharp1) eq 1 and sharp1 lt 1e5,ngdsharp)
  if ngdsharp gt 0 then begin
    wtsharp1 = fltarr(nobj)
    wtsharp1[gdsharp] = 1.0
    ;; Set wt to 1e-4 for ALLSTAR short exposures with S/N<5
    snr = fltarr(nobj)
    snr[measobj_index[mind]] = 1.087/meas[mind].err
    forced = intarr(nobj)-1
    forced[measobj_index[mind]] = meas[mind].forced
    lowsnrind = where(finite(sharp1) eq 1 and abs(sharp1) lt 1e5 and snr lt 5 and forced ne 1,nlowsnrind)
    if nlowsnrind gt 0 then wtsharp1[lowsnrind]=1e-4
    ;; weighted sum
    totsharp[gdsharp] += wtsharp1[gdsharp]*sharp1[gdsharp]
    totwtsharp[gdsharp] += wtsharp1[gdsharp]
  endif
Endfor
;; Make average CHI
gdchi = where(numchi gt 0,ngdchi)
avgchi = fltarr(nobj)+99.99
if ngdchi gt 0 then avgchi[gdchi]=totchi[gdchi]/numchi[gdchi]
obj.chi = avgchi
;; Make average SHARP
gdsharp = where(totwtsharp gt 0,ngdsharp)
avgsharp = fltarr(nobj)+99.99
if ngdsharp gt 0 then avgsharp[gdsharp]=totsharp[gdsharp]/totwtsharp[gdsharp]
obj.sharp = avgsharp


;; Weighted RA/DEC
totalwt = dblarr(nobj)
totalrawt = dblarr(nobj)
totaldecwt = dblarr(nobj)
for i=0,nexposure-1 do begin
  mind = eindex.index[eindex.lo[i]:eindex.hi[i]]  
  totalwt[measobj_index[mind]] += 1.0d0/meas[mind].err^2
  totalrawt[measobj_index[mind]] += meas[mind].ra * (1.0d0/meas[mind].err^2)
  totaldecwt[measobj_index[mind]] += meas[mind].dec * (1.0d0/meas[mind].err^2)
endfor
newra = totalrawt/totalwt
newdec = totaldecwt/totalwt
obj.ra = newra
obj.dec = newdec


; RA/DEC RMS
;  sqrt(mean(diff^2))
totalradiff = dblarr(nobj)
totaldecdiff = dblarr(nobj)
for i=0,nexposure-1 do begin
  mind = eindex.index[eindex.lo[i]:eindex.hi[i]]  
  totalradiff[measobj_index[mind]] += (obj[measobj_index[mind]].ra - meas[mind].ra)^2
  totaldecdiff[measobj_index[mind]] += (obj[measobj_index[mind]].dec - meas[mind].dec)^2
endfor
newrarms = sqrt( totalradiff/(obj.ndetall>1) ) * 3600 * cos(obj.dec/!radeg)
newdecrms = sqrt( totaldecdiff/(obj.ndetall>1) ) * 3600
bd = where(obj.ndetall eq 0,nbd)
if nbd gt 0 then newrarms[bd]=99.99
if nbd gt 0 then newdecrms[bd]=99.99
; Set rms=99.99 for numobs=1
oneobs = where(obj.ndetall eq 1,noneobs)
if noneobs gt 0 then newrarms[oneobs]=99.99
if noneobs gt 0 then newdecrms[oneobs]=99.99
obj.rarms = newrarms
obj.decrms = newdecrms

;; EBV
glactc,obj.ra,obj.dec,2000.0,glon,glat,1,/deg
obj.ebv = dust_getval(glon,glat,/noloop,/interp)

end
