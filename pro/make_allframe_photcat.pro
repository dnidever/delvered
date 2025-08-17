pro make_allframe_photcat,base

;; this uses the allframe output and calibrates the photometry
;; steps from allframe.pro and delvered_forcebrick.pro

logfile = -1

;; /d0/dnidever/delve/brk.6w45Tr

brick = '0280m685'

delvereddir = '/home/dnidever/projects/delvered/'
;; Load the brick information
brickstr = MRDFITS(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz',1,/silent)
;; Get the brick information
bind = where(brickstr.brickname eq brick,nbind)
brickstr1 = brickstr[bind[0]]
tilestr = MAKE_BRICK_WCS(brickstr1)

;; Step 1: Get the list of exposures/chips that overlap this brick
;;----------------------------------------------------------------
;printlog,logfile,'Step 1: Get list of chip files that overlap this brick'
;cenra = brickstr1.ra
;cendec = brickstr1.dec
;tempdir = '/d0/dnidever/delve/brk.6w45Tr'
;tmpfile = MKTEMP('tmp',/nodot,outdir=tempdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
;tmpfile += '.fits'
;spawn,[delvereddir+'bin/query_delvered_summary_table',strtrim(cenra,2),strtrim(cendec,2),tmpfile,'--lim','0.5'],out,errout,/noshell
;info = file_info(tmpfile)
tmpfile = '/d0/dnidever/delve/brk.6w45Tr/chstr.fits'
chstr = MRDFITS(tmpfile,1,/silent)
chstr.base = strtrim(chstr.base,2)
;file_delete,tmpfile,/allow
nchstr = n_elements(chstr)

;; Match to the mch files
mchfile = base+'.mch'
loadmch,mchfile,files,trans
match,file_basename(files,'.als'),chstr.base,ind1,ind2,/sort,count=nmatch
;; put them in the order in the mchfile
si = sort(ind1)
chstr0 = chstr
chstr = chstr[ind2[si]]
nchstr = n_elements(chstr)

MAKEMAG,base+'.tfr',base+'.makemag',error=magerror

; Prepend the ALF header to the makemag file
line1='' & line2='' & line3=''
openr,unit,/get_lun,file_basename(files[0],'.als')+'.alf'
readf,unit,line1
readf,unit,line2
readf,unit,line3
close,unit
free_lun,unit
head = [line1,line2,line3]
WRITELINE,base+'.makemag',head,/prepend
loadmakemag,base+'.makemag',mag

;; Add fake sextractor columns
newname = ['FLAG','PROB','MAG_AUTO','MAGERR_AUTO','BACKGROUND','THRESHOLD','ISOAREA',$
           'ASEMI','BSEMI','THETA','ELLIPTICITY','FWHM']
newvalues = ['0L','0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0','0.0']
temp = mag
undefine,mag
add_tags,temp,newname,newvalues,mag
undefine,temp

;; Write .mag file
mwrfits,mag,base+'.mag',/create


;; Step 4: Calculate coordinates
;;-------------------------------
magfile = base+'.mag'
printlog,logfile,'Step 4: Adding coordinates'
; Load the MCH file
LOADMCH,mchfile,alsfiles
nalsfiles = n_elements(alsfiles)
; Load the photometry file
instphot = PHOTRED_READFILE(magfile)
ninstphot = n_elements(instphot)
printlog,logfile,'Nstars = '+strtrim(ninstphot,2)
;; Converting to IDL X/Y convention, starting at (0,0)
;; DAOPHOT has X/Y start at (1,1)
HEAD_XYAD,tilestr.head,instphot.x-1.0,instphot.y-1.0,ra,dec,/degree

;; Step 5: Calibrating photometry with zero-points
;;-------------------------------------------------
printlog,logfile,'Step 5: Calibrating photometry with zero-points'
cmag = fltarr(ninstphot,nchstr)+99.99
cerr = fltarr(ninstphot,nchstr)+9.99
;; Chip loop
for i=0,nchstr-1 do begin
  ; id, x, y, ra, dec, unsolved magnitudes, chi, sharp
  imag = instphot.(2*i+3)
  ierr = instphot.(2*i+4)
  gdmag = where(imag lt 50,ngdmag)
  if ngdmag gt 0 then begin
    ;; exptime, aperture correction, zero-point
    ;; aperture correction is SUBTRACTIVE, makes it brighter
    ;; ZPTERM is a SUBTRACTIVE constant offset
    cmag[gdmag,i] = imag[gdmag] + 2.5*alog10(chstr[i].exptime) - chstr[i].apcor - chstr[i].calib_zpterm
    ;; Add zero-point error in quadrature
    cerr[gdmag,i] = sqrt(ierr[gdmag]^2+chstr[i].calib_zptermsig^2)
  endif
endfor
;; Calculate average photometry per filter
ufilt = chstr[uniq(chstr.filter,sort(chstr.filter))].filter
nufilt = n_elements(ufilt)
avgmag = fltarr(ninstphot,nufilt)
avgerr = fltarr(ninstphot,nufilt)
ndet = lonarr(ninstphot,nufilt)
for i=0,nufilt-1 do begin
  gdf = where(chstr.filter eq ufilt[i],ngdf)
  ;; Single exposure
  if ngdf eq 1 then begin
    avgmag[*,i] = cmag[*,gdf[0]]
    avgerr[*,i] = cerr[*,gdf[0]]
    ndet[*,i] = (cmag[*,gdf[0]] lt 50)
  ;; Multiple exposures
  endif else begin
    ;; Loop through all of the exposures and add up the flux,
    ;;  totalwt, etc.
    totalwt = dblarr(ninstphot)
    totalfluxwt = dblarr(ninstphot)
    for k=0,ngdf-1 do begin
      gdmag = where(cmag[*,gdf[k]] lt 50,ngdmag)
      if ngdmag gt 0 then begin
        totalwt[gdmag] += 1.0d0/cerr[gdmag,gdf[k]]^2
        totalfluxwt[gdmag] += 2.5118864d^cmag[gdmag,gdf[k]] * (1.0d0/cerr[gdmag,gdf[k]]^2)
	ndet[gdmag,i]++
      endif
    endfor
    newflux = totalfluxwt/totalwt
    newmag = 2.50*alog10(newflux)
    newerr = sqrt(1.0/totalwt)
    bdmag = where(finite(newmag) eq 0,nbdmag)
    if nbdmag gt 0 then begin
      newmag[bdmag] = 99.99
      newerr[bdmag] = 9.99
    endif
    avgmag[*,i] = newmag
    avgerr[*,i] = newerr
  endelse
endfor
;; Measure scatter
scatter = fltarr(ninstphot,nufilt)+99.99
for i=0,nufilt-1 do begin
  gdf = where(chstr.filter eq ufilt[i],ngdf)
  if ngdf gt 1 then begin
    totaldiff = fltarr(ninstphot)
    for k=0,ngdf-1 do begin
      gdmag = where(cmag[*,gdf[k]] lt 50,ngdmag)
      if ngdmag gt 0 then totaldiff[gdmag] += (avgmag[gdmag,i]-cmag[gdmag,gdf[k]])^2
    endfor
    scatter[*,i] = sqrt( totaldiff/(ndet[*,i]>1) )
    bd = where(ndet[*,i] le 1,nbd)
    if nbd gt 0 then scatter[bd,i]=99.99
  endif
endfor
;; Create final catalog schema
newschema = {objid:'',x:0.0,y:0.0,ra:0.0d0,dec:0.0d0}
;; Add columns for calibrated single-epoch photometry columns
cmagnames = strarr(nchstr)
cerrnames = strarr(nchstr)
for i=0,nufilt-1 do begin
  ind = where(chstr.filter eq ufilt[i],nind)
  cmagnames[ind] = strupcase(ufilt[i])+strtrim(lindgen(nind)+1,2)+'MAG'
  cerrnames[ind] = strupcase(ufilt[i])+strtrim(lindgen(nind)+1,2)+'ERR'
endfor
for i=0,nchstr-1 do newschema = create_struct(newschema,cmagnames[i],0.0,cerrnames[i],0.0)
;; Add columns for average photometry per filter
for i=0,nufilt-1 do newschema = create_struct(newschema,ufilt[i]+'MAG',0.0,ufilt[i]+'ERR',0.0,ufilt[i]+'SCATTER',0.0,'NDET'+ufilt[i],0L)
;; Extra columns
newschema = create_struct(newschema,'chi',0.0,'sharp',0.0,'prob',0.0,'ebv',0.0)
;; other SE columns
newschema = create_struct(newschema,'mag_auto',0.0,'magerr_auto',0.0,'asemi',0.0,'bsemi',0.0,'theta',0.0,'ellipticity',0.0,'fwhm',0.0)
;; in unique brick area
newschema = create_struct(newschema,'brickuniq',0B)
;; Create final catalog
phot = replicate(newschema,ninstphot)
struct_assign,instphot,phot,/nozero
phtags = tag_names(phot)
;; object IDs
phot.objid = brick+'.'+strtrim(instphot.id,2)
;; Stuff in the coordinates calculated above
phot.ra = ra
phot.dec = dec
;; Stuff in the calibrated single-epoch photometry columns
for i=0,nchstr-1 do begin
  magind = where(strupcase(phtags) eq cmagnames[i],nmagind)
  errind = where(strupcase(phtags) eq cerrnames[i],nerrind)
  phot.(magind) = cmag[*,i]
  phot.(errind) = cerr[*,i]
endfor
;; Stuff in the average photometry per filter
for i=0,nufilt-1 do begin
  magind = where(strupcase(phtags) eq strupcase(ufilt[i])+'MAG',nmagind)
  errind = where(strupcase(phtags) eq strupcase(ufilt[i])+'ERR',nerrind)
  scatind = where(strupcase(phtags) eq strupcase(ufilt[i])+'SCATTER',nscatind)
  detind = where(strupcase(phtags) eq 'NDET'+strupcase(ufilt[i]),ndetind)
  phot.(magind) = avgmag[*,i]
  phot.(errind) = avgerr[*,i]
  phot.(scatind) = scatter[*,i]
  phot.(detind) = ndet[*,i]
endfor

;; Calculate SFD E(B-V)
GLACTC,phot.ra,phot.dec,2000.0,glon,glat,1,/deg
phot.ebv = dust_getval(glon,glat,/noloop,/interp)

outfile = base+'_phot.fits'
print,'Writing final photometry to ',outfile
mwrfits,phot,outfile,/create


stop

end
