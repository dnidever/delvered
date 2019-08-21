;+
;
; DELVERED_FORCEBRICK
;
; Process a single DELVE brick and perform ALLFRAME FORCED photometry
;
; INPUTS:
;  brick   The DELVE brick name, e.g. 1234m045
;
; By D. Nidever  August 2019
;-

pro delvered_forcebrick,brick,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo,logfile=logfile

;; This bricks pre-processing script gets DELVE and community MC data ready
;; to run PHOTRED ALLFRAME on it.

t0 = systime(1)
CD,current=curdir

;; Not enough inputs
if n_elements(brick) eq 0 then begin
  print,'Syntax - delvered_forcebrick,brick,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo,logfile=logfile'
  return
endif

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/dl1/users/dnidever/delve/'
if FILE_TEST(delvedir,/directory) eq 0 then FILE_MKDIR,delvedir
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
;expfile = '/home/dnidever/projects/delvered/data/decam_mcs_20181009.fits.gz'
if n_elements(scriptsdir) gt 0 then scriptsdir=trailingslash(scriptsdir) else scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
if n_elements(irafdir) gt 0 then irafdir=trailingslash(irafdir) else irafdir='/home/dnidever/iraf/'
tempdir = '/tmp/'
if n_elements(workdir) eq 0 then begin
  undefine,workdir
  host = getenv('HOST')
  hostname = first_el(strsplit(host,'.',/extract))
  workdir = '/data0/dnidever/delve/'
  if n_elements(workdir) gt 0 then if FILE_TEST(workdir,/directory) eq 0 then FILE_MKDIR,workdir
endif
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Bricks directory
brickdir = trailingslash(delvedir)+'bricks/'
if file_test(brickdir,/directory) eq 0 then file_mkdir,brickdir
logsdir = brickdir+'logs/'
if file_test(logsdir,/directory) eq 0 then file_mkdir,logsdir
if n_elements(logfile) eq 0 then logfile=-1

; Start the logfile
;------------------
; format is delvered_forcebrick.brick.host.DATETIME.log
host = GETENV('HOST')
hostname = first_el(strsplit(host,'.',/extract))
jd = systime(/julian)
caldat,jd,month,day,year,hour,minute,second
smonth = strtrim(month,2)
if month lt 10 then smonth = '0'+smonth
sday = strtrim(day,2)
if day lt 10 then sday = '0'+sday
syear = strmid(strtrim(year,2),2,2)
shour = strtrim(hour,2)
if hour lt 10 then shour='0'+shour
sminute = strtrim(minute,2)
if minute lt 10 then sminute='0'+sminute
ssecond = strtrim(round(second),2)
if second lt 10 then ssecond='0'+ssecond
logtime = smonth+sday+syear+shour+sminute+ssecond
journalfile = logsdir+'delvered_brick.'+brick+'.'+hostname+'.'+logtime+'.log'
JOURNAL,journalfile

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)

;; Load the brick information
brickstr = MRDFITS(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz',1)

;; Get the brick information
bind = where(brickstr.brickname eq brick,nbind)
if nbind eq 0 then begin
  printlog,logfile,ibrick+' not in DELVE-MC brick list'
  return
endif
brickstr1 = brickstr[bind[0]]

;; Subdirectory is the first 4 digits of the brickname, e.g., 0952 of 0952m462, the RA portion of the name
subdir = brickdir+strmid(brick,0,4)+'/'
if file_test(subdir,/directory) eq 0 then file_mkdir,subdir
;; Create brick directory if it doesn't exist yet
bdir = subdir+brick+'/'
if file_test(bdir,/directory) eq 0 then file_mkdir,bdir
logfile = bdir+brick+'.'+logtime+'.log'

;; DECam imager
thisimager = {telescope:'BLANCO',instrument:'DECam',namps:62,separator:'_'}

;; Print information
printlog,logfile,'--------------------------------------------------------------------'
printlog,logfile,' Run ALLFRAME forced photometry on DELVE Brick = ',brick
printlog,logfile,'--------------------------------------------------------------------'
printlog,logfile,'DELVEDIR = ',delvedir
printlog,logfile,'SCRIPTSDIR = ',scriptsdir
printlog,logfile,'IRAFDIR = ',irafdir
printlog,logfile,'EXPOSUREDIR = ',expdir
printlog,logfile,'WORKINGDIR = ',workdir
printlog,logfile,'HOST = ',host



;; The photed.setup file
setup = ['##### REQUIRED #####',$
         'scriptsdir  '+scriptsdir,$
         'irafdir     '+irafdir,$
         'telescope   Blanco',$
         'instrument  DECAM',$
         'observatory CTIO',$
         'nmulti      10',$
         'nmulti_wcs       40',$
         'nmulti_daophot   30',$
         'nmulti_allframe  10',$
         'filtref     g,i,r,z,u',$
         'trans       delve.trans',$
         '##### OPTIONAL #####',$
         'sepfielddir  1',$
         'sepchipdir   1',$
         'keepmef      0',$
         'catformat    FITS',$
         'workdir      '+workdir,$
         'clean        1',$
         'skipcheck    1',$
         'redo         0',$
         'wcsrefname   GAIADR2',$
         'searchdist   20',$
         '#wcsrmslim   1.0',$
         'hyperthread   1',$
         'daopsfva      1',$
         'daofitradfwhm 1.0',$
         'psfcomsrc     0',$
         'psfcomglobal  0',$
         'psfcomgauss   0',$
         '#mchmaxshift  50.0',$
         'finditer      2',$
         'alfdetprog  sextractor',$
         '#alfnocmbimscale 0',$
         'alftrimcomb   0',$
         '#ddo51radoffset  1',$
         'cmbforce      1',$
         'keepinstr     1',$
         'avgmag        1',$
         'avgonlymag    0',$
         'todered       u,g,r,i,z,g-i',$
         '##### STAGES #####',$
         '#rename',$
         '#split',$
         '#wcs',$
         '#daophot',$
         '#zeropoint',$
         '#match',$
         ' allframe',$
         '#apcor',$
         '#astrom',$
         '#calib',$
         '#combine',$
         '#deredden',$
         '#save',$
         '#html']
WRITELINE,bdir+'photred.setup',setup
PHOTRED_LOADSETUP,setup

;; Print out some information about the brick
printlog,logfile,'RA = ',stringize(brickstr1.ra,ndec=5)
printlog,logfile,'DEC = ',stringize(brickstr1.dec,ndec=5)
printlog,logfile,'RA range  = [ ',stringize(brickstr1.ra1,ndec=5),',',stringize(brickstr1.ra2,nec=5),' ]'
printlog,logfile,'DEC range = [ ',stringize(brickstr1.dec1,ndec=5),',',stringize(brickstr1.dec2,nec=5),' ]'


;; Step 1: Get the list of exposures/chips that overlap this brick
;;----------------------------------------------------------------
printlog,logfile,'Step 1: Get list of chip files that overlap this brick'
;printlog,logfile,'Getting chips that overlap this brick'
cenra = brickstr1.ra
cendec = brickstr1.dec
tmpfile = MKTEMP('tmp',/nodot,outdir=tempdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
tmpfile += '.fits'
spawn,['query_delvered_table',strtrim(cenra,2),strtrim(cendec,2),tmpfile,'--lim','0.5'],out,errout,/noshell
info = file_info(tmpfile)
if info.size eq 0 then begin
  printlog,logfile,'No overlapping chips found'
  return
endif
chstr = MRDFITS(tmpfile,1,/silent)
file_delete,tmpfile,/allow
nchstr = n_elements(chstr)
printlog,logfile,'Found ',strtrim(nchstr,2),' overlapping chips within 0.5 deg of brick center'

;; Do more rigorous overlap checking
bvra = [brickstr1.ra1,brickstr1.ra2,brickstr1.ra2,brickstr1.ra1]
bvdec = [brickstr1.dec1,brickstr1.dec1,brickstr1.dec2,brickstr1.dec2]
olap = intarr(nchstr)
for i=0,nchstr-1 do begin
  hd1 = PHOTRED_READFILE(chstr[i].file,/header)
  nx = sxpar(hd1,'naxis1')
  ny = sxpar(hd1,'naxis2')
  head_xyad,hd1,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
  olap[i] = dopolygonsoverlap(bvra,bvdec,vra,vdec)
endfor
g = where(olap eq 1,ng)
if ng eq 0 then begin
  printlog,logfile,'No chips overlap this brick'
  return
endif
printlog,logfile,strtrim(ng,2),' chips overlap this brick'
chstr = chstr[g]
nchstr = ng

;; APPLY Quality cuts on the exposures
;;------------------------------------
printlog,logfile,'Applying quality cuts'
fwhmthresh = 2.0                ; seeing 2.0" threshold
printlog,logfile,'Applying quality cuts'
gdch = where(chstr.fwhm*chstr.pixscale le fwhmthresh,ngdch)
if ngdch eq 0 then begin
  printlog,logfile,'No chips passed the quality cuts'
  return
endif
printlog,logfile,strtrim(ngdch,2),' chips passed the quality cuts'
chstr = chstr[gdch]
nchstr = ngdch
chstr.file = strtrim(chstr.file,2)
chstr.base = strtrim(chstr.base,2)
add_tag,chstr,'newbase','',chstr

;; Create the symlinks
;;---------------------
printlog,logfile,'Copying over the necessary files to ',bdir
;; Chip loop
for i=0,nchstr-1 do begin
  printlog,logfile,'Copying over files for ',chstr[i].file
  ;; New filename
  chstr[i].newbase = 'chip'+string(chstr[i].chip,format='(i02)')+'/'+chstr[i].base
  ;; Make chip directory if necessary
  chdir = bdir+'chip'+string(chstr[i].chip,format='(i02)')+'/'
  if file_test(chdir,/directory) eq 0 then file_mkdir,chdir
  odir = file_dirname(chstr[i].file)+'/'
  obase = PHOTRED_GETFITSEXT(chstr[i].file,/basename)
  ;; Need symlinks to .psf, .als
  if file_test(odir+obase+'.psf') eq 0 then begin
    printlog,logfile,odir+obase+'.psf NOT FOUND.  Skipping this chip'
    goto,BOMB1
  endif
  if file_test(odir+obase+'.als') eq 0 then begin
    printlog,logfile,odir+obase+'.als NOT FOUND.  Skipping this chip'
    goto,BOMB1
  endif
  FILE_DELETE,chdir+chstr[i].base+['.psf','.als','.ap','.opt','.als.opt','.log'],/allow
  FILE_LINK,odir+obase+'.psf',chdir+chstr[i].base+'.psf'
  FILE_LINK,odir+obase+'.als',chdir+chstr[i].base+'.als'
  FILE_LINK,odir+obase+'.ap',chdir+chstr[i].base+'.ap'
  FILE_LINK,odir+obase+'.opt',chdir+chstr[i].base+'.opt'
  FILE_LINK,odir+obase+'.als.opt',chdir+chstr[i].base+'.als.opt'
  FILE_LINK,odir+obase+'.log',chdir+chstr[i].base+'.log'
  ;; Copy the fits, fits resource file and header files locally
  if file_test(odir+obase+'.fits') eq 0 then begin
    printlog,logfile,odir+obase+'.fits NOT FOUND.  Skipping this chip'
    goto,BOMB1
  endif
  FILE_COPY,odir+obase+'.fits',chdir,/over
  if file_test(odir+'.'+obase+'.fits') eq 1 then FILE_COPY,odir+'.'+obase+'.fits',chdir,/over
  if file_test(odir+obase+'.head') eq 1 then FILE_COPY,odir+obase+'.head',chdir,/over
  BOMB1:
endfor

; Make the tiling file
;---------------------
undefine,lines
; Lines with the tiling scheme first
nx = 3600
ny = 3600
step = 0.25
xref = nx/2
yref = ny/2
ntiles = 1
push,lines,'CENRA  = '+strtrim(brickstr1.ra,2)
push,lines,'NX     = '+strtrim(nx,2)
push,lines,'XSTEP  = '+strtrim(step,2)
push,lines,'XREF   = '+strtrim(xref+1,2)
push,lines,'CENDEC = '+strtrim(brickstr1.dec,2)
push,lines,'NY     = '+strtrim(ny,2)
push,lines,'YSTEP  = '+strtrim(step,2)
push,lines,'YREF   = '+strtrim(yref+1,2)
push,lines,'NTILES = '+strtrim(ntiles,2)
; Then one file with info for each tile
tilestr = {num:1,name:'tile1',x0:0,x1:nx-1,nx:nx,y0:0,y1:ny-1,ny:ny,nimages:nchstr}
ntiles = 1
for i=0,ntiles-1 do begin
  tilestr1 = tilestr[i]
  fmt = '(I-6,A8,6I8,I6)'
  ; Using IRAF indexing
  line = string(format=fmt,tilestr1.num,tilestr1.name,tilestr1.x0+1,tilestr1.x1+1,tilestr1.nx,$
                tilestr1.y0+1,tilestr1.y1+1,tilestr1.ny,tilestr1.nimages)
  push,lines,line
endfor
tilefile = bdir+brick+'.tiling'
printlog,logfile,'Writing tiling information to >>'+tilefile+'<<'
WRITELINE,tilefile,lines

;;  Make the header as well
MKHDR,tilehead,fltarr(5,5)
SXADDPAR,tilehead,'NAXIS1',nx2
SXADDPAR,tilehead,'CDELT1',step
SXADDPAR,tilehead,'CRPIX1',xref+1L
SXADDPAR,tilehead,'CRVAL1',brickstr1.ra
SXADDPAR,tilehead,'CTYPE1','RA---TAN'
SXADDPAR,tilehead,'NAXIS2',ny
SXADDPAR,tilehead,'CDELT2',step
SXADDPAR,tilehead,'CRPIX2',yref+1L
SXADDPAR,tilehead,'CRVAL2',brickstr1.dec
SXADDPAR,tilehead,'CTYPE2','DEC--TAN'
EXTAST,tilehead,tileast


;; Step 2: Run DAOMATCH_TILE.PRO on the files
;;--------------------------------------------
printlog,logfile,'Step 2: Matching up objects with DAOPHOT_TILE'
DAOMATCH_TILE,chstr.newbase+'.als',tilestr
mchfile = chstr[0].newbase+'.mch'

;; Step 3: Run ALLFRAME
;;----------------------
printlog,logfile,'Step 3: Run ALLFRAME'
ALLFRAME,mchfile,tile=tilestr,setupdir=curdir,scriptsdir=scriptsdir,irafdir=irafdir,$
         logfile=logfile,catformat='FITS',imager=thisimager,workdir=workdir


;; Step 4: Calculate coordinates
;;-------------------------------
printlog,logfile,'Step 4: Adding coordinates'
; Load the MCH file
LOADMCH,mchfile,alsfiles
nalsfiles = n_elements(alsfiles)
; Load the photometry file
phot = PHOTRED_READFILE(magfile)
nphot = n_elements(phot)
printlog,logfile,'Nstars = '+strtrim(nphot,2)
;; Converting to IDL X/Y convention, starting at (0,0)
;; DAOPHOT has X/Y start at (1,1)
HEAD_XYAD,tilehead,phot.x-1.0,phot.y-1.0,ra,dec,/degree


;; Step 5: Calibrating photometry with zero-points
;;-------------------------------------------------
printlog,logfile,'Step 5: Calibrating photometry with zero-points'
cmag = fltarr(nphot,nchstr)+99.99
cerr = fltarr(nphot,nchstr)+9.99
;; Chip loop
for i=0,nchstr-1 do begin
  ; id, x, y, unsolved magnitudes, chi, sharp
  imag = phot.(2*i+3)
  ierr = phot.(2*i+4)
  gdmag = where(imag lt 50,ngdmag)
  if ngdmag gt 0 then begin
    ;; exptime, aperture correction, zero-point
     ;; aperture correction is SUBTRACTIVE, makes it brighter
    ;; ZPTERM is a SUBTRACTIVE constant offset
    cmag[gdmag,i] = imag[gdmag] + 2.5*alog10(chstr[i].exptime) - chstr[i].apcor - chstr[i].zpterm
    ;; Add zero-point error in quadrature
    cerr[gdmag,i] = sqrt(ierr[gdmag]^2+chstr[i].zptermerr^2)
  endif
endfor
;; Calculate average photometry per filter
ufilt = chstr[uniq(chstr.filter,sort(chstr.filter))].filter
nufilt = n_elements(ufile)
avgmag = fltarr(nphot,nufilt)
avgerr = fltarr(nphot,nufilt)
for i=0,nufilt-1 do begin
  gdf = where(chstr.filter eq ufilt[i],ngdf)
  ;; Single exposure
  if ngdf eq 0 then begin
    avgmag[*,i] = cmag[*,gdf[0]]
    avgerr[*,i] = cerr[*,gdf[0]]
  ;; Multiple exposures
  endif else begin
    ;; Loop through all of the exposures and add up the flux, totalwt, etc.
    totalwt = dblarr(nphot)
    totalfluxwt = dblarr(nphot)
    for k=0,ngdf-1 do begin
      gdmag = where(cmag[*,gdf[k]] lt 50,ngdmag)
      if ngdmag gt 0 then begin
        totalwt[gdmag] += 1.0d0/cerr[gdmag,gdf[k]]^2
        totalfluxwt[gdag] += 2.5118864d^cmag[gdmag,gdf[k]] * (1.0d0/cerr[gdmag,gdf[k]]^2)
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
;; Create final catalog schema
newschema = {id:'',x:0.0,y:0.0,ra:0.0d0,dec:0.0d0}
;; Add columns for calibrated single-epoch photometry columns
cmagnames = strarr(nchstr)
cerrnames = strarr(nchstr)
for i=0,nufilt-1 do begin
  ind = where(chstr.filter eq ufilt[i],nind)
  cmagnames[ind] = strupcase(ufilt[i])+strtrim(lindgen(ind)+1,2)+'MAG'
  cerrnames[ind] = strupcase(ufilt[i])+strtrim(lindgen(ind)+1,2)+'ERR'
endfor
for i=0,nchstr-1 do newschema = create_struct(newschema,cmagnames[i],0.0,cerrnames[i],0.0)
;; Add columns for average photometry per filter
for i=0,nufilt-1 do newschema = create_struct(newschema,ufilt[i]+'MAG',0.0,ufilt[i]+'ERR')
;; Extra columns
newschema = create_struct(newschema,'chi',0.0,'sharp',0.0,'prob',0.0,'ebv',0.0)
;; Create final catalog
orig = phot
undefine,phot
phot = replicate(newschema,nphot)
struct_assign,orig,phot,/nozero
undefine,orig
phtags = tag_names(phot)
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
  phot.(magind) = avgmag[*,i]
  phot.(errind) = avgerr[*,i]
endfor

;; Calculate SFD E(B-V)
GLACTC,phot.ra,phot.dec,2000.0,glon,glat,1,/deg
phot.ebv = dust_getval(glon,glat,/noloop,/interp)

;; Saving final catalog
photfile = bdir+brick+'.fits'
printlog,logfile,'Writing photometry to ',photfile,'.gz'
MWRFITS,phot,photfile,/create
spawn,['gzip','-f',photfile],/noshell

;; Save metadata
metafile = bdir+brick+'_meta.fits'
printlog,logfile,'Writing meta-data to ',metafile
MWRFITS,chstr,metafile,/create

printlog,logfile,'DELVERED_FORCEBRICK done after '+strtrim(systime(1)-t0,2)+' sec.'

stop

JOURNAL

end
