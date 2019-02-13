pro delvered_prep,delvedir,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo

;; This pre-processing script gets DELVE and community MC data ready
;; from the NOAO mass store to be processed with PHOTRED.

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/dl1/users/dnidever/delve/'
if FILE_TEST(delvedir,/directory) eq 0 then FILE_MKDIR,delvedir
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
expfile = '/home/dnidever/projects/delvered/data/decam_mcs_20181009.fits.gz'
if n_elements(scriptsdir) gt 0 then scriptsdir=trailingslash(scriptsdir) else scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
if n_elements(irafdir) gt 0 then irafdir=trailingslash(irafdir) else irafdir='/home/dnidever/iraf/'
if n_elements(workdir) eq 0 then begin
  undefine,workdir
  host = getenv('HOST')
  hostname = first_el(strsplit(host,'.',/extract))
  workdir='/data0/dnidever/delve/'
  if n_elements(workdir) gt 0 then if FILE_TEST(workdir,/directory) eq 0 then FILE_MKDIR,workdir
endif
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'

;; Print information
print,'--------------------------------------------------------------------'
print,' PREPARING DELVE MC EXPOSURES FOR PHOTRED EXPOSURE-LEVEL PROCESSING'
print,'--------------------------------------------------------------------'
print,'DELVEDIR = ',delvedir
print,'SCRIPTSDIR = ',scriptsdir
print,'IRAFDIR = ',irafdir
print,'EXPOSUREDIR = ',expdir
print,'WORKINGDIR = ',workdir
print,'HOST = ',host

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)

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
         ' wcs',$
         ' daophot',$
         ' zeropoint',$
         ' match',$
         '#allframe',$
         ' apcor',$
         ' astrom',$
         ' calib',$
         ' combine',$
         ' deredden',$
         ' save',$
         '#html']


;; Get the catalog of DECam exposures
expfile = delvereddir+'data/decam_mcs_20181009.fits.gz'
print,'' & print,'Loading ',expfile
allexpstr = MRDFITS(expfile,1,/silent)
nallexp = n_elements(allexpstr)
print,strtrim(nallexp,2),' DECam exposures'
allexpstr.fluxfile = strtrim(allexpstr.fluxfile,2)
allexpstr.maskfile = strtrim(allexpstr.maskfile,2)
allexpstr.wtfile = strtrim(allexpstr.wtfile,2)
allexpstr.expnum = strtrim(allexpstr.expnum,2)
;; Fix names for hulk/thing
;;  /net/mss1 -> /mss1
allexpstr.fluxfile = strmid(allexpstr.fluxfile,4)
allexpstr.maskfile = strmid(allexpstr.maskfile,4)
allexpstr.wtfile = strmid(allexpstr.wtfile,4)

;; STEP 1) Pick the exposures that we want
;;-----------------------------------------
;; Which ones do we want, make sure they haven't been prepped yet
;; only g/r/i and exptime>=90s for now
filt = strmid(allexpstr.filter,0,1)
gdexp = where((filt eq 'g' or filt eq 'r' or filt eq 'i') and allexpstr.exposure ge 90,ngdexp)
print,strtrim(ngdexp,2),' EXPOSURES pass the cuts'
expstr = allexpstr[gdexp]
;; Check if they've been prepped yet
observatory,'ctio',obs
night = strarr(ngdexp)
for i=0,ngdexp-1 do begin
  expstr1 = expstr[i]
  dateobs = expstr1.date_obs
  dateobs = strjoin(strsplit(dateobs,' ',/extract),'T')
  ;; Get night number
  jd = date2jd(dateobs)
  jd -= obs.tz/24.0  ;; correct for time zone difference
  caldat,double(long(jd)),Month, Day, Year, Hour, Minute, Second  ;; get date for beginning of the night
  ;; Convert to YYYYMMDD format
  night1 = string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')
  night[i] = night1
endfor
;; Get unique nights
night_index = CREATE_INDEX(night)
nnights = n_elements(night_index.value)
print,strtrim(nnights,2),' unique nights of data'


;; KLUDGE
print,'REMOVING SMASH NIGHTS'
smashnights = file_search('/dl1/users/dnidever/smash/cp/red/photred/20??????/',/test_directory,count=nsmashnights)
smashnights = file_basename(smashnights)
match,night,smashnights,ind1,ind2,count=nmatch
if nmatch gt 0 then REMOVE,ind1,night
;night = night[ind1]
night_index = CREATE_INDEX(night)
nnights = n_elements(night_index.value)

;; Loop over the nights
;;---------------------------
;For n=17,nnights-1 do begin
;For n=61,nnights-1 do begin
For n=0,nnights-1 do begin
  inight = night_index.value[n]
  ind = night_index.index[night_index.lo[n]:night_index.hi[n]]
  nind = night_index.num[n]
  expstr1 = expstr[ind]
  print,'' & print,'----------------------------------'
  print,'Night ',strtrim(n+1,2),' ',strtrim(inight,2),' - ',strtrim(nind,2),' exposures'
  print,'==================================' & print,''

  ;; Get the WCS.inlist and WCS.success files and see if we already
  ;; have these exposures
  nightdir = expdir+strtrim(inight,2)+'/'
  setupfile = nightdir+'photred.setup'
  if file_test(nightdir,/directory) eq 1 then begin
    if file_test(setupfile) eq 0 then WRITELINE,setupfile,setup
    undefine,files
    READLIST,nightdir+'logs/WCS.inlist',inlist,/unique,setupdir=nightdir,count=ninlist,/silent
    if ninlist gt 0 then PUSH,files,inlist
    READLIST,nightdir+'logs/WCS.success',successlist,/unique,setupdir=nightdir,count=nsuccess,/silent
    if nsuccess gt 0 then PUSH,files,successlist
    ;; Some previous files found
    if n_elements(files) gt 0 then begin
      ;; Make these relative paths
      for k=0,n_elements(files)-1 do begin
        pos = strpos(files[k],inight)
        if pos ge 0 then begin
          len = strlen(inight)
          files[k] = strmid(files[k],pos+len+1)
        endif
        ;; I needed this for SMASH, but not for new data
        ;if strmid(files[k],4,5,/reverse_offset) eq '.fits' and file_test(files[k]) eq 0 then files[k]+='.fz'
      endfor
      files = files[uniq(files,sort(files))]  ;; only want unique ones
      ;; Make sure they exist
      gdf = where(file_test(nightdir+files) eq 1,ngdf,comp=bdf,ncomp=nbdf)
      if ngdf gt 0 then files=files[gdf] else undefine,files
    endif
    nfiles = n_elements(files)
    ;; Remove the existing exposures
    if n_elements(files) gt 0 then begin
      fbase = PHOTRED_GETFITSEXT(files,/basename)
      base = reform((strsplitter(fbase,'-',/extract))[0,*])
      MATCH,base,expstr1.expnum,ind1,ind2,/sort,count=nmatch
      if nmatch gt 0 then begin
        exptoadd = expstr1
        if nmatch lt nind then begin
          print,'Removing ',strtrim(nmatch,2),' exposures that already exist'
          REMOVE,ind2,exptoadd
        endif else undefine,exptoadd
      endif
    endif else exptoadd=expstr1
    ;; Load the existing "fields" file 
    if file_test(nightdir+'fields') eq 1 then begin
      oldfieldstr = IMPORTASCII(nightdir+'fields',fieldname=['shname','name'],fieldtypes=[7,7],/silent)
      add_tag,oldfieldstr,'fieldnum',0L,oldfieldstr
      oldfieldstr.fieldnum = long(strmid(oldfieldstr.shname,1))
    endif
  ;; No directory yet
  endif else begin
    print,strtrim(inight,2),' not started yet'
    FILE_MKDIR,[nightdir,nightdir+'logs']
    WRITELINE,setupfile,setup
    exptoadd = expstr1
  endelse
  ntoadd = n_elements(exptoadd)
  print,strtrim(ntoadd,2),' exposures to add'
  if ntoadd eq 0 then goto,NIGHTBOMB

  ;; Group the exposures
  fields = strarr(ntoadd)
  fieldcount = 1L
  if n_elements(oldfieldstr) gt 0 then fieldcount=max(oldfieldstr.fieldnum)+1
  undefine,newfieldstr
  if ntoadd gt 1 then begin
    matches = MATCHALL_SPH(exptoadd.ra,exptoadd.dec,exptoadd.ra,exptoadd.dec,0.1,nmatches)
    for j=0,ntoadd-1 do begin
      if matches[j+1] ne matches[j] then begin
        find = matches[matches[j]:matches[j+1]-1]
        ;; New field
        ;;  if min(ind)<j then it was found from a previous exposure
        if min(find) eq j then begin
          fieldname = 'F'+strtrim(fieldcount,2)
          fieldstr1 = {shname:fieldname,name:'',fieldnum:fieldcount}
          push,newfieldstr,fieldstr1
          fields[find] = fieldname
          fieldcount++
        endif
      endif
    endfor
  endif else begin
    fields = 'F'+strtrim(fieldcount,2)
    fieldstr1 = {shname:fields,name:'',fieldnum:fieldcount}
    push,newfieldstr,fieldstr1
  endelse
  field_index = CREATE_INDEX(fields)
  nfields = n_elements(field_index.num)
  print,strtrim(nfields,2),' fields found'
  print,''

  ;; Loop over the fields
  undefine,outfiles
  outfiles = strarr(62*ntoadd)
  ocount = 0LL
  For f=0,nfields-1 do begin
    find = field_index.index[field_index.lo[f]:field_index.hi[f]]
    nfind = n_elements(find)
    ifield = field_index.value[f]
    print,ifield,' - ',strtrim(nfind,2),' exposures'
    fexptoadd = exptoadd[find]
    ;; Create the field directory
    fielddir = nightdir+ifield+'/'
    if FILE_TEST(fielddir,/directory) eq 0 then FILE_MKDIR,fielddir

    ;; Get Gaia DR2 and other reference data for this field
    if FILE_TEST(nightdir+'refcat/',/directory) eq 0 then FILE_MKDIR,nightdir+'refcat/'
    savefile = nightdir+'refcat/'+ifield+'_refcat.fits'
    refcat = DELVERED_GETREFDATA(['c4d-g','c4d-r','c4d-i'],fexptoadd[0].ra,fexptoadd[0].dec,1.2,savefile=savefile)
    SPAWN,['gzip',savefile],/noshell

    ;; Loop over exposures
    For e=0,nfind-1 do begin
      filename = fexptoadd[e].fluxfile
      filebase = PHOTRED_GETFITSEXT(filename,/basename)
      ;; Get number of extensions
      ;;   use symlink to make fits_open think it's a normal FITS file
      tmpfile = MKTEMP('tmp',/nodot) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
      tmpfile += '.fits'
      FILE_LINK,filename,tmpfile
      FITS_OPEN,tmpfile,fcb & FITS_CLOSE,fcb
      ;; Get the field longname
      if e eq 0 then begin
        hd0 = HEADFITS(tmpfile,exten=0)
        object = sxpar(hd0,'OBJECT')
        if strtrim(object,2) eq '' then object='Exposure'+strtrim(fexptoadd[e].expnum,2)
        MATCH,newfieldstr.shname,ifield,indfield
        newfieldstr[indfield].name = strcompress(object,/remove_all)
      endif

      ;; Get the CCDNUM for the extensions
      MATCH,decam.name,fcb.extname,ind1,ind2,/sort,count=nmatch
      extnum = ind2
      ccdnum = decam[ind1].ccdnum
      ;; Loop over the extensions and create the resource files
      For c=0,fcb.nextend-1 do begin
        schip = string(ccdnum[c],format='(i02)')
        ;; Create chip directory if needed
        chipdir = fielddir+'chip'+schip+'/'
        if FILE_TEST(chipdir,/directory) eq 0 then FILE_MKDIR,chipdir
        outfile1 = chipdir+ifield+'-'+fexptoadd[e].expnum+'_'+schip+'.fits'
        WRITELINE,outfile1,''
        routfile1 = chipdir+'.'+ifield+'-'+fexptoadd[e].expnum+'_'+schip+'.fits'
        rlines = ['fluxfile = '+fexptoadd[e].fluxfile+'['+strtrim(extnum[c],2)+']',$
                  'wtfile = '+fexptoadd[e].wtfile+'['+strtrim(extnum[c],2)+']',$
                  'maskfile = '+fexptoadd[e].maskfile+'['+strtrim(extnum[c],2)+']']
        WRITELINE,routfile1,rlines
        outfiles[ocount] = ifield+'/chip'+schip+'/'+ifield+'-'+fexptoadd[e].expnum+'_'+schip+'.fits'  ; relative path
        ocount++

        ;; Save the reference catalog for this chip
        hd = HEADFITS(tmpfile,exten=extnum[c],/silent)
        ;; Temporarily fix NAXIS1/2 values
        sxaddpar,hd,'NAXIS1',sxpar(hd,'ZNAXIS1')
        sxaddpar,hd,'NAXIS2',sxpar(hd,'ZNAXIS2')
        nx = sxpar(hd,'naxis1')
        ny = sxpar(hd,'naxis2')
        HEAD_XYAD,hd,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
        gdrefcat = where(refcat.ra ge min(vra)-0.02 and refcat.ra le max(vra)+0.02 and $
                         refcat.dec ge min(vdec)-0.02 and refcat.dec le max(vdec)+0.02,ngdrefcat)
        refcat1 = refcat[gdrefcat]
        refcatfile = chipdir+ifield+'-'+fexptoadd[e].expnum+'_'+string(ccdnum[c],format='(i02)')+'_refcat.fits'
        MWRFITS,refcat1,refcatfile,/create
        SPAWN,['gzip',refcatfile],/noshell
      Endfor  ; chip loop
      FILE_DELETE,tmpfile,/allow  ; delete the temporary symlink
    Endfor  ; exposure loop
  Endfor  ; field loop
  outfiles = outfiles[0:ocount-1]  ;; trim outfiles

  ;; Create the WCS.inlist file
  if FILE_TEST(nightdir+'logs/',/directory) eq 0 then FILE_MKDIR,nightdir+'logs'
  WRITELINE,nightdir+'logs/WCS.inlist',outfiles

  ;; Write the "fields" file
  if n_elements(oldfieldstr) gt 0 then newfieldstr=[oldfieldstr,newfieldstr]
  WRITELINE,nightdir+'fields',newfieldstr.shname+'   '+newfieldstr.name

  ;; Copy scripts into the directory
  ;;  I DON'T THINK I NEED TO DO THIS
  NIGHTBOMB:
Endfor  ; night loop


;stop

end
