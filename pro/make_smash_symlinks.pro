pro make_smash_symlinks

;; Make symlinks to the SMASH data that has already been processedd

delvedir = '/dl1/users/dnidever/delve/exposures/'
smashdir = '/dl1/users/dnidever/smash/cp/red/photred/'
nights = FILE_SEARCH(smashdir+'20??????',/test_directory,count=nnights)
CD,current=origdir

scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
irafdir = '/home/dnidever/iraf/'
workdir = '/data0/dnidever/delve/'

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

print,'######################################'
print,strtrim(nnights,2),' SMASH nights to create symlinks for'
print,'######################################'

;; Night loop
;For i=2,nnights-1 do begin
For i=0,nnights-1 do begin
  inight = FILE_BASENAME(nights[i])
  CD,nights[i]
  CD,current=nightdir
  nightdir = trailingslash(nightdir)

  print,''
  print,'Making SMASH symlinks for night ',strtrim(i+1,2),' ',inight
  print,'-------------------------------------------'

  ;; Make sure the DELVE night directory exists
  if FILE_TEST(delvedir+inight,/directory) eq 0 then FILE_MKDIR,delvedir+inight

  ;; Load the fields file
  fieldstr = IMPORTASCII(nightdir+'fields',fieldname=['shname','name'],fieldtypes=[7,7],/silent)

  ;; Load from daophot.success
  READLIST,'logs/DAOPHOT.success',fitsfiles,setupdir='.',count=nfitsfiles,/silent
  ;; Convert fits to fits.fz
  allfield = strarr(nfitsfiles)
  for j=0,nfitsfiles-1 do begin
    ;; Make the names relative
    pos = strpos(fitsfiles[j],inight)
    if pos ge 0 then begin
      len = strlen(inight)
      fitsfiles[j] = strmid(fitsfiles[j],pos+len+1)
    endif
    if strmid(fitsfiles[j],4,5,/reverse_offset) eq '.fits' and file_test(fitsfiles[j]) eq 0 then fitsfiles[j]+='.fz'
    allfield[j] = first_el(strsplit(file_basename(fitsfiles[j]),'-',/extract))
  endfor
  ;; Get unique ones
  ui = uniq(fitsfiles,sort(fitsfiles))
  fitsfiles = fitsfiles[ui]
  allfield = allfield[ui]
  nfitsfiles = n_elements(fitsfiles)
  ;; Create index
  field_index = CREATE_INDEX(allfield)
  nfields = n_elements(field_index.value)

  ;; Field loop
  undefine,wcslines,daophotlines,matchlines
  wcslines = strarr(nfitsfiles)
  daophotlines = strarr(nfitsfiles)
  count = 0LL
  For f=0,nfields-1 do begin
    ind = field_index.index[field_index.lo[f]:field_index.hi[f]]
    nind = n_elements(ind)
    ifield = field_index.value[f]
    fitsfiles1 = fitsfiles[ind]
    MATCH,fieldstr.shname,ifield,fieldind
    fieldname = fieldstr[fieldind].name   ;; long field name
    print,strtrim(f+1,2),'/',strtrim(nfields,2),' ',ifield,' ',fieldname,' ',strtrim(nind,2),' files'

    ;; Make sure the DELVE night+field directory exists
    if FILE_TEST(delvedir+inight+'/'+ifield,/directory) eq 0 then FILE_MKDIR,delvedir+inight+'/'+ifield

    ;; File loop
    For k=0,nind-1 do begin
      ;; _cat.dat, _refcat.dat
      ;; als, ap, coo, plst, fits/fits.fz, psf
      base1 = first_el(PHOTRED_GETFITSEXT(fitsfiles1[k],/basename))
      FILE_LINK,nightdir+ifield+'/'+base1+['_cat.dat','_refcat.dat','.opt','.als.opt','.als','.ap','.coo','.plst','.psf','.psf.log','.log','a.als','a.ap'],$
                delvedir+inight+'/'+ifield
      FILE_LINK,nightdir+ifield+'/'+file_basename(fitsfiles1[k]),delvedir+inight+'/'+ifield

      ;; Update the lists, relative paths
      wcslines[count] = ifield+'/'+file_basename(fitsfiles1[k])
      daophotlines[count] = ifield+'/'+file_basename(fitsfiles1[k])
      count++
    Endfor  ; file loop

    ;; MATCH files, mch, raw
    mchfiles = FILE_SEARCH(nightdir+ifield+'/'+ifield+'-????????_??.mch',count=nmchfiles)
    rawfiles = FILE_DIRNAME(mchfiles)+'/'+FILE_BASENAME(mchfiles,'.mch')+'.raw'
    tfrfiles = FILE_DIRNAME(mchfiles)+'/'+FILE_BASENAME(mchfiles,'.mch')+'.tfr'
    FILE_LINK,mchfiles,delvedir+inight+'/'+ifield+'/'+file_basename(mchfiles)
    FILE_LINK,rawfiles,delvedir+inight+'/'+ifield+'/'+file_basename(rawfiles)
    FILE_LINK,tfrfiles,delvedir+inight+'/'+ifield+'/'+file_basename(tfrfiles)

    ;; Update the list, relative paths
    PUSH,matchlines,ifield+'/'+file_basename(mchfiles)

    sumfile = nightdir+fieldname+'_summary.fits'
    sumstr = MRDFITS(sumfile,1,/silent)
    cenra = median(sumstr.ra)
    cendec = median(sumstr.dec)

    ;; Get Gaia DR2 and other reference data for this field
;    if FILE_TEST(delvedir+inight+'/refcat/',/directory) eq 0 then FILE_MKDIR,delvedir+inight+'/refcat/'
;    savefile = delvedir+inight+'/refcat/'+ifield+'_refcat.fits'
;    refcat = DELVERED_GETREFDATA(['c4d-g','c4d-r','c4d-i'],cenra,cendec,1.2,savefile=savefile)
;    SPAWN,['gzip',savefile],/noshell
  Endfor  ; field loop

  ;; Copy apcor.lst and APCOR.success
  FILE_COPY,nightdir+'apcor.lst',delvedir+inight
  FILE_CHMOD,delvedir+inight+'/apcor.lst',/a_write
  READLINE,nightdir+'logs/APCOR.success',apcorlines,count=napcor
  bd = where(stregex(apcorlines,'.fits.fz',/boolean) eq 0,nbd)
  if nbd gt 0 then apcorlines[bd]+='.fz'
  for k=0,napcor-1 do apcorlines[k]=strmid(apcorlines[k],strpos(apcorlines[k],inight)+9)

  ;; Copy over the WCS.success, DAOPHOT.success and
  ;; MATCH.success/outlist files
  if FILE_TEST(delvedir+inight+'/logs',/directory) eq 0 then FILE_MKDIR,delvedir+inight+'/logs'
  WRITELINE,delvedir+inight+'/logs/WCS.success',wcslines
  WRITELINE,delvedir+inight+'/logs/DAOPHOT.success',daophotlines
  WRITELINE,delvedir+inight+'/logs/MATCH.outlist',matchlines
  WRITELINE,delvedir+inight+'/logs/APCOR.success',apcorlines

  ;; Copy the "fields" file
  FILE_COPY,nightdir+'fields',delvedir+inight
  FILE_CHMOD,delvedir+inight+'/fields',/a_write  

  ;; Make the setup file
  WRITELINE,delvedir+inight+'/photred.setup',setup
 
  ;stop

Endfor  ; night loop



stop


end
