pro make_smash_symlinks_night,inight

;; Make symlinks to the SMASH data that has already been processedd
;; For a SINGLE NIGHT

delvedir = '/dl1/users/dnidever/delve/exposures/'
smashdir = '/dl1/users/dnidever/smash/cp/red/photred/'
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


  CD,smashdir+inight
  CD,current=nightdir
  nightdir = trailingslash(nightdir)

  print,''
  print,'Making SMASH symlinks for night ',inight
  print,'-----------------------------------------'

  ;; Make sure the DELVE night directory exists
  if FILE_TEST(delvedir+inight,/directory) eq 0 then FILE_MKDIR,delvedir+inight

  ;; Load the fields file
  fieldstr = IMPORTASCII(nightdir+'fields',fieldname=['shname','name'],fieldtypes=[7,7],/silent)

  ;; Load from daophot.success
  undefine,fitsfiles
  READLIST,'logs/WCS.success',wcsfiles,setupdir='.',count=nwcsfiles,/silent
  if nwcsfiles gt 0 then push,fitsfiles,wcsfiles
  READLIST,'logs/DAOPHOT.success',daofiles,setupdir='.',count=ndaofiles,/silent
  if ndaofiles gt 0 then push,fitsfiles,daofiles
  nfitsfiles = n_elements(fitsfiles)
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
  undefine,matchlines
  wcslines = strarr(nfitsfiles)
  daophotlines = strarr(nfitsfiles)
  count = 0LL
  For f=0,nfields-1 do begin
    ind = field_index.index[field_index.lo[f]:field_index.hi[f]]
    nind = n_elements(ind)
    ifield = field_index.value[f]
    fitsfiles1 = fitsfiles[ind]
    MATCH,fieldstr.shname,ifield,fieldind,count=nmatch
    if nmatch eq 0 then begin
      print,ifield,' NOT FOUND in >>fields<< file'
      goto,FIELDBOMB
    endif
    fieldname = fieldstr[fieldind].name   ;; long field name
    print,strtrim(f+1,2),'/',strtrim(nfields,2),' ',ifield,' ',fieldname,' ',strtrim(nind,2),' files'

    ;; Make sure the DELVE night+field directory exists
    if FILE_TEST(delvedir+inight+'/'+ifield,/directory) eq 0 then FILE_MKDIR,delvedir+inight+'/'+ifield

    ;; File loop
    For k=0,nind-1 do begin
      ;; _cat.dat, _refcat.dat
      ;; als, ap, coo, plst, fits/fits.fz, psf
      base1 = first_el(PHOTRED_GETFITSEXT(fitsfiles1[k],/basename))
      chipnum = PHOTRED_GETCHIPNUM(base1,{namps:62,separator:'_'})
      chipdir = delvedir+inight+'/'+ifield+'/chip'+strtrim(chipnum,2)+'/'
      if file_test(chipdir,/directory) eq 0 then FILE_MKDIR,chipdir
      ;; Check that files exist
      if file_test(nightdir+ifield+'/'+base1+'.opt') eq 1 then begin
        FILE_LINK,nightdir+ifield+'/'+base1+['_cat.dat','_refcat.dat','.opt','.als.opt','.als','.ap','.coo','.plst','.psf','.psf.log','.log','a.als','a.ap'],$
                  chipdir
        FILE_LINK,nightdir+ifield+'/'+file_basename(fitsfiles1[k]),chipdir

        ;; Update the lists, relative paths
        wcslines[count] = ifield+'/chip'+strtrim(chipnum,2)+'/'+file_basename(fitsfiles1[k])
        daophotlines[count] = ifield+'/chip'+strtrim(chipnum,2)+'/'+file_basename(fitsfiles1[k])
        count++
      endif else print,base1+' NOT FOUND'
    Endfor  ; file loop

    ;; MATCH files, mch, raw
    mchfiles = FILE_SEARCH(nightdir+ifield+'/'+ifield+'-????????_??.mch',count=nmchfiles)
    if nmchfiles gt 0 then begin
      rawfiles = FILE_DIRNAME(mchfiles)+'/'+FILE_BASENAME(mchfiles,'.mch')+'.raw'
      tfrfiles = FILE_DIRNAME(mchfiles)+'/'+FILE_BASENAME(mchfiles,'.mch')+'.tfr'
      ;; Get chip subdirectories
      chipdirs = strarr(nmchfiles)
      for k=0,nmchfiles-1 do chipdirs[k]='chip'+strtrim(PHOTRED_GETCHIPNUM(file_basename(mchfiles[k],'.mch'),{namps:62,separator:'_'}),2)
      FILE_LINK,mchfiles,delvedir+inight+'/'+ifield+'/'+chipdirs+'/'+file_basename(mchfiles)
      FILE_LINK,rawfiles,delvedir+inight+'/'+ifield+'/'+chipdirs+'/'+file_basename(rawfiles)
      FILE_LINK,tfrfiles,delvedir+inight+'/'+ifield+'/'+chipdirs+'/'+file_basename(tfrfiles)

      ;; Update the list, relative paths
      PUSH,matchlines,ifield+'/'+chipdirs+'/'+file_basename(mchfiles)
    endif

    sumfile = nightdir+fieldname+'_summary.fits'
    sumstr = MRDFITS(sumfile,1,/silent)

    cenra = median([sumstr.ra])
    cendec = median([sumstr.dec])

    ;; Get Gaia DR2 and other reference data for this field
    if FILE_TEST(delvedir+inight+'/refcat/',/directory) eq 0 then FILE_MKDIR,delvedir+inight+'/refcat/'
    savefile = delvedir+inight+'/refcat/'+ifield+'_refcat.fits'
    if file_test(savefile) eq 0 and file_test(savefile+'.gz') eq 0 then begin
      refcat = DELVERED_GETREFDATA(['c4d-u','c4d-g','c4d-r','c4d-i','c4d-z','c4d-Y'],cenra,cendec,1.2,savefile=savefile)
      SPAWN,['gzip',savefile],/noshell
    endif
    FIELDBOMB:
  Endfor  ; field loop

  ;; Copy apcor.lst and APCOR.success
  FILE_COPY,nightdir+'apcor.lst',delvedir+inight
  FILE_CHMOD,delvedir+inight+'/apcor.lst',/a_write
  READLINE,nightdir+'logs/APCOR.success',apcorlines,count=napcor
  bd = where(stregex(apcorlines,'.fits.fz',/boolean) eq 0,nbd)
  if nbd gt 0 then apcorlines[bd]+='.fz'
  for k=0,napcor-1 do apcorlines[k]=strmid(apcorlines[k],strpos(apcorlines[k],inight)+9) ;; make relative
  ;; Get chip subdirectories
  chipdirs = strarr(napcor)
  for k=0,napcor-1 do chipdirs[k]='chip'+strtrim(PHOTRED_GETCHIPNUM(file_basename(apcorlines[k],'.fits.fz'),{namps:62,separator:'_'}),2)
  apcorlines0 = apcorlines
  apcorlines = file_dirname(apcorlines)+'/'+chipdirs+'/'+file_basename(apcorlines)  ;; F7/chip34/F7-00421651_34.fits.fz

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


end
