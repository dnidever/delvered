pro make_smash_symlinks

;; Make symlinks to the SMASH data that has already been processedd

delvedir = '/dl1/users/dnidever/delve/exposures/'
smashdir = '/dl1/users/dnidever/smash/cp/red/photred/'
nights = FILE_SEARCH(smashdir+'20??????',/test_directory,count=nnights)
CD,current=origdir

;; Night loop
For i=0,nnights-1 do begin
  inight = FILE_BASENAME(nights[i])
  CD,nights[i]
  CD,current=nightdir
  nightdir = trailingslash(nightdir)

  ;; Make sure the DELVE night directory exists
  if FILE_TEST(delvedir+inight,/directory) eq 0 then FILE_MKDIR,delvedir+inight

  ;; Load from daophot.success
  READLIST,'logs/DAOPHOT.success',fitsfiles,setupdir='.',count=nfitsfiles,/silent
  ;; Convert fits to fits.fz
  allfield = strarr(nfitsfiles)
  for j=0,nfitsfiles-1 do begin
    ;; Make the names relative
    pos = strpos(fitsfiles[j],inight)
    len = strlen(inight)
    fitsfiles[j] = strmid(fitsfiles[j],pos+len+1)
    if strmid(fitsfiles[j],4,5,/reverse_offset) eq '.fits' and file_test(fitsfiles[j]) eq 0 then fitsfiles[j]+='.fz'
    allfield[j] = first_el(strsplit(file_basename(fitsfiles[j]),'-',/extract))
  endfor
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

    ;; Make sure the DELVE night+field directory exists
    if FILE_TEST(delvedir+inight+'/'+ifield,/directory) eq 0 then FILE_MKDIR,delvedir+inight+'/'+ifield

    ;; File loop
    For k=0,nind-1 do begin
      ;; _cat.dat, _refcat.dat
      ;; als, ap, coo, plst, fits/fits.fz, psf
      base1 = first_el(PHOTRED_GETFITSEXT(fitsfiles1[k],/basename))
      FILE_LINK,nightdir+ifield+'/'+base1+['_cat.dat','_refcat.dat','.als','.ap','.coo','.plst','.psf','.psf.log','.log','a.als','a.ap'],$
                delvedir+inight+'/'+ifield
      FILE_LINK,nightdir+ifield+'/'+file_basename(fitsfiles1[k]),delvedir+inight+'/'+ifield

      ;; Update the lists
      wcslines[count] = delvedir+inight+'/'+ifield+'/'+file_basename(fitsfiles1[k])
      daophotlines[count] = delvedir+inight+'/'+ifield+'/'+file_basename(fitsfiles1[k])
      count++
    Endfor  ; file loop

    ;; MATCH files, mch, raw
    mchfiles = FILE_SEARCH(nightdir+ifield+'/'+ifield+'*.mch',count=nmchfiles)
    rawfiles = FILE_DIRNAME(mchfiles)+'/'+FILE_BASENAME(mchfiles,'.mch')+'.raw'
    tfrfiles = FILE_DIRNAME(mchfiles)+'/'+FILE_BASENAME(mchfiles,'.mch')+'.tfr'
    FILE_LINK,mchfiles,delvedir+inight+'/'+ifield+'/'+file_basename(mchfiles)
    FILE_LINK,rawfiles,delvedir+inight+'/'+ifield+'/'+file_basename(rawfiles)
    FILE_LINK,tfrfiles,delvedir+inight+'/'+ifield+'/'+file_basename(tfrfiles)

    ;; Update the list
    PUSH,matchlines,delvedir+inight+'/'+ifield+'/'+file_basename(mchfiles)

    ;; I need to download the reference data

    stop
  Endfor  ; field loop

  ;; Copy over the WCS.success, DAOPHOT.success and
  ;; MATCH.success/outlist files

 

Endfor  ; night loop



stop


end
