pro fix_smash_nights

;; chip0X directories missing for SMASH nights

delvedir = '/dl1/users/dnidever/delve/exposures/'
smashdir = '/dl1/users/dnidever/smash/cp/red/photred/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
CD,current=origdir

scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
irafdir = '/home/dnidever/iraf/'
workdir = '/data0/dnidever/delve/'
modeleqnfile = delvereddir+'params/modelmag_equations.txt'

fdir=['/dl1/users/dnidever/delve/exposures/20141218/F13',$
'/dl1/users/dnidever/delve/exposures/20141218/F14',$
'/dl1/users/dnidever/delve/exposures/20141218/F15',$
'/dl1/users/dnidever/delve/exposures/20141218/F16',$
'/dl1/users/dnidever/delve/exposures/20141218/F17',$
'/dl1/users/dnidever/delve/exposures/20141218/F8',$
'/dl1/users/dnidever/delve/exposures/20141218/F9',$
'/dl1/users/dnidever/delve/exposures/20160213/F3',$
'/dl1/users/dnidever/delve/exposures/20160213/F4',$
'/dl1/users/dnidever/delve/exposures/20160213/F5',$
'/dl1/users/dnidever/delve/exposures/20160213/F6',$
'/dl1/users/dnidever/delve/exposures/20160214/F1',$
'/dl1/users/dnidever/delve/exposures/20160214/F2',$
'/dl1/users/dnidever/delve/exposures/20160214/F3',$
'/dl1/users/dnidever/delve/exposures/20160214/F4',$
'/dl1/users/dnidever/delve/exposures/20160214/F5',$
'/dl1/users/dnidever/delve/exposures/20160214/F6',$
'/dl1/users/dnidever/delve/exposures/20160214/F7',$
'/dl1/users/dnidever/delve/exposures/20160215/F4',$
'/dl1/users/dnidever/delve/exposures/20160215/F5',$
'/dl1/users/dnidever/delve/exposures/20160215/F6',$
'/dl1/users/dnidever/delve/exposures/20160217/F18',$
'/dl1/users/dnidever/delve/exposures/20160217/F19',$
'/dl1/users/dnidever/delve/exposures/20160217/F2',$
'/dl1/users/dnidever/delve/exposures/20160217/F20',$
'/dl1/users/dnidever/delve/exposures/20160217/F21',$
'/dl1/users/dnidever/delve/exposures/20160217/F22',$
'/dl1/users/dnidever/delve/exposures/20160217/F23',$
'/dl1/users/dnidever/delve/exposures/20160217/F24',$
'/dl1/users/dnidever/delve/exposures/20160217/F3',$
'/dl1/users/dnidever/delve/exposures/20160217/F4',$
'/dl1/users/dnidever/delve/exposures/20160217/F5',$
'/dl1/users/dnidever/delve/exposures/20160217/F6',$
'/dl1/users/dnidever/delve/exposures/20160217/F7',$
'/dl1/users/dnidever/delve/exposures/20160217/F8',$
'/dl1/users/dnidever/delve/exposures/20160217/F9']

nights=['20141218', '20160213', '20160214', '20160215', '20160217']
nnights = n_elements(nights)

  ;; Load the list of DECam images
  ;expstr = mrdfits('/dl1/users/dnidever/nsc/instcal/v3/lists/decam_instcal_list.fits.gz',1)
  expstr = mrdfits('/dl1/users/dnidever/nsc/instcal/v3/lists/decam_instcal_list_full.fits.gz',1)
  expstr.expnum = strtrim(expstr.expnum,2)
  expstr.plver = strtrim(expstr.plver,2)
  expstr.fluxfile = strtrim(expstr.fluxfile,2)
  expstr.maskfile = strtrim(expstr.maskfile,2)
  expstr.wtfile = strtrim(expstr.wtfile,2)
  expstr_expnumplver = expstr.expnum+'-'+expstr.plver

  ;; Load the DECam extension name to ccdnum conversion file
  decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)

;; Night loop
for n=0,nnights-1 do begin
  inight = nights[n]

  CD,smashdir+inight
  CD,current=nightdir
  nightdir = trailingslash(nightdir)

  print,''
  print,'Making SMASH symlinks for night ',inight
  print,'-----------------------------------------'

  ;; Make sure the DELVE night directory exists
  ;if FILE_TEST(delvedir+inight,/directory) eq 0 then FILE_MKDIR,delvedir+inight

  ;; Fix the DELVE log files
  logfiles = file_search(delvedir+inight+'/logs/*',count=nlogfiles)
  for j=0,nlogfiles-1 do begin
    print,'Fixing ',logfiles[j]
    readline,logfiles[j],lines
    lines2 = lines
    b = where(stregex(lines,'/chip?[1-9]/',/boolean) eq 1,nb)
    if nb gt 0 then begin
      for k=0,nb-1 do begin
        line1 = lines[b[k]]
        dum = strsplit(line1,'/',/extract)
        schip1 = dum[1]
        schip2 = 'chip0'+strmid(schip1,4,1)
        dum2 = dum
        dum2[1] = schip2
        line2 = strjoin(dum2,'/')
        lines2[b[k]] = line2
      endfor
      writeline,logfiles[j],lines2
    endif
  endfor


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
    print,ifield
    fitsfiles1 = fitsfiles[ind]
    MATCH,fieldstr.shname,ifield,fieldind,count=nmatch
    if nmatch eq 0 then begin
      print,ifield,' NOT FOUND in >>fields<< file'
      goto,FIELDBOMB
    endif
    fieldname = fieldstr[fieldind].name   ;; long field name

    ;; Make sure the DELVE night+field directory exists
    ;if FILE_TEST(delvedir+inight+'/'+ifield,/directory) eq 0 then FILE_MKDIR,delvedir+inight+'/'+ifield

    ;; Check the number of chip directories
    chdir = file_search(delvedir+inight+'/'+ifield+'/chip*',/test_directory,count=nchdir)
    if nchdir ge 55 then goto,FIELDBOMB


    ;; Gentral position of field
    sumfile = nightdir+fieldname+'_summary.fits'
    sumstr = MRDFITS(sumfile,1,/silent)
    cenra = median([sumstr.ra])
    cendec = median([sumstr.dec])

    ;; Get Gaia DR2 and other reference data for this field
    if FILE_TEST(delvedir+inight+'/refcat/',/directory) eq 0 then FILE_MKDIR,delvedir+inight+'/refcat/'
    savefile = delvedir+inight+'/refcat/'+ifield+'_refcat.fits'
    ;if (file_test(savefile) eq 0 and file_test(savefile+'.gz') eq 0) or keyword_set(redo) then begin
    ;  refcat = DELVERED_GETREFDATA(['c4d-u','c4d-g','c4d-r','c4d-i','c4d-z','c4d-Y','c4d-VR'],cenra,cendec,1.5,savefile=savefile)
    ;  SPAWN,['gzip','-f',savefile],/noshell
    ;endif
    refcat = MRDFITS(savefile+'.gz',1,/silent)

    ;;; Match FITS files with the original CP c4d files
    ;fbase = PHOTRED_GETFITSEXT(fitsfiles1,/basename)
    ;dum = strsplitter(fbase,'-',/extract)
    ;fexpnum = strmid(reform(dum[1,*]),0,8)
    ;fui = uniq(fexpnum,sort(fexpnum))
    ;fexpnum = fexpnum[fui]
    ;nexpnum = n_elements(fexpnum)
    ;plver = strarr(nexpnum)
    ;for k=0,nexpnum-1 do plver[k]=sxpar(PHOTRED_READFILE(fitsfiles1[fui[k]],exten=1,/header),'plver')
    ;plver = strtrim(plver,2)
    ;MATCH,expstr_expnumplver,fexpnum+'-'+plver,ind1,ind2,/sort,count=nmatch
    ;if nmatch ne nexpnum then stop,'Not all exposures have matches in DECam master list'
    ;fexpnum = fexpnum[ind2]
    ;expstr1 = expstr[ind1]

    ;print,strtrim(f+1,2),'/',strtrim(nfields,2),' ',ifield,' ',fieldname,' ',strtrim(nexpnum,2),' exposures'

    ;; Use chip10/ files as a template
    files10 = file_search(delvedir+inight+'/'+ifield+'/chip10/F*',count=nfiles10)
    print,strtrim(nfiles10,2),' files'

    ;; chip loop
    for c=1,9 do begin
      schip = string(c,format='(i02)')
      ;; file loop
      for k=0,nfiles10-1 do begin 
        file1 = files10[k]
        base = file1
        ext = first_el(strsplit(base,'.',/extract),/last)
        if stregex(base,'.fits',/boolean) eq 0 then begin
          ;; /dl1/users/dnidever/delve/exposures/20141218/F13/chip10/F13-00389253_10.als
          file2 = repstr(file1,'/chip10/','/chip'+schip+'/')
          file2 = repstr(file2,'_10','_'+schip)
          lfile1 = file_readlink(file1)
          ;; /dl1/users/dnidever/smash/cp/red/photred/20141218/F13/F13-00389253_10.als
          lfile2 = repstr(lfile1,'_10','_'+schip)
        endif else begin
        ;; FITS file

          ;; Resource file
          stop
        endelse

stop
      endfor
    endfor

    ;; Loop over exposures for this field
    For e=0,nexpnum-1 do begin
      fexpnum1 = fexpnum[e]
      print,'  ',strtrim(e+1,2),' ',fexpnum1
      fluxfile = repstr(expstr1[e].fluxfile,'/net/mss1/','/mss1/')
      maskfile = repstr(expstr1[e].maskfile,'/net/mss1/','/mss1/')
      wtfile = repstr(expstr1[e].wtfile,'/net/mss1/','/mss1/')
      ;; Get number of extensions
      ;;   use symlink to make fits_open think it's a normal FITS file
      tmpfluxfile = MKTEMP('tmp',/nodot,outdir=workdir) & TOUCHZERO,tmpfluxfile+'.fits' & FILE_DELETE,[tmpfluxfile,tmpfluxfile+'.fits'],/allow
      tmpfluxfile += '.fits'
      FILE_LINK,fluxfile,tmpfluxfile
      FITS_OPEN,tmpfluxfile,fcb & FITS_CLOSE,fcb

      ;; Get the CCDNUM for the extensions
      MATCH,decam.name,fcb.extname,ind1,ind2,/sort,count=nmatch
      extnum = ind2
      ccdnum = decam[ind1].ccdnum
      nccdnum = nmatch

      ;; Chip loop
      For c=0,nccdnum-1 do begin
        chipnum1 = ccdnum[c]
stop
        extnum1 = extnum[c]  ; extension for the CP files
        chbase1 = ifield+'-'+fexpnum1+'_'+string(chipnum1,format='(i02)')
        ;; _cat.dat, _refcat.dat
        ;; als, ap, coo, plst, fits/fits.fz, psf
        chipdir1 = delvedir+inight+'/'+ifield+'/chip'+string(chipnum1,format='(i02)')+'/'
        if file_test(chipdir1,/directory) eq 0 then FILE_MKDIR,chipdir1
        ;; Check that files exist
        if file_test(nightdir+ifield+'/'+chbase1+'.opt') eq 1 then begin
          FILE_DELETE,chipdir1+chbase1+['_cat.dat','.opt','.als.opt','.als','.ap','.coo','.plst','.psf','.psf.log','.log','a.als','a.ap'],/allow
          FILE_LINK,nightdir+ifield+'/'+chbase1+['_cat.dat','.opt','.als.opt','.als','.ap','.coo','.plst','.psf','.psf.log','.log','a.als','a.ap'],chipdir1

          ;; Create FITS resource file 
          FILE_DELETE,chipdir1+chbase1+['.fits','.fits.fz'],/allow
          outfile1 = chipdir1+chbase1+'.fits'
          WRITELINE,outfile1,''
          routfile1 = chipdir1+'.'+chbase1+'.fits'
          rlines = ['fluxfile = '+strtrim(fluxfile,2)+'['+strtrim(extnum1,2)+']',$
                    'wtfile = '+strtrim(wtfile,2)+'['+strtrim(extnum1,2)+']',$
                    'maskfile = '+strtrim(maskfile,2)+'['+strtrim(extnum1,2)+']']
          WRITELINE,routfile1,rlines
          
          ;; Save the reference catalog for this chip
          hd = HEADFITS(tmpfluxfile,exten=extnum1,/silent)
          ;; Temporarily fix NAXIS1/2 values
          sxaddpar,hd,'NAXIS1',sxpar(hd,'ZNAXIS1')
          sxaddpar,hd,'NAXIS2',sxpar(hd,'ZNAXIS2')
          nx = sxpar(hd,'naxis1')
          ny = sxpar(hd,'naxis2')
          HEAD_XYAD,hd,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
          gdrefcat = where(refcat.ra ge min(vra)-0.02 and refcat.ra le max(vra)+0.02 and $
                           refcat.dec ge min(vdec)-0.02 and refcat.dec le max(vdec)+0.02,ngdrefcat)
          refcat1 = refcat[gdrefcat]
          FILE_DELETE,chipdir1+chbase1+'_refcat.dat',/allow
          refcatfile = chipdir1+chbase1+'_refcat.fits'
          MWRFITS,refcat1,refcatfile,/create
          SPAWN,['gzip','-f',refcatfile],/noshell

          ;; Update the lists, relative paths
          wcslines[count] = chipdir1+'/'+chbase1+'.fits'
          daophotlines[count] = chipdir1+'/'+chbase1+'.fits'
          count++
        endif else print,chbase+' NOT FOUND'
        CHIPBOMB:
      Endfor  ; chip loop
      ;; Delete the temporary file link
      FILE_DELETE,tmpfluxfile
    Endfor  ; exposure loop

stop

    ;; MATCH files, mch, raw
    mchfiles = FILE_SEARCH(nightdir+ifield+'/'+ifield+'-????????_??.mch',count=nmchfiles)
    if nmchfiles gt 0 then begin
      rawfiles = FILE_DIRNAME(mchfiles)+'/'+FILE_BASENAME(mchfiles,'.mch')+'.raw'
      tfrfiles = FILE_DIRNAME(mchfiles)+'/'+FILE_BASENAME(mchfiles,'.mch')+'.tfr'
      ;; Get chip subdirectories
      chipdirs = strarr(nmchfiles)
      for k=0,nmchfiles-1 do chipdirs[k]='chip'+string(PHOTRED_GETCHIPNUM(file_basename(mchfiles[k],'.mch'),{namps:62,separator:'_'}),format='(i02)')
      FILE_DELETE,delvedir+inight+'/'+ifield+'/'+chipdirs+'/'+file_basename(mchfiles),/allow
      FILE_DELETE,delvedir+inight+'/'+ifield+'/'+chipdirs+'/'+file_basename(rawfiles),/allow
      FILE_DELETE,delvedir+inight+'/'+ifield+'/'+chipdirs+'/'+file_basename(tfrfiles),/allow
      FILE_LINK,mchfiles,delvedir+inight+'/'+ifield+'/'+chipdirs+'/'+file_basename(mchfiles)
      FILE_LINK,rawfiles,delvedir+inight+'/'+ifield+'/'+chipdirs+'/'+file_basename(rawfiles)
      FILE_LINK,tfrfiles,delvedir+inight+'/'+ifield+'/'+chipdirs+'/'+file_basename(tfrfiles)

      ;; Update the list, relative paths
      PUSH,matchlines,ifield+'/'+chipdirs+'/'+file_basename(mchfiles)
    endif


    FIELDBOMB:
  Endfor  ; field loop

  ;; Copy apcor.lst and APCOR.success
  FILE_COPY,nightdir+'apcor.lst',delvedir+inight,/over
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

  ;;; Copy over the WCS.success, DAOPHOT.success and
  ;;; MATCH.success/outlist files
  ;if FILE_TEST(delvedir+inight+'/logs',/directory) eq 0 then FILE_MKDIR,delvedir+inight+'/logs'
  ;;;WRITELINE,delvedir+inight+'/logs/WCS.success',wcslines
  ;WRITELINE,delvedir+inight+'/logs/WCS.inlist',wcslines   ;; we want to redo WCS with GaiaDR2
  ;WRITELINE,delvedir+inight+'/logs/DAOPHOT.success',daophotlines
  ;WRITELINE,delvedir+inight+'/logs/MATCH.outlist',matchlines
  ;WRITELINE,delvedir+inight+'/logs/APCOR.success',apcorlines

  ;;; Copy the "fields" file
  ;FILE_COPY,nightdir+'fields',delvedir+inight,/over,/allow
  ;FILE_CHMOD,delvedir+inight+'/fields',/a_write  

  ;; Make the setup file
  ;WRITELINE,delvedir+inight+'/photred.setup',setup

endfor  ; night loop

;stop


end
