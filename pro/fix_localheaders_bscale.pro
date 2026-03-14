pro fix_localheaders_bscale,nights

;; fix local headers for _d2 exposures (gain units) that were replaced
;; with CP versions (ADU units)
;; add BSCALE of the gain to the local headers

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl2/dnidever/delve/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)

;; go through the nights
nnights = n_elements(nights)
if nnights eq 0 then begin
  nights = file_search(delvedir+'exposures/20??????',/test_directory,count=nnights)
  nights = file_basename(nights)
endif
print,strtrim(n_elements(nights),2),' nights'

undefine,oldfiles

;; Night loop
For n=0,nnights-1 do begin
  inight = nights[n]

  ;; check all local header files and see if bunit='electrons'

  files = file_search(delvedir+'exposures/'+inight+'/F*/chip01/F*-*_01.fits.head',count=nfiles)
  print,strtrim(n+1,2),' ',strtrim(inight,2),' ',strtrim(nfiles,2),' files'

  undefine,oldfiles1
  for i=0,nfiles-1 do begin
    readline,files[i],hlines
    bunit = sxpar(hlines,'bunit',count=nbunit)
    if nbunit gt 0 and strlowcase(bunit) eq 'electrons' then begin
      push,oldfiles1,files[i]

      ;; Loop over all of the chips for this exposure
      dir = file_dirname(file_dirname(files[i]))
      base = file_basename(files[i],'.fits.head')
      base = strmid(base,0,strlen(base)-3)
      print,'  Fixing ',base
      for j=1,62 do begin
        chip = j
        schip = string(j,format='(I02)')
        ;; Loading resource file and checking large fluxfile name
        rchipfile = dir+'/chip'+schip+'/.'+base+'_'+schip+'.fits'
        if j eq 1 then print,'    ',rchipfile
        if file_test(rchipfile) eq 0 then continue
        readline,rchipfile,rlines
        ind = where(stregex(rlines,'fluxfile = ',/boolean) eq 1,nind)
        fluxfile = (strsplit(rlines[ind[0]],' ',/extract))[2]
        lo = strpos(fluxfile,'[')
        hi = strpos(fluxfile,']')
        fext = strmid(fluxfile,lo+1,hi-lo-1)
        fluxfile = strmid(fluxfile,0,lo)
        fluxbase = file_basename(fluxfile,'.fits.fz')
        if j eq 1 then print,'    ',fluxfile
        if strmid(fluxbase,0,3) eq 'c4d' then begin
          fluxtag = first_el(strsplit(fluxbase,'_',/extract),/last)
          if fluxtag eq 'd1' or fluxtag eq 'd2' then stop,'still pointing to old DES file'
        endif else begin
          stop,'check if this is the old des version'
        endelse
        ;; Loading local header file
        chipfile = dir+'/chip'+schip+'/'+base+'_'+schip+'.fits.head'
        if file_test(chipfile) eq 0 then continue
        readline,chipfile,hlines

        ;; Load the actual data
        ;;   and measure the new median of the CP file
        tmpdir = MKTEMP('rsrc',/directory,/nodot)
        FILE_CHMOD,tmpdir,/a_execute
        tfluxfile = tmpdir+'/flux.fits'
        FILE_WAIT,fluxfile
        SPAWN,['funpack','-E',fext,'-O',tfluxfile,fluxfile],/noshell
        FITS_READ,tfluxfile,newim1,newhead1
        newmed = median(newim1)
        file_delete,[tfluxfile,tmpdir],/allow

        ;; old des sky value
        oldsky = sxpar(hlines,'skybrite',count=nskybrite)
        if nskybrite eq 0 then stop,'no skybrite found in '+chipfile
        bscale = oldsky/newmed

        ;;dum = sxpar(hlines,'bscale',count=nbscale)
        ;;if nbscale gt 0 then begin
        ;;  if j eq 1 then print,'    already have bscale in '+chipfile
        ;;  continue
        ;;endif
        ;;gain = sxpar(hlines,'arawgain',count=ngain)
        ;;if ngain eq 0 then stop,'no gain found in '+chipfile

        sxaddpar,hlines,'bscale',bscale,' scale to des-like electrons'
        ;; save a copy of the original one
        jd = systime(/julian)
        caldat,jd,month,day,year,hour,minute,second
        ;;timestamp = string(year,month,day,hour,minute,second,format='(I04,I02,I02,I02,I02,I02)')
        ;;bakchipfile = chipfile+'.bak'+timestamp
        ;;file_move,chipfile,bakchipfile,/allow,/over
        writeline,chipfile,hlines
      endfor  ; chip loop
    endif  ; electrons
  endfor  ; exposure loop

  print,' '
  print,strtrim(n,2),' ',inight,' ',strtrim(n_elements(oldfiles1),2)
  print,' '
  push,oldfiles,oldfiles1

  nightbomb:
Endfor ; night loop


stop

end
