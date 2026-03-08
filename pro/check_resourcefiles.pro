pro check_resourcefiles,nights

;; There are files in the mass store that don't exist anymore
;; fix resource files and replace with new versions

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

undefine,fixablefiles
undefine,missingfiles
undefine,missingexp

allexp = mrdfits('/net/dl2/dnidever/delve/exposures/all_exposures.fits',1)
expdt = allexp[0]
struct_assign,{dum:''},expdt

;; Night loop
For n=0,nnights-1 do begin
  inight = nights[n]
  ;;print,strtrim(n+1,2),' ',inight

  ;; Get the exposures file
  expfile = delvedir+'exposures/'+inight+'/'+inight+'_exposures.fits'
  if file_test(expfile) eq 0 then begin
    print,expfile,' NOT FOUND'
    goto,NIGHTBOMB
  endif
  expstr = mrdfits(expfile,1,/silent)
  expstr.expnum = strtrim(expstr.expnum,2)
  expstr.fluxfile = strtrim(expstr.fluxfile,2)
  expstr.maskfile = strtrim(expstr.maskfile,2)
  expstr.wtfile = strtrim(expstr.wtfile,2)

  undefine,missingfiles1

  fdir = file_search(delvedir+'exposures/'+inight+'/F*',/test_directory,count=nfdir)

  ;; Loop over the fields
  For f=0,nfdir-1 do begin
    ifield = file_basename(fdir[f])

    expfiles = file_search(delvedir+'exposures/'+inight+'/'+ifield+'/chip01/'+ifield+'-????????_01.fits',count=nexpfiles)
    expnum = file_basename(expfiles,'_01.fits')
    len = strlen(ifield)
    expnum = strmid(expnum,len+1)
    nexpnum = n_elements(expnum)

    ;; Loop over exposures
    For e=0,nexpnum-1 do begin
      iexpnum = expnum[e]

      ;; Get the associated mass store file
      MATCH,expstr.expnum,iexpnum,ind1,ind2,/sort,count=nmatch
      ;;if nmatch eq 0 then stop,'no match to exposure list'
      if nmatch eq 0 then print,'no match to exposure list'
      if nmatch eq 0 then continue
      expstr1 = expstr[ind1]
      push,allexposures1,expstr1.fluxfile

      ;; Check if the file exists in the mass store
      ;;if file_test(expstr1.fluxfile) and file_test(expstr1.maskfile) and file_test(expstr1.wtfile) then continue
      ;;if file_test(expstr1.fluxfile) then continue
      base = file_basename(expstr1.fluxfile,'.fits.fz')
      basehead = strmid(base,0,17)

      ;; Loop over the chips
      For c=0,61 do begin
        ccdnum = c+1

        ;if ccdnum eq 61 then goto,CHIPBOMB
        schip = string(ccdnum,format='(i02)')
        ;chipdir = fielddir+'chip'+schip+'/'
        ind = where(decam.ccdnum eq ccdnum,nind)
        extname = decam[ind[0]].name

        file = delvedir+'exposures/'+inight+'/'+ifield+'/chip'+schip+'/'+ifield+'-'+string(iexpnum,format='(i08)')+'_'+schip+'.fits'
        rfile = delvedir+'exposures/'+inight+'/'+ifield+'/chip'+schip+'/.'+ifield+'-'+string(iexpnum,format='(i08)')+'_'+schip+'.fits'

        if file_test(rfile) eq 0 then begin
          if ccdnum ne 2 and ccdnum ne 61 then print,file,' NOT FOUND'
          goto,CHIPBOMB
        endif


        ;; Load resource file
        READLINE,rfile,rlines
        newrlines = rlines
        ;; Remove blank spaces
        arr = strsplitter(rlines,' ',/extract)
        sz = size(arr)
        tags = ['fluxfile','wtfile','maskfile']
        ;;for k=0,n_elements(tags)-1 do begin
        ftind = where(arr[0,*] eq 'fluxfile',nftind)
        if nftind eq 0 then stop,'cannot find '+tags[k]
        farr = reform(arr[*,ftind])
        if n_elements(farr) gt 3 then stop,'more than 3 elements'
        hi = strpos(farr[2],'[')
        fext = strmid(farr[2],hi)
        oldfluxfile = strmid(farr[2],0,hi)
        base = file_basename(oldfluxfile,'.fits.fz')
        basehead = strmid(base,0,17)
        testoldfluxfile = file_test(oldfluxfile)
        ;;if file_test(oldfluxfile) eq 0 then print,oldfluxfile,' not found'
        wtind = where(arr[0,*] eq 'wtfile',nwtind)
        warr = reform(arr[*,wtind])
        hi = strpos(warr[2],'[')
        oldwtfile = strmid(warr[2],0,hi)
        testoldwtfile = file_test(oldwtfile)
        mtind = where(arr[0,*] eq 'maskfile',nmtind)
        marr = reform(arr[*,mtind])
        hi = strpos(marr[2],'[')
        oldmaskfile = strmid(marr[2],0,hi)
        testoldmaskfile = file_test(oldmaskfile)

        ;; All the files exist in the mass store, SKIP
        if testoldfluxfile eq 1 and testoldwtfile eq 1 and testoldmaskfile eq 1 then continue

        push,missingfiles1,oldfluxfile

        CHIPBOMB:
     Endfor  ; chip loop      


    Endfor   ; exposure loop
    BOMB:
  Endfor ; field loop

  print,strtrim(n+1,2),' ',inight,n_elements(allexposures1),n_elements(missingfiles1)
  if n_elements(missingfiles1) gt 0 then push,missingfiles,missingfiles1

  ;;if nmissing gt 0 then stop

  nightbomb:
Endfor ; night loop

print,strtrim(n_elements(missingfiles),2),' missing files'

stop

end
