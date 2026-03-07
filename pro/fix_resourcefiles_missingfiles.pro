pro fix_resourcefiles_missingfiles,nights

;; There are files in the mass store that don't exist anymore
;; fix resource files and replace with new versions

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl2/dnidever/delve/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)

;; new mass store filenames
newstr = mrdfits('/net/dl2/dnidever/nsc/instcal/v4/lists/decam_instcal_01302026.fits',1)
add_tag,newstr,'basehead','',newstr
for i=0,n_elements(newstr)-1 do newstr[i].basehead=strmid(newstr[i].base,0,17)

;; Frank's reprocessed images
frankfiles = importascii('/home/dnidever/delve_archive.txt',fieldnames=['filename'])
frankfiles.filename = repstr(frankfiles.filename,'/net/archive/','/net/mss1/archive/')
add_tag,frankfiles,'base','',frankfiles
frankfiles.base = file_basename(frankfiles.filename)
add_tag,frankfiles,'basehead','',frankfiles
for i=0,n_elements(frankfiles)-1 do frankfiles[i].basehead = strmid(frankfiles[i].base,0,17)
gd = where(stregex(frankfiles.filename,'ooi',/boolean) eq 1,ngd)
frankfiles = frankfiles[gd]

;; Frank reprocessed these as well
frankfiles2 = importascii('/home/dnidever/stillmissing_030326.txt',fieldnames=['basehead','expnum'])
add_tag,frankfiles2,'filename','',frankfiles2
daysofmonth = [31,28,31,30,31,30,31,31,30,31,30,31]
for i=0,n_elements(frankfiles2)-1 do begin
  base = frankfiles2[i].basehead
  night = '20'+strmid(base,4,6)
  syear = strmid(night,0,4)
  smonth = strmid(night,4,2)
  sday = strmid(night,6,2)
  year = long(syear)
  month = long(smonth)
  day = long(sday)
  newfilename = '/net/mss1/archive/pipe/'+night+'/ct4m/2012B-0001/'+base+'_ooi_*_delve.fits.fz'
  newfiles = file_search(newfilename,count=nnewfiles)
  if nnewfiles gt 0 then begin
    frankfiles2[i].filename = newfiles[0]
    continue
  endif
  night2 = syear+smonth+string(day+1,format='(I02)')
  newfilename2 = '/net/mss1/archive/pipe/'+night2+'/ct4m/2012B-0001/'+base+'_ooi_*_delve.fits.fz'
  newfiles2 = file_search(newfilename2,count=nnewfiles2)
  if nnewfiles2 gt 0 then begin
    frankfiles2[i].filename = newfiles2[0]
    continue
  endif
  night3 = syear+smonth+string(day-1,format='(I02)')
  if day-1 le 0 then begin
    night3 = syear+string(month-1,format='(I02)')+string(daysofmonth[month-2],format='(I02)')
  endif
  if day-1 le 0 and month-1 le 0 then begin
    night3 = string(year-1,format='(I04)')+'1231'
  endif
  newfilename3 = '/net/mss1/archive/pipe/'+night3+'/ct4m/2012B-0001/'+base+'_ooi_*_delve.fits.fz'
  newfiles3 = file_search(newfilename3,count=nnewfiles3)
  if nnewfiles3 gt 0 then begin
    frankfiles2[i].filename = newfiles3[0]
    continue
  endif
  newfilename4 = '/net/mss1/archive/pipe/'+night+'/ct4m/*/'+base+'_ooi_*_delve.fits.fz'
  newfiles4 = file_search(newfilename4,count=nnewfiles4)
  if nnewfiles4 gt 0 then begin
    frankfiles2[i].filename = newfiles4[0]
    continue
  endif
  newfilename5 = '/net/mss1/archive/pipe/'+night2+'/ct4m/*/'+base+'_ooi_*_delve.fits.fz'
  newfiles5 = file_search(newfilename5,count=nnewfiles5)
  if nnewfiles5 gt 0 then begin
    frankfiles2[i].filename = newfiles5[0]
    continue
  endif
  newfilename6 = '/net/mss1/archive/pipe/'+night3+'/ct4m/*/'+base+'_ooi_*_delve.fits.fz'
  newfiles6 = file_search(newfilename6,count=nnewfiles6)
  if nnewfiles6 gt 0 then begin
    frankfiles2[i].filename = newfiles6[0]
    continue
  endif
  print,'no new filename for '+base
  stop
endfor
add_tag,frankfiles2,'base','',frankfiles2
frankfiles2.base = file_basename(frankfiles2.filename)
;; now combine the two frank lists
temp = replicate({filename:'',base:'',basehead:''},n_elements(frankfiles2))
struct_assign,frankfiles2,temp
push,frankfiles,temp


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

  undefine,allexposures1
  undefine,fixablefiles1
  undefine,missingfiles1

  ;; Check that all of these files exist in the mass store
  exists = file_test(expstr.fluxfile)
  ;;mexists = file_test(expstr.maskfile)
  ;;wexists = file_test(expstr.wtfile)
  missing = where(exists eq 0,nmissing)
  if nmissing eq 0 then begin
    print,strtrim(n+1,2),' ',inight,n_elements(expstr),nmissing
    continue
  endif
  ;;if nmissing gt 0 then begin
  ;;  push,missingfiles,expstr[missing].fluxfile
  ;;  ;;push,allmissing,expstr[missing]
  ;;endif else begin
  ;;  continue
  ;;endelse

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
      if file_test(expstr1.fluxfile) then continue
      base = file_basename(expstr1.fluxfile,'.fits.fz')
      basehead = strmid(base,0,17)


      ;; Get new frank files
      find = where(frankfiles.basehead eq basehead,nfind)
      if nfind gt 0 then push,fixablefiles,base
      ;;if nfind gt 0 then continue

      ;; Get the new filename
      ind = where(newstr.basehead eq basehead,nind)
      if nind gt 0 then push,fixablefiles,base
      ;;if nind gt 0 then continue

      ;; missing
      push,missingfiles,base
      temp = replicate(expdt,1)
      struct_assign,expstr1,temp
      if n_elements(missingexp) eq 0 then missingexp=temp else push,missingexp,temp
      push,missingfiles1,base
      ;;continue

;;stop

      ;;if nind eq 0 then push,missingfiles1,base
      ;;if nind gt 0 then push,fixablefiles1,base
      ;;continue

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


        ;; Get the new filename
        ;;ind = where(newstr.basehead eq basehead,nind)
        frankind = where(frankfiles.basehead eq basehead,nfrankind)
        if nfrankind gt 0 then begin
          newfluxfile = frankfiles[frankind].filename
        endif else begin
          newind = where(newstr.basehead eq basehead,nnewind)
          if nnewind eq 0 then stop,'no new filename'
          newfluxfile = newstr[newind].filename
        endelse

        if c eq 0 then print,newfluxfile

        ;continue

        if file_test(newfluxfile) eq 0 then stop,newfluxfile+' NOT FOUND'
        newrlines[ftind] = 'fluxfile = '+newfluxfile+fext
        newwtfile = repstr(newfluxfile,'ooi','oow')      
        if file_test(newwtfile) eq 0 then stop,newwtfile,' NOT FOUND'
        newrlines[wtind] = 'wtfile = '+newwtfile+fext
        newmaskfile = repstr(newfluxfile,'ooi','ood')
        if file_test(newmaskfile) eq 0 then stop,newmaskfile,' NOT FOUND'
        newrlines[mtind] = 'maskfile = '+newmaskfile+fext
        ;; Move old resource to backup
        caldat,systime(/julian),mt,dy,yr,hr,mn,sc
        datestamp = string(yr,mt,dy,hr,mn,sc,format='(I04,I02,I02,I02,I02,I02)')
        backuprfile = rfile+'.bak'+datestamp
        FILE_MOVE,rfile,backuprfile
        ;; Write new resource file
        WRITELINE,rfile,newrlines

        ;;stop


        CHIPBOMB:
     Endfor  ; chip loop      


    Endfor   ; exposure loop
    BOMB:
  Endfor ; field loop

  print,strtrim(n+1,2),' ',inight,n_elements(allexposures1),n_elements(fixablefiles1),n_elements(missingfiles1)
  if n_elements(missingfiles1) gt 0 then push,missingfiles,missingfiles1

  ;;if nmissing gt 0 then stop

  nightbomb:
Endfor ; night loop


stop

end
