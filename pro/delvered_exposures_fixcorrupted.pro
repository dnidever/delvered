pro delvered_exposures_deletefiles,nightdir
  ;; delete files

  CD,current=origdir
  CD,nightdir  

  ;; Starting fresh
  print,'Starting fresh.  Deleting old files'
  ;; Delete old files
  ;;--------------------
  ;; zero-out fits files
  fitsfiles = file_search('F*/chip*/F*_??.fits',count=nfitsfiles)
  if nfitsfiles gt 0 then begin
    print,'Deleting image-level files'
    touchzero,fitsfiles
    base = repstr(fitsfiles,'.fits','')
    ext = ['_cat.dat','.fits.head','.opt','.als.opt','.coo','.ap','.grp','.lst','.cmn.lst',$
           '.psfini.ap','.lst1','.lst1.chi','.lst2','.lst2.chi','.plst','.psf','.psf.log',$
           '.als','.als.inp','.log','.nei','.nst','a.ap','a.als','a.als.inp','a.log',$
           's.fits.fz','.mch','.tfr','.raw','.ast','.input','.phot','.ast']
    for j=0,n_elements(ext)-1 do file_delete,base+ext[j],/allow
    ;; remove HEADER from resource files
    print,'Updating resource files'
    for j=0,nfitsfiles-1 do begin
      dir1 = file_dirname(fitsfiles[j])
      base1 = file_basename(fitsfiles[j])
      rfile1 = dir1+'/.'+base1
      if file_test(rfile1) eq 1 then begin
        READLINE,rfile1,rlines
        bd = where(stregex(rlines,'HEADER',/boolean) eq 1,nbd)
        if nbd gt 0 then begin
          REMOVE,bd,rlines
          WRITELINE,rfile1,rlines
        endif
      endif
    endfor

    ;; chip directory-level files
    print,'Deleting chip directory-level files'
    fitsdirs = file_dirname(fitsfiles)
    udirs = fitsdirs[uniq(fitsdirs,sort(fitsdirs))]
    ;; scripts: daophot.sh, lstfilter.py, srcfilter.pro, goodpsf.pro, apcor.opt
    for j=0,n_elements(udirs)-1 do file_delete,udirs[j]+'/'+['daophot.sh','lstfilter.py','srcfilter.pro','goodpsf.pro','photo.opt','apcor.opt'],/allow
    ;; batch files, wfit, dopt, dao, initpsf, match,
    for j=0,n_elements(udirs)-1 do begin
      bfiles = file_search(udirs[j]+'/'+['wfit*.batch*','dopt*.batch*','dao*.sh*','initpsf*batch*','match*batch*'],count=nbfiles)
      if nbfiles gt 0 then file_delete,bfiles,/allow
    endfor
    ;; filters
    file_delete,udirs+'/filters',/allow
    ;; dangling rsrcXXXXXX/ directories
    for j=0,n_elements(udirs)-1 do begin
      rdirs = file_search(udirs[j]+'/rsrc*',/test_directory,count=nrdirs)
      if nrdirs gt 0 then begin
        rfiles = file_search(rdirs+'/*',count=nrfiles)
        if nrfiles gt 0 then file_delete,rfiles,/allow
        file_delete,rdirs,/allow
      endif
    endfor

    ;; field-level files: cmb, dered, final
    fdirs = file_dirname(udirs)  ; field directories
    fdirs = fdirs[uniq(fdirs,sort(fdirs))]
    for j=0,n_elements(fdirs)-1 do begin
      ffiles = file_search(fdirs[j]+'/'+['*.cmb','*.dered','*.final','*.dat','*.fits.gz','*_summary.fits','filters'],count=nffiles)
      if nffiles gt 0 then file_delete,ffiles,/allow
    endfor
  endif

  ;; daogrow files
  dfiles = file_search('daogrow*/*',count=ndfiles)
  if ndfiles gt 0 then begin
    print,'Deleting daogrow files'
    file_delete,dfiles,/allow
    ddirs = file_dirname(dfiles)
    ddirs = ddirs[uniq(ddirs,sort(ddirs))]
    file_delete,ddirs,/allow
  endif

  ;; summary files, delvered night summary file
  sfiles = file_search(['apcor*','delve.trans','*.dat','*.final','*.fits.gz','*_summary.fits','extinction','filters','idlbatch','runbatch'],count=nsfiles)
  if nsfiles gt 0 then file_delete,sfiles,/allow

  ;; zero-out lists and log files
  lfiles = file_search('logs/*',count=nlfiles)
  if nlfiles gt 0 then file_delete,lfiles,/allow
  ;; Put all fits files in WCS.inlist
  WRITELINE,'logs/WCS.inlist',fitsfiles

  CD,origdir

end


pro delvered_exposures_fixfiles,corruptedstr,fixedfiles
  ;; Fix the resource files for the corrupted exposures

  ;; Chip name and CCDNUM relations
  decam = IMPORTASCII('/home/dnidever/projects/delvered/data/decam.txt',/header,/silent)

  undefine,fixedfiles
  
  ;; Exposure loop
  for i=0,n_elements(corruptedstr)-1 do begin
    expnum = corruptedstr[i].expnum
    fluxfile = corruptedstr[i].filename
    wtfile = repstr(fluxfile,'ooi','oow')
    maskfile = repstr(fluxfile,'ooi','ood')
    if file_test(fluxfile) eq 0 then stop,fluxfile,' not found'
    if file_test(wtfile) eq 0 then stop,wtfile,' not found'
    if file_test(maskfile) eq 0 then stop,maskfile,' not found'
    chip1file = file_search('F*/chip01/F*-'+expnum+'_01.fits',count=nchip1file)
    field = (strsplit(chip1file,'/',/extract))[0]
    ;; Chip loop
    for j=1,62 do begin
      chip = string(j,format='(I02)')
      chind = where(decam.ccdnum eq j)
      extname = decam[chind[0]].name
      chipfile = field+'/chip'+chip+'/'+field+'-'+expnum+'_'+chip+'.fits'
      if file_test(chipfile) eq 0 then continue
      resourcefile = file_dirname(chipfile)+'/.'+file_basename(chipfile)
      print,'fixing ',resourcefile
      if file_test(resourcefile) eq 1 then begin
        readline,resourcefile,lines
        newlines = lines
        find = where(stregex(lines,'^fluxfile',/boolean) eq 1,nfind)
        newlines[find[0]] = 'fluxfile = '+fluxfile+'['+extname+']'
        wind = where(stregex(lines,'^wtfile',/boolean) eq 1,nfind)
        newlines[wind[0]] = 'wtfile = '+wtfile+'['+extname+']'
        mind = where(stregex(lines,'^maskfile',/boolean) eq 1,nfind)
        newlines[mind[0]] = 'maskfile = '+maskfile+'['+extname+']'
        file_move,resourcefile,resourcefile+'.orig',/over
        writeline,resourcefile,newlines
      endif else begin
        ;; make resource file fresh
        print,'  making fresh'
        newlines = 'fluxfile = '+fluxfile+'['+extname+']'
        push,newlines,'wtfile = '+wtfile+'['+extname+']'
        push,newlines,'maskfile = '+maskfile+'['+extname+']'
        writeline,resourcefile,newlines
      endelse

      push,fixedfiles,chipfile

    endfor  ;; chip loop
  endfor  ;; exposure loop

end


;+
;
; DELVERED_EXPOSURES_FIXCORRUPTED
;
; This is a script to run PHOTRED on DELVE MC exposures to create
; single-image level catalogs and PSFs.
;
; FIX corrupted images and rerun them.
;
; INPUTS:
;  input         What nights to run PHOTRED on.  Either an array or
;                  a range such as 20160101-20160506.
;  /uselocal     Do all processing in a local drive.  Off by default.
;  /startfresh   Start completely fresh.  Remove old files and start
;                  from the very beginning.
;  =doapcor      Control whether APCOR is run.
;  =dozeropoint  Control whether ZEROPOINT is run.
;
; OUTPUTS:
;  PHOTRED will be run on each night and a final summary file
;  will be created for each night.
;
; USAGE:
;  IDL>delvered_exposures,'20160101'
;
; By D. Nidever  Feb 2019
;-

pro delvered_exposures_fixcorrupted,input,delvedir=delvedir,redo=redo,dozeropoint=dozeropoint,$
                       doapcor=doapcor,uselocal=uselocal,startfresh=startfresh,stp=stp,$
                       nmulti=nmulti,workdir=workdir

;; Defaults
;if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl2/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Logs directory
logsdir = expdir+'logs/'
if file_test(logsdir,/directory) eq 0 then file_mkdir,logsdir
if n_elements(workdir) eq 0 then $
  workdir = '/data0/dnidever/delve/'  ;; temporary work directory
if n_elements(nmulti) eq 0 then nmulti=20
if n_elements(uselocal) eq 0 then uselocal=0

COMMON photred,setup

;; Not enough inputs
if n_elements(input) eq 0 then begin
  print,'Syntax - delvered_exposures,input,delvedir=delvedir,redo=redo,stp=stp'
  return
endif

; Start the logfile
;------------------
; format is delvered_exposures.host.DATETIME.log
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
logfile = logsdir+'delvered_exposures.'+hostname+'.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

;; Get all of the nights
dirs = FILE_SEARCH(expdir+'20??????',/test_directory,count=ndirs)
dirs = FILE_BASENAME(dirs)
numdirs = long(dirs)

;; Parse the input nights
for i=0,n_elements(input)-1 do begin
  input1 = input[i]
  ;; See if there is a -
  if strpos(input1,'-') ne -1 then begin
    nightrange = strsplit(input1,'-',/extract)
    ind_dirs = where(numdirs ge long(nightrange[0]) and numdirs le long(nightrange[1]),nind_dirs)
    if nind_dirs gt 0 then push,nights,dirs[ind_dirs] else print,'No directories found matching ',input1
  endif else begin
    MATCH,dirs,input1,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then push,nights,input1 else print,input1,' directory not found'
  endelse
endfor
nnights = n_elements(nights)
if nnights eq 0 then begin
  print,'No nights to process'
  return
endif

;; Which email program is available
mailprog = ''
spawn,'which sendmail',out,errout
if file_test(out[0]) eq 1 then mailprog='sendmail'
if mailprog eq '' then begin
  spawn,'which mail',out,errout
  if file_test(out[0]) eq 1 then mailprog='mail'
endif

;; Load the list of corrupted images
corrupted = MRDFITS('/net/dl2/dnidever/delve/bricks/db/delvemc_corrupted_exposures.fits',1)

; Print info
;-----------
print,''
print,'#######################################################'
print,'Starting DELVERED_EXPOSURES   ',systime(0)
print,'Running on ',host
print,'#######################################################'
print,''


; Make sure we have the right printlog.pro, not Markwardt's version
tempprogs = strsplit(!path,':',/extract)+'/printlog.pro'
test = file_test(tempprogs)
ind = where(test eq 1,nind)
bd = where(stregex(tempprogs[ind],'markwardt',/boolean) eq 1,nbd)
if nbd gt 0 then begin
  baddir = file_dirname(tempprogs[ind[bd]])
  print,"There is a version of Markwardt's PRINTLOG.PRO in "+baddir
  print,'Please rename this program (i.e. printlog.pro.orig)'
  return
endif

print,'Processing ',strtrim(nnights,2),' nights of data'

;#########################################
;#  STARTING THE PROCESSING
;#########################################
CD,current=origdir

;; Night loop
FOR i=0,nnights-1 do begin
  t0 = systime(1)
  inight = nights[i]
  CD,expdir+inight
  print,'' & print,'================================='
  print,' RUNNING PHOTRED on ',inight
  print,'=================================' & print,''

  if file_test('photred.setup') eq 0 then begin
    print,'NO photred.setup'
    goto,nightbomb
  endif
  ; LOAD THE SETUP FILE
  undefine,setup
  PHOTRED_LOADSETUP,setup,count=count
  if (count lt 1) then begin
    print,'Problem with photred.setup file'
    goto,nightbomb
  endif

  ;; Delete old files
  if keyword_set(startfresh) then $
       DELVERED_EXPOSURES_DELETE,expdir+inight

  ;; Check for corrupted files
  ;;--------------------------
  ;; Get the list of all the exposures
  chipfiles = FILE_SEARCH('F*/chip01/F*_01.fits',count=nchipfiles)
  expnum = file_basename(chipfiles,'_01.fits')
  expnum = (strsplitter(expnum,'-',/extract))[1,*]
  uexpnum = expnum[uniq(expnum,sort(expnum))]
  print,strtrim(n_elements(uexpnum),2),' unique exposures'
  ;; match to corrupted exposures list
  MATCH,uexpnum,corrupted.expnum,ind1,ind2,/sort
  dum = where(ind1 gt -1,nmatch)
  print,strtrim(nmatch,2),' corrupted images'
  if nmatch eq 0 then continue
  corrupted_expnum = uexpnum[ind1]
  corruptedstr = corrupted[ind2]
  print,strjoin(corrupted_expnum,' ')

  ;; Fix the resource files
  DELVERED_EXPOSURES_FIXFILES,corruptedstr,fixedfiles
  fixedfields = (strsplitter(fixedfiles,'/',/extract))[0,*]
  fixedfields = fixedfields[uniq(fixedfields,sort(fixedfields))]
  print,strtrim(n_elements(fixedfields),2),' fields have corrupted files'
  print,strjoin(fixedfields,' ')

  ;; Write fixed files list to a file
  WRITELINE,'fixedcorruptedfiles.lst',fixedfiles

  ;; gunzip any compressed ap, raw files
  print,'Uncompressing gzipped files.  This will take several minutes'
  ;; the "yes n" will repeatedly print "n" with a carriage return to any
  ;; questions/prompts that gunzip gives for an existing uncompressed file
  SPAWN,['find | grep .ap.gz | xargs -L 1 -I % gunzip -f %'],out,errout
  SPAWN,['find | grep .raw.gz | xargs -L 1 -I % gunzip -f %'],out,errout

  ;; Copy everything to the work directory
  if keyword_set(uselocal) then begin
    if file_test(workdir+inight) eq 0 then FILE_MKDIR,workdir+inight
    print,''
    print,'Copying all files to temporary directory >>'+workdir+inight+'<<'
    SPAWN,['rsync','-av',expdir+inight+'/',workdir+inight+'/'],out1,errout1,/noshell
    if n_elements(errout1) gt 1 or errout1[0] ne '' then begin
      print,'There was a problem with the rsync '
      printline,errout1
      return
    endif
    print,strtrim(n_elements(out1)-4,2)+' files/directories copied'
    printline,out1
    ;; Go to the temporary directory
    CD,workdir+inight
  endif

  ;; Adding files to WCS.inlist
  print,'Adding files to WCS.inlist'
  READLINE,'logs/WCS.inlist',wcslist,count=nwcslist
  push,wcslist,fixedfiles
  WRITELINE,'logs/WCS.inlist',wcslist

  ;; Make sure the WCS.inlist files are relative
  READLINE,'logs/WCS.inlist',wcslist,count=nwcslist
  if nwcslist gt 0 then begin
    wcslist = repstr(wcslist,expdir+inight+'/','')
    WRITELINE,'logs/WCS.inlist',wcslist
  endif

  DELVERED_WCS,/redo,nmulti=nmulti

  ;; DAOPHOT will automatically pick up the output from WCS
  DELVERED_DAOPHOT,/redo,nmulti=nmulti,workdir=workdir

  ;; Put all als files of affected fields in MATCH.inlist
  READLINE,'logs/DAOPHOT.success',daolist,count=ndaolist
  alslist = repstr(daolist,'.fits','.als')
  alsfields = (strsplitter(alslist,'/',/extract))[0,*]
  undefine,gd
  for f=0,n_elements(fixedfields)-1 do begin
    gd1 = where(alsfields eq fixedfields[f],ngd1)
    if ngd1 gt 0 then push,gd,gd1
  endfor
  READLINE,'logs/MATCH.inlist',matchlist,count=nmatchlist
  push,matchlist,alslist[gd]
  matchlist = matchlist[uniq(matchlist,sort(matchlist))]
  WRITELINE,'logs/MATCH.inlist',matchlist

  DELVERED_MATCH,/redo,nmulti=nmulti

  ;; Run APCOR on all exposures
  READLINE,'logs/DAOPHOT.success',daolist
  READLINE,'logs/APCOR.inlist',apcorlist
  push,apcorlist,daolist
  apcorlist = apcorlist[uniq(apcorlist,sort(apcorlist))]
  WRITELINE,'logs/APCOR.inlist',apcorlist
  print,systime()
  PHOTRED_APCOR,/redo

  ;; ZEROPOINT, rerun fixed exposures
  READLINE,'logs/ZEROPOINT.inlist',zerolist
  push,zerolist,fixedfiles
  zerolist = zerolist[uniq(zerolist,sort(zerolist))]
  WRITELINE,'logs/ZEROPOINT.inlist',zerolist
  print,systime()
  DELVERED_ZEROPOINT,/redo,nmulti=nmulti


  ;; Rerun astrom-save on all files/fields
  READLINE,'logs/ASTROM.inlist',astromlist
  READLINE,'logs/ASTROM.success',astromsuccess
  push,astromlist,astromsuccess
  READLINE,'logs/MATCH.outlist',matchoutlist
  push,astromlist,matchoutlist
  astromlist = astromlist[uniq(astromlist,sort(astromlist))]
  WRITELINE,'logs/ASTROM.inlist',astromlist

  PHOTRED_ASTROM,/redo,nmulti=nmulti
  PHOTRED_CALIB,/redo,nmulti=nmulti
  PHOTRED_COMBINE,/redo,nmulti=nmulti
  PHOTRED_DEREDDEN,/redo,nmulti=nmulti

  PHOTRED_SAVE,/redo,/sumquick,nmulti=nmulti

  print,'DELVERED FINISHED'
  stages = ['WCS','DAOPHOT','MATCH','APCOR','ASTROM','ZEROPOINT','CALIB','COMBINE','DEREDDEN','SAVE']
  PHOTRED_SUMMARY,outlines=outlines,stages=stages,/quick

  ;; Copy everything back to permanent directory
  if keyword_set(uselocal) then begin
    print,''
    print,'Copying all files back to permanent directory >>'+expdir+inight+'<<'
    print,''
    SPAWN,['rsync','-av',workdir+inight+'/',expdir+inight+'/'],out2,errout2,/noshell
    if n_elements(errout2) gt 1 or errout2[0] ne '' then begin
      print,'There was a problem with the rsync.'
      printline,errout2
      return
    endif
    printline,out2

    ;; SHOULD I CHECK HERE THAT EVERYTHING WAS COPIED CORRECTLY???

    CD,origdir

    ;; Delete temporary directory
    print,''
    print,'Deleting temporary directory >>'+workdir+inight+'<<'
    print,''
    SPAWN,['rm','-R',workdir+inight],/noshell

    ;; Fix the absolute paths in the summary files
    print,'Fixing absolute paths in summary files'
    sumfiles = file_search(expdir+inight+'/*_summary.fits',count=nsumfiles)
    for j=0,nsumfiles-1 do begin
      print,strtrim(j+1,2),' ',sumfiles[j]
      expstr = mrdfits(sumfiles[j],1,/silent)
      chstr = mrdfits(sumfiles[j],2,/silent)
      chstr.file = repstr(chstr.file,workdir,expdir)
      MWRFITS,expstr,sumfiles[j],/create
      MWRFITS,chstr,sumfiles[j],/silent
    endfor
  endif  ; /uselocal

  ;; Create the nightly summary file
  DELVERED_NIGHTSUMMARY,inight,delvedir=delvedir,redo=redo

  dt = systime(1)-t0
  print,'dt = ',strtrim(dt,2),' sec.'

  ;; Send email that this night is done
  CASE mailprog of
  'sendmail': begin
    undefine,elines
    push,elines,'From: noreply.delvered@noirlab.edu'
    push,elines,'To: dnidever@montana.edu'
    push,elines,'Subject: delvered_exposures night='+inight+' FINISHED on '+hostname
    push,elines,'Content-Type: text/html'
    push,elines,'MIME-Version: 1.0'
    push,elines,'<pre>'
    push,elines,'delvered_exposures night='+inight+' FINISHED'
    push,elines,systime(0)
    push,elines,'HOST='+hostname
    push,elines,'dt='+strtrim(dt,2)+' sec.'
    push,elines,''
    push,elines,outlines
    push,elines,'</pre>'
    tempfile = mktemp('mail',/nodot)
    writeline,tempfile,elines
    cmd = 'cat '+tempfile+' | sendmail -t'
    spawn,cmd,out,errout
    file_delete,tempfile
  end
  'mail': begin
    undefine,elines
    push,elines,'delvered_exposures night='+inight+' FINISHED'
    push,elines,systime(0)
    push,elines,'HOST='+hostname
    push,elines,'dt='+strtrim(dt,2)+' sec.'
    push,elines,''
    push,elines,outlines
    tempfile = mktemp('mail',/nodot)
    writeline,tempfile,elines
    cmd = 'cat '+tempfile+' | mail  -s "delvered_exposures night='+inight+' FINISHED on '+hostname+'" dnidever@noao.edu'
    spawn,cmd,out,errout
    file_delete,tempfile
  end
  else: print,'No mail program available'
  ENDCASE

  NIGHTBOMB:
ENDFOR

; End logfile
;------------
JOURNAL

;; Send out final email that we are done
sinput = strtrim(input,2)
if n_elements(sinput) gt 1 then sinput='['+strjoin(sinput,',')+']'
body = 'delvered_exposures '+sinput+' FINISHED at '+systime(0)+' on HOST='+hostname
cmd = 'echo "'+body+'" |  mail -s "delvered_exposures FINISHED on '+hostname+'" dnidever@montana.edu'
spawn,cmd,out,errout

if keyword_set(stp) then stop

end
