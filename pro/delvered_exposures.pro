;+
;
; DELVERED_EXPOSURES
;
; This is a script to run PHOTRED on DELVE MC exposures to create
; single-image level catalogs and PSFs.
;
; INPUTS:
;  input    What nights to run PHOTRED on.  Either an array or
;            a range such as 20160101-20160506.
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

pro delvered_exposures,input,delvedir=delvedir,redo=redo,stp=stp

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/net/dl1/users/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Logs directory
logsdir = expdir+'logs/'
if file_test(logsdir,/directory) eq 0 then file_mkdir,logsdir
workdir = '/data0/dnidever/delve/'  ;; temporary work directory

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

  ;; Copy everything to the work directory
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

  ;; Make sure the WCS.inlist files are relative
  READLINE,'logs/WCS.inlist',wcslist,count=nwcslist
  if nwcslist gt 0 then begin
    wcslist = repstr(wcslist,expdir+inight+'/','')
    WRITELINE,'logs/WCS.inlist',wcslist
  endif

  if READPAR(setup,'WCS') ne '0' then DELVERED_WCS,redo=redo
  if READPAR(setup,'DAOPHOT') ne '0' then PHOTRED_DAOPHOT,redo=redo
  if READPAR(setup,'MATCH') ne '0' then DELVERED_MATCH,redo=redo
  if READPAR(setup,'APCOR') ne '0' then PHOTRED_APCOR,redo=redo
  if READPAR(setup,'ASTROM') ne '0' then PHOTRED_ASTROM,redo=redo
  if READPAR(setup,'ZEROPOINT') ne '0' then DELVERED_ZEROPOINT,redo=redo
  if READPAR(setup,'CALIB') ne '0' then PHOTRED_CALIB,redo=redo
  if READPAR(setup,'COMBINE') ne '0' then PHOTRED_COMBINE,redo=redo
  if READPAR(setup,'DEREDDEN') ne '0' then PHOTRED_DEREDDEN,redo=redo
  if READPAR(setup,'SAVE') ne '0' then PHOTRED_SAVE,redo=redo,/sumquick

  print,'DELVERED FINISHED'
  stages = ['WCS','DAOPHOT','MATCH','APCOR','ASTROM','ZEROPOINT','CALIB','COMBINE','DEREDDEN','SAVE']
  PHOTRED_SUMMARY,outlines=outlines,stages=stages,/quick

  ;; Create the nightly summary file
  DELVERED_NIGHTSUMMARY,inight,delvedir=delvedir,redo=redo

  ;; Copy everything back to permanent directory
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
    print,strtrim(j+1,2),' ',sumfiles[i]
    expstr = mrdfits(sumfiles[i],1,/silent)
    chstr = mrdfits(sumfiles[i],2,/silent)
    chstr.file = repstr(chstr.file,workdir,expdir)
    MWRFITS,expstr,sumfiles[i],/create
    MWRFITS,chstr,sumfiles[i]
  endfor

  dt = systime(1)-t0
  print,'dt = ',strtrim(dt,2),' sec.'

  ;; Send email that this night is done
  undefine,elines
  push,elines,'From: dnidever@noao.edu'
  push,elines,'To: dnidever@noao.edu'
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

  NIGHTBOMB:
ENDFOR

; End logfile
;------------
JOURNAL

;; Send out final email that we are done
sinput = strtrim(input,2)
if n_elements(sinput) gt 1 then sinput='['+strjoin(sinput,',')+']'
body = 'delvered_exposures '+sinput+' FINISHED at '+systime(0)+' on HOST='+hostname
cmd = 'echo "'+body+'" |  mail -s "delvered_exposures FINISHED on '+hostname+'" dnidever@noao.edu'
spawn,cmd,out,errout

if keyword_set(stp) then stop

end
