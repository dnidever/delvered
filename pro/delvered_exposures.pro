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
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/dl1/users/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Logs directory
logsdir = expdir+'logs/'
if file_test(logsdir,/directory) eq 0 then file_mkdir,logsdir

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
print,'############################################'
print,'Starting DELVERED_EXPOSURES   ',systime(0)
print,'Running on ',host
print,'############################################'
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
  inight = nights[i]
  CD,expdir+inight
  print,'' & print,'================================='
  print,' RUNNING PHOTRED on ',inight
  print,'=================================' & print,''

  PHOTRED_WCS,redo=redo
  PHOTRED_DAOPHOT,redo=redo
  PHOTRED_MATCH,redo=redo
  PHOTRED_APCOR,redo=redo
  PHOTRED_ASTROM,redo=redo
  DELVERED_ZEROPOINT,redo=redo
  PHOTRED_CALIB,redo=redo
  PHOTRED_COMBINE,redo=redo
  PHOTRED_DEREDDEN,redo=redo
  PHOTRED_SAVE,redo=redo

  print,'PHOTRED FINISHED'
  PHOTRED_SUMMARY

  ;; Create the nightly summary file
  DELVERED_NIGHTSUMMARY,inight,delvedir=delvedir,redo=redo
ENDFOR

; End logfile
;------------
JOURNAL

if keyword_set(stp) then stop

end
