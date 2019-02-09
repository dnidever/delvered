;+
;
; DELVERED_EXPOSURES
;
; This is a script to run PHOTRED on DELVE MC exposures to create
; single-image level catalogs and PSFs.
;
; INPUTS:
;  nights   What nights to run PHOTRED on.  Either an array or
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

pro delvered_exposures,nights,delvedir=delvedir,redo=redo,stp=stp

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/dl1/users/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Logs directory
logsdir = expdir+'logs/'

;; Not enough inputs
if n_elements(nights) eq 0 then begin
  print,'Syntax - delvered_exposures,nights,delvedir=delvedir,redo=redo,stp=stp'
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

;; Parse the input nights
for i=0,n_elements(nights)-1 do begin
  inight = nights[i]
  ;; See if there is a -
  if strpos(inight,'-') ne -1 then begin
    nightrange = strsplit(inight,'-',/extract)
    MATCH,dirs,nightrange,ind1,ind2,/sort,count=nmatch
stop
  endif else begin
    MATCH,dirs,inight,ind1,ind2,/sort,count=nmatch
    stop
  endelse
endfor


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
  DELVERED_ZERPOINT,redo=redo
  PHOTRED_CALIB,redo=redo
  PHOTRED_COMBINE,redo=redo
  PHOTRED_DEREDDEN,redo=redo
  PHOTRED_SAVE,redo=redo

  print,'PHOTRED FINISHED'

  PHOTRED_SUMMARY

  ;; Create the nightly summary file
  sumfiles = FILE_SEARCH('*_summary.fits',count=nsumfiles)
  undefine,expstr,chipstr
  fields = IMPORTASCII('fields',fieldnames=['shname','name'],fieldtypes=[7,7],/silent)
  For j=0,nsumfiles-1 do begin
    base = FILE_BASENAME(sumfiles[j],'_summary.fits')
    MATCH,fields.name,base,ind1,ind2,/sort,count=nmatch
    shname = fields[ind1[0]].shname
    expstr1 = MRDFITS(sumfiles[j],1,/silent)
    add_tag,expstr1,'fieldname',base,expstr1
    add_tag,expstr1,'field',shname,expstr1
    PUSH,expstr,expstr1
    chipstr1 = MRDFITS(sumfiles[j],2,/silent)
    add_tag,chipstr1,'fieldname',base,chipstr1
    PUSH,chipstr,chipstr1
  Endfor
  nightsumfile = expdir+inight+'/'+inight+'_summary.fits'
  print,'Writing nightly summary file to ',nightsumfile
  MWRFITS,expstr,nightsumfile,/create
  MWRFITS,chipstr,nightsumfile,/silent

ENDFOR

; End logfile
;------------
JOURNAL

if keyword_set(stp) then stop

end
