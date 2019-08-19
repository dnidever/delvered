;+
;
; DELVERED_BRICKS
;
; This is a script to run PHOTRED ALLFRAME on DELVE MC exposures to create
; forced-photometry catalogs.
;
; INPUTS:
;  input    What bricks to run PHOTRED/ALLFRAME on.  Either an array or
;            a range such as 100-110.
;
; OUTPUTS:
;  PHOTRED_ALLFRAME will be run on each brick and a final catalog and
;  summary file created for each brick.
;
; USAGE:
;  IDL>delvered_bricks,100
;
; By D. Nidever  Aug 2019
;-

pro delvered_bricks,input,delvedir=delvedir,redo=redo,stp=stp

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/dl1/users/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Bricks directory
brickdir = trailingslash(delvedir)+'bricks/'
;; Logs directory
logsdir = expdir+'logs/'
if file_test(logsdir,/directory) eq 0 then file_mkdir,logsdir

;; Not enough inputs
if n_elements(input) eq 0 then begin
  print,'Syntax - delvered_bricks,input,delvedir=delvedir,redo=redo,stp=stp'
  return
endif

; Start the logfile
;------------------
; format is delvered_bricks.host.DATETIME.log
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
logfile = logsdir+'delvered_bricks.'+hostname+'.'+smonth+sday+syear+shour+sminute+ssecond+'.log'
JOURNAL,logfile

;; Get all of the bricks
dirs = FILE_SEARCH(brickdir+'*',/test_directory,count=ndirs)
dirs = FILE_BASENAME(dirs)
numdirs = long(dirs)

;; Parse the input nights
for i=0,n_elements(input)-1 do begin
  input1 = input[i]
  ;; See if there is a -
  if strpos(input1,'-') ne -1 then begin
    brickrange = strsplit(input1,'-',/extract)
    ind_dirs = where(numdirs ge long(brickrange[0]) and numdirs le long(brickrange[1]),nind_dirs)
    if nind_dirs gt 0 then push,bricks,dirs[ind_dirs] else print,'No directories found matching ',input1
  endif else begin
    MATCH,dirs,input1,ind1,ind2,/sort,count=nmatch
    if nmatch gt 0 then push,bricks,input1 else print,input1,' directory not found'
  endelse
endfor
nbricks = n_elements(bricks)
if nbricks eq 0 then begin
  print,'No bricks to process'
  return
endif

;; Load the brick information
brickstr = MRDFITS(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz',1)

; Print info
;-----------
print,''
print,'############################################'
print,'Starting DELVERED_BRICKS   ',systime(0)
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

print,'Processing ',strtrim(nbricks,2),' bricks'


;#########################################
;#  STARTING THE PROCESSING
;#########################################
CD,current=origdir

;; Brick loop
FOR i=0,nbricks-1 do begin
  ibrick = bricks[i]
  CD,brickdir+ibrick
  print,'' & print,'================================='
  print,' RUNNING PHOTRED on BRICK ',ibrick
  print,'=================================' & print,''

  ;; Get the brick information
  bind = where(brickstr.brickname eq ibrick,nbind)
  if nbind eq 0 then begin
    printlog,logfile,ibrick+' not in DELVE-MC brick list'
    goto,BOMB
  endif
  brickstr1 = brickstr[bind[0]]

  ;; Prepare the input for ALLFRAME??
  ;; Need to setup symlinks to exposure-level files
  PHOTRED_MATCH,redo=redo
  PHOTRED_ALLFRAME,redo=redo
  PHOTRED_ASTROM,redo=redo
  PHOTRED_CALIB,redo=redo
  PHOTRED_COMBINE,redo=redo  ; only 1 brick/tile, nothing to combine
  PHOTRED_DEREDDEN,redo=redo
  PHOTRED_SAVE,redo=redo

  print,'PHOTRED FINISHED'

  ;;PHOTRED_SUMMARY

  ;;; Create the brick summary file
  ;sumfiles = FILE_SEARCH('*_summary.fits',count=nsumfiles)
  ;undefine,brickstr
  ;fields = IMPORTASCII('fields',fieldnames=['shname','name'],fieldtypes=[7,7],/silent)
  ;For j=0,nsumfiles-1 do begin
  ;  base = FILE_BASENAME(sumfiles[j],'_summary.fits')
  ;  MATCH,fields.name,base,ind1,ind2,/sort,count=nmatch
  ;  shname = fields[ind1[0]].shname
  ;  expstr1 = MRDFITS(sumfiles[j],1,/silent)
  ;  add_tag,expstr1,'fieldname',base,expstr1
  ;  add_tag,expstr1,'field',shname,expstr1
  ;  PUSH,brickstr,expstr1
  ;  chipstr1 = MRDFITS(sumfiles[j],2,/silent)
  ;  add_tag,chipstr1,'fieldname',base,chipstr1
  ;  PUSH,chipstr,chipstr1
  ;Endfor
  ;bricksumfile = brickdir+ibrick+'/'+ibrick+'_summary.fits'
  ;print,'Writing brick summary file to ',bricksumfile
  ;MWRFITS,expstr,bricksumfile,/create
  ;MWRFITS,chipstr,bricksumfile,/silent

  BOMB:
ENDFOR

; End logfile
;------------
JOURNAL

if keyword_set(stp) then stop

end
