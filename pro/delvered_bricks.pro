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
;  =nmulti  Number of simultaneous jobs to run.
;  =redo    Redo bricks that were previously done.
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

pro delvered_bricks,input,nmulti=nmulti,redo=redo,stp=stp

t0 = systime(1)
  
;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/dl1/users/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Bricks directory
brickdir = trailingslash(delvedir)+'bricks/'
;; Logs directory
logsdir = expdir+'logs/'
if file_test(logsdir,/directory) eq 0 then file_mkdir,logsdir
tempdir = '/tmp/'
if n_elements(workdir) eq 0 then begin
  undefine,workdir
  host = getenv('HOST')
  hostname = first_el(strsplit(host,'.',/extract))
  workdir = '/data0/dnidever/delve/'
  if n_elements(workdir) gt 0 then if FILE_TEST(workdir,/directory) eq 0 then FILE_MKDIR,workdir
endif
if n_elements(nmulti) eq 0 then nmulti=10

;; Not enough inputs
if n_elements(input) eq 0 then begin
  print,'Syntax - delvered_bricks,input,nmulti=nmulti,redo=redo,stp=stp'
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
logtime = smonth+sday+syear+shour+sminute+ssecond
logfile = logsdir+'delvered_bricks.'+hostname+'.'+logtime+'.log'
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

;; Check if bricks were previously done
if not keyword_set(redo) then begin
  outfiles = brickdir+bricks+'/'+bricks+'.fits.gz'
  bd = where(file_test(outfiles) eq 1,nbd,comp=gd,ncomp=ngd)
  if ngd eq 0 then begin
    print,'All bricks were previously processed.  Nothing to do.'
    return
  endif
  if nbd gt 0 then begin
    print,'REDO not set and ',strtrim(nbd,2),' bricks were previously processed. Removing them from current list to process.'
    bricks = bricks[gd]
    nbricks = ngd
  endif
endif
  
print,'Processing ',strtrim(nbricks,2),' bricks'

;#########################################
;#  STARTING THE PROCESSING
;#########################################
cmd = 'delvered_forcebrick,"'+bricks+'"'
if keyword_set(redo) then cmd += ',/redo'
cmddirs = strarr(nbricks)+workdir
PBS_DAEMON,cmd,cmddirs,jobs=jobs,prefix='dlvbrcks',/idle,/hyperthread,nmulti=nmulti,wait=5

print,'DELVERED_BRICKS done after ',strtrim(systime(1)-t0,2),' sec'

; End logfile
;------------
JOURNAL

if keyword_set(stp) then stop

end
