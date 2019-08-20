pro delvered_bricks_prep,brick,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo

;; This bricks pre-processing script gets DELVE and community MC data ready
;; to run PHOTRED ALLFRAME on it.

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/dl1/users/dnidever/delve/'
if FILE_TEST(delvedir,/directory) eq 0 then FILE_MKDIR,delvedir
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
;expfile = '/home/dnidever/projects/delvered/data/decam_mcs_20181009.fits.gz'
if n_elements(scriptsdir) gt 0 then scriptsdir=trailingslash(scriptsdir) else scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
if n_elements(irafdir) gt 0 then irafdir=trailingslash(irafdir) else irafdir='/home/dnidever/iraf/'
tempdir = '/tmp/'
if n_elements(workdir) eq 0 then begin
  undefine,workdir
  host = getenv('HOST')
  hostname = first_el(strsplit(host,'.',/extract))
  workdir='/data0/dnidever/delve/'
  if n_elements(workdir) gt 0 then if FILE_TEST(workdir,/directory) eq 0 then FILE_MKDIR,workdir
endif
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'
;; Bricks directory
brickdir = trailingslash(delvedir)+'bricks/'
if file_test(brickdir,/directory) eq 0 then file_mkdir,brickdir

;; Print information
print,'--------------------------------------------------------------------'
print,' PREPARING DELVE MC EXPOSURES FOR PHOTRED BRICK-LEVEL PROCESSING'
print,'--------------------------------------------------------------------'
print,'DELVEDIR = ',delvedir
print,'SCRIPTSDIR = ',scriptsdir
print,'IRAFDIR = ',irafdir
print,'EXPOSUREDIR = ',expdir
print,'WORKINGDIR = ',workdir
print,'HOST = ',host

;; Chip name and CCDNUM relations
decam = IMPORTASCII(delvereddir+'data/decam.txt',/header,/silent)

;; The photed.setup file
setup = ['##### REQUIRED #####',$
         'scriptsdir  '+scriptsdir,$
         'irafdir     '+irafdir,$
         'telescope   Blanco',$
         'instrument  DECAM',$
         'observatory CTIO',$
         'nmulti      10',$
         'nmulti_wcs       40',$
         'nmulti_daophot   30',$
         'nmulti_allframe  10',$
         'filtref     g,i,r,z,u',$
         'trans       delve.trans',$
         '##### OPTIONAL #####',$
         'sepfielddir  1',$
         'sepchipdir   1',$
         'keepmef      0',$
         'catformat    FITS',$
         'workdir      '+workdir,$
         'clean        1',$
         'skipcheck    1',$
         'redo         0',$
         'wcsrefname   GAIADR2',$
         'searchdist   20',$
         '#wcsrmslim   1.0',$
         'hyperthread   1',$
         'daopsfva      1',$
         'daofitradfwhm 1.0',$
         'psfcomsrc     0',$
         'psfcomglobal  0',$
         'psfcomgauss   0',$
         '#mchmaxshift  50.0',$
         'finditer      2',$
         'alfdetprog  sextractor',$
         '#alfnocmbimscale 0',$
         'alftrimcomb   0',$
         '#ddo51radoffset  1',$
         'cmbforce      1',$
         'keepinstr     1',$
         'avgmag        1',$
         'avgonlymag    0',$
         'todered       u,g,r,i,z,g-i',$
         '##### STAGES #####',$
         '#rename',$
         '#split',$
         ' wcs',$
         ' daophot',$
         ' zeropoint',$
         ' match',$
         '#allframe',$
         ' apcor',$
         ' astrom',$
         ' calib',$
         ' combine',$
         ' deredden',$
         ' save',$
         '#html']

;; Load the brick information
brickstr = MRDFITS(delvereddir+'data/delvemc_bricks_0.25deg.fits.gz',1)

;; Get the brick information
bind = where(brickstr.brickname eq brick,nbind)
if nbind eq 0 then begin
  printlog,logfile,ibrick+' not in DELVE-MC brick list'
  return
endif
brickstr1 = brickstr[bind[0]]

;; Get the list of exposures/chips that overlap this brick
print,'Getting chips that overlap this brick'
cenra = brickstr1.ra
cendec = brickstr1.dec
tmpfile = MKTEMP('tmp',/nodot,outdir=tempdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
tmpfile += '.fits'
spawn,['query_delvered_table',strtrim(cenra,2),strtrim(cendec,2),tmpfile,'--lim 0.5'],out,errout,/noshell
chstr = MRDFITS(tmpfile,1,/silent)
stop
file_delete,tmpfile,/allow

;; Get the catalog of DECam exposures
;;  Pick the most recent file
;expfile = delvereddir+'data/decam_mcs_20181009.fits.gz'
expfile = delvereddir+'data/decam_mcs_20190817.fits.gz'
print,'' & print,'Loading ',expfile
allexpstr = MRDFITS(expfile,1,/silent)
nallexp = n_elements(allexpstr)
allexpstr.fluxfile = strtrim(allexpstr.fluxfile,2)
allexpstr.maskfile = strtrim(allexpstr.maskfile,2)
allexpstr.wtfile = strtrim(allexpstr.wtfile,2)
allexpstr.expnum = strtrim(allexpstr.expnum,2)
;; Fix names for hulk/thing
;;  /net/mss1 -> /mss1
allexpstr.fluxfile = strmid(allexpstr.fluxfile,4)
allexpstr.maskfile = strmid(allexpstr.maskfile,4)
allexpstr.wtfile = strmid(allexpstr.wtfile,4)

;; Create summary file of ALL exposures and ALL chips

;; Get all of the nights
dirs = FILE_SEARCH(expdir+'20??????',/test_directory,count=ndirs)
dirs = FILE_BASENAME(dirs)
numdirs = long(dirs)
;; Nightly summary files
allnightsumfiles = dirs+'/'+file_basename(dirs)+'_summary.fits.gz'
gdnightsumfiles = where(file_test(allnightsumfiles) eq 1,nnightsumfiles)
print,strtrim(nnightsumfiles,2),' nightly summary files found'
nightsumfiles = allnightsumfiles[gdnightsumfiles]
;; Load in the data
expcount = 0L
chcount = 0L
for i=0,nnightsumfiles-1 do begin
  expstr1 = MRDFITS(nightsumfiles[i],1,/silent)
  nexpstr1 = n_elements(expstr1)
  chstr1 = MRDFITS(nightsumfiles[i],1,/silent)
  nchstr1 = n_elements(chstr1)
  if i eq 0 then begin
    expschema = expstr1[0]
    struct_assign,{dum:''},expschema
    expstr = replicate(expschema,50000L)
    nexpstr = n_elements(expstr)
    chschema = chstr1[0]
    struct_assign,{dum:''},chschema
    chstr = replicate(chschema,2500000L)
    nchstr = n_elements(chstr)
  endif
  ;; add more elements
  if (expcount+nexpstr1) gt nexpstr then begin
    orig = expstr
    expstr = replicate(expschema,nexpstr+20000L)
    expstr[0] = orig
    nexpstr = n_elements(expstr)
    undefine,orig
  endif
  if (chcount+nchstr1) gt nchstr then begin
    orig = chstr
    chstr = replicate(chschema,nchstr+200000L)
    chstr[0] = orig
    nchstr = n_elements(chstr)
    undefine,orig
  endif
  ;; Stuff in new exposures
  expstr[expcount:expcount+nexpstr1-1] = expstr1
  expcount += nexpstr1
  chstr[chcount:chcount+nchstr1-1] = chstr1
  chcount += nchstr1
endfor
;; Trim extra elements
if nexpstr gt expcount then expstr=expstr[0:expcount-1]
if nchstr gt chcount then chstr=chstr[0:chcount-1]

;; APPLY Quality cuts on the exposures

stop


;; STEP 1) Pick the exposures that we want
;;-----------------------------------------
;; Which ones do we want, make sure they haven't been prepped yet
;; only g/r/i and exptime>=90s for now
filt = strmid(allexpstr.filter,0,1)
;gdexp = where((filt eq 'g' or filt eq 'r' or filt eq 'i') and allexpstr.exposure ge 90,ngdexp)
gdexp = where((filt eq 'u' or filt eq 'g' or filt eq 'r' or filt eq 'i' or filt eq 'z' or filt eq 'Y') and allexpstr.exposure ge 90,ngdexp)
print,strtrim(ngdexp,2),' EXPOSURES pass the cuts'
expstr = allexpstr[gdexp]
;; Check if they've been prepped yet
observatory,'ctio',obs
night = strarr(ngdexp)
for i=0,ngdexp-1 do begin
  expstr1 = expstr[i]
  dateobs = expstr1.date_obs
  dateobs = strjoin(strsplit(dateobs,' ',/extract),'T')
  ;; Get night number
  jd = date2jd(dateobs)
  jd -= obs.tz/24.0  ;; correct for time zone difference
  caldat,double(long(jd)),Month, Day, Year, Hour, Minute, Second  ;; get date for beginning of the night
  ;; Convert to YYYYMMDD format
  night1 = string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')
  night[i] = night1
endfor
;; Get unique nights
night_index = CREATE_INDEX(night)
nnights = n_elements(night_index.value)
print,strtrim(nnights,2),' unique nights of data'

;; SMASH directories
smashnights = file_search('/dl1/users/dnidever/smash/cp/red/photred/20??????/',/test_directory,count=nsmashnights)
smashnights = file_basename(smashnights)

;; Loop over the nights
;;---------------------------
undefine,cmd,cmddir
For n=0,nnights-1 do begin
  inight = night_index.value[n]
  ind = night_index.index[night_index.lo[n]:night_index.hi[n]]
  nind = night_index.num[n]
  expstr1 = expstr[ind]
  nightdir = expdir+strtrim(inight,2)+'/'

  ;; Is this a SMASH night?
  MATCH,inight,smashnights,ind1,ind2,count=nmatch
  if nmatch gt 0 then begin
    print,'SMASH night.  Skipping.'
    goto,NIGHTBOMB
  endif

  if file_test(nightdir,/directory) eq 0 then FILE_MKDIR,nightdir
  expfile = nightdir+inight+'_exposures.fits'
  MWRFITS,expstr1,expfile,/create
  print,'Saving exposure structure to ',expfile
  push,cmd,'delvered_prep_night,"'+expfile+'"'
  NIGHTBOMB:
Endfor

ncmd = n_elements(cmd)
cmddir = strarr(ncmd)+workdir

stop

;; Run PBS_DAEMON
print,'Starting the JOBS'
nmulti = 10
PBS_DAEMON,cmd,cmddir,jobs=jobs,nmulti=nmulti,/hyperthread,/idle,prefix='dlvntprep',wait=5

;stop

end
