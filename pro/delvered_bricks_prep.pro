pro delvered_bricks_prep,brick,scriptsdir=scriptsdirs,irafdir=irafdir,workdir=workdir,redo=redo,logfile=logfile

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
if n_elements(logfile) eq 0 then logfile=-1

;; Print information
printlog,logfile,'--------------------------------------------------------------------'
printlog,logfile,' PREPARING DELVE MC EXPOSURES FOR PHOTRED BRICK-LEVEL PROCESSING'
printlog,logfile,'--------------------------------------------------------------------'
printlog,logfile,'DELVEDIR = ',delvedir
printlog,logfile,'SCRIPTSDIR = ',scriptsdir
printlog,logfile,'IRAFDIR = ',irafdir
printlog,logfile,'EXPOSUREDIR = ',expdir
printlog,logfile,'WORKINGDIR = ',workdir
printlog,logfile,'HOST = ',host

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
printlog,logfile,'Getting chips that overlap this brick'
cenra = brickstr1.ra
cendec = brickstr1.dec
tmpfile = MKTEMP('tmp',/nodot,outdir=tempdir) & TOUCHZERO,tmpfile+'.fits' & FILE_DELETE,[tmpfile,tmpfile+'.fits'],/allow
tmpfile += '.fits'
spawn,['query_delvered_table',strtrim(cenra,2),strtrim(cendec,2),tmpfile,'--lim','0.5'],out,errout,/noshell
info = file_info(tmpfile)
if info.size eq 0 then begin
  printlog,logfile,'No overlapping chips found'
  return
endif
chstr = MRDFITS(tmpfile,1,/silent)
file_delete,tmpfile,/allow
nchstr = n_elements(chstr)
printlog,logfile,'Found ',strtrim(nchstr,2),' overlapping chips within 0.5 deg of brick center'

;; Do more rigorous overlap checking
bvra = [brickstr1.ra1,brickstr1.ra2,brickstr1.ra2,brickstr1.ra1]
bvdec = [brickstr1.dec1,brickstr1.dec1,brickstr1.dec2,brickstr1.dec2]
olap = intarr(nchstr)
for i=0,nchstr-1 do begin
  hd1 = PHOTRED_READFILE(chstr[i].file,/header)
  nx = sxpar(hd1,'naxis1')
  ny = sxpar(hd1,'naxis2')
  head_xyad,hd1,[0,nx-1,nx-1,0],[0,0,ny-1,ny-1],vra,vdec,/degree
  olap[i] = dopolygonsoverlap(bvra,bvdec,vra,vdec)
endfor
g = where(olap eq 1,ng)
if ng eq 0 then begin
  printlog,logfile,'No chips overlap this brick'
  return
endif
printlog,logfile,strtrim(ng,2),' chips overlap this brick'
chstr = chstr[g]
nchstr = ng

;; APPLY Quality cuts on the exposures
;;------------------------------------
fwhmthresh = 2.0    ; seeing 2.0" threshold
printlog,logfile,'Applying quality cuts'
gdch = where(chstr.fwhm*chstr.pixscale le fwhmthresh,ngdch)
if ngdch eq 0 then begin
  printlog,logfile,'No chips passed the quality cuts'
  return
endif
printlog,logfile,strtrim(ngdch,2),' chips passed the quality cuts'
chstr = chstr[gdch]
nchstr = ngdch
chstr.file = strtrim(chstr.file,2)
chstr.base = strtrim(chstr.base,2)

;; Create the symlinks
;;---------------------
;; Subdirectory is the first 4 digits of the brickname, e.g., 0952 of 0952m462, the RA portion of the name
subdir = brickdir+strmid(brick,0,4)+'/'
if file_test(subdir,/directory) eq 0 then file_mkdir,subdir
;; Create brick directory if it doesn't exist yet
bdir = subdir+brick+'/'
if file_test(bdir,/directory) eq 0 then file_mkdir,bdir
printlog,logfile,'Copying over the necessary files to ',bdir
;; Chip loop
for i=0,nchstr-1 do begin
  printlog,logfile,'Copying over files for ',chstr[i].file
  ;; Make chip directory if necessary
  chdir = bdir+'chip'+string(chstr[i].chip,format='(i02)')+'/'
  if file_test(chdir,/directory) eq 0 then file_mkdir,chdir
  odir = file_dirname(chstr[i].file)+'/'
  obase = PHOTRED_GETFITSEXT(chstr[i].file,/basename)
  ;; Need symlinks to .psf, .als
  if file_test(odir+obase+'.psf') eq 0 then begin
    printlog,logfile,odir+obase+'.psf NOT FOUND.  Skipping this chip'
    goto,BOMB1
  endif
  if file_test(odir+obase+'.als') eq 0 then begin
    printlog,logfile,odir+obase+'.als NOT FOUND.  Skipping this chip'
    goto,BOMB1
  endif
  FILE_DELETE,chdir+chstr[i].base+['.psf','.als'],/allow
  FILE_LINK,odir+obase+'.psf',chdir+chstr[i].base+'.psf'
  FILE_LINK,odir+obase+'.als',chdir+chstr[i].base+'.als'
  ;; Copy the fits, fits resource file and header files locally
  if file_test(odir+obase+'.fits') eq 0 then begin
    printlog,logfile,odir+obase+'.fits NOT FOUND.  Skipping this chip'
    goto,BOMB1
  endif
  FILE_COPY,odir+obase+'.fits',chdir,/over
  if file_test(odir+'.'+obase+'.fits') eq 1 then FILE_COPY,odir+'.'+obase+'.fits',chdir,/over
  if file_test(odir+obase+'.head') eq 1 then FILE_COPY,odir+obase+'.head',chdir,/over
  BOMB1:
endfor

; Make the tiling file
;---------------------
undefine,lines
; Lines with the tiling scheme first
nx = 3600
ny = 3600
step = 0.25
xref = nx/2
yref = ny/2
ntiles = 1
push,lines,'CENRA  = '+strtrim(brickstr1.ra,2)
push,lines,'NX     = '+strtrim(nx,2)
push,lines,'XSTEP  = '+strtrim(step,2)
push,lines,'XREF   = '+strtrim(xref+1,2)
push,lines,'CENDEC = '+strtrim(brickstr1.dec,2)
push,lines,'NY     = '+strtrim(ny,2)
push,lines,'YSTEP  = '+strtrim(step,2)
push,lines,'YREF   = '+strtrim(yref+1,2)
push,lines,'NTILES = '+strtrim(ntiles,2)
; Then one file with info for each tile
tilestr = {num:1,name:'tile1',x0:0,x1:nx-1,nx:nx,y0:0,y1:ny-1,ny:ny,nimages:nchstr}
ntiles = 1
for i=0,ntiles-1 do begin
  tilestr1 = tilestr[i]
  fmt = '(I-6,A8,6I8,I6)'
  ; Using IRAF indexing
  line = string(format=fmt,tilestr1.num,tilestr1.name,tilestr1.x0+1,tilestr1.x1+1,tilestr1.nx,$
                tilestr1.y0+1,tilestr1.y1+1,tilestr1.ny,tilestr1.nimages)
  push,lines,line
endfor
tilefile = bdir+brick+'.tiling'
printlog,logfile,'Writing tiling information to >>'+tilefile+'<<'
WRITELINE,tilefile,lines

;; Put files into logs/MATCH.inlist file


stop


;; STEP 1) Pick the exposures that we want
;;-----------------------------------------
;; Which ones do we want, make sure they haven't been prepped yet
;; only g/r/i and exptime>=90s for now
filt = strmid(allexpstr.filter,0,1)
;gdexp = where((filt eq 'g' or filt eq 'r' or filt eq 'i') and allexpstr.exposure ge 90,ngdexp)
gdexp = where((filt eq 'u' or filt eq 'g' or filt eq 'r' or filt eq 'i' or filt eq 'z' or filt eq 'Y') and allexpstr.exposure ge 90,ngdexp)
printlog,logfile,strtrim(ngdexp,2),' EXPOSURES pass the cuts'
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
printlog,logfile,strtrim(nnights,2),' unique nights of data'

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
    printlog,logfile,'SMASH night.  Skipping.'
    goto,NIGHTBOMB
  endif

  if file_test(nightdir,/directory) eq 0 then FILE_MKDIR,nightdir
  expfile = nightdir+inight+'_exposures.fits'
  MWRFITS,expstr1,expfile,/create
  printlog,logfile,'Saving exposure structure to ',expfile
  push,cmd,'delvered_prep_night,"'+expfile+'"'
  NIGHTBOMB:
Endfor

ncmd = n_elements(cmd)
cmddir = strarr(ncmd)+workdir

stop

;; Run PBS_DAEMON
printlog,logfile,'Starting the JOBS'
nmulti = 10
PBS_DAEMON,cmd,cmddir,jobs=jobs,nmulti=nmulti,/hyperthread,/idle,prefix='dlvntprep',wait=5

;stop

end
