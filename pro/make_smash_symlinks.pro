pro make_smash_symlinks

;; Make symlinks to the SMASH data that has already been processedd

delvedir = '/dl1/users/dnidever/delve/exposures/'
smashdir = '/dl1/users/dnidever/smash/cp/red/photred/'
nmulti=5

nights = FILE_SEARCH(smashdir+'20??????',/test_directory,count=nnights)
CD,current=origdir

scriptsdir = '/home/dnidever/projects/PHOTRED/scripts/'
irafdir = '/home/dnidever/iraf/'
workdir = '/data0/dnidever/delve/'

print,'######################################'
print,strtrim(nnights,2),' SMASH nights to create symlinks for'
print,'######################################'

;; Night loop
undefine,cmd,cmddir
For i=0,nnights-1 do begin
  inight = FILE_BASENAME(nights[i])
  ;; Skip if the directory exists
  if file_test(delvedir+inight,/directory) eq 0 then begin
    cmd1 = "make_smash_symlinks_night,'"+inight+"'"
    push,cmd,cmd1
  endif else print,nights[i],' exists already.  Skipping.'
Endfor
ncmd = n_elements(cmd)
cmddir = strarr(ncmd)+'/data0/dnidever/delve/'

;; Run PBS_DAEMON
print,'Starting the JOBS'
PBS_DAEMON,cmd,cmddir,jobs=jobs,nmulti=nmulti,/hyperthread,/idle,prefix='smash',wait=5

stop


end
