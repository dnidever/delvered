;+
;
; DELVERED_NIGHTSUMMARY
;
; This creates a summary file of all the taken on one night.
;
; INPUTS:
;  night    The night name, e.g. 20160101.
;
; OUTPUTS:
;  A night summary file will be created in the night directory
;  with name NIGHT_summary.fits.
;
; USAGE:
;  IDL>delvered_nightsummary,'20160101'
;
; By D. Nidever  Aug 2019
;-

pro delvered_nightsummary,night,delvedir=delvedir,redo=redo,stp=stp

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/dl1/users/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'

;; Not enough inputs
if n_elements(night) eq 0 then begin
  print,'Syntax - delvered_nightsummary,night,delvedir=delvedir,redo=redo,stp=stp'
  return
endif
inight = strtrim(night,2)

CD,current=origdir

;; Check if directory exists
if file_test(expdir+inight) eq 0 then begin
  print,'Directory '+expdir+inight+' NOT FOUND'
  return
endif
;; Check if output file already exists
nightsumfile = expdir+inight+'/'+inight+'_summary.fits'
if file_test(nightsumfile) eq 1 and not keyword_set(redo) then begin
  print,nightsumfile+' EXISTS and /redo NOT set'
  return
endif

CD,expdir+inight

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
print,'Writing nightly summary file to ',nightsumfile
MWRFITS,expstr,nightsumfile,/create
MWRFITS,chipstr,nightsumfile,/silent

if keyword_set(stp) then stop

end
