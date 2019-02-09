;+
;
; DELVERED_ZEROPOINT
;
; This program calculates NSC-like exposure-level zero-points
; for DELVERED_EXPOSURES.
;
; By D. Nidever  Feb 2019
;-

pro delvered_zeropoint,redo=redo

COMMON photred,setup

print,''
print,'############################'
print,'RUNNING DELVERED_ZEROPOINT'
print,'############################'
print,''

CD,current=curdir

; Does the logs/ directory exist?
testlogs = file_test('logs',/directory)
if testlogs eq 0 then FILE_MKDIR,'logs'

; Print date and time to the logfile
printlog,logfile,''
printlog,logfile,'Starting DELVERED_ZEROPOINT  ',systime(0)

; LOAD THE SETUP FILE if not passed
;-----------------------------------
; This is a 2xN array.  First colume are the keywords
; and the second column are the values.
; Use READPAR.PRO to read it
if n_elements(setup) eq 0 then begin
  PHOTRED_LOADSETUP,setup,count=count
  if count lt 1 then return
endif

; Are we redoing?
doredo = READPAR(setup,'REDO')
if keyword_set(redo) or (doredo ne '-1' and doredo ne '0') then redo=1

; Get the scripts directory from setup
scriptsdir = READPAR(setup,'SCRIPTSDIR')
if scriptsdir eq '' then begin
  printlog,logfile,'NO SCRIPTS DIRECTORY'
  return
endif

; LOAD THE "imagers" FILE
;----------------------------
printlog,logfile,'Loading imager information'
imagerstest = FILE_TEST(scriptsdir+'/imagers')
if (imagerstest eq 0) then begin
  printlog,logfile,'NO >>imagers<< file in '+scriptsdir+'  PLEASE CREATE ONE!'
  return
endif
; The columns need to be: Telescope, Instrument, Naps, separator
imagers_fieldnames = ['telescope','instrument','observatory','namps','separator']
imagers_fieldtpes = [7,7,7,3,7]
imagers = IMPORTASCII(scriptsdir+'/imagers',fieldnames=imagers_fieldnames,$
                      fieldtypes=imagers_fieldtypes,comment='#')
imagers.telescope = strupcase(strtrim(imagers.telescope,2))
imagers.instrument = strupcase(strtrim(imagers.instrument,2))
imagers.observatory = strupcase(strtrim(imagers.observatory,2))
singleind = where(imagers.namps eq 1,nsingle)
if nsingle gt 0 then imagers[singleind].separator = ''
if (n_tags(imagers) eq 0) then begin
  printlog,logfile,'NO imagers in '+scriptsdir+'/imagers'
  return
endif

; What IMAGER are we using??
;---------------------------
ind_imager = where(imagers.telescope eq telescope and imagers.instrument eq instrument,nind_imager)
if nind_imager eq 0 then begin
  printlog,logfile,'TELESCOPE='+telescope+' INSTRUMENT='+instrument+' NOT FOUND in >>imagers<< file'
  return
endif
thisimager = imagers[ind_imager[0]]
; print out imager info
printlog,logfile,''
printlog,logfile,'USING IMAGER:'
printlog,logfile,'Telescope = '+thisimager.telescope
printlog,logfile,'Instrument = '+thisimager.instrument
printlog,logfile,'Namps = '+strtrim(thisimager.namps,2)
printlog,logfile,"Separator = '"+thisimager.separator+"'"
printlog,logfile,''


;###################
; GETTING INPUTLIST
;###################

;; Get all of the files from DAOPHOT.success
READLIST,curdir+'/logs/DAOPHOT.success',fitsfiles,/unique,/fully,setupdir=curdir,count=nfitsfiles,logfile=logfile,/silent
;; Get the unique exposures
allbase = strarr(nfitsfiles)
for i=0,nfitsfiles-1 do begin
  if strmid(fitsfile,6,7,/reverse_offset) eq 'fits.fz' then allbase[i]=FILE_BASENAME(fitsfiles[i],'.fits.fz') else $
    allbase[i]=FILE_BASENAME(fitsfiles[i],'.fits')
endfor
;; Remove the ccdnum suffix
expbase = allbase
if thisimager.namps gt 1 then $
  for i=0,nfitsfiles-1 do expbase[i] = first_el(strsplit(expbase[i],thisimager.separator,/extract))
expindex = CREATE_INDEX(expbase)
nexp = n_elements(expindex.value)
printlog,logfile,strtrim(nexp,2),' exposures to process'

;########################################
;#  PROCESSING THE FILES
;########################################
printlog,logfile,''
printlog,logfile,'-----------------------'
printlog,logfile,'PROCESSING THE FILES'
printlog,logfile,'-----------------------'
printlog,logfile,''
printlog,logfile,systime(0)

;; Exposure loop
FOR i=0,nfields-1 do begin
  ind = expindex.index[expindex.lo[i]:expindex.hi[i]]
  nind = n_elements(ind)

  ;; Load all of the ALS files and add coordinates
  printlog,logfile,'Loading data for ',expindex.value[i]
  undefine,phot
  For j=0,nind-1 do begin
    dir1 = FILE_DIRNAME(fitsfiles[ind[j]])
    base1 = allbase[ind[j]]
    alsfile1 = dir1+'/'+base1+'.als'
    ccdnum = PHOTRED_GETCHIPNUM(base1,thisimager)
    ifield = first_el(strsplit(base1,'-',/extract))
    als0 = PHOTRED_READFILE(alsfile1,count=nals)
    ;; Change ID to a string
    schema = {id:'',x:0.0,y:0.0,ra:0.0d0,dec:0.0d0,mag:0.0,err:0.0,sky:0.0,iter:0.0,chi:0.0,sharp:0.0}
    als = REPLICATE(schema,nals)
    STRUCT_ASSIGN,als0,als,/nozero
    ;; Get the coordinates
    if strpos(fitsfiles[ind[j]],'.fits.fz') ne 1 then exten=1 else exten=0
    head1 = PHOTRED_READFILE(fitsfiles[ind[j]],exten=exten,/header)
    ;; Converting to IDL X/Y convention, starting at (0,0)
    ;; DAOPHOT has X/Y start at (1,1)
    HEAD_XYAD,head1,als.x-1.0,als.y-1.0,ra,dec,/degree
    als.ra = ra
    als.dec = dec
    ;; Modify the IDs
    if thisimager.namps gt 1 then begin
      ; FIELD_EXT.IDNUMBER, i.e. 190L182a_5.17366
      ;---------------------------------------------
      id2 = ifield+'_'+strtrim(ccdnum,2)+'.'+strtrim(als.id,2)
      als.id = id2
    endif else begin
      ;; Updating the IDs
      ;; FIELD.IDNUMBER, i.e. 190L182a.17366
      ;;---------------------------------------------
      id2 = ifield+'.'+strtrim(als.id,2)
      als.id = id2
    endelse
    PUSH,phot,als
  Endfor

  ;; Load the reference data for this field

  ;; Get the model magnitudes

  ;; Calculate the zeropoint

  stop

ENDFOR

;; Write the night.trans file
;; DO WE HAVE PREVIOUS DATA TO INCLUDE???

stop

end
