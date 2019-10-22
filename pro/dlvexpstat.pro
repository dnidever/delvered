pro dlvexpstat,delvedir=delvedir

;; DELVERED exposure status

;; Defaults
if n_elements(delvedir) gt 0 then delvedir=trailingslash(delvedir) else delvedir = '/dl1/users/dnidever/delve/'
;; Exposures directory
expdir = trailingslash(delvedir)+'exposures/'

;; Get all of the nights
dirs = FILE_SEARCH(expdir+'20??????',/test_directory,count=ndirs)
dirs = FILE_BASENAME(dirs)
numdirs = long(dirs)
print,strtrim(ndirs,2),' DELVE nights'

;; Night loop
str = replicate({night:'',fieldsfile:-1,setupfile:-1,exposurefile:-1,nightsumfile:-1,logsdir:-1,wcs_success:-1L,wcs_failure:0L,$
                 daophot_success:-1L,daophot_failure:0L,match_success:-1L,match_failure:0L,apcor_success:-1L,apcor_failure:0L,$
                 astrom_success:-1L,astrom_failure:0L,zeropoint_success:-1L,zeropoint_failure:0L,calib_success:-1L,calib_failure:0L,$
                 combine_success:-1L,combine_failure:0L,deredden_success:-1L,deredden_failure:0L,save_success:-1L,save_failure:0L,$
                 success:lonarr(10)-1,failure:lonarr(10)},ndirs)
tags = tag_names(str)
stages = ['WCS','DAOPHOT','MATCH','APCOR','ASTROM','ZEROPOINT','CALIB','COMBINE','DEREDDEN','SAVE']
nstages = n_elements(stages)
ngivesum = 0
For i=0,ndirs-1 do begin
  inight = dirs[i]
  idir = expdir+inight+'/'
  str[i].night = inight
  str[i].fieldsfile = file_test(idir+'fields')
  str[i].setupfile = file_test(idir+'photred.setup')
  str[i].exposurefile = file_test(idir+inight+'_exposures.fits')
  str[i].nightsumfile = file_test(idir+inight+'_summary.fits')
  str[i].logsdir = file_test(idir+'logs/',/directory)
  if str[i].logsdir eq 1 then begin
    for j=0,nstages-1 do begin
      sind = where(tags eq stages[j]+'_SUCCESS',nsind)
      if file_test(idir+'/logs/'+stages[j]+'.success') then str[i].(sind[0])=file_lines(idir+'logs/'+stages[j]+'.success')
      find = where(tags eq stages[j]+'_FAILURE',nfind)
      if file_test(idir+'/logs/'+stages[j]+'.failure') then str[i].(find[0])=file_lines(idir+'logs/'+stages[j]+'.failure')
      str[i].success[j] = str[i].(sind[0])
      str[i].failure[j] = str[i].(find[0])
    endfor
  endif

  if str[i].logsdir eq 1 and min(str[i].success) gt 0 then begin
    if ngivesum mod 40 eq 0 then print,'NIGHT       WCS    DAOPHOT   MATCH    APCOR   ASTROM   ZEROPT  CALIB   COMBINE  DERED   SAVE  NIGHTSUM    DONE'

    comment = ''
    done = 0
    if min(str[i].failure) eq 0 and max(str[i].failure) eq 0 and min(str[i].success) gt 0 and $
       str[i].fieldsfile eq 1 and str[i].setupfile eq 1 and str[i].exposurefile eq 1 and str[i].nightsumfile eq 1 then done=1
    if done eq 1 then comment = '       FINISHED'

    format = '(A-10, I4,A1,I-4, I4,A1,I-4, I4,A1,I-4, I4,A1,I-4, I4,A1,I-4, I3,A1,I-3, I4,A1,I-4, I4,A1,I-4, I3,A1,I-3, I3,A1,I-3, I5,A-15)'
    print,inight,str[i].wcs_success,'/',str[i].wcs_failure,str[i].daophot_success,'/',str[i].daophot_failure,$
          str[i].match_success,'/',str[i].match_failure,str[i].apcor_success,'/',str[i].apcor_failure,$
          str[i].astrom_success,'/',str[i].astrom_failure,str[i].zeropoint_success,'/',str[i].zeropoint_failure,$
          str[i].calib_success,'/',str[i].calib_failure,str[i].combine_success,'/',str[i].combine_failure,$
          str[i].deredden_success,'/',str[i].deredden_failure,str[i].save_success,'/',str[i].save_failure,$
          str[i].nightsumfile,comment,format=format
    ngivesum++
  endif  
Endfor

end
