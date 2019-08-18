pro make_exposures_list,version

;; Make list of DECAM MC exposures of DELVE+community data
;; using the NSC master decam list

if n_elements(version) eq 0 then version = 'v3'
nscdir = '/dl1/users/dnidever/nsc/instcal/'+version+'/lists/'
if n_elements(delvereddir) gt 0 then delvereddir=trailingslash(delvereddir) else delvereddir = '/home/dnidever/projects/delvered/'
str = mrdfits(nscdir+'decam_instcal_list.fits.gz',1)
nstr = n_elements(str)
print,strtrim(nstr,2),' total public DECam exposures'
glactc,str.ra,str.dec,2000.0,glon,glat,1,/deg
cel2lmc,str.ra,str.dec,lmcpa,lmcrad
cel2smc,str.ra,str.dec,smcpa,smcrad
gd = where(abs(glat) gt 10 and (lmcrad lt 25 or smcrad lt 15),ngd)
mc = str[gd]
print,strtrim(ngd,2),' MC exposures'
jd = systime(/julian)
caldat,jd,month,day,year,hour
outfile = delvereddir+'data/decam_mcs_'+string(year,format='(i04)')+string(month,format='(i02)')+string(day,format='(i02)')+'.fits'
print,'Writing to ',outfile
MWRFITS,mc,outfile,/create
if file_test(outfile+'.gz') eq 1 then file_delete,outfile+'.gz'
spawn,['gzip',outfile],/noshell

;stop

end
