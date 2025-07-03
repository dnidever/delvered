pro delvered_coverage,pix,version,redo=redo

; Make the coverage map for a single nside=128 NSC healpix

; Combine all of the data
;NSC_ROOTDIRS,dldir,mssdir,localdir
;dir = dldir+'users/dnidever/nsc/instcal/'+version+'/'
;dir = '/net/dl2/dnidever/nsc/instcal/'+version+'/'
dir = '/net/dl2/dnidever/delve/'
nside = 128
nside2 = 4096
radeg = 180.0d0 / !dpi


; Not enough inputs
if n_elements(pix) eq 0 then begin
  print,'Syntax - delvered_coverage,pix,redo=redo,version=version'
  return
endif

; Does the coverage map already exist
covfile = dir+'coverage/'+strtrim(long(pix)/1000,2)+'/'+strtrim(pix,2)+'_coverage.fits'
;covfile = dir+'combine/coverage/'+strtrim(pix,2)+'_coverage.fits'
if file_test(covfile) eq 1 and not keyword_set(redo) then begin
  print,covfile,' EXISTS and /redo NOT set.'
  return
endif

; Get healpix boundary coordinates
PIX2VEC_RING,nside,pix,vec,vertex
vertex = transpose(reform(vertex))  ; [1,3,4] -> [4,3]
VEC2ANG,vec,hcendec,hcenra,/astro
VEC2ANG,vertex,hdec,hra,/astro

; Rotate to tangent plane
ROTSPHCEN,hra,hdec,hcenra,hcendec,hlon,hlat,/gnomic
mmhlon = minmax(hlon)
mmhlat = minmax(hlat)

; Get the pixel numbers for nside=4096 healpix that are within
; this larger pixel
QUERY_POLYGON,nside2,vertex,listpix,nlistpix

; Check that they belong inside this healpix
PIX2ANG_RING,nside2,listpix,listtheta,listphi
ANG2PIX_RING,nside,listtheta,listphi,listpix1
bdpix = where(listpix1 ne pix,nbdpix)

; Not 1024 pixels, use python
if (nlistpix ne 1024) or (nbdpix gt 0) then begin
  print,'Not 1024 healpix.  Using python.'
  listpix0 = listpix
  ;; Use python code to get 4096 pixels
  tempfile = MKTEMP('hp')
  file_delete,tempfile+'.fits',/allow
  step = 100
  pylines = 'import healpy as hp; from astropy.io import fits; import numpy as np;'+$
            'v=hp.pix2vec('+strtrim(nside,2)+','+strtrim(pix,2)+');'+$
            'radius=hp.nside2resol('+strtrim(nside,2)+');'+$
            'pix2=hp.query_disc('+strtrim(nside2,2)+',v,radius=2*radius);'+$
            'theta,phi=hp.pix2ang('+strtrim(nside2,2)+',pix2);'+$
            'pix1=hp.ang2pix('+strtrim(nside,2)+',theta,phi);'+$
            'gd = (pix1 == '+strtrim(pix,2)+');'+$
            'pix = pix2[gd];'+$
            "fits.writeto('"+tempfile+".fits'"+',pix)'
  writeline,tempfile,pylines
  file_chmod,tempfile,'755'o
  spawn,['python',tempfile],out,errout,/noshell
  listpix = MRDFITS(tempfile+'.fits',0,/silent)
  file_delete,[tempfile,tempfile+'.fits'],/allow
  nlistpix = n_elements(listpix)
endif

; Initialize the coverage structure
print,'Creating coverage structure for pixel ',strtrim(pix,2)
covstr = replicate({pix:0L,pix128:long(pix),ra:0.0d0,dec:0.0d0,nobj:0L,coverage:0.0,nexp:0,ucoverage:0.0,unexp:0,udepth:-9999.0,gcoverage:0.0,gnexp:0,gdepth:-9999.0,$
                    rcoverage:0.0,rnexp:0,rdepth:0.0,icoverage:0.0,inexp:0,idepth:-9999.0,zcoverage:0.0,$
                    znexp:0,zdepth:-9999.0,ycoverage:0.0,ynexp:0,ydepth:-9999.0,vrcoverage:0.0,vrnexp:0,vrdepth:-9999.0},nlistpix)
covstr.pix = listpix
PIX2ANG_RING,nside2,covstr.pix,theta,phi
covstr.ra = phi*radeg
covstr.dec = 90-theta*radeg
covtags = tag_names(covstr)

; Find the overlapping DELVE-MC bricks
bstr = MRDFITS('/home/dnidever/projects/delvered/data/delvemc_bricks_0.25deg.fits.gz',1)
dist = sphdist(bstr.ra,bstr.dec,covstr.ra,covstr.dec,/deg)

stop

; Does the combined object file exist?
objfile = dir+'combine/'+strtrim(long(pix)/1000,2)+'/'+strtrim(pix,2)+'.fits.gz'
;objfile = dir+'combine/'+strtrim(pix,2)+'.fits.gz'
if file_test(objfile) eq 0 then begin
  print,objfile,' NOT FOUND'
  goto,SAVEFILE
endif

;; Check that it's not blank
head0 = headfits(objfile,exten=1,errmsg=errmsg)
if errmsg ne '' then begin
  print,objfile,' IS BLANK'
  goto,SAVEFILE
endif

; Load the list of exposures
expstr = MRDFITS(objfile,1,/silent)
nexpstr = n_elements(expstr)
expstr.file = strtrim(expstr.file,2)
expstr.base = strtrim(expstr.base,2)
expstr.dateobs = strtrim(expstr.dateobs,2)
expstr.filter = strtrim(expstr.filter,2)
add_tag,expstr,'success',0,expstr
add_tag,expstr,'hoverlap',0,expstr
add_tag,expstr,'chipindx',-1LL,expstr

; Load the object table
obj = MRDFITS(objfile,2,/silent)
ophi = obj.ra/radeg
otheta = (90-obj.dec)/radeg 
ANG2PIX_RING,nside2,otheta,ophi,obj_pix4096

; Load the chip summary structure for each exposure
print,'Loading the chip summary information'
print,' Number                Exposure          Nchips  Noverlap'
nchstr = 0LL
cnt = 0LL
for i=0,nexpstr-1 do begin
  ; Construct metadata filename
  base = expstr[i].base
  instrument = 'c4d'
  if stregex(base,'k4m',/boolean) eq 1 then instrument='k4m'
  if stregex(base,'ksb',/boolean) eq 1 then instrument='ksb'
  dateobs = expstr[i].dateobs
  night = strmid(dateobs,0,4)+strmid(dateobs,5,2)+strmid(dateobs,8,2)
  metafile = dir+instrument+'/'+night+'/'+base+'/'+base+'_meta.fits'

  undefine,chstr1
  if file_test(metafile) eq 1 then begin
    ; Load the chip summary structure
    chstr1 = MRDFITS(metafile,2,/silent)
    nchstr1 = n_elements(chstr1)
    chstr1.expdir = strtrim(chstr1.expdir,2)
    chstr1.filename = strtrim(chstr1.filename,2)
    add_tag,chstr1,'hoverlap',0,chstr1   ; does it overlap healpix
    add_tag,chstr1,'filter',strtrim(expstr[i].filter,2),chstr1
    expstr[i].success = 1

    ; Check if the chips overlap the healpix
    for j=0,nchstr1-1 do begin
      ROTSPHCEN,chstr1[j].vra,chstr1[j].vdec,hcenra,hcendec,lon,lat,/gnomic      
      chstr1[j].hoverlap = DOPOLYGONSOVERLAP(hlon,hlat,lon,lat)      
    endfor
    expstr[i].hoverlap = max(chstr1.hoverlap)

    ; Only keep overlapping chips
    gdchstr1 = where(chstr1.hoverlap eq 1,ngdchstr1,comp=bdchstr1,ncomp=nbdchstr1)
    print,i+1,base,nchstr1,ngdchstr1,format='(I5,A35,I6,I6)'
    if ngdchstr1 eq 0 then begin
      print,'No chips overlap this healpix'
      expstr[i].hoverlap = 0
      expstr[i].nchips = 0
      expstr[i].chipindx = -1
      goto,BOMB
    endif
    ; Remove non-overlapping chips
    if nbdchstr1 gt 0 then begin
      REMOVE,bdchstr1,chstr1
      nchstr1 = n_elements(chstr1)
    endif
    expstr[i].chipindx = cnt
    expstr[i].nchips = nchstr1

    ; Start CHSTR structure
    if nchstr eq 0 then begin
      schema_chstr = chstr1[0]
      struct_assign,{dum:''},schema_chstr
      chstr = REPLICATE(schema_chstr,10000L)
      nchstr = n_elements(chstr)
    endif

    ; Add new elements
    if cnt+nchstr1 gt nchstr then begin
      old = chstr
      chstr = REPLICATE(schema_chstr,nchstr+10000L)
      chstr[0:nchstr-1] = old
      nchstr = n_elements(chstr)
      undefine,old
    endif

    ; Stuff into the CHSTR structure
    newchstr1 = REPLICATE(schema_chstr,nchstr1)
    struct_assign,chstr1,newchstr1
    chstr[cnt:cnt+nchstr1-1] = newchstr1  
    cnt += nchstr1

  ; Metadata file not found
  endif else begin
    print,metafile,' NOT FOUND'
  endelse
  BOMB:
endfor
; Trim off the extra elements
chstr = chstr[0:cnt-1]
nchstr = n_elements(chstr)

; Trim off any exposures with no overlapping chips
expstr0 = expstr
bdexp = where(expstr.hoverlap eq 0,nbdexp)
if nbdexp gt 0 then REMOVE,bdexp,expstr


filters = ['u','g','r','i','z','Y','VR']
nfilters = n_elements(filters)

; Loop over the small pixels and figure out the coverage and depth
; and number of exposures
for i=0,nlistpix-1 do begin

  ; Get the objects for this pixel
  MATCH,listpix[i],obj_pix4096,pind1,pind2,/sort,count=nobjmatch
  covstr[i].nobj = nobjmatch

  ; Get healpix boundary coordinates
  PIX2VEC_RING,nside2,listpix[i],vec1,vertex1
  vertex1 = transpose(reform(vertex1))  ; [1,3,4] -> [4,3]
  VEC2ANG,vec1,hcendec1,hcenra1,/astro
  VEC2ANG,vertex1,hdec1,hra1,/astro

  ; Rotate to tangent plane
  ROTSPHCEN,hra1,hdec1,hcenra1,hcendec1,hlon1,hlat1,/gnomic

  dx = 1d-3  ; gives ~20x20 pixel
  rlon = minmax(hlon1)
  rlat = minmax(hlat1)
  lon0 = rlon[0]
  lat0 = rlat[0]
  hx = (hlon1-lon0)/dx
  hy = (hlat1-lat0)/dx
  nx = ceil((rlon[1]-lon0)/dx)+1
  ny = ceil((rlat[1]-lat0)/dx)+1

  ; Mask image for which pixels are in the healpix region
  mask = bytarr(nx,ny)
  in = polyfillv(hx,hy,nx,ny)
  mask[in] = 1
  maskpix = total(mask)

  ; Loop over each filter
  for f=0,nfilters-1 do begin
    filtind = where(chstr.filter eq filters[f],nfiltind)

    ; Get columns for this filter
    covind = where(covtags eq strupcase(filters[f])+'COVERAGE',ncovind)
    numind = where(covtags eq strupcase(filters[f])+'NEXP',nnumind)
    depind = where(covtags eq strupcase(filters[f])+'DEPTH',ndepind)

    ; Coverage map for this filter
    numim = bytarr(nx,ny)
    depthim = fltarr(nx,ny)

    ; Now loop over each chip
    alloverlap = 0
    for c=0,nfiltind-1 do begin
      j = filtind[c]
      ROTSPHCEN,chstr[j].vra,chstr[j].vdec,hcenra1,hcendec1,lon1,lat1,/gnomic
      ; Check if they overlap
      overlap = DOPOLYGONSOVERLAP(hlon1,hlat1,lon1,lat1)
      if overlap eq 1 then begin
        ; Get the chip overlap image
        ; Transform to pixel-based tangent plane system
        vx = (lon1-lon0)/dx
        vy = (lat1-lat0)/dx
        ; Pixels in the chip
        cin = POLYFILLV(vx,vy,nx,ny)
        cmask = bytarr(nx,ny)
        if n_elements(cin) gt 1 or cin[0] ne -1 then cmask[cin]=1  ; some good overlap
        cmask = cmask * mask   ; pixels in healpix region

        numim[*,*] += cmask               ; number of chips that overlap
        depthim[*,*] += cmask*chstr[j].depth95  ; sum of depth of all pixels
      endif  ; overlap
    endfor  ; chip loop

    ; Calculate coverage
    gdpix = where(numim gt 0,ngdpix)
    overlapfrac = float(ngdpix) / maskpix
    ; Calculate average depth image
    mndepthim = depthim/(numim > 1)
    if ngdpix gt 0 then depth = median(mndepthim[gdpix]) else depth=-9999.0
    ; Number of chips
    nchipoverlap = max(numim)

    ; Stuff into the coverage structure
    covstr[i].(covind) = overlapfrac
    covstr[i].(depind) = depth
    covstr[i].(numind) = nchipoverlap

    ; Add to total coverage and exposure for this pixel
    covstr[i].coverage >= overlapfrac
    covstr[i].nexp += nchipoverlap
  endfor  ; filter loop

endfor  ; pixel loop

; Save the coverage map
SAVEFILE:
print,'Writing coverage information to ',covfile
if file_test(file_dirname(covfile),/directory) eq 0 then file_mkdir,file_dirname(covfile)
MWRFITS,covstr,covfile,/create

;cat = mrdfits(objfile,2,/silent)
;plotc,covstr.ra,covstr.dec,covstr.coverage,ps=1,/ysty
;for j=0,nchstr-1 do oplot,[chstr[j].vra,chstr[j].vra[0]],[chstr[j].vdec,chstr[j].vdec[0]],ps=-1   
;oplot,cat.ra,cat.dec,ps=3,co=200
;
;g = where(covstr.zdepth gt 0,ng)
;mn = min(covstr[g].zdepth)-0.2
;mx = max(covstr[g].zdepth)+0.2
;plotc,covstr.ra,covstr.dec,covstr.zdepth,ps=1,/ysty,min=mn,max=mx
;for j=0,nchstr-1 do oplot,[chstr[j].vra,chstr[j].vra[0]],[chstr[j].vdec,chstr[j].vdec[0]],ps=-1   
;oplot,cat.ra,cat.dec,ps=3,co=200

;stop

end
