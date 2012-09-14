PRO JWST_TRAINING
;======================================================
;+
; NAME: 					JWST_TRAINING
;
; DESCRIPTION:    IDL procedures for completing JWST exercises.
;
; INPUTS:
;
; KEYWORDS:
;
; AUTHOR: 				Derek Hammer (STScI)
;
; HISTORY:				Written Sep 5, 2012
;-
;======================================================

imarr=readfits('MIRI_VM2T00003582_1_IM_S_2008-09-14T11h02m41.fits')
imarr=double(imarr)
hdr=headfits('MIRI_VM2T00003582_1_IM_S_2008-09-14T11h02m41.fits',exten=0)


; PROBLEM 2A: Plot values of 3 pixels over each frame.
xvals=[2,376,282]
yvals=[657,693,441]
tframe = fxpar(hdr,'TFRAME')    ;get exposure time of each frame
time=(findgen(n_elements(imarr(0,0,*)))+1.0)*float(tframe)

; print results to file
radec2file, time, imarr(xvals(0),yvals(0),*), ofile='pixel1_vs_time.dat'
radec2file, time, imarr(xvals(1),yvals(1),*), ofile='pixel2_vs_time.dat'
radec2file, time, imarr(xvals(2),yvals(2),*), ofile='pixel3_vs_time.dat'


; PROBLEM 2B: Measure the slope (cnts vs. time) for every pixel in image. Create slope image.
;--Subtract a reference frame (1st frame) before measuring slope.
base=imarr(*,*,0)   ;init base image
slopim=fltarr(n_elements(imarr(*,0,0)),n_elements(imarr(0,*,0)))  ;init slope image

for x=0L,n_elements(imarr(*,0,0))-1 do begin
  for y=0L, n_elements(imarr(0,*,0))-1 do begin
    tmpy=float(reform(imarr(x,y,*)-base(x,y)))
    ;tmpfit=linfit(time(1:n_elements(time)-2),tmpy(1:n_elements(tmpy)-2),yfit=yfit,/double)
    tmpfit=linfit(findgen(40),tmpy(0:n_elements(tmpy)-1),yfit=yfit,/double)
    slopim(x,y) = tmpfit(1)
  endfor
endfor

;stats w/o clipping
ave0=avg(slopim)
sig0=sigma(slopim)
med0=median(slopim)

;3-sigma clipped values
meanclip, slopim, ave, sig, subs=subs
med=median(slopim(subs))


print, 'Median = '+strtrim(string(med0),2)
print, 'Median (3-sigma clipped) = '+strtrim(string(med),2)
print
print, 'Average = '+strtrim(string(ave0),2)+'  +/- '+strtrim(string(sig0),2)
print, 'Average (3-sigma clipped) = '+strtrim(string(ave),2)+'  +/- '+strtrim(string(sig),2)
print


cghistoplot, slopim, bin=1, xr=[-10, 25], background=cgcolor('white'), axiscolor=cgcolor('black'), $
  datacolor=cgcolor('red'), xtitle='Slope [cnts/frame]', ytitle='Number of pixels', output='slope_histo.pdf'

; /LSCALE keyword is necessary to save image w/o weird 32-bit number limitations
mwrfits, slopim, 'slope.fits', hdr, /create, /LSCALE


END
