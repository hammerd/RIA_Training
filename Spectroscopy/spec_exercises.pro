;======================================================
;+
; NAME: 					SpecTraining
;
; DESCRIPTION:    Routines to perform exercises in training manual.
;
; INPUTS:
;
; KEYWORDS:
;
; AUTHOR: 				Derek Hammer (STScI)
;
; HISTORY:				Written Aug 17, 2012
;-
;======================================================
PRO Spec_Exercises


path='./data/exercises/'

;======================================================
; Exercise 1: Open "x1d" files from lbgu17qnq and o8k401010. What kind of data is stored,
; and what are the general similarities/differences between both instruments?
;======================================================
scos=mrdfits(path+'lbgu17qnq_x1d.fits',1,hc)
stis=mrdfits(path+'o8k401010_x1d.fits',1,hs)
splot, stis.wavelength, stis.flux


;======================================================
; Exercise 2: Open the raw files. What kind of data does each contain?
; Do the COS and STIS files contain the same type of raw data? Why?
; STIS=2-d images; COS=tables
;======================================================
scosa_raw=mrdfits(path+'lbgu17qnq_rawtag_a.fits',1,hc)
scosb_raw=mrdfits(path+'lbgu17qnq_rawtag_b.fits',1,hc)
rstis=mrdfits(path+'o8k401010_raw.fits',1,hs)


;======================================================
; Exercise 3:Compare reprocessed COS/STIS spectra to archive. Are they identical?
;---RESULT: STIS is matched perfectly. COS is mismatched probably b/c I'm using a 
;           more recent calcos version than the archive.
;======================================================
rscos=mrdfits(path+'cos_reprocess/lbgu17qnq_x1d.fits',1,hrc)
rstis=mrdfits(path+'stis_reprocess/o8k401010_x1d.fits',1,hrs)

;---STIS comparison
diffstis=double(stis.flux - rstis.flux)             ;STIS datasets are identical
bad=where(abs(diffstis)/double(stis.flux) gt 1d-8)  ;none found

;---COS comparison
diffcosa=double(scos(0).flux) - double(rscos(0).flux)
diffcosb=double(scos(1).flux) - double(rscos(1).flux)
bada=where(abs(diffcosa)/double(scos(0).flux) gt 1d-8,nbada)
badb=where(abs(diffcosb)/double(scos(1).flux) gt 1d-8,nbadb)  ;MAJORITY ARE UNMATCHED!!

;---DIFFERENCES BETWEEN HEADER FILES
;randomseed=1345217950L             ; this is listed in archive header, 1345487487L in reprocessed.
;gsagtab='lref$w4h183811_gsag.fits' ; this tag is listed in reprocessed header, 'N/A' in archive.
tt=mrdfits(path+'cos_reprocess/lbgu17qnq_x1d.fits',0,htt)
ss=mrdfits(path+'lbgu17qnq_x1d.fits',0,hss)
tmp=lonarr(n_elements(hss))
for i=0, n_elements(hss)-1 do begin & tmp(i)=strcmp(strtrim(htt(i),2), strtrim(hss(i),2)) & endfor
bad=where(tmp ne 1,nbad)
for i=0, nbad-1 do begin & print & print, strtrim(htt(bad(i)),2), '   ', strtrim(hss(bad(i)),2) & endfor


;======================================================
; Exercise 4: Turn off FLATCORR calibration switch and reprocess images. Plot the difference.
;======================================================
rscos1=mrdfits(path+'cos_reprocess1/lbgu17qnq_x1d.fits',1,hrc1)
plot, rscos(0).wavelength, (rscos1(0).flux - rscos(0).flux)/rscos(0).flux, yr=[-0.01,0.01],ys=1, $
 xtitle='wavelength ('+tex2IDL('\AA')+')', ytitle='Flux ratio [ (Flux_Noflat-Flux_Flat)/Flux_Flat]'


;======================================================
; Exercise 5: Convolve the COS spectra with a boxcar kernel. Try different sizes to find the best match
; to the STIS spectra.
;======================================================
plot, stis.wavelength, stis.flux, xr=[1300,1430],xs=1, yr=[0,7d-14]
oplot, scos(0).wavelength, smooth(scos(0).flux,100)*1.00, color=cgcolor('green')  
oplot, scos(0).wavelength, smooth(scos(0).flux,200)*1.25, color=cgcolor('blue')   ; data are shifted slightly in y-direction for clarity
oplot, scos(0).wavelength, smooth(scos(0).flux,300)*1.50, color=cgcolor('red')


;======================================================
; Exercise 6: Estimate the S/N of the COS and STIS continuum (use low-order polynomial fit).
;---first identify regions of the continuum that lack absorption/emission lines.
;======================================================
plot, stis.wavelength, stis.flux, xr=[1300,1430],xs=1, yr=[0,7d-14]
oplot, scos(0).wavelength, smooth(scos(0).flux,200),color=cgcolor('blue')
gds=where(stis.wavelength ge 1372. and stis.wavelength le 1408.,ngds)
gdc=where(scos(0).wavelength ge 1372. and scos(0).wavelength le 1408., ngdc)

;---fit polynomials to the spectra
sfit=poly_fit(stis.wavelength(gds), stis.flux(gds),3, yfit=syfit)
oplot, stis.wavelength(gds), syfit, color=cgcolor('cyan')
cfit=poly_fit(scos(0).wavelength(gdc), smooth(scos(0).flux(gdc),51),3, yfit=cyfit)
oplot, scos(0).wavelength(gdc), cyfit, color=cgcolor('magenta')

;---calculate residuals and 3-sigma clipped errors
meanclip, stis.flux(gds) - syfit, smean, ssig
print, avg(syfit/ssig)  ; S/N=30
meanclip, scos(0).flux(gdc) - cyfit, cmean, csig
print, avg(cyfit/csig)  ; S/N=4  (also tried a very small region and got the same results: lam=1372-1377)


;======================================================
; Exercise 7: Perform linear interpolation of STIS spectra onto COS scale. 
;---Take the flux ratio, and see where they differ.
;======================================================
lstis=interpol(stis.flux,stis.wavelength, scos(0).wavelength)
meanclip, lstis/scos(0).flux, amean, asig ; RESULT:  mean=1.013  sig=0.592   ~1.3% increase.

;--plot COS spectra plus interpolated STIS spectra
op1=Obj_New('cgOverPlot', [1290,1430],[1,1],Color=cgcolor('red'),LINESTYLE=2, thick=3)
cgplot, scos(0).wavelength, lstis/scos(0).flux, yr=[0,14],xr=[1290,1430],xs=1,ys=1,color=cgcolor('black'), $
  output='STISinterp_vsCOS.pdf', xtitle='wavelength ('+tex2IDL('$\AA$')+')', ytitle='flux ratio (STIS_interp/COS)', $
  oplots=[op1]


;======================================================
; Exercise 8: Perform fit to Ly-alpha emission line in spectra (use COS B).
;======================================================
gdp=where(scos(1).wavelength ge 1218. and scos(1).wavelength le 1252., ngdc)  ;wavelen range for plot
radec2file, scos(1).wavelength(gdp), scos(1).flux(gdp), ofile='Lya_spec.dat'

;--perform Gaussian fit and record results to file
foff=0.5d-13
gdf=where(scos(1).wavelength ge 1218. and scos(1).wavelength le 1252., ngdc)              ;define wavelength range for fit
dd=gaussfit(scos(1).wavelength(gdf), smooth(scos(1).flux(gdf),3)-foff, coeffs, nterms=3)  ;perform gaussian fit
radec2file, scos(1).wavelength(gdf), dd+foff, ofile='Lya_fit.dat'                         ;save results to file

;--plot the spectrum+fit
op2=Obj_New('cgOverPlot',scos(1).wavelength(gdf), dd+foff, color='red',LINESTYLE=2, thick=2)        ;IDL object file (for oplot)
cgplot, scos(1).wavelength, scos(1).flux,axiscolor=cgcolor('black'),background=cgcolor('white'), $
xr=[1218,1252],xs=1, xtitle='wavelength ('+tex2IDL('$\AA$')+')', oplots=op2, $
ytitle='Flux (erg sec!e-1!n cm!e-2!n !3 '+string(197b)+'!x!e-1!n)'



;======================================================
; ASSIGNMENT 3.1
;======================================================
; SEE PYTHON CODE FOR ASSIGNMENT #1





END
