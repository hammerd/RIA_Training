;======================================================
;+
; NAME: 					IDL_TRAINING
;
; DESCRIPTION:    IDL routines written to address training exercises.
;
; INPUTS:         NA
;
; KEYWORDS:
;     EXT_APPS    flag to call external functions inside this program (e.g., ds9, fv).
;                 The apps are rather annoying and require that you close them in order to continue,
;                 so the default is to NOT allow calls to external apps.
;                  
; WARNINGS:       Assumes all external input files are in the working IDL directory.
; 
;                 THIS PROCEDURE WAS NOT WRITTEN TO BE RUN FROM BEGINNING TO END AND
;                 PROVIDE AESTHETIC FORMATTED OUTPUT FOR EACH EXERCISE. THE PRIMARY PURPOSE 
;                 WAS TO DOCUMENT THE PROCEDURE FOR EACH EXERCISE, WHICH MAY THEN BE COPIED 
;                 AND PASTED INTO AN IDL SESSION. THAT BEING SAID, IT WILL RUN FROM BEGINNING
;                 TO END WITHOUT ERRORS, BUT THE OUTPUT FORMAT WILL BE FUNKY AND PLOTTING
;                 WINDOWS MAY HAVE WEIRD COLOR TABLE ISSUES.
;
; AUTHOR: 				Derek Hammer (STScI)
;
; HISTORY:				Written Aug 13, 2012
;-
;======================================================
PRO IDL_TRAINING, EXT_APPS=ext_apps

;close all IDL windows
nwin=!d.window
while(nwin ge 0) do begin 
  wdelete, nwin
  nwin=!d.window
endwhile
  

; Exercise 2.0 - What's the difference between a/b when a=3, and b=2 and 2.0?
; Answer: integer division gives integer quotient (floor), while float gives float quotient. 
print, 3/2, 3/2.


; Exercise 2.1  - Create the sequence 0.1, 0.2, 0.3...1.4 using findgen.
print, findgen(14)/10. + 0.1



; Exercise 2.2:  Create the sequence -3.2, -3.0, -2.8…-1.0 using findgen.
print, findgen(12)/5. - 3.2


; Exercise 2.3:  Create a 2x3 array with 10 in all places.
print, fltarr(2,3) > 10.0


; Exercise 2.4: Create two 3x3 different arrays (A and B). Calculate A+B,A/B,A*B
A=[[5.0,-1.0,3.0], [2.0,0.0,1.0],[3.0,2.0,1.0]]
B=fltarr(3,3) > 2.0
print, A+B
print, A/B
print, A*B



; Exercise 2.5: Using array A, investigate how to obtain its inverse. Verify that AxA^-1 = I.
print, A#invert(A)

 
 
; Exercise 2.6: Create a random real 1000 element array. Use WHERE to create another array
; with those elements with counts lower than 0.5 and higher than -0.5.
rando=randomn(12.,1000)
good=where(abs(rando) lt 0.5,ngood)


; Exercise 2.7: Given the string "j12345678_flt.fits", use STRMID to extract rootname of file.
fname='j12345678_flt.fits'
print, strmid(fname, 0, strpos(fname, '.'))


; Exercise 2.8: Use the command SORT to create array [7,4,3,2,1] from array "a" in example.
a=[4,3,7,1,2]
print, a(reverse(sort(a)))
       

; Exercise 2.9: Given the array a=randomn(10,1000), what are the moments of a?
a=randomn(10,1000)
print, moment(a)
   
 
; Exercise 2.10: Create your own input file with 3 columns and at least 4 rows.
; Develop a code using READCOL to import the columns into IDL vectors.
readcol, 'nuvlf_binsz0.5.dat', mag, lf, lferru, lferrl
print, mag, lf, lferru, lferrl


; Exercise 2.11: Use FORPRINT to write a file with the columns from previous exercise
; but inverting the column order.
forprint, lferrl, lferru, lf, mag, textout='nuvlf_binsz0.5_invertedcol.dat'


; Exercise 2.12: Use SPAWN in an IDL script that opens ds9 and loads 2 FITS files in different frames.
if keyword_set(ext_apps) then spawn, 'ds9 -tile -frame 1 m101_blue.fits -frame 2 j91c12biq_flt.fits'


; Exercise 2.13: Open all 6 extensions of j91c12biq_flt.fits using -multiframe flag from IDL.
; Can you make the images blink?
if keyword_set(ext_apps) then spawn, 'ds9 -multiframe j91c12biq_flt.fits'
if keyword_set(ext_apps) then spawn, 'ds9 -multiframe j91c12biq_flt.fits -blink'


; Exercise 2.14(a): What is the structure of the NICMOS file n8ws10siq_ima.fits?
; Answer: 76 extensions. One primary header plus 5 extensions [SCI,ERR,DQ,SAMP,TIME]
;         for 15 observed fields.
if keyword_set(ext_apps) then spawn, 'fv n8ws10siq_ima.fits&'


; Exercise 2.14(b): What is the size in pixels of the SCI extension?  A: 256
im=mrdfits('n8ws10siq_ima.fits',1,h)
print, size(im)


; Exercise 2.15: Save a 700x700 subarray for extension 4 of "ib6w71lxq_flt.fits."
im=mrdfits('ib6w71lxq_flt.fits',4,h)
subim=im(2000:2699,700:1399)
data = 255b - bytscl(subim,min=0.0,max=100.)
window, 0, xs=700,ys=700
cgimage, data,/axes, axkeywords={xtickname:strarr(8)+' ', ytickname: strarr(8)+' '}, output='ib6w71lxq_flt.ps'
wdelete

; Exercise 2.16: Read the primary header of the WFC3 image ib6w71xq_flt.fits
hdr=headfits('ib6w71lxq_flt.fits',exten=0)


; Exercise 2.17: Extract the filter, image bias correction(?), name of dark image, & gain in Amp D for above fits.
print, sxpar(hdr, 'FILTER'), sxpar(hdr, 'BIASCORR'), '  ', sxpar(hdr, 'DARKFILE'), sxpar(hdr, 'CCDOFSTD')


; Exercise 2.18: Create a vector with 10 values equal to the square of the index.
print, findgen(10)^2


; Exercise 2.19: Create the following sequence using a FOR...DO statement:2001-01-01, 2001-02-01, 2001-03-01, 2001-04-01
for i=0, 3 do print, '2001-0'+strtrim(string(i+1),2)+'-01'


; Exercise 2.20: Create IDL sequence that loops through all 4 WFPC2 images in uba3010em_c0f.fits and displays them to 
; the screen in order. Use CGIMAGE.
imarr=mrdfits('uba3010em_c0f.fits',0,h)
for i=0, n_elements(imarr(0,0,*))-1 do cgimage, bytscl(imarr(*,*,i),min=0.,max=30.)



; Exercise 2.21: Using MAST, download the calibrated assoc product ib3p11010 (ims ib3p11p7q,ib3p11p8q,ib3p11phq,ib3p11q9q).
; Loop through all images and prints the following info to screen: root name/target/filter/exposure time/obsdate.

fnames=file_search('ib3p11*flt*',count=nfiles)
print, ' ROOT NAME   TARGET NAME  FILTER    EXPOSURE TIME     DATE-OBS'
print, '--------------------------------------------------------------'
for i=0, nfiles-1 do begin
  hdr=headfits(fnames(i),exten=0)
  print, strtrim(sxpar(hdr,'ROOTNAME'),2), '      ', strtrim(sxpar(hdr,'TARGNAME'),2), '      ', $
    strtrim(sxpar(hdr,'FILTER'),2),'      ', strtrim(sxpar(hdr,'EXPTIME'),2), '      ', strtrim(sxpar(hdr,'DATE-OBS'),2) &
endfor



; Exercise 2.22: Given a string array of rootnames, identify HST camera according to rootname.
root=['j12345678', 'u12345678','n12345678','i12345678','l12345678']
for i=0, n_elements(root)-1 do $
  case strmid(root(i),0,1) OF $
     'j': print, root(i)+' = ACS' & $
     'u': print, root(i)+' = WFPC2' & $
     'n': print, root(i)+' = NICMOS' & $
     'i': print, root(i)+' = WFC3' & $
     'l': print, root(i)+' = COS' & $
     ELSE: print, root(i)+' = Not a valid identifier' & $
  endcase
  
  
;Exercise 2.23: Define a function that returns the area of a circle given the radius.
;  *** Function "Circle_Area" is defined at the end of this procedure ***
print, circle_area(2.)          ; area of circle of radius=2


; Exercise: 2.24: Write IDL function that returns the slope/intercept of a line given
; two data points in array format (x0,y0) and (x1,y1). Beware of vertical line case.
; *** FUNCTION "DATA2Line" IS DEFINED AT THE END OF THE FILE***
data1=[1.,1.]
data2=[1.,2.]
print, data2line(data1, data2)


; Exercise 3.1: Plot y=f(x)=exp(x)*x^-2 and create a PDF file. Select the appropriate range in x/y.
x=findgen(1000)/100. + 0.1
y=exp(x) * x^(-2)
cgplot, x, y, xr=[0,12],yr=[0,250],xs=1,ys=1,color=cgcolor('red'), background=cgcolor('white'), axiscolor=cgcolor('black'), xtitle='X', ytitle='Y', output='exercise3.1.pdf'


; Exercise 3.2: Plot mathematical functions of your choice in contigous plots.
wdelete
x=findgen(1000)/1000.*2.*!Pi
y1=sin(x)
y2=cos(x)
y3=y1*y2
cgplot, x, y1, xr=[0,!Pi*2.],yr=[-1.2,1.2],xs=1, ys=1, axiscolor=cgcolor('black'), background=cgcolor('white'), xtitle='X (radians)', ytitle='SIN (x)',/noerase, layout=[3,1,0]
oplot, x, fltarr(n_elements(x)), color=cgcolor('black'), linestyle=2
cgplot, x, y2, xr=[0,!Pi*2.],yr=[-1.2,1.2],xs=1, ys=1, axiscolor=cgcolor('black'), background=cgcolor('white'), ytitle='COS (x)', xtitle='X (radians)',/noerase, layout=[3,1,2]
oplot, x, fltarr(n_elements(x)), color=cgcolor('black'), linestyle=2
cgplot, x, y3, xr=[0,!Pi*2.],yr=[-1.2,1.2],xs=1, ys=1, axiscolor=cgcolor('black'), background=cgcolor('white'), ytitle='SIN (x) * COS (x)', xtitle='X (radians)',/noerase, layout=[3,1,3]
oplot, x, fltarr(n_elements(x)), color=cgcolor('black'), linestyle=2


; Exercise 3.3: Read the 1st column in ngc4214_336.dat and create a histogram.
; Choose appropriate binsize and do NOT fill polygons.
readcol, 'ngc4214_336.dat', f336
xmin=floor(min(f336))
xmax=ceil(max(f336))
binsz=0.2
cghistoplot, f336, backcolorname=cgcolor('white'), datacolorname=cgcolor('blue'), axiscolor=cgcolor('black'), binsize=binsz,xr=[floor(xmin),ceil(xmax)],xs=1,mininput=floor(xmin), maxinput=ceil(xmax), xtitle='F336 mag'


; Exercise 3.4: Create a plot that displays the same contour as in the example, but without the background fits image.
im = mrdfits('m101_blue.fits', 0, hdr)
im=im[300:800,300:800]
nlevels=4.
levelvals=(findgen(nlevels)+1.0)*(max(im)-min(im))/(nlevels+1.) + min(im)
cgcontour, im, findgen(n_elements(im(*,0))), findgen(n_elements(im(0,*))), levels=levelvals, axiscolor=cgcolor('black'),c_colors=[cgcolor('blue'), cgcolor('green'),cgcolor('red'), cgcolor('firebrick')], background=cgcolor('white')


; Exercise 4.1: Show how many extensions are in fits file.
im = mrdfits('j8hm01xaq_flt.fits', 0, hdr)
print, sxpar(hdr, 'NEXTEND')


; Exercise 4.2: When was the image created? Wo is the PI? Physical units of chip? Exposure time?
print, '  Image creation = '+strtrim(sxpar(hdr,'DATE'),2)
print, '  PI = '+strtrim(sxpar(hdr,'PR_INV_F'),2)+' '+strtrim(sxpar(hdr,'PR_INV_L'),2)
im = mrdfits('j8hm01xaq_flt.fits', 4, hdr1)
image500=im[900:1399,800:1299]
print, '  Physical Units of Chip = '+strtrim(sxpar(hdr1, 'BUNIT'),2)
print, '  Exposure Time = '+strtrim(sxpar(hdr, 'EXPTIME'),2)


; Exercise 4.3: Construct image/plot and overplot detected sources via FIND.
find, image500,xcoo,ycoo,flux,sharp, roundness, 40., 2.0, [-1.0,1.0], [0.2,1.0]
set_plot, 'ps'
page_height=27.94
page_width=21.59
plot_left=5.
plot_bottom=5.
xsize=14.
ysize=14.
device, filename='photometry_02.ps', xsize=page_width, ysize=page_height,xoffset=0.,yoffset=0.,scale_factor=1.0,/portrait
data = 255b - bytscl( image500, min = 1, max = 40)
ss = size(data, /dimensions)
cgloadct, 0
tvlct, redvector, greenvector, bluevector, /get
cgimage, data, position = [plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width,(plot_bottom + ysize) / page_height]
cgplot, [0], [0],xcharsize = 1, ycharsize = 1, thick = 2, xrange= [0, ss[0]],yrange= [0, ss[1]],xtitle = "x pixels",ytitle = "y pixels", xstyle = 1, ystyle = 1, /nodata,$
  /normal, /noerase, position = [plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height]
oplot, xcoo,ycoo,color=cgcolor('red'),psym=symcat(9)
device, /close
set_plot, 'X'


; Exercise 4.4: What is the meaning of the flux flag using APER? What happens if you omit /flux?
; Answer: the flux flag reports photometry in flux units instead of magnitude (default).

; Exercise 4.5: Write a code that outputs the flux for each star in the format: xcoo,ycoo,flux,flux error, sky,sky error
aper, image500,xcoo,ycoo,flux,eflux,sky,skyerr,1, [2., 5., 10., 15., 20.], [20,30],[-32767d,80000d], /flux, /silent
forprint, xcoo, ycoo, flux(1,*), eflux(1,*), sky, skyerr, textout=1

; Exercise 4.6: Plot the fluxerror as a function of flux. Use blue filled circles as symbols, and label plot.
cgplot, flux(1,*)/30., eflux(1,*)/30., background=cgcolor('white'), axiscolor=cgcolor('black'), symcolor=cgcolor('blue'), psym=symcat(16), xtitle='Flux (e-/sec)', ytitle='Flux Error (e-/sec)'



; Exercise 5.1: Can you show how many extensions the fits file has? What is the structure of the fits file?
im=mrdfits('ngc4151_hband.fits',0,hdr)
print, sxpar(hdr, 'NEXTEND')              ; answer=1+primary header
im=mrdfits('ngc4151_hband.fits',1,hdr1)   ; extension1 is a 3-D array: there are 2040 images of size 59x59 pixels.


; Exercise 5.2:  Repeat the same procedure with the other two dimensions in order to build 2 arrays: one for x, one for y.
;   In what units are these axes given?  Answer: RA/DEC (degrees)
crpix1=sxpar(hdr1, 'CRPIX1')
cdelt1=sxpar(hdr1, 'CDELT1')
crval1=sxpar(hdr1, 'CRVAL1')
cunit1=sxpar(hdr1, 'CUNIT1')
xval=dindgen(59) + 1.0
xpix=crval1 + (xval-crpix1)*cdelt1

crpix2=sxpar(hdr1, 'CRPIX2')
cdelt2=sxpar(hdr1, 'CDELT2')
crval2=sxpar(hdr1, 'CRVAL2')
cunit2=sxpar(hdr1, 'CUNIT2')
yval=dindgen(59) + 1.0
ypix=crval2 + (yval-crpix2)*cdelt2

crpix3=sxpar(hdr1, 'CRPIX3')
cdelt3=sxpar(hdr1, 'CDELT3')
crval3=sxpar(hdr1, 'CRVAL3')
cunit3=sxpar(hdr1, 'CUNIT3')
lval=dindgen(2040) + 1.0
lambda=crval3 + (lval-crpix3)*cdelt3


; Exercise 5.3: What is the wavelength range of this cube? 
; Repeat the procedure for another plane w/an interesting IR wavelength.
print, minmax(lambda)/1d4   ; 1.4767 -- 1.8029 um

ss = size(im, /dimensions)
contour_x=findgen(ss[0])
contour_y=findgen(ss[1])
plot_left=3.
plot_bottom=4.
page_height=27.94
page_width=21.59
xsize=15.
ysize=15.
cgloadct, 33, ncolors = 256, bottom = 0, clip = [0, 256], /reverse
tvlct, redvector, greenvector, bluevector, /get
set_plot, 'ps'

device, filename='ngc4151_1.65um.ps', xsize=page_width, ysize=page_height,xoffset=0.,yoffset=0.,scale_factor=1.0,/portrait
bytim=255b - bytscl(im(*,*,i), min=0., max=100.)
cgimage, bytim, position = [plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width, (plot_bottom + ysize) / page_height]
cgcontour, im[*,*,i], contour_x, contour_y,levels = [5, 10, 20, 30, 40, 50, 100, 200, 600],xs=1,ys=1,axiscolor='black',xcharsize=0.8, ycharsize=0.8, $
c_colors=cgcolor('white'), /noerase, position=[plot_left / page_width, plot_bottom / page_height, (plot_left + xsize) / page_width, (plot_bottom+ysize) / page_height]
device, /close
set_plot, 'X'


; Exercise 6.1: What values of the Hubble constant do you get? Can you determine error in H0?
readcol, 'hubble.dat', distance,vel
hfit=linfit(distance, vel, sigma=sig, yfit=yfit)
print, 'H0 = '+strtrim(string(hfit(1)),2)+' +/- '+strtrim(string(sig(1)),2)+' [km/s/Mpc]'


; Exercise 6.2: Plot empirical data and the fit
cgplot, distance, vel,color=cgcolor('blue'), background=cgcolor('white'),axiscolor=cgcolor('black'),psym=symcat(15),symsize=2, xtitle='Distance (Mpc)', ytitle='Velocity (km/sec)',xrange=[0,18],xs=1
oplot, distance, yfit, color=cgcolor('red'), thick=1.5


; Exercise 6.3: Calculate the FWHM for this line. 
spec = mrdfits('o5bn02010_x1d.fits', 1, hdr)
wave=spec.wavelength
flux=spec.flux
i=sort(wave)
swave=wave(i)
sflux=flux(i)
scaled_sflux=sflux / (1d-13)
wline = swave[34610: 34900]
fline = scaled_sflux[34610: 34900]
fit = gaussfit(wline, fline, coeff, nterms = 3)

fwhm=2.*sqrt(2.*alog(2))*coeff(2) ;answer=0.31 Ang


; Exercise 6.4: Produce IDL code to recreate Fig 6.2. Include vertical line to 
; show position of the spectral line shown in Figure 3.4.
cgplot, wline, fline, color=cgcolor('cadetblue'), background=cgcolor('white'), axiscolor=cgcolor('black'),xtitle='Wavelength (!3'+string(197b)+'!n)', $
  ytitle='Flux (10!e-13!n erg sec!e-1!n cm!e-2!n !3 '+string(197b)+'!x!e-1!n)',xrange=[1547.,1549.],xs=1,yrange=[-1,12],ys=1
oplot, wline, coeff(0)*exp(-1.*(wline-coeff(1))^2/(2.*coeff(2)^2)), color=cgcolor('red'),linestyle=5, thick=1.2
oplot, coeff(1)*[1,1], [0,2],thick=3, color=cgcolor('black')
xyouts, coeff(1)+0.02, 2.2,strtrim(coeff(1),2), orientation=90.,charsize=1.5

END



;Exercise 2.23: Define a function that returns the area of a circle given the radius.
FUNCTION CIRCLE_AREA, radius
; PURPOSE:
;     Return the area of a circle given the radius.
;
; INPUT:
;    radius   -scalar radius of the circle.
; 
; OUTPUT:
;     Floating pt scalar giving the area of the circle.
;
  return, !Pi*radius^2
END



; Exercise: 2.24: Write IDL function that returns the slope/intercept of a line given
; two data points in array format (x0,y0) and (x1,y1). Beware of vertical line case.
FUNCTION data2Line, xy1,xy2
;
; PURPOSE:
;    Computes slope and y-intercept given two sets of x & y positions.
;
; INPUT:
;    xy1    -2-element array giving x/y position of 1st data point [x1,y1].
;    xy2    -2-element array giving x/y position of 2nd data point [x2,y2].
;
; OUTPUT:
;     2-element floating pt array giving [slope, y-intercept]. Vertical line returns [-9999,-9999].
;
  xy1=float(xy1)
  xy2=float(xy2)
  deltax=xy2(0)-xy1(0)
  deltay=xy2(1)-xy1(1)
  slope=deltay/deltax
  y0=xy2(1) - slope * xy2(0)
  CASE finite(slope) OF
    0: line=[-9999.,-9999.]
    1: line=[slope, y0]
  ENDCASE
  return, line
END


