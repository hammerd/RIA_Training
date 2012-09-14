#! /usr/bin/env python

'''
ABOUT:
This program performs point-source photometry using the IRAF tasks daofind and daoimage on a single HST image or a list of images.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER for STScI, 2012

HISTORY:
Sept. 2012: Original script (v0.1).


FUTURE IMPROVEMENTS:
-Add ds9 & tvmark command to inspect sources detections from daofind.
-Add automated procedure to measure the PSF.
-Add procedure to remove "INDEF" vlaues in external file (replace with -9999.0)
-Add procedure to annotate graphs usign current plot limits, instead of hard-wired plot positions.
-Possible change method by whcih we perform source detection. Should probably be performed on e-/sec
 images, not e-, and weight detections using rms map.


USE:
ptsrc_photom.py
'''

__author__='D.M. HAMMER'
__version__= 0.1


import pyraf, os, glob, argparse, pdb, pyfits, pylab, fileinput
import numpy as np
from pyraf import iraf
from iraf import noao, digiphot, daophot
from iraf import images, imutil


def mk_counts_image(image, outfile='default'):

	'''FUNCTION TO CONVERT CNTS/SEC IMAGE TO COUNTS & RESUPPLY THE BACKGROUND'''

	# Parse input parameters
	if outfile == 'default':
		osplit=image.split('.fits')
		outfile=osplit[0]+'_cts.fits'

	# Read in fits image and header
	f = pyfits.open(image)
	fdata = f[0].data
	fheader = f[0].header
	ftable = f[1].data
	f.close()

	# Extract relevant info from the header
	exptime = fheader['texptime']
	back = ftable.field('mdrizsky')
	ipxscl = fheader['D001ISCL']
	opxscl = fheader['D001SCAL']

	# Add background and correct for different pixel scale (original backgrd is measured in raw images)
	fdata_cnts = np.copy(fdata) * exptime + np.sum(back)/2.0*(opxscl/ipxscl)**2

	# Save cnts image to disk
	file_query = os.access(outfile, os.R_OK)
	if file_query == True: os.remove(outfile)		# remove output file if it already exists
	pyfits.writeto(outfile, fdata_cnts,header=fheader)

	return outfile



def run_daofind(image, extension=0,outfile='default',dthreshold=3.0, fwhmpsf=2.5, backsigma=-1.0,rdnoise=5.2):

	'''THIS PROCEDURE RUNS DAOFIND ON INPUT IMAGE'''

	# Parse input parameters
	if outfile == 'default': outfile = image+'0.coo.1'

	# Read in fits header
	f = pyfits.open(image)
	fheader = f[0].header
	f.close()

	# Extract relevant info from the header
	exptime = fheader['texptime']
	instr = fheader['INSTRUME']
	if instr == 'WFC3':
		filter = fheader['FILTER']
	else: #assuming ACS
		filter = fheader['FILTER1']
		if filter[0] == 'C': filter == fheader['FILTER2']
	ipxscl = fheader['D001ISCL']
	opxscl = fheader['D001SCAL']
	num_flts = float(fheader['NDRIZIM'])/2.0
	#df_max = 10000000.	# upper limit for "good" pixels (sometimes used to ID bad pixels)---Not used here.

	# Perform read noise correction
	rdnoise_corr = np.sqrt(num_flts * (rdnoise * opxscl/ipxscl)**2)
		
	# Perform background noise calculation
	if backsigma < 0.0:
		backrms=iraf.imstatistics(image+'[0]', fields='stddev', nclip=10, lsigma=3.0, usigma=3.0, cache='yes', format='no',Stdout=1)
		backsigma=float(backrms[0])


	'''Run daofind'''
	# remove old daofind files
	file_query = os.access(outfile, os.R_OK)	
	if file_query == True: os.remove(outfile)
	iraf.daofind.unlearn()
	iraf.daofind(image=image+'[0]', interactive='no', verify='no',output=outfile, fwhmpsf=fwhmpsf, sigma=backsigma, \
	readnoise=rdnoise_corr, itime=exptime, threshold=dthreshold)

	# Display results of daofind (***WORK IN PROGRESS***)
	#!ds9 &
	#iraf.display(file+'['+exten+']',1)
	#iraf.tvmark(1,ofile+'0.coo.1',mark = 'circle', radii = 8, color = 205)

	return outfile		# return name of coordinate file



def run_daophot(image, outfile='default', coordfile='NA', apertures='5.0,16.66', annulus=17.0, dannulus=3.0, calgorithm='centroid', salgorithm='median', fwhmpsf=2.5, backsigma=-1.0,rdnoise=5.2):

	'''THIS PROCEDURE RUNS DAOPHOT ON INPUT IMAGE'''

	# Parse input parameters
	if outfile == 'default': outfile = image + '0.mag.1'
	if coordfile == 'NA': coordfile = image + '0.coo.1'

	# Read in fits header
	f = pyfits.open(image)
	fheader = f[0].header
	f.close()

	# Extract relevant info from the header
	exptime = fheader['texptime']
	filter = fheader['FILTER2']
	ipxscl = fheader['D001ISCL']
	opxscl = fheader['D001SCAL']
	num_flts = float(fheader['NDRIZIM'])/2.0
	dp_zmag = -2.5 * np.log10(float(fheader['PHOTFLAM'])) + float(fheader['PHOTZPT'])       # zero-pt for each filter
	#df_max = 10000000.      # upper limit for "good" pixels (sometimes used to ID bad pixels). Not using this parameter right now.

	# Perform read noise correction
	rdnoise_corr = np.sqrt(num_flts * (rdnoise * opxscl/ipxscl)**2)

	# Perform background noise calculation
	if backsigma < 0.0:
		backrms=iraf.imstatistics(image+'[0]', fields='stddev', nclip=10, lsigma=3.0, usigma=3.0, cache='yes', format='no',Stdout=1)
		backsigma=float(backrms[0])


	'''Run daophot'''
	# Remove old phot output files
	file_query = os.access(outfile, os.R_OK)      
	if file_query == True: os.remove(outfile)

	# Run phot
	iraf.phot.unlearn()         # reset daophot parameters to default values
	iraf.phot(image=image+'[0]', interactive='no', verify='no', coords=coordfile, output=outfile, fwhmpsf=fwhmpsf, \
      		sigma=backsigma, readnoise=rdnoise_corr, itime=exptime, calgorithm=calgorithm, salgorithm=salgorithm, \
      		annulus=annulus, dannulus=dannulus, apertures=apertures,zmag=dp_zmag)


	return outfile 		# return name of output catalog



def replace_filevalue(file, orgval, newval):

	''' REPLACE UNWANTED VALUES IN EXTERNAL FILE '''

	for line in fileinput.input(file, inplace = 1):
		print line.replace(str(orgval), str(newval)),
	fileinput.close()



def calc_apcorr(file):

	'''  MEASURE APPERTURE CORRECTION GIVEN TWO MAGNITUDES'''
	'''    file -- daophot filename (expects only 2 circular aperture magnitude estimates in file)'''

	# Save a trimmed daophot catalog with columns relevant to this study (xc,yc,mag1,mag2)
        trimfile = file+'.trimmed'
        file_query = os.access(trimfile, os.R_OK)
        if file_query == True: os.remove(trimfile)      # remove previous trimmed files
        iraf.txdump(file,'xcenter,ycenter,mag', 'yes', Stdout=trimfile)

        # Measure the aperture correction assuming a pt source
        xc, yc, mag1,mag2 = np.loadtxt(trimfile, unpack=True)
        deltamag=mag1-mag2

	# Use somewhat hacky method to extract box in mag.vs.deltamag plot where we will estimate median deltamag
	xgd=[21.8,25.0]		# fixed--works well for both filters
	#estimate y-range from a rough median estimate
	tmp_median = np.median(deltamag[(mag1 >= xgd[0])&(mag1 <= xgd[1])&(np.abs(deltamag) < 1)])
	ygd=[tmp_median-0.1, tmp_median+0.1]
        f_deltamag=deltamag[(mag1 >= xgd[0])&(mag1 <= xgd[1])&(deltamag > ygd[0])&(deltamag <ygd[1])]     # avoid saturated and faint sources
	ac05 = np.median(f_deltamag)


	# Plot the aperture correction on delta_mag vs. mag scatter diagram
        pylab.scatter(mag1, deltamag,s=1.0,c='r',marker='o')
        pylab.ylim(-0.5,1.0)
        pylab.xlim(19.0,29.0)
        pylab.xlabel('m3')
        pylab.ylabel('m3 - m16')
        pylab.axhline(0.0, c='k', ls='--', linewidth=3)
        pylab.axhline(ac05)
	pylab.vlines(xgd,ygd[0],ygd[1], colors='red',linestyles='dotted')
	pylab.hlines(ygd,xgd[0],xgd[1], colors='red',linestyles='dotted')
        pylab.annotate('ac05='+str(ac05),[19.3,ac05-0.07])
        pylab.annotate('m_corr = m1 - ac05 - AC05',[19.3,-0.4])
        fsplit=file.split('_')
        pylab.annotate(fsplit[0], [20.0,-0.2],color='r')
        pylab.savefig('deltamag_vs_mag_'+fsplit[0]+'.pdf')
        pylab.clf()

	return ac05



def mk_cmd(mag, deltamag, outfile = 'CMD.pdf'):
	pdb.set_trace()
	pylab.scatter(deltamag, mag, s=1.0,c='r',marker='o')
	pylab.ylim(32,17.0)
	pylab.xlim(-2,4)
	pylab.xlabel('F606W - F814W')
	pylab.ylabel('F606W [STMag]')
	pylab.savefig(outfile)
	pylab.clf()



if __name__=='__main__':

	'''THIS PROCEDURE ADDRESSES PHOTOMETRY TRAINING'''

        # Parse input parameters
        parser = argparse.ArgumentParser(description='Create stellar CMD using DAO find/phot procedures.')
        parser.add_argument('-im', '--images',default='f*drc_sci.fits', type=str, help='Input fits file(s). Default is all drizzled science images in working directory.')
        parser.add_argument('-ext', '--extension', default=1, type=int, help='Fits extension for photometry (default=1).')
        options = parser.parse_args()

	# Initialize filename and aperture correction variables
	file_list = glob.glob(options.images)
	cnts_name = ['f606w_drc_sci_cts.fits', 'f814w_drc_sci_cts.fits']	#originally left the name variables as blank arrays (i.e., "=[]")
        find_name = [cnts_name[0]+'0.coo.1', cnts_name[1]+'0.coo.1']
	phot_name = [cnts_name[0]+'0.mag.1', cnts_name[1]+'0.mag.1']
	col_name = [cnts_name[0]+'0.f814colmag.1', cnts_name[1]+'0.f606colmag.1']
	tcol_name = [cnts_name[0]+'0.f814colmag.1.trimmed', cnts_name[1]+'0.f606colmag.1.trimmed']

	ac05 = np.array([-9999.0,-9999.0])	# 0.09"-->0.5" aperture correction variable
	AC05 = 0.088	# aperture correction from 0.5" --> INF (Sirianni et al. 2005).


	'''Generate source catalogs'''
	# catalogs for src detection & photometry within same image
        for file, ff in zip(file_list,xrange(len(file_list))):
	#	mk_counts_image(file, outfile=cnts_name[ff])					
	#	run_daofind(cnts_name[ff], outfile=find_name[ff], dthreshold=4.0)			
	#	run_daophot(cnts_name[ff], coordfile=find_name[ff], outfile=phot_name[ff])
	#	replace_filevalue(phot_name[ff], 'INDEF',-9999.0)
		ac05[ff] = calc_apcorr(phot_name[ff])

	# color catalogs--src detection in one image, photometry in other (***KLUGE*** ASSUMES ONLY 2 INPUT IMAGES)
	# f606w color cats
	#run_daophot(cnts_name[0], outfile=col_name[0], coordfile=find_name[1])
        #replace_filevalue(col_name[0], 'INDEF', -9999.0)
	#file_query = os.access(tcol_name[0], os.R_OK)
	#if file_query == True: os.remove(tcol_name[0])      # remove previous trimmed files
	#iraf.txdump(col_name[0],'xcenter,ycenter,mag', 'yes', Stdout=tcol_name[0])

	#f814w color cats
	#run_daophot(cnts_name[1], outfile=col_name[1], coordfile=find_name[0])
        #replace_filevalue(col_name[1], 'INDEF', -9999.0)
	#file_query = os.access(tcol_name[1], os.R_OK)
	#if file_query == True: os.remove(tcol_name[1])      # remove previous trimmed files
	#iraf.txdump(col_name[1],'xcenter,ycenter,mag', 'yes', Stdout=tcol_name[1])


        '''make CMD'''
	# extract color mags (apply aperture corrections)
        xca, yca, mag1a,mag2a = np.loadtxt(phot_name[0]+'.trimmed', unpack=True)
	xcb, ycb, mag1b,mag2b = np.loadtxt(tcol_name[1], unpack=True)
	tmaga = mag1a - ac05[0] - AC05
	tmagb = mag1b - ac05[1] - AC05
	mk_cmd(tmaga, tmaga-tmagb, outfile='CMD_color.pdf')


	# extract individual mags -- must match catalogs by position (apply aperture corrections)
	#iraf.tmatch()

