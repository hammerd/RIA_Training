#! /usr/bin/env python


'''
ABOUT:
This program addresses assignment#1 of the spectroscopy training program. Specifically, we will split the original exposures
into smaller durations in order to identify a flare that occurred during observations.  We use the IRAF COS task splittag & x1dspec
to accomplish these tasks. We also use IRAF splot to view and measure line indices that trace the flare event.


DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER for STScI, 2012

HISTORY:
Sept. 2012: Original script (v0.1).


FUTURE IMPROVEMENTS:

USE:
spec_ass1.py
'''


__author__='D.M. HAMMER'
__version__= 0.1


import pyraf
from pyraf import iraf
from iraf import stsdas, hst_calib, hstcos
import os, glob, argparse, pdb, pylab


if __name__=='__main__':

    '''THIS FUNCTION ADDRESSES SPECTROSCOPY TRAINING ASSIGNMENT #1'''

	# parse input parameters
        parser = argparse.ArgumentParser(description='Identify flare event captured during COS observations.')
        parser.add_argument('-ev', '--events',default='*corrtag*.fits', type=str, help='Input fits file(s). Default is all corrtag files in working directory.')
	parser.add_argument('-odir', '--outdir', default='./splittag15/', type=str, help='Directory to output split files and analysis.')
        options = parser.parse_args()
	outdir=options.outdir

	# Initialize filename and aperture correction variables
        file_list = glob.glob(options.events)
	obshdr = [x[0:9] for x in file_list[::2]]	# the [::2] selects every other file (there are 2 corrtag files for each spectra)
	ctag_list = [x[0:18]+'?.fits' for x in file_list[::2]]
	prefix_ctag_list = [x[0:9] for x in ctag_list]
	prefix_stag_list = [x[0:9]+'_stag' for x in file_list[::2]]

	# remove all previous stag files in the "outdir" directory
	oldf = glob.glob(outdir+'/*stag*')
	for x in oldf: os.remove(x)

	# Create 1d spectra in smaller time intervals
	for file, ff in zip(ctag_list, xrange(len(ctag_list))):
		# split COS obs into spectra with smaller time duration
		iraf.splittag(file, outdir+prefix_stag_list[ff], increment=15.0)

		# assemble 1d spectra for each split observation.
		stag_list = glob.glob(outdir+obshdr[ff]+'_stag*corrtag*.fits')
		for ss in xrange(len(stag_list)):
			iraf.x1dcorr(outdir+prefix_stag_list[ff]+'_'+str(ss)+'_corrtag_?.fits')
			iraf.tomultispec(outdir+prefix_stag_list[ff]+'_'+str(ss)+'_x1d.fits', outdir+prefix_stag_list[ff]+'_'+str(ss)+'_x1d.imh', flux_col='FLUX', wave_col='WAVELENGTH')


	# Load 1d spectra into viewer
	pdb.set_trace()
	iraf.splot('*.imh')
