#! /usr/bin/env python

'''
ABOUT:
This program calls IRAF task daofind on a single image or a list of images.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER for STScI, 2012

HISTORY:
2012: Trial program.

USE:
python py_daofind.py
'''

__author__='D.M. HAMMER'
__version__= 0.1

import pyraf, os, glob, argparse, pdb
from pyraf import iraf
from iraf import noao, digiphot, daophot



if __name__=='__main__':
	parser = argparse.ArgumentParser(description='Run IRAF task daofind.')
        parser.add_argument('-f', '--file',default='*_ima.fits', type=str, help='Input filename. Default is all *_imh.fits files.')
        parser.add_argument('-ext', '--extension', default=1, type=int, help='Fits extension for photometry (only one).')
        options = parser.parse_args()

        #pdb.set_trace()
	exten = str(options.extension)
	infile=options.file

	# Generate a list of all fits files
	file_list = glob.glob(infile)
	print file_list

	# Loop through all fits files and run daofind
	for file in file_list:
		# Test for old files daofind, and remove them if they exist
		file_query = os.access(file + '1.coo.1', os.R_OK)
		if file_query == True: os.remove(file + '1.coo.1')

		# Run daofind on one image
		iraf.daofind(image=infile+'['+exten+']', interactive='no', verify='no')
	
