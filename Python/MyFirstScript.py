#! /usr/bin/env python

'''
ABOUT:
This program reads data from an external file and makes line plots using matplotlib.

DEPENDS:
Python 2.5.4

AUTHOR:
D. HAMMER for STScI, 2012

HISTORY:
2012: Trial program.

USE:
python MyFirstScript.py
'''

__author__='D.M. HAMMER'
__version__= 0.2


import numpy as np
import pylab as pl
import argparse
import pdb


def mkplot(outfile,xx,yy1,yy2,yy3,yy4,xlab='Slope [e-/s]',ylab='Slope Unc. [e-/s]'):
	pl.loglog(xx,yy1,'b--',xx,yy2,'r:',xx,yy3,'g-',xx,yy4,'m-.', linewidth=3)
	pl.legend(('Random Uncertainty', 'Correlated Uncertainty', 'Both', 'Equation'), 'best')
	pl.ylabel(ylab)
	pl.xlabel(xlab)
	pl.savefig(outfile)
	pl.clf()
	print 'Saved file to: ', outfile
	return


if __name__=='__main__':

	parser = argparse.ArgumentParser(description='Make a plot.')
	parser.add_argument('-f', '--file',default='Gordon2005_Fig16.txt', type=str, help='Input file.')
	options = parser.parse_args()

	pdb.set_trace()

	infile=options.file
	slope_outfile='fig16_slope.pdf'	
	yint_outfile='fig16_yint.pdf'

        slope, ran_slope_unc, corr_slope_unc, both_slope_unc, eqn_slope_unc, \
	ran_yint_unc, corr_yint_unc, both_yint_unc, eqn_yint_unc = \
        np.loadtxt(infile,unpack=True)

	mkplot(slope_outfile, slope, ran_slope_unc, corr_slope_unc, both_slope_unc, eqn_slope_unc)
	mkplot(yint_outfile, slope, ran_yint_unc, corr_yint_unc, both_yint_unc, eqn_yint_unc, ylab='Y-Intercept Unc. [e-]')
