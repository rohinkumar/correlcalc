# from fileios import *
# msg = 'Enter Absolute Path to file: '
# f_name = raw_input(msg).strip()
#
# path = file_data_and_path(f_name)
# if path != None:
#        print 'Path:',path
# from Tkinter import Tk
# from tkFileDialog import askopenfilename
#
# Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
# filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
# print(filename)
# def get_filename(file_type):
#     while True:
#         print('enter ' + file_type + ' filename: ')
#         filename = input()
#         print(filename)
#         try:
#             with open(filename, 'r') as f:
#                 my_file = f.read()
#             return my_file
#         except FileNotFoundError:
#             print('No such file.  Check file name and path and try again.')
#
#
# x = get_filename('TEMPLATE')
# print(x)
# -*- coding: utf-8 -*-
"""To add test methods.
"""
# from time import sleep
# from halo import Halo
# from time import time
#
# def rocket_launch():
#     #spinner = Halo({'spinner': 'shark'})
#     spinner = Halo({
#     'spinner': {
#         'interval': 100,
#         'frames': ['-', '\\', '|', '/', '-']
#     }
# })
#     spinner.start()
#     while(1):
#         spinner.text = 'Running... Time Elapsed: {} seconds'.format(time())
#         sleep(10)
#         break
#     spinner.succeed('Rocket launched')
#
# rocket_launch()

from antpcf import *

# bins = np.arange(0.01, 0.201, 0.01)
# atpcf('/Users/rohin/Downloads/DR7-Full.ascii', bins, randfile='/Users/rohin/Downloads/random-DR7-Ful.ascii',permetric='apzdth', parmetric='apdz', weights=True)
# tpcf('/Users/rohin/Downloads/DR3-ns.ascii',bins,randfile='/Users/rohin/Downloads/random-DR3-ns.ascii',weights=True)
# def pmethod():
bins = np.arange(0.002, 0.06, 0.002)
correl = tpcf('./testw.dat',bins,randfile='./testw.dat',weights=True)
#    return correl
# pool = multiprocessing.Pool(processes=ncount)
# correl = pool.map(pmethod)
# print correl
# atpcf('./testw.dat',bins,randfile='./testw.dat',permetric='apzdth',parmetric='apdz',method='ls',weights=True)
# blha=readfitsfile('/Users/rohin/Documents/ipy_notebooks/galsurveystudy/input/galaxy_DR12v5_CMASS_North.fits','data')
# dr12gcmn, weights = datprep('/Users/rohin/Documents/ipy_notebooks/galsurveystudy/input/galaxy_DR12v5_LOWZ_South.fits','data','lcdm')
# dat = readfitsfile('/Users/rohin/Documents/ipy_notebooks/galsurveystudy/input/galaxy_DR12v5_LOWZ_South.fits','data')
# weights = dat['WEIGHT_SYSTOT']
# import pyfits
# dpy = pyfits.open('/Users/rohin/Documents/ipy_notebooks/galsurveystudy/input/galaxy_DR12v5_LOWZ_South.fits')
# dpyd = dpy[1].data
# wts = dpyd['WEIGHT_SYSTOT']
# print(wts)
# print(min(wts))
# print(max(wts))
# print(dr12gcmn)
# print(weights)
# print (min(weights))
# print(max(weights))
# dr12gls=tpcf('/Users/rohin/Documents/ipy_notebooks/galsurveystudy/input/galaxy_DR12v5_LOWZ_South.fits',bins,randfile='/Users/rohin/Downloads/random0_DR12v5_LOWZ_South.fits',weights=True)
# Planck run with changed parameters in param.py
# corrdr3milne = tpcf('/Users/rohin/Downloads/DR3-ns.ascii', bins, randfile='/Users/rohin/Downloads/random-DR3-ns.ascii', weights=True, geometry='open', cosmology='lc')
# corrdr3milne = tpcf('/Users/rohin/Downloads/DR3-ns.ascii', bins, weights=True, mask='/Users/rohin/Documents/ipy_notebooks/galsurveystudy/masks/window.dr72safe0.ply')
# corrdr3milne = tpcf('./testw.dat', bins, weights=True, mask='/Users/rohin/Documents/ipy_notebooks/galsurveystudy/masks/window.dr72safe0.ply')
