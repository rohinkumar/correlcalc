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

from tpcf import *
bins=np.arange(0.002,0.082,0.002)
#atpcf('/Users/rohin/Downloads/DR7-Full.ascii',bins,randfile='/Users/rohin/Documents/ipy_notebooks/correlcalc-nb/randcat_DR72-2x.dat',permetric='apzdth',parmetric='apdz',method='ls')
tpcf('/Users/rohin/Downloads/DR3-ns.ascii',bins,randfile='/Users/rohin/Downloads/random-DR3-ns.ascii',weights=True)
#tpcf('./testw.dat',bins,randfile='./testw.dat',weights=True)
