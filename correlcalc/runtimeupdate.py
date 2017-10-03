# -*- coding: utf-8 -*-
"""Update runtime graphic.
"""
from time import sleep
from halo import Halo
from time import time

#def runupdate(text,fintext):
    #spinner = Halo({'spinner': 'shark'})
spinner = Halo({
    'spinner': {
        'interval': 100,
        'frames': ['-', '\\', '|', '/', '-']
    }    })

#spinner.text = text  #'Running...'# Time Elapsed: {} seconds'.format(time())
#spinner.start()
#spinner.succeed(fintext)


