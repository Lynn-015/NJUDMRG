#!/usr/bin/python
'''
Job selector and Runner.

Copy this file to `./run.py`(Do not modify this file!), uncomment specific task to run it!

    $ cp run_sample.py run.py
    $ vim run.py   ###do some modification.
    $ python run.py
'''

############ model_bc3 ###############
from model_bc3.views import *
#show the lattice structure of BC3
#bc3lattice()

#plot the band structure of BC3
#bc3band()

#show the fermi surface of BC3
#bc3fs()

#plot the dispersion of BC3 using 3D plot.
#bc3disper()

############# model_hex6 ################
from model_hex6.views import *
#calculate the ground state energy of a hexagon structure.
#hexE0()
#compare the ground state energy of a hexagon structure with the exact result(without interaction).
hexcompare()
