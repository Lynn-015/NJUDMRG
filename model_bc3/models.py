#!/usr/bin/python
#-*-coding:utf-8-*-
'''
Author: K. Leo
Date: 2015-11
Description: this is a simple honeycomb lattice model(BC3) as a tutorial.
'''

from tba.hgen.generator import KHGenerator
from tba.hgen.oplib import *
from tba.lattice.latticelib import Honeycomb_Lattice
from tba.hgen.spaceconfig import SpaceConfig
from tba.lattice.bond import show_bonds
from tba.lattice.group import TranslationGroup
import time

class BC3(object):
    '''
    This is a tight-binding model for bc3.

    Construct
    ------------------------
    BC3(N=(100,100)), N is the size of honeycomb lattice.

    Attributes
    -----------------
    hgen:
        Hamiltonian generator, a <KHGenerator> instance.
    '''
    def __init__(self,N=(100,100)):
        #define a Hilber space Configuration with 1-nnambu(not superconducting), 2-spin, 2-atom, 1-orbital.
        spaceconfig=SpaceConfig([1,2,2,1],kspace=True)
        #define a hamiltonian generator in KSpace
        hgen=KHGenerator(spaceconfig=spaceconfig,propergauge=False)

        #define a lattice with size N, and will be used by hgen.
        rlattice=Honeycomb_Lattice(N=N)
        #set up the periodic boundary condition in both x and y direction.
        tg=TranslationGroup(Rs=[ai*ni for ai,ni in zip(rlattice.a,rlattice.N)],per=[True,True])
        rlattice.usegroup(tg)
        hgen.uselattice(rlattice)

        #define the parameters for this model.
        hgen.register_params({
            '-mu':0.5,
            't1':0.62,
            't2':0.,
            't3':-0.38,
            })

        #get bonds for a unit cell, which can be extracted from rlattice.cbonds
        b1s,b2s,b3s=hgen.rlattice.cbonds[1:4]

        #define the operators form hopping terms defined on specific bonds.
        #bind it to the specific parameters.
        op_t1=op_simple_hopping(label='hop1',spaceconfig=spaceconfig,bonds=b1s)
        hgen.register_operator(op_t1,param='t1')

        op_t2=op_simple_hopping(label='hop2',spaceconfig=spaceconfig,bonds=b2s)
        hgen.register_operator(op_t2,param='t2')

        op_t3=op_simple_hopping(label='hop3',spaceconfig=spaceconfig,bonds=b3s)
        hgen.register_operator(op_t3,param='t3')

        op_n=op_simple_onsite(label='n',spaceconfig=spaceconfig)
        hgen.register_operator(op_n,param='-mu')

        self.hgen=hgen


