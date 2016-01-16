#!/usr/bin/python
#-*-coding:utf-8-*-
#Author: K. Leo
#Date: 2014-06
#Description: this is a mean field program

from numpy import *
from tba.hgen.generator import RHGenerator
from tba.hgen.oplib import op_simple_hopping,op_U
from tba.lattice.structure import Structure
from tba.hgen.spaceconfig import SuperSpaceConfig,SpaceConfig
import time,pdb,warnings

class Hex6(object):
    '''This is a tight-binding model for a hexagon.'''
    def __init__(self,t,U=0,occ=True):
        self.t,self.U=t,U
        self.occ=occ

        #occupation representation will use <SuperSpaceConfig>, otherwise <SpaceConfig>.
        if self.occ:
            spaceconfig=SuperSpaceConfig([1,2,6,1])
        else:
            spaceconfig=SpaceConfig([1,2,6,1],kspace=False)
            if abs(U)>0: warnings.warn('U is ignored in non-occupation representation.')
        hgen=RHGenerator(spaceconfig=spaceconfig)

        #define the operator of the system
        hgen.register_params({
            't1':self.t,
            'U':self.U,
            })

        #define a structure and initialize bonds.
        rlattice=Structure(sites=[(0.,0),(0,sqrt(3.)/3),(0.5,sqrt(3.)/2),(0.5,-sqrt(3.)/6),(1.,0),(1.,sqrt(3.)/3)])
        hgen.uselattice(rlattice)

        b1s=rlattice.getbonds(1)  #the nearest neighbor

        #add the hopping term.
        op_t1=op_simple_hopping(label='hop1',spaceconfig=spaceconfig,bonds=b1s)
        hgen.register_operator(op_t1,param='t1')

        #add the hubbard interaction term if it is in the occupation number representation.
        if self.occ:
            op_ninj=op_U(label='ninj',spaceconfig=spaceconfig)
            hgen.register_operator(op_ninj,param='U')

        self.hgen=hgen
