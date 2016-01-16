#!/usr/bin/python
'''
Matrix Product State.
'''
from numpy import *
from matplotlib.pyplot import *
from matplotlib import patches
from matplotlib.collections import LineCollection
from scipy.linalg import svd,qr,rq
from scipy import sparse as sps
from utils import bcast_dot
import pdb,time

class OpString(object):
    '''
    Operator String.
    '''
    def __init__(self,nsite):
        self.__opdict__={}
        self.nsite=nsite

    def __getitem__(self,l):
        return self.__opdict__.get(l)

    def __setitem__(self,l,op):
        self.__opdict__[l]=op

    def __iter__(self):
        for i in xrange(self.nsite):
            yield self.__opdict__.get(l)

    @property
    def oplist(self):
        '''A list of operators defined on sites.'''
        opl=[None]*self.nsite
        for l in self.__opdict__:
            opl[l]=self.__opdict__[l]

    @property
    def siteindices(self):
        '''The site indices with valid data.'''
        return self.__opdict__.keys()

class MPO(object):
    '''
    Matrix product operator.

    WL:
        The Matrix product operator datas.
    '''
    def __init__(self,WL):
        self.WL=WL

    def __str__(self):
        return self.WL.__str__()

    def serialize(self):
        '''
        Return The serialized form of operator.
        '''
        O=w[0]
        for w in self.WL[1:]:
            O=O.dot(w)
        return O[0,0]

    @property
    def nsite(self):
        '''Number of sites.'''
        return len(self.WL)
