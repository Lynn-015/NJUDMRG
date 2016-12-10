#!/usr/bin/python
'''
Matrix Product State.
'''
from numpy import *
from matplotlib.pyplot import *
from matplotlib import patches
from matplotlib.collections import LineCollection
from scipy.linalg import svd,qr,rq,norm,block_diag
from mps import MPS,VidalMPS,KMPS,BMPS
from scipy import sparse as sps
import pdb,time

def state2VMPS(state,sitedim,tol=1e-8):
    '''
    Parse a normal state into a Vidal Matrix produdct state.

    Parameters
    --------------
    state:
        The target state, 1D array.
    sitedim:
        The dimension of a single site, integer.
    tol:
        The tolerence of singular value, float.

    *return*:
        A <VidalMPS> instance.

    Note
    ---------------
    `svd` method is used in decomposition.
    '''
    nsite=int(round(log(len(state))/log(sitedim)))
    GL,LL=[],[]
    ri=1

    for i in xrange(nsite):
        state=state.reshape([sitedim*ri,-1])
        U,S,V=svd(state,full_matrices=False)
        #remove zeros from v
        kpmask=abs(S)>tol
        ri=kpmask.sum()
        S=S[kpmask]
        LL.append(S)
        state=S[:,newaxis]*V[kpmask]
        U=U[:,kpmask]
        U=U.reshape([-1,sitedim,ri])
        ai=swapaxes(U,0,1)
        if i==0:
            gi=ai
        else:
            gi=ai/LL[-2][:,newaxis]
        GL.append(gi)
    S=LL[-1][0]
    return VidalMPS(GL,LL[:-1],S)

def state2MPS(state,sitedim,l,method='qr',tol=1e-8):
    '''
    Parse a normal state into a Matrix produdct state.

    state:
        The target state, 1D array.
    sitedim:
        The dimension of a single site, integer.
    l:
        The division point of left and right canonical scanning, integer between 0 and number of site.
    method:
        The method to extract A,B matrices.
        * 'qr'  -> get A,B matrices by the method of QR decomposition, faster, rank revealing in a non-straight-forward way.
        * 'svd'  -> get A,B matrices by the method of SVD decomposition, slow, rank revealing.
    tol:
        The tolerence of singular value, float.

    *return*:
        A <MPS> instance.
    '''
    nsite=int(round(log(len(state))/log(sitedim)))
    AL,BL=[],[]
    ri=1
    assert(method=='svd' or method=='qr')
    assert(l>=0 and l<=nsite)

    for i in xrange(l):
        state=state.reshape([sitedim*ri,-1])
        if method=='svd':
            U,S,V=svd(state,full_matrices=False)
            #remove zeros from v
            kpmask=abs(S)>tol
            ri=kpmask.sum()
            state=S[kpmask,newaxis]*V[kpmask]
            U=U[:,kpmask]
        else:
            U,state=qr(state,mode='economic')
            kpmask=sum(abs(state),axis=1)>tol
            ri=kpmask.sum()
            state=state[kpmask]
            U=U[:,kpmask]
        ai=swapaxes(U.reshape([-1,sitedim,ri]),0,1)
        AL.append(ai)

    ri=1
    for i in xrange(nsite-l):
        state=state.reshape([-1,sitedim*ri])
        if method=='svd':
            U,S,V=svd(state,full_matrices=False)
            #remove zeros from v
            kpmask=abs(S)>tol
            ri=kpmask.sum()
            state=S[kpmask]*U[:,kpmask]
            V=V[kpmask,:]
        else:
            state,V=rq(state,mode='economic')
            kpmask=sum(abs(state),axis=0)>tol
            ri=kpmask.sum()
            state=state[:,kpmask]
            V=V[kpmask]
        bi=swapaxes(V.reshape([ri,sitedim,-1]),0,1)
        BL.append(bi)
    BL=BL[::-1]
    S=state.diagonal()
    return KMPS(AL,BL,S=S)

def mps_add(*args):
    '''
    Add <KMPS>.

    Parameters
    -----------
    args:
        <MPS> instances to be added.
    '''
    if len(args)<=1:
        raise ValueError('At least 2 args is required.')
    AL=[]
    BL=[]
    hndim=args[0].hndim
    na=len(args[0].AL)
    nb=len(args[0].BL)
    nsite=na+nb
    for i in xrange(na):
        if i==0:
            ai=[concatenate([mps.AL[i][j] for mps in args],axis=1) for j in xrange(hndim)]
        elif i==nsite-1:
            ai=[concatenate([mps.AL[i][j] for mps in args],axis=0) for j in xrange(hndim)]
        else:
            ai=[block_diag(*[mps.AL[i][j] for mps in args]) for j in xrange(hndim)]
        AL.append(ai)
    for i in xrange(nb):
        if i+na==0:
            bi=[concatenate([mps.BL[i][j] for mps in args],axis=1) for j in xrange(hndim)]
        elif i+na==nsite-1:
            bi=[concatenate([mps.BL[i][j] for mps in args],axis=0) for j in xrange(hndim)]
        else:
            bi=[block_diag(*[mps.BL[i][j] for mps in args]) for j in xrange(hndim)]
        BL.append(bi)
    S=concatenate([mps.S for mps in args])
    return args[0].__class__(AL=AL,BL=BL,S=S)

def mps_unique(mps):
    '''
    View <MPS> as unique.
    '''
    mps.token=id(l)

def mps_viewsame(*mpses):
    '''
    View a set of <MPS> instances as same.
    '''
    mps0=mpses[0]
    for mps in mpses[1:]:
        mps.token=mps0.token

def compress(self,mps):
    '''
    Compress a MPS to compact form.
    '''
    raise Exception('Not Implemented')

