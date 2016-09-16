#!/usr/bin/env python
"""Tight-binding chain"""

from __future__ import division
from scipy.sparse import kron,identity
from scipy.sparse.linalg import eigsh
import numpy as np

e0=0.2;t=-1.0
"""
e0:on-site energy
t:hopping t
"""
N=np.array([[0,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,2]])
Z=np.array([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,1]])
Cm_up=np.array([[0,0,1,0],[0,0,0,1],[0,0,0,0],[0,0,0,0]])
Cm_dn=np.array([[0,1,0,0],[0,0,0,0],[0,0,0,-1],[0,0,0,0]])
Cp_up=Cm_up.conjugate().transpose()
Cp_dn=Cm_dn.conjugate().transpose()

def fermions_chain_1(L):
    """Method 1"""
    cp_u={}
    cp_d={}
    cm_u={}
    cm_d={}
    n={}

    for i in range(L):
        cm_u[i]=Cm_up
        cm_d[i]=Cm_dn
        cp_u[i]=Cp_up
        cp_d[i]=Cp_dn
        n[i]=N
        for l in range(i):
            cm_u[i]=kron(Z,cm_u[i])
            cm_d[i]=kron(Z,cm_d[i])
            cp_u[i]=kron(Z,cp_u[i])
            cp_d[i]=kron(Z,cp_d[i])
            n[i]=kron(identity(4),n[i])
        for l in range(i+1,L):
            cm_u[i]=kron(cm_u[i],identity(4))
            cm_d[i]=kron(cm_d[i],identity(4))
            cp_u[i]=kron(cp_u[i],identity(4))
            cp_d[i]=kron(cp_d[i],identity(4))
            n[i]=kron(n[i],identity(4))

    H=e0*n[L-1]+t*cp_u[L-1].dot(cm_u[0])+t*cp_d[L-1].dot(cm_d[0])\
               +t*cp_u[0].dot(cm_u[L-1])+t*cp_d[0].dot(cm_d[L-1])
    for i in range(L-1):
        H=H+e0*n[i]+t*cp_u[i].dot(cm_u[i+1])+t*cp_d[i].dot(cm_d[i+1])\
                   +t*cp_u[i+1].dot(cm_u[i])+t*cp_d[i+1].dot(cm_d[i])

    E0,psi0=eigsh(H,k=1,which='SA')
    print "E/L(Method 1):",E0/L

def fermions_chain_2(L):
    """Method 2"""
    H=e0*N
    cp_up=Cp_up
    cm_up=Cm_up
    cp_dn=Cp_dn
    cm_dn=Cm_dn
    d=4

    for i in range(L-1):
        H=kron(H,identity(4))+kron(identity(d),e0*N)\
         +t*kron(cp_up,identity(4)).dot(kron(Z,cm_up))\
         +t*kron(cp_dn,identity(4)).dot(kron(Z,cm_dn))\
         +t*kron(Z,cp_up).dot(kron(cm_up,identity(4)))\
         +t*kron(Z,cp_dn).dot(kron(cm_dn,identity(4)))
        cp_up=kron(Z,cp_up)
        cm_up=kron(Z,cm_up)
        cp_dn=kron(Z,cp_dn)
        cm_dn=kron(Z,cm_dn)
        d=d*4
        
    H=H+t*kron(Cp_up,identity(4**(L-1))).dot(cm_up)\
       +t*kron(Cp_dn,identity(4**(L-1))).dot(cm_dn)\
       +t*cp_up.dot(kron(Cm_up,identity(4**(L-1))))\
       +t*cp_dn.dot(kron(Cm_dn,identity(4**(L-1))))

    E0,psi0=eigsh(H,k=1,which='SA')
    print "E/L(Method 2):",E0/L

def exact(L):
    """Exact solution"""
    E0=0.0
    for i in range(0,L,1):
        Ek=e0+t*np.exp(1j*2*np.pi*i/L)+t*np.exp(-1j*2*np.pi*i/L)
        if Ek<0:
            E0=E0+e0+t*np.exp(1j*2*np.pi*i/L)+t*np.exp(-1j*2*np.pi*i/L)
    print "E/L(Exact):",2*E0/L

if __name__=="__main__":    

    fermions_chain_1(6)

    fermions_chain_2(6)

    exact(6)

