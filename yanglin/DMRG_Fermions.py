#!/usr/bin/env python
"""Tight-binding chain"""

from __future__ import division
from scipy.sparse import kron,identity
from scipy.sparse.linalg import eigsh
from copy import copy
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

class Block(object):  
    def __init__(self):
        """Initialize the block"""
        self.H=e0*N
        self.cp_up=Cp_up
        self.cm_up=Cm_up
        self.cp_dn=Cp_dn
        self.cm_dn=Cm_dn
        self.l=1
        self.d=4
        self.z=Z

    def enlarge(self):
        """Enlarge the block"""
        self.H=kron(self.H,identity(4))+kron(identity(self.d),e0*N)\
              +t*kron(self.cp_up,identity(4)).dot(kron(Z,self.cm_up))\
              +t*kron(self.cp_dn,identity(4)).dot(kron(Z,self.cm_dn))\
              +t*kron(Z,self.cp_up).dot(kron(identity(4),self.cm_up))\
              +t*kron(Z,self.cp_dn).dot(kron(identity(4),self.cm_dn))
        self.cp_up=kron(self.z,Cp_up)
        self.cm_up=kron(self.z,Cm_up)
        self.cp_dn=kron(self.z,Cp_dn)
        self.cm_dn=kron(self.z,Cm_dn)
        self.l=self.l+1
        self.d=self.d*4
        self.z=kron(Z,self.z)

class Super_Block(object):
    def __init__(self,l_block,r_block):
        """Form super-block"""
        self.l_block=l_block
        self.r_block=r_block
        self.H=kron(self.l_block.H,identity(self.r_block.d))\
              +kron(identity(self.l_block.d),self.r_block.H)\
              +t*kron(self.l_block.cp_up,identity(self.r_block.d)).dot(kron(self.l_block.z,self.r_block.cm_up))\
              +t*kron(self.l_block.cp_dn,identity(self.r_block.d)).dot(kron(self.l_block.z,self.r_block.cm_dn))\
              +t*kron(self.l_block.z,self.r_block.cp_up).dot(kron(self.l_block.cm_up,identity(self.r_block.d)))\
              +t*kron(self.l_block.z,self.r_block.cp_dn).dot(kron(self.l_block.cm_dn,identity(self.r_block.d)))
        self.l=self.l_block.l+self.r_block.l
        self.d=self.l_block.d*self.r_block.d

    def eigen(self):
        """Get ground state energy"""
        self.E0,self.psi0=eigsh(self.H,k=1,which='SA')
        return self.E0

    def truncate(self,m):
        """Truncate the left block"""
        self.psi0=self.psi0.reshape([self.l_block.d,-1],order='C')
        rho=np.dot(self.psi0,self.psi0.conjugate().transpose())
        e_vals,e_vecs=np.linalg.eigh(rho)

        eigens=[]
        for e_val,e_vec in zip(e_vals,e_vecs.transpose()):
            eigens.append((e_val,e_vec))
        eigens.sort(reverse=True,key=lambda x:x[0])

        t_matrix=np.zeros((self.l_block.d,m),dtype='d',order='F')
        for i,(e_val,e_vec) in enumerate(eigens[:m]):
            t_matrix[:,i]=e_vec
        self.error=1-sum([x[0] for x in eigens[:m]])

        self.l_block.H=t_matrix.conjugate().transpose().dot(self.l_block.H.dot(t_matrix))
        self.l_block.cp_up=t_matrix.conjugate().transpose().dot(self.l_block.cp_up.dot(t_matrix))
        self.l_block.cm_up=t_matrix.conjugate().transpose().dot(self.l_block.cm_up.dot(t_matrix))
        self.l_block.cp_dn=t_matrix.conjugate().transpose().dot(self.l_block.cp_dn.dot(t_matrix))
        self.l_block.cm_dn=t_matrix.conjugate().transpose().dot(self.l_block.cm_dn.dot(t_matrix))
        self.l_block.d=m
        l_block=self.l_block

        return self.error

class Infinite_System(object):
    """Infinite chain dmrg algorithm"""
    def __init__(self,L):
        self.block=Block()
        self.L=L
        
    def dmrg(self,m):
        while 2*self.block.l<self.L:
            self.block.enlarge()
            self.super_block=Super_Block(self.block,self.block)
            self.E0=self.super_block.eigen()
            #self.error=self.super_block.truncate(m)

class Finite_System(object):
    """Finite chain dmrg algorithm"""
    def __init__(self,L):
        self.block=Block()
        self.L=L
    
    def dmrg(self,m_start,m_list):
        blocks={}    
        blocks["l",self.block.l]=copy(self.block)
        blocks["r",self.block.l]=copy(self.block)

        while 2*self.block.l<self.L:
            self.block.enlarge()
            self.super_block=Super_Block(self.block,self.block)
            self.E0=self.super_block.eigen()
            #self.error=self.super_block.truncate(m_start)     
            blocks["l",self.block.l]=copy(self.block)
            blocks["r",self.block.l]=copy(self.block)

        l_label,r_label="l","r"
        l_block=self.block
        
        for m in m_list:
            while True:                
                r_block=blocks[r_label,self.L-l_block.l-2]
                if r_block.l==1:
                    rblock=copy(l_block)
                    l_block=copy(r_block)
                    r_block=copy(rblock) 
                    l_label,r_label=r_label,l_label
                l_block.enlarge()
                r_block.enlarge()
                self.super_block=Super_Block(l_block,r_block)
                self.E0=self.super_block.eigen()
                #self.error=self.super_block.truncate(m)

                blocks[l_label,l_block.l]=copy(l_block)

                if l_label=="l" and 2*l_block.l==self.L:
                    break    
    
def finite_dmrg(L,m_start,m_list):
    finite_chain=Finite_System(L)
    finite_chain.dmrg(m_start,m_list)
    print "Finite System(L=%d,dmrg):\n"%L,\
          "E/L:",finite_chain.E0/L;\
          #"Truncation Error:",finite_chain.error
    
def infinite_dmrg(L,m):
    infinite_chain=Infinite_System(L)
    infinite_chain.dmrg(m)
    print "Infinite System(L=%d,dmrg):\n"%L,\
          "E/L:",infinite_chain.E0/L;\
          #"Truncation Error:",infinite_chain.error

def exact(L):
    """Exact solution"""
    E0=0.0
    for i in range(0,L,1):
        Ek=e0+t*np.exp(1j*2*np.pi*i/L)+t*np.exp(-1j*2*np.pi*i/L)
        if Ek<0:
            E0=E0+e0+t*np.exp(1j*2*np.pi*i/L)+t*np.exp(-1j*2*np.pi*i/L)
    print "L=%d(exact):\n"%L,\
          "E/L:",2*E0/L

if __name__=="__main__":
    m_list=[10,20,30,40,40]
    finite_dmrg(6,10,m_list)
    
    infinite_dmrg(6,40)
    
    exact(6)



        
        


        
