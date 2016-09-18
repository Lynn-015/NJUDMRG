from __future__ import division

import numpy as np
from scipy.sparse.linalg import eigsh

from hgen import NRGHGen,LHGen,RHGen,Term
from nrg import NRGEngine
from dmrg import DMRGEngine

J=0.1
Sp=np.array([[0,1],[0,0]])
Sm=np.array([[0,0],[1,0]])
Sz=np.array([[1,0],[0,-1]])*0.5

terms=[Term(op1=Sp,op2=Sm,param=J/2,dist=1,label=['Sp','Sm']),Term(op1=Sm,op2=Sp,param=J/2,dist=1,label=['Sm','Sp']),Term(op1=Sz,op2=Sz,param=J,dist=1,label=['Sz','Sz'])]

hgen=NRGHGen(terms)
lhgen=LHGen(terms)
rhgen=RHGen(terms,5)

def test():
	for i in range(5):
		hgen.enlarge()

	val,vec=eigsh(hgen.H,1,which='SA')
	print val

def testNRG():
	nrg=NRGEngine(hgen)
	nrg.full_sweep(N=5,m=10)

def testDMRG():
	dmrg=DMRGEngine(lhgen,rhgen)
	dmrg.infinite(m=10)
	dmrg.finite(mwarmup=10,mlist=[10,20,30,40])

if __name__=='__main__':
	#test()
	#testNRG()
	testDMRG()
