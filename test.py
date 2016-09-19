from __future__ import division

import numpy as np
from scipy.sparse.linalg import eigsh

from hgen import LHGen,RHGen,SpinTerm
from dmrg import DMRGEngine

J=0.1
Sp=np.array([[0,1],[0,0]])
Sm=np.array([[0,0],[1,0]])
Sz=np.array([[1,0],[0,-1]])*0.5

terms=[SpinTerm(op1=Sp,op2=Sm,param=J/2,label=['Sp','Sm']),SpinTerm(op1=Sm,op2=Sp,param=J/2,label=['Sm','Sp']),SpinTerm(op1=Sz,op2=Sz,param=J,label=['Sz','Sz'])]

lhgen=LHGen(terms,2)
rhgen=RHGen(terms,2,10)

def testidmrg():
	dmrg=DMRGEngine(lhgen,rhgen)
	dmrg.infinite(m=10)

def testfdmrg():
	dmrg=DMRGEngine(lhgen,rhgen)
	dmrg.finite(mwarmup=10,mlist=[10,20,30,40])
if __name__=='__main__':
	#testidmrg()
	testfdmrg()
