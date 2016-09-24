from __future__ import division

import numpy as np
from scipy.sparse.linalg import eigsh

from hgen import LHGen,RHGen
from dmrg import DMRGEngine
from ops import SpinOp,BoseOp,Op,SpinTerm,BoseTerm

J=1.
t=-0.1
U=1

cdag=BoseOp(label='c+',mat=np.array([[0,0],[1,0]]))
c=BoseOp(label='c',mat=np.array([[0,1],[0,0]]))
n=BoseOp(label='n',mat=np.array([[0,0],[0,1]]))

Sp=SpinOp(label='Sp',mat=np.array([[0,1],[0,0]]))
Sm=SpinOp(label='Sm',mat=np.array([[0,0],[1,0]]))
Sz=SpinOp(label='Sz',mat=np.array([[1,0],[0,-1]])*0.5)

I=Op(label='I',mat=np.array([[1,0],[0,1]]))

sterms=[SpinTerm(op1=Sp,op2=Sm,param=J/2),SpinTerm(op1=Sm,op2=Sp,param=J/2),SpinTerm(op1=Sz,op2=Sz,param=J)]
bterms=[BoseTerm(op1=cdag,op2=c,param=t),BoseTerm(op1=c,op2=cdag,param=t)]

lhgen=LHGen(sterms,2)
rhgen=RHGen(sterms,2,10)

blhgen=LHGen(bterms,2)
brhgen=RHGen(bterms,2,20)

def test_spin_idmrg():
	dmrg=DMRGEngine(lhgen,rhgen)
	dmrg.infinite(m=20)

def test_spin_fdmrg():
	dmrg=DMRGEngine(lhgen,rhgen)
	dmrg.finite(mwarmup=10,mlist=[10,20,30,40,40])

def test_bose_idmrg():
	dmrg=DMRGEngine(blhgen,brhgen)
	dmrg.infinite(m=5)

def test_bose_fdmrg():
	dmrg=DMRGEngine(blhgen,brhgen)
	dmrg.finite(mwarmup=5,mlist=[5,10,15,20])

#def test_bosespin_idmrg():

if __name__=='__main__':
	test_spin_idmrg()
	#test_spin_fdmrg()
	#test_bose_idmrg()
	#test_bose_fdmrg()
