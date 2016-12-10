import numpy as np
from scipy.linalg import kron
from scipy.sparse import identity
from copy import deepcopy

from mps import MPS,ket2mps,overlap,expect
from dmrg import DMRGEngine
from testdmrg import lhgen,rhgen

ket=np.random.rand(2**4)
mps=ket2mps(ket,2,4,cano='mixed',div=2)
ket2=np.random.rand(2**4)
mps2=ket2mps(ket2,2,4,cano='mixed',div=2)

opl=[np.array([[0.5,0.],[0.,-0.5]]),np.array([[1.,0.],[0.,1.]]),np.array([[1.,0.],[0.,1.]]),np.array([[1.,0.],[0.,1.]])]
OP=deepcopy(opl[0])
for i in range(3):
	OP=kron(OP,np.array([[1.,0.],[0.,1.]]))

def test_tomps():
	dmrg=DMRGEngine(lhgen,rhgen)  
	dmrg.finite(mwarmup=10,mlist=[10])
	dmps=dmrg.tomps()
	return dmps

def test_shape(mps):
	for M in mps.Ms:
		print M.shape

def test_shift(mps):
	mps.shift(site=mps.L/2-1)
	test_shape(mps)
	test_cano(mps)
	mps.shift(direct='l',site=mps.L-2)

def test_cano(mps):
	for A in mps.As:
		print np.tensordot(A.conjugate().transpose(),A,axes=([1,2],[1,0]))-identity(A.shape[-1])
	for B in mps.Bs:                                                                                                                                                         
		print np.tensordot(B,B.conjugate().transpose(),axes=([1,2],[1,0]))-identity(B.shape[0])

def test_toket():
	tket=mps.toket()
	print tket-ket	

def test_overlap():
	mps.contract_s()
	mps2.contract_s()
	ovlap=overlap(mps,mps2)
	ov=ket.conjugate().transpose().dot(ket2)
	print ovlap-ov

def test_expect():
	mps.contract_s()
	print expect(mps,opl)
	print ket.conjugate().transpose().dot(OP).dot(ket)


if __name__=='__main__':
	#dmps=test_tomps()
	#test_shape(dmps)
	#test_cano(dmps)
	#test_toket()
	#test_overlap()
	#test_shift(mps)
	test_expect(mps,opl)
                                                                                                                                                                                           
