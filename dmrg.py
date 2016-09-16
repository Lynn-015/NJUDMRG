'''
dmrg
'''
import numpy as np
from scipy.sparse.linalg import eigsh

from hgen import LHGen,SuperBlock,Term

class DMRGEngine(object):
	'''
	dmrg engine
	'''
	def __init__(self,lhgen):
		self.lhgen=lhgen

	def single_step(self,m): #only for infinite algorithm
		self.lhgen.enlarge()
		self.sblock=SuperBlock(self.lhgen,self.lhgen)
		val,vec=eigsh(self.sblock.H,1,which='SA')
		print val
		if self.lhgen.H.shape[0]>m:
			psi=vec.reshape([self.lhgen.D,self.lhgen.D])
			phoA=psi.dot(psi.conjugate().transpose())
			vals,vecs=eigsh(phoA,m,which='SA')
			U=vecs
			self.lhgen.truncate(U)

	def infinite(self,N,m):
		for i in range(N-1):
			self.single_step(m)
