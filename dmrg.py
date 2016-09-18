'''
dmrg
'''
import numpy as np
from scipy.sparse.linalg import eigsh
from copy import deepcopy

from hgen import LHGen,SuperBlock,Term

class DMRGEngine(object):
	'''
	dmrg engine
	construct:DMRGEngine(lhgen,rhgen)
	
	attributes:
	lhgen:left hamiltonian generator
	rhgen:right hamiltonian generator
	N:length of whole chain
	lblocks:a list to store left generators
	rblocks:a list to store right generators
	sblock:super block	
	'''
	def __init__(self,lhgen,rhgen):
		self.lhgen=lhgen
		self.rhgen=rhgen
		self.N=rhgen.N
		self.lblocks=[deepcopy(self.lhgen)]
		self.rblocks=[deepcopy(self.rhgen)]

	def single_step(self,m): #only for infinite algorithm
		'''single dmrg step.m:number of kept states'''
		self.lhgen.enlarge()
		self.rhgen.enlarge()
		self.sblock=SuperBlock(self.lhgen,self.rhgen)
		val,vec=eigsh(self.sblock.H,1,which='SA')
		print val
		if self.lhgen.H.shape[0]>m:
			psi=vec.reshape([self.lhgen.D,self.rhgen.D])
			phoA=psi.dot(psi.conjugate().transpose())
			vals,vecs=eigsh(phoA,m,which='SA')
			U=vecs
			self.lhgen.truncate(U)
		self.lblocks.append(deepcopy(self.lhgen))
		self.rblocks.append(deepcopy(self.rhgen))

	def infinite(self,m):
		'''infinite algorithm'''
		for i in range(self.N-1):
			self.single_step(m)
	
	def right_sweep(self,m):
		'''sweep one site towards right'''
		self.lhgen.enlarge()
		self.rblocks.pop(-1)
		self.sblock=SuperBlock(self.lhgen,self.rblocks[-1])
		val,vec=eigsh(self.sblock.H,1,which='SA')
		print val
		if self.lhgen.D>m:
			psi=vec.reshape([self.lhgen.D,self.sblock.rhgen.D])
			phoA=psi.dot(psi.conjugate().transpose())
			vals,vecs=eigsh(phoA,m,which='SA')
			U=vecs
			self.lhgen.truncate(U)
		self.lblocks.append(deepcopy(self.lhgen))

	def left_sweep(self,m):
		'''sweep one site towards left'''
		self.rhgen.enlarge()
		self.lblocks.pop(-1)
		self.sblock=SuperBlock(self.lblocks[-1],self.rhgen)
		val,vec=eigsh(self.sblock.H,1,which='SA')
		print val
		if self.rhgen.D>m:
			psi=vec.reshape([self.sblock.lhgen.D,self.rhgen.D])
			phoA=psi.transpose().dot(psi.conjugate())
			vals,vecs=eigsh(phoA,m,which='SA')
			V=vecs
			self.rhgen.truncate(V)
		self.rblocks.append(deepcopy(self.rhgen))

	def finite(self,mwarmup,mlist):
		'''finite algorithm'''
		self.infinite(mwarmup) #mind the initialize problem
		for m in mlist:
			for i in range(self.N/2,self.N-1):
				self.right_sweep(m)
			for i in range(1,self.N-1):
				self.left_sweep(m)
			for i in range(1,self.N/2):
				self.right_sweep(m)

