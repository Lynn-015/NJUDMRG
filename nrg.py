'''
one dimensional nrg,or block decimation
'''
import numpy as np
from scipy.sparse.linalg import eigsh

from hgen import NRGHGen

class NRGEngine(object):
	'''
	one dimensional NRG Engine
	'''
	def __init__(self,hgen):
		self.hgen=hgen

	def single_step(self,m):
		self.hgen.enlarge()
		if self.hgen.H.shape[0]>m:
			vals,vecs=eigsh(self.hgen.H,m,which='SA')
			U=vecs
			self.hgen.truncate(U)
			print vals[0]

	def full_sweep(self,N,m):
		for i in range(N-1):
			self.single_step(m)

			
