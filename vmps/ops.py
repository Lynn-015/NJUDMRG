import numpy as np
from scipy.sparse import kron
from mpo import MPO,op2mpo,add

class OpUnit(object):
	def __init__(self,mat,site=None):
		self.mat=mat
		self.d=self.mat.shape[0]
		self.site=site

class OpString(object):
	def __init__(self,ops,L,param=1.):
		self.ops=ops
		self.L=L
		self.d=self.ops[0].d
		self.param=param
		
		self.mat=np.array([[1.]])
		for i in range(L):
			siteop=False
			for op in self.ops:
				if op.site==i:
					self.mat=kron(self.mat,op.mat)
					siteop=True
			if siteop==False:
				self.mat=kron(self.mat,identity(self.d)
		self.mat*=self.param
	
	def tompo(self):
		self.mpo=op2mpo(self.mat)
		
class OpCollection(object):
	def __init__(self,opstrs):
		self.opstrs=opstrs
		self.mpo=self.opstrs[0].mpo
		for i in range(1,len(self.opstrs)):
			self.mpo=add(self.mpo,self.opstrs[i].mpo)
		

