import numpy as np
from copy import deepcopy
from scipy.sparse import kron

class Op(object):
	'''
	operator
	in the future,different representation may be considered
	
	attributes:
	label:str,label of the operator
	mat:operator matrix
	site:site of the operator
	'''
	def __init__(self,mat,label='',site=None):
		self.label=label
		self.mat=mat
		self.site=site


class SpinOp(Op):
	'''
	spin operator
	'''
	def __init__(self,mat,label='',site=None):
		Op.__init__(self,mat,label,site=None)

class BoseOp(Op):
	'''
	bose operator

	attributes:
	spin:up or dn
	'''
	def __init__(self,mat,label='',site=None,spin=None):
		Op.__init__(self,mat,label,site=None)
		self.spin=spin

class SpinTerm(object):
	'''
	two site spin-spin term
	support operators on every site and on specific sites

	attributes:
	op1/2:left and right operator 
	dist:site distance between op1 and op2
	param:interaction parameter
	'''
	def __init__(self,param=1.,op1=None,op2=None,dist=1):
		self.param=param
		self.op1=deepcopy(op1)
		self.op2=deepcopy(op2)
		self.dist=dist

class BoseTerm(object): 
	'''
	two site bose term in occupation representation
	'''
	def __init__(self,param=1.,op1=None,op2=None,dist=1):
		self.param=param
		self.op1=deepcopy(op1)
		self.op2=deepcopy(op2)
		self.dist=dist

		if self.dist!=0:			
			if self.op1.spin=='up':
				self.op1.mat=kron(self.op1.mat,identity(self.op1.mat.shape[0]))
			elif self.op1.spin=='dn':
				self.op1.mat=kron(identity(self.op1.mat.shape[0]),self.op1.mat)
		
			if self.op2.spin=='up':
				self.op2.mat=kron(self.op2.mat,identity(self.op2.mat.shape[0]))
			elif self.op2.spin=='dn':
				self.op2.mat=kron(identity(self.op2.mat.shape[0]),self.op2.mat)

		elif self.dist==0 and self.op2 is None:
			if self.op1.spin=='up':
				self.op1.mat=kron(self.op1.mat,identity(self.op1.mat.shape[0]))
			elif self.op1.spin=='dn':
				self.op1.mat=kron(identity(self.op1.mat.shape[0]),self.op1.mat)

		elif self.dist==0 and self.op1 is None:
			if self.op2.spin=='up':
				self.op2.mat=kron(self.op2.mat,identity(self.op2.mat.shape[0]))
			elif self.op2.spin=='dn':
				self.op2.mat=kron(identity(self.op2.mat.shape[0]),self.op2.mat)			
		else:
			self.op1.mat=kron(self.op1.mat,self.op2.mat)
			self.op2=None













