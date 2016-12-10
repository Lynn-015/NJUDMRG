'''operator and hamiltonian terms'''
import numpy as np
from copy import deepcopy
from scipy.sparse import kron,identity

Z=np.array([[-1,0],[0,1]])
Zs=np.array([[1,0,0,0],[0,-1,0,0],[0,0,-1,0],[0,0,0,1]])

class Op(object):
	'''
	single operator,can be used for spin and bose operator
	in the future,fermi operator and different representation may be considered
	
	attributes:
	mat:matrix of the operator
	label:label of the operator,str
	site:site of the operator
	spin:spin of the operator,'up'/'dn' or None
	'''
	def __init__(self,mat,label=None,site=None): #should site and spin be attributes of Op or Onsite/Twosite?
		self.mat=mat
		self.label=label
		self.site=site

class Term(object): #spin and bose term
	def __init__(self,ops=[],param=1.,dists=[1],label=None):
		self.ops=deepcopy(ops)
		self.param=param
		self.dists=dists
		self.label=label

class FTerm(object):
	def __init__(self,ops=[],param=1.,dists=[1],label=None):
		self.ops=deepcopy(ops)
		self.param=param
		self.dists=dists
		self.label=label
		for i in range(len(self.ops)):
			if (len(self.ops)-i)%2==0:
				self.ops[i].mat=self.ops[i].mat.dot(Z)	

class SFTerm(object):
	def __init__(self,ops=[],param=1.,dists=[1],label=None):
		self.ops=deepcopy(ops)
		self.param=param
		self.dists=dists
		self.label=label
		for i in range(len(self.ops)):
			if (len(self.ops)-i)%2==0:
				self.ops[i].mat=self.ops[i].mat.dot(Zs)	

oplib={}

oplib['s+']=Op(label='S+',mat=np.array([[0,1],[0,0]]))
oplib['s-']=Op(label='S-',mat=np.array([[0,0],[1,0]]))
oplib['sz']=Op(label='Sz',mat=np.array([[1,0],[0,-1]])*0.5)

oplib['c+']=Op(label='c+',mat=np.array([[0,1],[0,0]]))
oplib['c']=Op(label='c',mat=np.array([[0,0],[1,0]]))
oplib['n']=Op(label='n',mat=np.array([[1,0],[0,0]]))

#spin fermi
oplib['Cup+']=Op(label='Cup+',mat=np.array([[0,0,1,0],[0,0,0,1],[0,0,0,0],[0,0,0,0]]))
oplib['Cup']=Op(label='Cup',mat=np.array([[0,0,0,0],[0,0,0,0],[1,0,0,0],[0,1,0,0]]))
oplib['Cdn+']=Op(label='Cdn+',mat=np.array([[0,-1,0,0],[0,0,0,0],[0,0,0,1],[0,0,0,0]]))
oplib['Cdn']=Op(label='Cdn',mat=np.array([[0,0,0,0],[-1,0,0,0],[0,0,0,0],[0,0,1,0]]))

#spin bose



