'''
one-dimensional hamiltonian generator for dmrg
'''
import numpy as np
from scipy.sparse import kron,identity
from copy import deepcopy

class HGen(object):
	'''
	father class of LHGen and RHGen
	
	attributes:
	l:length of chain
	d:dimension of single site hilbert space
	D:dimension of total hibert space
	H:hamiltonian matrix
	pterms:terms need to be handled
	terms:hamiltonian terms
	'''
	def __init__(self,terms,d):
		self.l=1
		self.d=d
		self.D=self.d
		self.H=np.zeros([self.d,self.d])
		self.terms=terms
		self.pterms=[]

class LHGen(HGen):
	'''
	left hamiltonian generator
	construct:LHGen(terms,d)

	attributes:
	same as HGen
	'''
	def __init__(self,terms,d):
		HGen.__init__(self,terms,d)
		for term in self.terms:
			if term.op2 is None and (term.op1.site is None or term.op1.site==1): #onsite term
				self.H=self.H+term.op1.mat*term.param
			elif term.op2 is not None and (term.op1.site is None or term.op1.site==1):
				pterm=deepcopy(term)
				pterm.op1.site=1
				pterm.op2.site=1+pterm.dist
				self.pterms.append(pterm)

	def enlarge(self):
		'''enlarge one site towards right'''
		self.l=self.l+1
		self.H=kron(self.H,identity(self.d))
		pts=[]
		for pterm in self.pterms:
			if pterm.op2.site==self.l:
				self.H=self.H+kron(pterm.op1.mat,pterm.op2.mat)*pterm.param
			else:
				pterm.op1.mat=kron(pterm.op1.mat,identity(self.d))
				pts.append(pterm)
		self.pterms=deepcopy(pts)
		for term in self.terms:
			if term.op2 is None and (term.op1.site is None or term.op1.site==self.l):
				self.H=self.H+kron(identity(D),term.op1.mat)*term.param
			elif term.op2 is not None and (term.op1.site is None or term.op1.site==self.l):
				pterm=deepcopy(term)
				pterm.op1.mat=kron(identity(self.D),pterm.op1.mat)
				pterm.op1.site=self.l
				pterm.op2.site=self.l+pterm.dist
				self.pterms.append(pterm)
		self.D=self.D*self.d

	def truncate(self,U):
		'''truncate H and part_terms with U'''
		self.H=U.conjugate().transpose().dot(self.H.dot(U))
		for pterm in self.pterms:
			pterm.op1.mat=U.conjugate().transpose().dot(pterm.op1.mat.dot(U))
		self.D=self.H.shape[0] 

class RHGen(HGen):
	'''
	right hamiltonian generator
	construct:RHGen(terms,d,N)

	attributes:
	N:length of whole chain
	same as HGen 
	'''
	def __init__(self,terms,d,N): #mirror image not used
		HGen.__init__(self,terms,d)
		self.N=N
		for term in self.terms:
			if term.op1 is None and (term.op2.site is None or term.op2.site==N):
				self.H=self.H+term.op2.mat*term.param
			elif term.op1 is not None and (term.op2.site is None or term.op2.site==N):
				pterm=deepcopy(term)
				pterm.op2.site=N
				pterm.op1.site=N-pterm.dist
				self.pterms.append(pterm)

	def enlarge(self):
		self.l=self.l+1
		self.H=kron(identity(self.d),self.H)
		pts=[]
		for pterm in self.pterms:
			if pterm.op1.site==self.N-self.l+1:
				self.H=self.H+kron(pterm.op1.mat,pterm.op2.mat)*pterm.param
			else:
				pterm.op2.mat=kron(identity(self.d),pterm.op2.mat)
				pts.append(pterm)
		self.pterms=deepcopy(pts)
		for term in self.terms:
			if term.op1 is None and (term.op2.site is None or term.op2.site==self.N-self.l+1):
				self.H=self.H+kron(term.op2.mat,identity(self.D))*term.param
			elif term.op1 is not None and (term.op2.site is None or term.op2.site==self.N-self.l+1):
				pterm=deepcopy(term)
				pterm.op2.mat=kron(pterm.op2.mat,identity(self.D))
				pterm.op2.site=self.N-self.l+1
				pterm.op1.site=self.N-self.l+1-pterm.dist
				self.pterms.append(pterm)
		self.D=self.D*self.d
	
	def truncate(self,V):
		self.H=V.conjugate().transpose().dot(self.H.dot(V))
		for pterm in self.pterms:
			pterm.op2.mat=V.conjugate().transpose().dot(pterm.op2.mat.dot(V))
		self.D=self.H.shape[0]

class SuperBlock(object):
	'''
	super block hamiltonian generator

	attributes:
	lhgen:left hamiltonian generator
	rhgen:right hamiltonian generator
	H:super block hamiltonian matrix
	'''
	def __init__(self,lhgen,rhgen):
		self.lhgen=deepcopy(lhgen)
		self.rhgen=deepcopy(rhgen)
		self.L=self.lhgen.l+self.rhgen.l
		self.H=kron(self.lhgen.H,identity(self.rhgen.D))+kron(identity(self.lhgen.D),self.rhgen.H)
		for lpterm in self.lhgen.pterms:
			for rpterm in self.rhgen.pterms:
				if (lpterm.op1.label,lpterm.op2.label)==(rpterm.op1.label,rpterm.op2.label) and (lpterm.op1.site+self.rhgen.N-self.L,lpterm.op2.site+self.rhgen.N-self.L)==(rpterm.op1.site,rpterm.op2.site): 
					#this doesn't work in mirror image
					self.H=self.H+kron(lpterm.op1.mat,rpterm.op2.mat)*lpterm.param
