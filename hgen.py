'''
one-dimensional hamiltonian generator for dmrg
'''
import numpy as np
from scipy.sparse import kron,identity
from copy import deepcopy

class LHGen(object):
	'''
	left hamiltonian generator
	construct:LHGen(terms)

	attributes:
	l:length of left chain
	d:dimension of single site hilbert space
	D:dimension of total hilbert space
	H:hamiltonian matrix
	part_terms:terms need to be handled
	terms:hamiltonian terms
	'''
	def __init__(self,terms):
		self.l=1
		self.d=terms[0].op1.shape[0]
		self.D=self.d
		self.H=np.zeros([self.d,self.d])
		self.part_terms=[]
		self.terms=terms
		for term in self.terms:
			if term.op2 is None and (term.site1 is None or term.site1==1):
				self.H=self.H+term.op1*term.param
			elif term.op2 is not None and (term.site1 is None or term.site1==1):
				pterm=deepcopy(term)
				pterm.site1=1
				pterm.site2=1+pterm.dist
				self.part_terms.append(pterm)

	def enlarge(self):
		'''enlarge one site towards right'''
		self.l=self.l+1
		self.H=kron(self.H,identity(self.d))
		pts=[]
		for pterm in self.part_terms:
			if pterm.site2==self.l:
				self.H=self.H+kron(pterm.op1,pterm.op2)*pterm.param
			else:
				pterm.op1=kron(pterm.op1,identity(self.d))
				pts.append(pterm)
		self.part_terms=deepcopy(pts)
		for term in self.terms:
			if term.op2 is None and (term.site1 is None or term.site1==self.l):
				self.H=self.H+kron(identity(D),term.op1)*term.param
			elif term.op2 is not None and (term.site1 is None or term.site1==self.l):
				pterm=deepcopy(term)
				pterm.op1=kron(identity(self.D),pterm.op1)
				pterm.site1=self.l
				pterm.site2=self.l+pterm.dist
				self.part_terms.append(pterm)
		self.D=self.D*self.d

	def truncate(self,U):
		'''truncate H and part_terms with U'''
		self.H=U.conjugate().transpose().dot(self.H.dot(U))
		for pterm in self.part_terms:
			pterm.op1=U.conjugate().transpose().dot(pterm.op1.dot(U))
		self.D=self.H.shape[0] 

class RHGen(object):
	'''
	right hamiltonian generator
	construct:RHGen(terms,N)

	attributes:
	l:length of right chain
	N:length of whole chain
	d,D,H,part_terms,terms:same as LHGen

	methods:
	enlarge():enlarge one site towards left
	truncate(V):truncate H and part_terms with V 
	'''
	def __init__(self,terms,N): #mirror image not used
		self.l=1
		self.N=N
		self.d=terms[0].op1.shape[0]
		self.D=self.d
		self.H=np.zeros([self.d,self.d])
		self.part_terms=[]
		self.terms=terms
		for term in self.terms:
			if term.op1 is None and (term.site2 is None or term.site2==N):
				self.H=self.H+term.op2*term.param
			elif term.op1 is not None and (term.site2 is None or term.site2==1):
				pterm=deepcopy(term)
				pterm.site2=N
				pterm.site1=N-pterm.dist
				self.part_terms.append(pterm)

	def enlarge(self):
		self.l=self.l+1
		self.H=kron(identity(self.d),self.H)
		pts=[]
		for pterm in self.part_terms:
			if pterm.site1==self.N-self.l+1:
				self.H=self.H+kron(pterm.op2,pterm.op1)*pterm.param
			else:
				pterm.op2=kron(identity(self.d),pterm.op1)
				pts.append(pterm)
		self.part_terms=deepcopy(pts)
		for term in self.terms:
			if term.op1 is None and (term.site2 is None or term.site2==self.N-self.l+1):
				self.H=self.H+kron(term.op2,identity(self.D))*term.param
			elif term.op1 is not None and (term.site2 is None or term.site2==self.N-self.l+1):
				pterm=deepcopy(term)
				pterm.op2=kron(pterm.op2,identity(self.D))
				pterm.site2=self.N-self.l+1
				pterm.site1=self.N-self.l+1-pterm.dist
				self.part_terms.append(pterm)
		self.D=self.D*self.d
	
	def truncate(self,V):
		self.H=V.conjugate().transpose().dot(self.H.dot(V))
		for pterm in self.part_terms:
			pterm.op2=V.conjugate().transpose().dot(pterm.op2.dot(V))
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
		self.H=kron(self.lhgen.H,identity(self.rhgen.D))+kron(identity(self.lhgen.D),self.rhgen.H)
		for lpterm in self.lhgen.part_terms:
			for rpterm in self.rhgen.part_terms:
				if all(lpterm.label)==all(rpterm.label) and (lpterm.site1+self.rhgen.N-self.lhgen.l-self.rhgen.l,lpterm.site2+self.rhgen.N-self.lhgen.l-self.rhgen.l)==(rpterm.site1,rpterm.site2): 
					#this doesn't work in mirror image
					self.H=self.H+kron(lpterm.op1,rpterm.op2)*lpterm.param

class Term(object):
	'''
	hamiltonian terms.just support simple one and two site operator at now
	besides,it support operators on every site and on specific sites

	attributes:
	op1/2:left and right operator 
	dist:site distance between op1 and op2
	param:interaction parameter
	site1/2:site of op1/2
	label:['op1','op2'],used to form superblock
	'''
	def __init__(self,op1=None,op2=None,site1=None,site2=None,dist=None,param=1.,label=['','']):
		self.op1=op1
		self.op2=op2 #if op2==None,it's an onsite term
		self.dist=dist
		self.param=param
		self.site1=site1 #if site1 and site2 are all None,it means that this term is on every site
		self.site2=site2
		self.label=label
