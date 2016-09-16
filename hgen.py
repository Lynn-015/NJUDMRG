'''
one-dimensional hamiltonian generator used for nrg and dmrg
'''
import numpy as np
from scipy.sparse import kron,identity
from copy import deepcopy

class NRGHGen(object):
	'''
	one dimensional hamiltonian generator for nrg,without considering fermion sign problem 
	
	attributes:
	l:length of chain
	d:dimension of single site hilbert space
	D:dimension of total hilbert space
	H:hamiltonian matrix
	part_terms:terms need to be handled
	terms:all terms	
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
			if term.op2 is None and(term.site1 is None or term.site1==self.l):
				self.H=self.H+kron(identity(D),term.op1)*term.param
			elif term.op2 is not None and((term.site1 is None) or (term.site1==self.l)):
				pterm=copy(term)
				pterm.op1=kron(identity(self.D),pterm.op1)
				pterm.site1=self.l
				pterm.site2=self.l+pterm.dist
				self.part_terms.append(pterm)
		self.D=self.D*self.d
		
	def truncate(self,U):
		self.H=U.conjugate().transpose().dot(self.H.dot(U))
		for term in self.part_terms:
			term.op1=U.conjugate().transpose().dot(term.op1.dot(U))
		self.D=self.H.shape[0]

'''
class RHgen(Hgen):
	def __init__(self,terms,N):
		self.l=1
		self.d=terms[0].op1.shape[0]
		self.D=self.d
		self.H=np.zeros([self.d,self.d])
		self.part_terms=[]
		for terms in self.terms:
			if term.op1==None and(term.site2==None or term.site2==N):
				self.H=self.H+term.op2*term.param
			elif term.op1!=None and(term.site2==None or term.site2==N):
				self.part_terms.append(term)
				self.part_terms[-1].site1=2*l-term.dist
				self.part_terms[-1].site2=2*l
	
	def enlarge(self):
		self.l=self.l+1
		self.H=kron(identity(self.d),self.H)
		for term in self.part_terms:
			if term.site1==1:
				term.op2=kron(term.op1,term.op2)
				self.H=self.H+term.op2*term.param
				self.part_terms.remove(term)
			else:
				term.op2=kron(identity(self.d),term.op2)
		for term in self.terms:
			if term.op1==None and(term.site2==None or term.site1==

class SuperBlock(object):
	def __init__(LHgen,RHgen):
		self.lpterms=self.LHgen.part_terms
		self.rpterms=self.RHgen.part_terms
		self.H=kron(self.LHgen.H,identity(self.RHgen.D))+kron(identity(self.LHgen.D),self.RHgen.H)
		for lterm in self.lpterms:
			for rterm in self.rpterms:
				if lterm.label==rterm.label and...:
					self.H=self.H+kron(self.lterm.op1,self.rterm.op2)*self.lterm.param

				else raise() 
'''

class Term(object):
	'''
	a class to store hamiltonian terms.just support simple one and two site operator at now
	besides,it support operators on every site and on specific sites
	for now,site1 is the current site,site2 is next site

	attributes:
	op1/2:operator on site1/2 
	dist:site distance between op1 and op2
	param:interaction parameter
	site1/2:site of op1/2
	'''
	def __init__(self,op1=None,op2=None,site1=None,site2=None,dist=None,param=1.):
		self.op1=op1
		self.op2=op2 #if op2==None,it's an onsite term
		self.dist=dist
		self.param=param
		self.site1=site1 #if site1 and site2 are all None,it means that this term is on every site
		self.site2=site2
