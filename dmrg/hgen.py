import numpy as np
from scipy.sparse import kron,identity
from copy import copy,deepcopy

from ops import Z,Zs
from utils import index_map

class HGen(object):
	def __init__(self,terms,L,d=2,part='left',fermi=False,sfermi=False,sectors=np.array([0.5,-0.5])):
		self.l=1;self.d=d;self.D=self.d
		self.H=np.zeros([self.d,self.d])
		self.terms=terms;self.pterms=[]
		self.L=L;self.part=part
		self.single_site_sectors=sectors
		self.basis_sector_array=copy(self.single_site_sectors)
		self.basis_by_sector=index_map(self.basis_sector_array)
		
		if fermi==True:
			self.I=Z
		elif sfermi==True:
			self.I=Zs
		else:
			self.I=identity(self.d)

		if self.part=='left':
			for term in self.terms:
				if len(term.ops)==1 and (term.ops[0].site is None or term.ops[0].site==1):
					self.H+=term.ops[0].mat*term.param
				elif len(term.ops)>1 and (term.ops[0].site is None or term.ops[0].site==1):
					pterm=deepcopy(term)
					pterm.ops[0].site=1
					pterm.current_index=0
					pterm.current_op=deepcopy(pterm.ops[0])
					for i in range(len(pterm.dists)): #if sites are given,this step can be skipped
						pterm.ops[i+1].site=pterm.dists[i]+pterm.ops[i].site
					self.pterms.append(pterm)
		else:
			for term in self.terms:
				if len(term.ops)==1 and (term.ops[-1].site is None or term.ops[-1].site==self.L):
					self.H+=term.ops[-1].mat*term.param
				elif len(term.ops)>1 and (term.ops[-1].site is None or term.ops[-1].site==self.L):
					pterm=deepcopy(term)
					pterm.ops[-1].site=self.L
					pterm.current_index=len(term.ops)-1
					pterm.current_op=deepcopy(pterm.ops[-1])
					for i in range(len(pterm.dists)):
						pterm.ops[-i-2].site=pterm.ops[-i-1].site-pterm.dists[-i-1]
					self.pterms.append(pterm)
		
	def enlarge(self):
		self.l+=1
		if self.part=='left':
			self.H=kron(self.H,identity(self.d))
			pts=[]
			for pterm in self.pterms:
				if pterm.ops[pterm.current_index+1].site==self.l:
					pterm.current_index+=1					
					pterm.current_op.mat=kron(pterm.current_op.mat,pterm.ops[pterm.current_index].mat) #other attribute?
				else:
					pterm.current_op.mat=kron(pterm.current_op.mat,self.I)
				if pterm.current_index<len(pterm.ops)-1:
					pts.append(pterm)
				else:
					self.H+=pterm.current_op.mat*pterm.param
			
			self.pterms=deepcopy(pts)
			for term in self.terms:
				if len(term.ops)==1 and (term.ops[0].site is None or term.ops[0].site==self.l):
					self.H+=kron(identity(self.D),term.ops[0].mat)*term.param
				elif len(term.ops)>1 and (term.ops[0].site is None or term.ops[0].site==self.l):
					pterm=deepcopy(term)
					pterm.current_index=0
					pterm.current_op=deepcopy(pterm.ops[0])
					pterm.current_op.mat=kron(identity(self.D),pterm.current_op.mat)
					pterm.ops[0].site=self.l
					for i in range(len(pterm.dists)):
						pterm.ops[i+1].site=pterm.dists[i]+pterm.ops[i].site
					self.pterms.append(pterm) 

			self.basis_sector_array=np.add.outer(self.basis_sector_array,self.single_site_sectors).flatten()

		else:
			self.H=kron(identity(self.d),self.H)
			pts=[]
			for pterm in self.pterms:
				if pterm.ops[pterm.current_index-1].site==self.L-self.l+1:
					pterm.current_index-=1
					pterm.current_op.mat=kron(pterm.ops[pterm.current_index].mat,pterm.current_op.mat)
				else:
					pterm.current_op.mat=kron(self.I,pterm.current_op.mat)
				if pterm.current_index>0:
					pts.append(pterm)
				else:
					self.H+=pterm.current_op.mat*pterm.param
			
			self.pterms=deepcopy(pts)
			for term in self.terms:
				if len(term.ops)==1 and (term.ops[-1].site is None or term.ops[-1].site==self.L-self.l+1):
					self.H+=kron(term.ops[-1].mat,identity(self.D))*term.param
				elif len(term.ops)>1 and (term.ops[-1].site is None or term.ops[-1].site==self.L-self.l+1):
					pterm=deepcopy(term)
					pterm.current_index=len(pterm.ops)-1
					pterm.current_op=deepcopy(pterm.ops[-1])
					pterm.current_op.mat=kron(pterm.current_op.mat,identity(self.D))
					pterm.ops[-1].site=self.L-self.l+1
					for i in range(len(pterm.dists)):
						pterm.ops[-i-2].site=pterm.ops[-i-1].site-pterm.dists[-i-1]
					self.pterms.append(pterm)

			self.basis_sector_array=np.add.outer(self.single_site_sectors,self.basis_sector_array).flatten()

		self.basis_by_sector=index_map(self.basis_sector_array)
		self.D*=self.d

	def transform(self,T):
		self.H=T.conjugate().transpose().dot(self.H.dot(T))
		for pterm in self.pterms:
			pterm.current_op.mat=T.conjugate().transpose().dot(pterm.current_op.mat.dot(T))
		self.D=self.H.shape[0]
