import numpy as np
from copy import copy,deepcopy
from scipy.sparse import identity

'''all mps tensors has physical indices in the middle'''

class MPS(object):
	def __init__(self,d,L,As=[],Bs=[],S=None): 
		self.As=As
		self.Bs=Bs
		self.S=S #1d array
		self.L=L #length of the chain
		self.d=d #physical dimension
		self.Ms=self.As+self.Bs
	
	def contract_s(self):
		if len(self.As)<len(self.Bs):
			self.Psi=np.tensordot(np.diag(self.S),self.Bs[0],1)
			self.Bs=self.Bs[1:]
		else:
			self.Psi=np.tensordot(self.As[-1],np.diag(self.S),1)
			self.As=self.As[:-1]
		self.Ms=[]
		self.Ms.extend(self.As)
		self.Ms.append(self.Psi)
		self.Ms.extend(self.Bs)
                                                                                                                                                                                       
	def toket(self):
		ket=np.array([[1.]])
		for M in self.Ms:
			ket=np.tensordot(ket,M,1)
		return ket.flatten()		
		
	def shift(self,direct='r',site=1): #shift S 
		if direct=='r':
			for i in range(site):
				self.Bs[0]=np.tensordot(np.diag(self.S),self.Bs[0],1)
				self.Bs[0]=self.Bs[0].reshape(-1,self.Bs[0].shape[-1])
				U,S,Vdag=np.linalg.svd(self.Bs[0],full_matrices=False)
				A=U.reshape(-1,self.d,U.shape[-1])
				self.As.append(A)
				self.S=S
				self.Bs=self.Bs[1:]
				self.Bs[0]=np.tensordot(Vdag,self.Bs[0],1)
		else:
			for i in range(site):
				self.As[-1]=np.tensordot(self.As[-1],np.diag(self.S),1)
				self.As[-1]=self.As[-1].reshape(self.As[-1].shape[0],-1)
				U,S,Vdag=np.linalg.svd(self.As[-1],full_matrices=False)
				B=Vdag.reshape(Vdag.shape[0],self.d,-1)
				self.S=S
				self.As=self.As[:-1]
				self.As[-1]=np.tensordot(self.As[-1],U,1)	

def ket2mps(ket,d,L,cano='mixed',div=None): #dimension of state should be d^L
	mps=MPS(d,L) #L is the length of the chain,not the list
	c=copy(ket)

	if cano=='left':
		mps.cano='left'
		div=L
	elif cano=='right':
		mps.cano='right'
		div=0 #div is the length of As
	
	As=[]
	a=1
	for i in range(div):
		psi=c.reshape(a*d,-1)
		U,S,Vdag=np.linalg.svd(psi,full_matrices=False)
		A=U.reshape(-1,d,U.shape[-1])
		As.append(A)
		c=np.tensordot(np.diag(S),Vdag,1)
		a=S.shape[0]

	Bs=[]
	a=1
	for i in range(L-div):
		psi=c.reshape(-1,a*d)
		U,S,Vdag=np.linalg.svd(psi,full_matrices=False)
		B=Vdag.reshape(Vdag.shape[0],d,-1)
		Bs.append(B)
		c=np.tensordot(U,np.diag(S),1)
		a=S.shape[0]

	if cano=='mixed':
		mps.cano='mixed'
		As[-1]=np.tensordot(As[-1],U,1)		
	mps.As=As
	mps.Bs=Bs
	mps.Bs.reverse()
	mps.S=S
	return mps

#to check the validity of dmrg-tomps,we must calculate some expectation
def overlap(mps1,mps2): #contract s first	
	overlap=np.array([[1.]])
	for i in range(len(mps1.Ms)):
		overlap=np.tensordot(mps1.Ms[i].conjugate().transpose(),np.tensordot(overlap,mps2.Ms[i],1),axes=([1,2],[1,0])) #some index problem
	return float(overlap)

def expect(mps,opl): #opl is a list of matrices including identities.but this is not necessary
	#contract s first
	expect=np.array([[1.]])
	for i in range(len(mps.Ms)):
		overlap=np.tensordot(mps.Ms[i].conjugate().transpose(),np.tensordot(expect,mps.Ms[i],1),1)
		expect=np.tensordot(opl[i],overlap,axes=([0,1],[1,2]))
	return float(expect)

