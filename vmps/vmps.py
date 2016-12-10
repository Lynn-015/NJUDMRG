import numpy as np
from scipy.sparse.linalg import eigsh

class VMPSEngine(object):
	def __init__(self,Hmpo,gmps):
		self.Hmpo=deepcopy(Hmpo)
		self.mps=deepcopy(gmps)
		self.Ls=[]
		self.Rs=[]

	def contract(self): #calculate Rs
		self.mps.contract_s() #right canonicalized
		R=np.array([[[1.]]])
		for l in range(self.Hmpo.L-2):
			r=np.tensordot(self.mps.Ms[-l-1].conjugate().transpose(),np.tensordot(self.Hmpo.Ws[-l-1],self.mps.Ms[-l-1],axes=(2,1)),axes=(1,1))
			R=np.tensordot(r,R,axes=([1,3,5],[0,1,2]))
			self.Rs.append(R)
		self.Rs.reverse()

	def right_sweep(self,m):
		L=np.array([[[1.]]])
		for l in range(self.Hmpo.L-1):
			self.H=np.tensordot(L,tensordot(self.Hmpo.Ws[l],np.tensordot(self.Hmpo.Ws[l+1],self.Rs[0],axes=(3,1)),axes=(3,0)),axes=(1,0))
			self.H.transpose(0,2,4,6,1,3,5,7)
			self.H.reshape(self.H.shape(0)*self.Hmpo.d**2*self.H.shape[3])
			lamda,v=eigsh(self.H,1,which='SA') #quantum number 0 block
			M=v.reshape(-1,self.Hmpo.d*self.Rs[0].shape[2]) #add zero block?
			self.Rs=self.Rs[1:]
			U,S,Vdag=np.linalg.svd(M,full_matrices=False)
			self.mps.Ms[i]=U.reshape(-1,self.Hmpo.d,U.shape[-1]) #blockize?
			S=S[:m] #renormalize?
			B=np.tensordot(np.diag(S),Vdag,1)
			U,S,Vdag=np.linalg.svd(B,fullmatrices=False) #blockize?
			self.mps.Ms[i+1]=U.reshape(U.shape[0],self.Hmpo.d,-1)
			self.mps.Ms[i+2]=np.tensordot(np.tensordot(np.diag(S),Vdag,1),self.mps.Ms[i+2],1)
			l=np.tensordot(self.mps.Ms[i].conjugate().transpose(),np.tensordot(self.Hmpo.Ws[i],self.mps.Ms[i],axes=(2,1),axes=(1,1)))
			L=np.tensordot(L,l,axes=([0,1,2],[0,2,4]))
			self.Ls.append(L)
	
	def left_sweep(self):
		R=np.array([[[1.]]])
		for l in range(self.Hmpo.L-1):
			

