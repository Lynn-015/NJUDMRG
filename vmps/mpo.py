import numpy as np
from copy import deepcopy

from mps import MPS

class MPO(object):
	def __init__(self,d,L,Ws=[]):
		self.Ws=Ws 
		self.d=d
		self.D=self.d**2
		self.L=L		

	def compress(self,m): #remain m largest value
		o=np.array([[[1.]]])
		for l in range(self.L):
			self.Ws[l]=np.tensordot(o,self.Ws[l],1)
			self.Ws[l].reshape(-1,self.Ws[l][-1]) #blockize
			U,S,Vdag=np.svd(self.Ws[l],full_matrices=False)
			W=U.reshape(-1,self.d,self.d,U.shape[-1])
			S.sort()
			S=S[:m] #renormalize?		
			o=np.tensordot(np.diag(S),Vdag,1)

def op2mpo(op,d,L): #opstring to mpo,opsting is an array
	o=deepcopy(op)
	D=d**2

	a=1;Ws=[]
	for i in range(L):
		O=o.reshape((a,d,d**(L-i-1),d,d**(L-i-1)))
		O.transpose((0,1,3,2,4)).reshape(-1,d**(2*(L-i-1)))
		#blockize? seems not
		U,S,Vdag=np.linalg.svd(O,full_matrices=False)
		W=U.reshape(-1,d,d,U.shape[-1]) #blockize seem no need since add and compress
		Ws.append(W)
		o=np.tensordot(np.diag(S),Vdag,1)
		a=S.shape[0]
	Ws[-1]*=float(o)
	return MPO(d,L,Ws)

def add(mpo1,mpo2):
	Ws=[]
	for l in range(mpo1.L):
		for i in range(mpo1.d):
			for j in range(mpo1.d):
				if l==0:
					W[:,i,j,:]=np.hstack(mpo1.Ws[l][:,i,j,:],mpo2.Ws[l][:,i,j,:])
				elif l==mpo1.L-1:
					W[:,i,j,:]=np.vstack(mpo1.Ws[l][:,i,j,:],mpo2.Ws[l][:,i,j,:])
				else:
					W[:,i,j,:]=block_diag(mpo1.Ws[l][:,i,j,:],mpo2.Ws[l][:,i,j,:])
				Ws.append(W)
	mpo=MPO(mpo1.d,mpo1.L,Ws)
	return mpo
'''
def act(mpo,mps):
	Ns=[]
	tmps=deepcopy(mps) #target mps
	tmps.contract_s() #contract S to A or B
	for i in range(mpo.L):
		N=contract(mpo.Ws[i],tmps.Ms[i],axes=([2,],[1,])) #details later
		N.transpose((0,3,1,2)).reshape(N.shape[0]*N.shape[1],mpo.d,-1)
		Ns.append(N)
		rmps.Ms=Ns #result mps
	return rmps
'''
