import numpy as np

class MPS(object):
	def __init__(self,matrices=[]):
		self.matrices=matrices
		self.L=len(self.matrices)
			
	@property
	def state(self):
		state=self.matrices[0]
		for in in range(1,self.L):
			state.dot(self.matrices[i])
		return state.reshape(-1,1)
		

def state2mps(state,d,L):
	'''
	state is a state ket
	assume it is already in the 'right' base
	'''
	Ms=[]
	for i in range(1,L):
		psi=state.reshape(d**i,d**(L-i))
		U,S,Vdag=np.linalg.svd(psi,full_matrices=False)
		M=U.reshape(-1,d,U.shape[1])).transpose(1,0)
		Ms.append(M)
	Ms.append(S.dot(Vdag))
	return Ms
