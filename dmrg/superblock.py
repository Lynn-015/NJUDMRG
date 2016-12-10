import numpy as np
from scipy.sparse import kron,identity,lil_matrix
from scipy.sparse.linalg import eigsh
from copy import copy,deepcopy

from utils import index_map

class SuperBlock(object):
	def __init__(self,lhgen,rhgen,target_sector=0.,joint=True):
		self.lhgen=deepcopy(lhgen)
		self.rhgen=deepcopy(rhgen)
		self.L=self.lhgen.l+self.rhgen.l
		self.H=kron(self.lhgen.H,identity(self.rhgen.D))+kron(identity(self.lhgen.D),self.rhgen.H)
		if joint==True:
			for lpterm in self.lhgen.pterms:
				for rpterm in self.rhgen.pterms:
					if lpterm.label==rpterm.label: #label must include all the important imformation
						self.H=self.H+kron(lpterm.current_op.mat,rpterm.current_op.mat)*lpterm.param

		self.target_sector=target_sector
		self.sector_indices={}
		self.rsector_indices={}
		self.restricted_basis_indices=[]
		for sys_sec,sys_basis_states in self.lhgen.basis_by_sector.items():
			self.sector_indices[sys_sec]=[]
			env_sec=target_sector-sys_sec
			if env_sec in self.rhgen.basis_by_sector:
				self.rsector_indices[env_sec]=[]
				for i in sys_basis_states:
					i_offset=self.rhgen.D*i
					for j in self.rhgen.basis_by_sector[env_sec]:
						current_index=len(self.restricted_basis_indices)
						self.sector_indices[sys_sec].append(current_index)
						self.rsector_indices[env_sec].append(current_index)
						self.restricted_basis_indices.append(i_offset+j)
		self.restricted_superblock_hamiltonian=self.H.todense()[:,self.restricted_basis_indices][self.restricted_basis_indices,:]
		#self.restricted_superblock_hamiltonian.to? 
					
	def eigen(self,psi0_guess=None):
		if psi0_guess is not None:
			restricted_psi0_guess=psi0_guess[self.restricted_basis_indices]			
		else:
			restricted_psi0_guess=None
			
		energys,restricted_psi0s=eigsh(self.restricted_superblock_hamiltonian,k=1,which="SA",v0=restricted_psi0_guess)
		self.restricted_psi0=restricted_psi0s[:,0]

		self.full_psi0=np.zeros([self.lhgen.D*self.rhgen.D,1],dtype='d')
		for i,z in enumerate(self.restricted_basis_indices):
			self.full_psi0[z,0]=self.restricted_psi0[i]

		if psi0_guess is not None:
			overlap=np.absolute(np.dot(psi0_guess.conjugate().transpose(),self.full_psi0).item())
			overlap/=np.linalg.norm(psi0_guess)*np.linalg.norm(self.full_psi0)
			print "overlap =",overlap

		return energys[0],self.restricted_psi0
	
	def transmat(self,m,use_qn=True):
		rho_block_dict={}
		for sys_sec,indices in self.sector_indices.items():
			if indices:
				psi0_sector=self.restricted_psi0[indices]
				psi0_sector=psi0_sector.reshape([len(self.lhgen.basis_by_sector[sys_sec]),-1],order="C")
				rho_block_dict[sys_sec]=np.dot(psi0_sector,psi0_sector.conjugate().transpose())
				
		possible_eigenstates=[]
		for sector,rho_block in rho_block_dict.items():
			evals,evecs=np.linalg.eigh(rho_block)
			current_sector_basis=self.lhgen.basis_by_sector[sector]
			for eval,evec in zip(evals,evecs.transpose()):
				possible_eigenstates.append((eval,evec,sector,current_sector_basis))
		possible_eigenstates.sort(reverse=True,key=lambda x:x[0])
	
		my_m=min(len(possible_eigenstates),m)
		transformation_matrix=lil_matrix((self.lhgen.D,my_m),dtype='d')
		self.new_sector_array=np.zeros((my_m,),dtype='d')

		possible_eigenstates=possible_eigenstates[:my_m]
		possible_eigenstates.sort(reverse=False,key=lambda x:x[3])
		for i,(eval,evec,sector,current_sector_basis) in enumerate(possible_eigenstates):
			for j,v in zip(current_sector_basis,evec):
				transformation_matrix[j,i]=v
			self.new_sector_array[i]=sector
		self.new_basis_by_sector=index_map(self.new_sector_array)

		transformation_matrix=transformation_matrix.tocsr()

		self.s=[]
		for i in range(my_m):
			self.s.append((possible_eigenstates[i][0])**0.5)

		return transformation_matrix
		
	def rtransmat(self,m,use_qn=True): #should be discarded
		rho_block_dict={}
		for env_sec,indices in self.rsector_indices.items():
			if indices: #
				psi0_sector=self.restricted_psi0[indices]
				psi0_sector=psi0_sector.reshape([-1,len(self.rhgen.basis_by_sector[env_sec])],order="C")
				rho_block_dict[env_sec]=np.dot(psi0_sector.transpose(),psi0_sector.conjugate())

		possible_eigenstates=[]
		for sector,rho_block in rho_block_dict.items():
			evals,evecs=np.linalg.eigh(rho_block)
			current_sector_basis=self.rhgen.basis_by_sector[sector]
			for eval,evec in zip(evals,evecs.transpose()):
				possible_eigenstates.append((eval,evec,sector,current_sector_basis))
		possible_eigenstates.sort(reverse=True,key=lambda x:x[0])

		my_m=min(len(possible_eigenstates),m)
		transformation_matrix=lil_matrix((self.rhgen.D,my_m),dtype='d') #
		self.rnew_sector_array=np.zeros((my_m,),dtype='d')
		
		possible_eigenstates=possible_eigenstates[:my_m]
		possible_eigenstates.sort(reverse=False,key=lambda x:x[3])
		for i,(eval,evec,sector,current_sector_basis) in enumerate(possible_eigenstates):
			for j,v in zip(current_sector_basis,evec):
				transformation_matrix[j,i]=v #
			self.rnew_sector_array[i]=sector
		self.rnew_basis_by_sector=index_map(self.rnew_sector_array)
		
		transformation_matrix=transformation_matrix.tocsr()
		return transformation_matrix
