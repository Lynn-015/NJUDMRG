import numpy as np

from marker import bcadd

class Tensor():
	def __init__(self,array,markers=None):
		self.array=array
		self.shape=self.array.shape
		self.markers=markers #a list of markers for each axes
		self.dim=len(self.shape)		

	def reshape(self,shape,markers=None): #split axes
		self.shape=shape
		self.array.reshape(shape)
		self.markers=markers

	def merge(self,axe1,axe2): #can only merge neighbor axes
		self.markers[axe1]=bcadd(self.markers[axe1],self.markers[axe2])
		del(self.markers[axe2])
		self.shape[axe1]+=self.shape[axe2]
		del(self.shape[axe2])
		self.array.reshape(self.shape)
		
	'''
	def split(self,axe,shape,markers):	#can only split to neighbor axes	
		nshape=()
		nmarkers=[]
		for i in range(self.rank):
			if self.shape[i]!=axe:	
				nshape.append(self.shape[i])
				nmarkers.append(self.markers[i])				
			else:
				nshape.extend(shape)
				nmarkers.extend(markers)
		self.shape=nshape
		self.nmarkers=nmarkers	
	'''

	def transpose(self,order)	 
		self.array.transpose(order)
		markers=[]
		for i in len(order):
			markers.append(deepcopy(self.markers[order[i]]))
		self.markers=markers
		self.shape=self.array.shape

	'''
	def blockize(self): #only blockize a two dimensional matrix
		perm0=self.markers[0].reorder
		perm1=self.markers[1].reorder
		a=np.array(self.shape)
		for i in self.shape[0]:
			for j in self.shape[1]:
				a[i][j]=self.array[perm0[i]][perm1[j]]
		self.array=a

	def take(self,i): #take some subspace(block)
		a=np.array(self.markers[0].sizes[i],self.markers[1].sizes[i])
		a=self.array[self.markers[0].divs[i]:self.markers[0].divs[i]+self.makers[0].sizes[i]][...] #wrong
		return a

def tensordot(t1,t2,axes=None):
	array=np.tensordot(t1.array,t2.array,axes)
	markers=t1.markers.del()+t2.markers.del()
	return Tensor(array,markers)

def add(t1,t2):

def svd(psi):
	
def solve(H): #solve eigen problem
'''	
