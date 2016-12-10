import numpy as np
from scipy.sparse import kron,identity

from mpo import MPO,op2mpo

class TestMPO(object):
	def __init__(self):
		self.op=np.array([[0.5,0.],[0.,-0.5]])
		for i in range(5):
			self.op=kron(self.op,identity(2))
		self.mpo=op2mpo(self.op,2,6)

	def shape(self):
		for i in range(6):
			print self.mpo.Ws[i].shape

if __name__=='__main__':
	tmpo=TestMPO()
	tmpo.shape()
