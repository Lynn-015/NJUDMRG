import numpy as np
from scipy.sparse import kron,identity
from scipy.sparse import kron

from ops import OpUnit,OpString,OpCollection
'''heisenberg model'''

L=6
J=1.
sp=np.array([[0,1],[0,0]])
sm=np.array([[0,0],[1,0]])
sz=np.array([[1,0],[0,-1]])*0.5

opstrs=[]
for i in range(L-1):
	opstr=OpString([OpUnit(sp,site=i),OpUnit(sm,site=i+1)])
	opstr2=
	opstr3=
	opstrs.extend([opstr,opstr2,opstr3])

opcol=OpCollenction(opstrs)	

up=np.array([1.,0])
dn=np.array([0,1.])

state=up
for i in range(L-1):
	state=kron(state,up)
	
gmps=ket2mps(state)
