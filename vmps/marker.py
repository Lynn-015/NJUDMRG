import numpy as np

#change to sparse form later

class Marker(object):
	def __init__(self,qns,sizes=[],divs=[],sign=None,Q='M'):
		self.qns=qns #quantum numbers of each block
		self.sizes=sizes #sizes of each block.meaningless before blockize
		self.divs=divs #divs of each block.meaningless before blockize
		self.sign=sign #to be determined later for each axes
		self.Q=Q

	def reorder(self):
		perm=list(np.range(len(self.qns)))
		qns_perm=[(qn,per) for qn in self.qns for per in perm]
		qns_perm=sorted(qns_perm key=lamda x:x[0])
		self.qns=[]
		perm=[]
		for qn,per in qns_perm:
			self.qns.append(qn)
			perm.append(per)
		return perm
		
	def compact(self): #after reorder
		q=self.qns[0]
		qs=[];qs.append(q)
		self.sizes=[1]
		for qn in self.qns[1:]:
			if qn!=q:
				q=qn;qs.append(qn)
				self.sizes.append(1)
			else:
				self.sizes[-1]+=1
		self.qns=qs	
		self.divs=[]
		size=0
		for i in range(len(self.sizes)-1):
			self.divs.append(size+self.sizes[i])
			size+=self.sizes[i]

def bcadd(mk1,mk2):
	marker=Marker([])
	marker.qns=list((np.array(mk1.qns)[:,np.newaxis]+np.array(marker.qns)).flatten()) #sign?
	return marker

sz_in=Marker([1,-1],[1,1],[1],1)
sz_out=Marker([1,-1],[1,1],[1],-1)

class MarkerGen(object): #marker generator of a  mps/mpo chain
	def __init__(self,L,M=0): #quantum numbers need to be predicted and selected
		self.marker=Marker([0]) #begin marker
		self.M=M
		self.L=L
		self.l=0
		
	def enlarge(self,markers): #enlarge one site
		for marker in markers:
			self.marker.qns=list((np.array(self.marker.qns)[:,np.newaxis]+np.array(marker.qns*marker.sign)).flatten())
		self.l+=1
			
	def filt(self):
		qs=[]
		for qn in self.marker.qns:
			if qn+(self.L-self.l)>=self.M and qn-(self.L-self.l)<=self.M:
				qs.append(qn)
			self.marker.qns=qs
				
