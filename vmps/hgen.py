import numpy as np
from scipy.sparse import kron,identity
from copy import deepcopy

class Hgen(object):
	def __init__(self,L,terms=[]):
		self.L=L
		self.terms=terms
