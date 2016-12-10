import numpy as np
from copy import deepcopy

from marker import Marker,MarkerGen,sz_in,sz_out

class TestMarker(object):
	def __init__(self):
		mgen=MarkerGen(6)
		markers=[]
		for i in range(6):
			mgen.enlarge([sz_in])
			mgen.filt()
			marker=deepcopy(mgen.marker)
			marker.reorder()
			marker.compact()
			markers.append(marker)
			print marker.qns,marker.sizes,marker.divs

if __name__=='__main__':
	TestMarker()
