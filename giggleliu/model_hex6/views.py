from models import *
from scipy.sparse.linalg import eigsh
from scipy.linalg import eigvalsh
from matplotlib.pyplot import *
import pdb


def hexE0():
    '''lowest energy of hexagonal structure.'''
    U=1.
    t=1.
    model=Hex6(t=t,U=U,occ=True)
    rlattice=model.hgen.rlattice
    #The Hamiltonian is given by HGnerator.H(), it is a sparse matrix!
    t0=time.time()
    h=model.hgen.H()
    t1=time.time()
    Emin=eigsh(h,which='SA',k=1)[0]
    t2=time.time()
    print 'The Ground State Energy for hexagon(t = %s,U = %s) is %s.'%(t,U,Emin)
    print 'Time Elapse: get H -> %s, get Emin -> %s'%(t1-t0,t2-t1)
    pdb.set_trace()

def hexcompare():
    '''compare the result with exact result in the non-interaction limit.'''
    U=0.
    t=1.
    model_exact=Hex6(t=t,U=U,occ=False)
    model_occ=Hex6(t=t,U=U,occ=True)
    rlattice=model_occ.hgen.rlattice
    h_occ=model_occ.hgen.H()
    h_exact=model_exact.hgen.H()
    Emin=eigsh(h_occ,which='SA',k=1)[0]
    E_excit=eigvalsh(h_exact)
    Emin_exact=sum(E_excit[E_excit<0])
    print 'The Ground State Energy for hexagon(t = %s) is %s, tolerence %s.'%(t,Emin,Emin-Emin_exact)
    pdb.set_trace()
