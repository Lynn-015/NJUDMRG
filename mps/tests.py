#!/usr/bin/python
'''
Tests for MPS and MPO
'''
from mpslib import *

def test_mps(hndim):
    '''
    Test function for constructing matrix product

    hndim:
        The number of states on each site.
    '''
    nsite=10   #number of sites
    vec=random.random(hndim**nsite)  #a random state in form of 1D array.
    t0=time.time()
    mps=state2MPS(vec,sitedim=hndim,l=nsite/2,method='svd')     #parse the state into a <MPS> instance.
    t1=time.time()
    print 'Get MPS: %s, Elapse -> %s'%(mps,t1-t0)

    t0=time.time()
    qu=mps.query(zeros(nsite))   #query the state of a specified site config.
    t1=time.time()
    print 'Query the first site: %s(true: %s), Elapse -> %s'%(qu,vec[0],t1-t0)

    t0=time.time()
    nstate=mps.state            #recover the 1D array state representation.
    t1=time.time()
    print 'State tolerence %s, Elapse -> %s'%(sum(abs(nstate-vec)),t1-t0)

    #show it graphically.
    ion()
    mps.show()

    pdb.set_trace()

def test_vmps(hndim):
    '''
    Test function for constructing matrix product

    hndim:
        The number of states on each site.
    '''
    nsite=10   #number of sites
    vec=random.random(hndim**nsite)/sqrt(hndim**nsite/2.)  #a random state in form of 1D array.
    t0=time.time()
    vmps=state2VMPS(vec,sitedim=hndim)     #parse the state into a <MPS> instance.

    t0=time.time()
    nstate=vmps.state            #recover the 1D array state representation.
    t1=time.time()
    print 'State tolerence %s, Elapse -> %s'%(sum(abs(nstate-vec)),t1-t0)

    vmps.check_canonical()
    ion()
    vmps.show()
    print '\nChanging to canonical form!'
    print '###########################'
    mps=vmps.canonical(6)
    t0=time.time()
    nstate=mps.state            #recover the 1D array state representation.
    t1=time.time()
    print 'State tolerence %s, Elapse -> %s'%(sum(abs(nstate-vec)),t1-t0)
    t1=time.time()
    print 'Get MPS: %s, Elapse -> %s'%(mps,t1-t0)
    #mps.show()
    pdb.set_trace()


def test_mps(hndim):
    '''
    Test function for constructing matrix product

    hndim:
        The number of states on each site.
    '''
    nsite=10   #number of sites
    l=3
    vec=random.random(hndim**nsite)  #a random state in form of 1D array.
    vec2=random.random(hndim**nsite)  #a random state in form of 1D array.
    t0=time.time()
    mps=state2MPS(vec,sitedim=hndim,l=l,method='svd')     #parse the state into a <MPS> instance.
    mps2=state2MPS(vec2,sitedim=hndim,l=l,method='svd')     #parse the state into a <MPS> instance.
    t1=time.time()
    print 'Get MPS: %s, Elapse -> %s'%(mps,t1-t0)

    t0=time.time()
    qu=mps.query(zeros(nsite))   #query the state of a specified site config.
    t1=time.time()
    print 'Query the first site: %s(true: %s), Elapse -> %s'%(qu,vec[0],t1-t0)

    t0=time.time()
    nstate=mps.state            #recover the 1D array state representation.
    t1=time.time()
    print 'State tolerence %s, Elapse -> %s'%(sum(abs(nstate-vec)),t1-t0)

    print 'Checking for unitary for M-matrices.\n',mps.check_canonical()

    #test for multiplication
    time.time()
    mpsbra=mps.hconj
    overlap=mpsbra*mps
    t1=time.time()
    print 'Overlap %s(true %s), Elapse -> %s.'%(overlap,(vec.conj()*vec).sum(),t1-t0)

    #show it graphically.
    ion()
    mps.show()
    mpsbra.show(offset=(0,2))

    pdb.set_trace()

def test_compress(hndim):
    '''
    Test for addition of two <MPS> instances.
    '''
    nsite=10   #number of sites
    l=3
    vec=random.random(hndim**nsite)  #a random state in form of 1D array.
    vec2=random.random(hndim**nsite)  #a random state in form of 1D array.
    #vec2=vec
    vecadded=vec+vec2

    mps=state2MPS(vec,sitedim=hndim,l=l,method='svd')     #parse the state into a <MPS> instance.
    mps2=state2MPS(vec2,sitedim=hndim,l=l,method='svd')     #parse the state into a <MPS> instance.
    mpsadded=mps_add(mps,mps2)

    print mpsadded
    t0=time.time()
    nstate=mps.state            #recover the 1D array state representation.
    nstateadded=mpsadded.state            #recover the 1D array state representation.
    t1=time.time()
    print 'State tolerence %s Elapse -> %s'%(sum(abs(nstateadded-vecadded)),t1-t0)

    t0=time.time()
    mpsadded.compress()
    t1=time.time()
    nstateadded=mpsadded.state            #recover the 1D array state representation.
    print mpsadded
    print 'State tolerence(after compress) %s Elapse -> %s'%(sum(abs(nstateadded-vecadded)),t1-t0)
    pdb.set_trace()

def test_add(hndim):
    '''
    Test for addition of two <MPS> instances.
    '''
    nsite=10   #number of sites
    l=3
    vec=random.random(hndim**nsite)   #a random state in form of 1D array.
    vec2=random.random(hndim**nsite)  #a random state in form of 1D array.
    vecadded=vec+vec2

    mps=state2MPS(vec,sitedim=hndim,l=l,method='svd')     #parse the state into a <MPS> instance.
    mps2=state2MPS(vec2,sitedim=hndim,l=l,method='svd')     #parse the state into a <MPS> instance.
    mpsadded=mps_add(mps,mps2)
    print mps
    print mpsadded

    t0=time.time()
    nstate=mps.state            #recover the 1D array state representation.
    nstate2=mps2.state            #recover the 1D array state representation.
    nstateadded=mpsadded.state            #recover the 1D array state representation.
    t1=time.time()
    print 'State tolerence(first,second,added) %s, %s, %s Elapse -> %s'%(sum(abs(nstate-vec)),sum(abs(nstate2-vec2)),sum(abs(nstateadded-vecadded)),t1-t0)

    pdb.set_trace()

def test_move(hndim):
    '''
    Test for canonical move of <MPS>.
    '''
    ion()
    nsite=10   #number of sites
    l=3
    vec=random.random(hndim**nsite)  #a random state in form of 1D array.
    mps=state2MPS(vec,sitedim=hndim,l=l,method='svd')     #parse the state into a <MPS> instance.
    t0=time.time()
    nstate=mps.state            #recover the 1D array state representation.
    t1=time.time()
    mps.show()
    print 'State tolerence %s, Elapse -> %s'%(sum(abs(nstate-vec)),t1-t0)
    pdb.set_trace()

    t0=time.time()
    mps>>(1,1e-5)
    t1=time.time()
    nstate=mps.state            #recover the 1D array state representation.
    cla()
    mps.show()
    print 'State tolerence %s, Elapse -> %s'%(sum(abs(nstate-vec)),t1-t0)
    pdb.set_trace()
 
    t0=time.time()
    mps<<(4,1e-6)
    t1=time.time()
    nstate=mps.state            #recover the 1D array state representation.
    cla()
    mps.show()
    print 'State tolerence %s, Elapse -> %s'%(sum(abs(nstate-vec)),t1-t0)
    pdb.set_trace()

if __name__=='__main__':
    test_compress(3)
