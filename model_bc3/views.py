'''
Viewes for this project, jobs can be defined here.
'''
from numpy import *
from models import *
from tba.lattice.path import path_k
from tba.hgen.mesh import Hmesh,Emesh,Gmesh
from tba.lattice.fs import FSHandler
from matplotlib.pyplot import *
import pdb

def bc3lattice():
    '''show lattice structures of bc3 model'''
    model=BC3((20,30))  #lattice size (20 x 30)
    ion()
    rlattice=model.hgen.rlattice  #the lattice of this model
    rlattice.show_sites(color='r')    #show sites in red
    rlattice.initbonds()
    print len(rlattice.getbonds(1))
    rlattice.show_bonds((1,),color='k')  #show nearest neighbor bonds(indicated by (1,)) in black.
    axis('equal')
    pdb.set_trace()
    cla()
    rlattice.groups.get('translation').per=(True,False)
    rlattice.initbonds()
    print len(rlattice.getbonds(1))
    rlattice.show_sites(color='r')    #show sites in red
    rlattice.show_bonds((1,),color='k')  #show nearest neighbor bonds(indicated by (1,)) in black.
    axis('equal')
    pdb.set_trace()

def bc3band():
    '''plot the band.'''
    #first, create a <BC3> instance, using the default lattice size (100 x 100).
    model=BC3()
    hgen=model.hgen   #the hamiltonian generator of this model.
    kspace=hgen.kspace #the <KSpace> used by this hamiltonian generator.

    #then generate a path(KPath instance), vertices and number of k-samples as parameters.
    #path_k will generate a route along given verices.
    path=path_k(vertices=[kspace.G,kspace.K[0],kspace.M[0],kspace.G],N=300)

    #evaluate energy(through hgen.Ek function) for k defined on this path.
    #a list of energy is returned by path.eval
    ekl=path.eval(hgen.Ek)

    #finally, display the energies using path.plot.
    ion()
    path.plot(ekl)
    pdb.set_trace()

def bc3disper():
    '''plot the dispresion of BC3'''
    model=BC3((70,80))
    hgen=model.hgen

    #hgen.gethkmesh will generator a mesh of hamiltonian on the kspace.kmesh,
    #alternatively, you may specify the target kmesh by calling hgen.gethkmesh(kmesh)
    hkmesh=hgen.gethkmesh()

    #Create a <Hmesh> instance and call the method Hmesh.getemesh to genrator a mesh of energies.
    #multi-processing is supported in this method.
    ekmesh=Hmesh(hkmesh).getemesh()
    ion()
    kspace=hgen.kspace

    #this time, we will get the Brillouing zone of this kspace, which is identified by the vertices in k-space.
    bzone=kspace.get_bzone()
    kmesh=kspace.kmesh

    ############## Homework 3, uncomment and complete the following statement #################
    #kmesh=bzone.k2bzone(kspace.kmesh,b=kspace.b)
    Emesh(ekmesh).show(kmesh)   #Create an Emesh instance, note this Emesh supports method of surf the dispersion(show).
    pdb.set_trace()

def bc3fs():
    '''
    Show the fermi surface of BC3.
    '''
    ion()
    filling=0.25
    geta=2e-2
    model=BC3((100,100))
    hgen=model.hgen
    kmesh=hgen.kspace.kmesh
    hkmesh=hgen.gethkmesh()
    ek=sort(Hmesh(hkmesh).getemesh().ravel())
    #get the chemical potential
    mu=ek[int(filling*len(ek))]-hgen.params['-mu']
    hgen.params['-mu']=-mu

    #generate a new set of hamiltonian
    hkmesh=hgen.gethkmesh()
    #get the retarded Green's function with smearing factor geta=0.02.
    gkmesh=Hmesh(hkmesh).getgmesh(w=0,geta=geta,tp='r')
    #get the spectrum function of this Green's function, trace over it to get electron density.
    Amesh=trace(Gmesh(gkmesh,geta=geta,tp='r').Amesh,axis1=-1,axis2=-2).real
    #show the spectrum function.
    pcolormesh(kmesh[...,0],kmesh[...,1],Amesh)

    #search for fermi surface - advance version
    for i in [0]:
        print 'Searching for %s-th Band!'%i
        #get the function of energy for the i-th band
        efunc=lambda k:hgen.Ek(k)[i]
        #create a handler instance for fermi surface. with resolution ability 2% and tolerence 1e-5
        fsh=FSHandler(efunc,resolution=0.02,tol=1e-5)
        #get the surface by name, 'GKM' for `G` pocket, 'M' pockets, and 'K' pockets will be searched.
        fs=fsh.get_fs_byname('GKM',hgen.kspace,nseg=50)
        #display the fermi surface.
        fs.show(color='r')

    #show the Brillouin zone.
    hgen.kspace.get_bzone().show(edgecolor='r')
    pdb.set_trace()


