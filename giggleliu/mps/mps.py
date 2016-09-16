#!/usr/bin/python
'''
Matrix Product State.
'''
from numpy import *
from matplotlib.pyplot import *
from matplotlib import patches
from matplotlib.collections import LineCollection
from scipy.linalg import svd,qr,rq
from scipy import sparse as sps
from utils import bcast_dot
import pdb,time

LEFT_CANONICAL='-'
RIGHT_CANONICAL='|'
NO_CANONICAL='x'

class MPSBase(object):
    '''
    The Base class of Matrix Product state.

    Attributes
    -----------
    hndim:
        The number of channels on each site.
    nsite:
        The number of sites.
    '''
    @property
    def hndim(self):
        '''The number of state in a single site.'''
        raise Exception('Not Implemented!')

    @property
    def nsite(self):
        '''Number of sites.'''
        raise Exception('Not Implemented!')

    @property
    def state(self):
        '''The state representation of this MPS'''
        raise Exception('Not Implemented!')

    def show(self,*args,**kwargs):
        '''Show this MPS graphically'''
        raise Exception('Not Implemented!')

class VidalMPS(MPSBase):
    '''
    Matrix Product state in the standard Vidal form.

    Construct
    -----------
        VidalMPS(GL,LL)

    Attributes
    -----------
    GL:
        A list of Gamma Matrices.
    LL:
        A list of Lambda Matrices.
    S:
        The remaining factor for Vidal representation.
    '''
    def __init__(self,GL,LL,S):
        assert(len(GL)==len(LL)+1)
        self.LL=LL
        self.GL=GL
        self.S=S

    @property
    def hndim(self):
        '''The number of state in a single site.'''
        return len(self.GL[0])

    @property
    def nsite(self):
        '''Number of sites.'''
        return len(self.GL)

    @property
    def state(self):
        '''The state representation of this MPS'''
        ULG=[swapaxes(ai,0,1).reshape([-1,ai.shape[-1]]) for ai in self.GL]
        LL=[identity(1)*self.S]+self.LL
        state=identity(1)
        for ui,li in zip(ULG[::-1],LL[::-1]):
            state=reshape(state,(ui.shape[1],-1))
            state=ui.dot(state)
            state=state.reshape([li.shape[0],-1])*li[:,newaxis]
        return state.ravel()

    def show(self,offset=(0,0)):
        '''
        Display this MPS instance graphically.

        offset:
            The global displace of the MPS chain.
        '''
        r=0.25
        rs=0.35
        barend=(r+rs)*1.3
        textoffset=barend+0.1
        color_A='k'
        color_S=color_A
        color_B=color_A
        edgecolor='k'
        facecolor='none'
        ax=gca()
        lc=[]
        for i in xrange(self.nsite):
            xi,yi=2*i+offset[0],offset[1]
            ax.add_patch(patches.Circle(xy=(xi,yi),radius=r,facecolor=facecolor,edgecolor=edgecolor,hatch='xxxx'))
            lc.append([(xi,yi+r),(xi,yi+barend)])
            text(xi,yi+textoffset,r'$\Gamma^{\sigma_%s}$'%i,horizontalalignment='center')
            if i!=self.nsite-1:
                lc.append([(xi+r,yi),(xi+1-rs,yi)])
                ax.add_patch(patches.Polygon(xy=[(xi+1,yi+rs),(xi+1+rs,yi),(xi+1,yi-rs),(xi+1-rs,yi)],facecolor=facecolor,hatch='xxxx'))
                lc.append([(xi+1,yi+rs),(xi+1,yi+barend)])
                lc.append([(xi+1+rs,yi),(xi+2-r,yi)])
                text(xi+1,yi+textoffset,r'$\Lambda^{[%s]}$'%i,horizontalalignment='center')
        lc=LineCollection(lc,color=edgecolor,lw=2)
        ax.add_collection(lc)
        axis('equal')
        ax.autoscale()

    def get_rho(self,l=None):
        '''
        Density matrix between site l and l+1.
        
        l:
            The site index.
        '''
        if l is None:
            rl=[diag(li**2) for li in self.LL]
            return rl
        else:
            rl=diag(self.LL[l]**2)
            return rl

    def von_Neumann_entropy(self,l=None):
        '''
        The von-Neumann entropy at the link between A|B split at l.

        l:
            The index of division point.
        '''
        print 'Warning , not tested.'
        if l is None:
            return [-sum((li*log(li)/log(2.))**2) for li in self.LL]
        else:
            return -sum((li*log(self.LL[l])/log(2.))**2)

    def check_canonical(self,tol=1e-5):
        '''
        Check the canonicallity of this VidalMPS instance.

        tol:
            The rolerence of the result from orthogonality.
        '''
        rl=self.get_rho()
        nsite=self.nsite
        #check for left canonical
        i_l=[sum([gij.T.conj().dot(ri).dot(gij) for gij in gi],axis=0) for ri,gi in zip([1]+rl,self.GL)]
        diff_l=array([sum(abs(ii-identity(ii.shape[-1]))) for ii in i_l])
        #check for right canonical
        i_r=[sum([gij.dot(ri).dot(gij.T.conj())*(1. if i!=0 else 1./self.S**2) for gij in gi],axis=0) for i,ri,gi in zip(arange(nsite),rl+[self.S**2],self.GL)]
        diff_r=array([sum(abs(ii-identity(ii.shape[-1]))) for ii in i_r])
        cl=diff_l<tol
        cr=diff_r<tol
        print 'Checking canonicallity for left sweep ->',cl
        print 'Checking canonicallity for right sweep ->',cr
        return all(cl) and all(cr)

    def canonical(self,l='L'):
        '''
        Get the canonical form for this MPS.

        Parameters
        -----------
        l:
            The specific canonical form.

            * `L` -> left canonical form.
            * `R` -> right canonical form.
            * <integer> -> the mixed canonical form with division index l(0 to self.nsite).
        '''
        nsite=self.nsite
        hndim=self.hndim
        if l=='L':
            sindex=nsite
        elif l=='R':
            sindex=0
        else:
            sindex=l
        assert(0<=sindex and sindex<=nsite)
        GL=self.GL
        LL=[ones(1)]+self.LL+[ones(1)]
        AL=[LL[i][:,newaxis]*GL[i].reshape([hndim,LL[i].shape[0],-1]) for i in xrange(sindex)]
        BL=[GL[i].reshape([hndim,-1,LL[i+1].shape[0]])*LL[i+1] for i in xrange(sindex,nsite)]
        S=LL[sindex]*self.S
        return KMPS(AL,BL,S,sindex)

class MPS(MPSBase):
    '''
    Matrix product states.

    Attributes
    -----------------
    AL/BL:
        The sequence of A/B-matrices.
    l/S:
        The division point of left and right scan, and the singular value matrix at the division point.

        Also, it is the index of the non-unitary M-matrix(for the special case of l==N, S is attached to the right side of N-1-th matrix)
    is_ket:
        It is a ket if True else bra.
    token:
        The unique token to identify this MPS.
    '''
    def __init__(self,AL,BL,S,is_ket,token):
        assert(ndim(S)==1)
        self.AL=AL
        self.BL=BL
        self.S=S
        #constructing ML
        N=self.nsite
        self.is_ket=is_ket
        self.token=token
        self.canonical_mask=[LEFT_CANONICAL]*len(AL)+[RIGHT_CANONICAL]*len(BL)

    def __str__(self):
        string='<MPS,%s>\n'%(self.nsite)
        string+='\n'.join(['  A[s=%s] (%s x %s)'%(len(a),a[0].shape[0],a[0].shape[1]) for a in self.AL])+'\n'
        string+='  S      %s\n'%(self.S.shape,)
        string+='\n'.join(['  B[s=%s] (%s x %s)'%(len(a),a[0].shape[0],a[0].shape[1]) for a in self.BL])
        return string

    def __mul__(self,target):
        hndim=self.hndim
        if self.is_ket or not target.is_ket:
            raise Exception('Not implemented for multipling ket on the left side.')
        else:
            S=identity(1)
            for mi,tmi in zip(self.ML,target.ML):
                S=sum([mi[j].dot(S).dot(tmi[j]) for j in xrange(hndim)],axis=0)
            return S

    def __rmul__(self,target):
        return target.__mul__(self)

    def __lshift__(self,k):
        '''Left move l-index by k.'''
        tol=1e-8
        if isinstance(k,tuple):
            k,tol=k
        for i in xrange(k):
            self.__canonical_move__(False,tol=tol)

    def __rshift__(self,k):
        '''Right move l-index by k.'''
        tol=1e-8
        if isinstance(k,tuple):
            k,tol=k
        for i in xrange(k):
            self.__canonical_move__(True,tol=tol)

    def __canonical_move__(self,right,tol):
        '''
        Move l-index by one with specific direction.
        
        Parameters
        --------------
        right:
            Move l to right if True.
        tol:
            The tolerence for compression.
        '''
        nsite=self.nsite
        l=self.l
        hndim=self.hndim
        if (l>=nsite and right) or (l<=0 and not right):
            raise ValueError('Can not move to %s for l = %s with %s sites in total.'%('right' if right else 'left',l,nsite))
        if right:
            B0=self.BL.pop(0)
            B0=swapaxes(self.S[:,newaxis]*B0,0,1)
            M=reshape(B0,(-1,B0.shape[-1]))
            U,S,V=svd(M,full_matrices=False)
            kpmask=abs(S)>tol
            U,S,V=U[:,kpmask],S[kpmask],V[kpmask]
            self.AL.append(swapaxes(reshape(U,(-1,hndim,U.shape[1])),0,1))
            self.S=S
            if len(self.BL)>0:
                self.BL[0]=bcast_dot(V,self.BL[0])
            else:
                self.S*=V[0,0]
        else:
            A0=self.AL.pop(-1)
            A0=swapaxes(A0*self.S,0,1)
            M=reshape(A0,(A0.shape[0],-1))
            U,S,V=svd(M,full_matrices=False)
            kpmask=abs(S)>tol
            U,S,V=U[:,kpmask],S[kpmask],V[kpmask]
            self.BL.insert(0,swapaxes(reshape(V,(V.shape[0],hndim,-1)),0,1))
            self.S=S
            if len(self.AL)>0:
                self.AL[-1]=bcast_dot(self.AL[-1],U)
            else:
                self.S*=U[0,0]

    @property
    def hndim(self):
        '''The number of state in a single site.'''
        if len(self.AL)>0:
            return len(self.AL[0])
        else:
            return len(self.BL[0])

    @property
    def nsite(self):
        '''Number of sites.'''
        return len(self.AL)+len(self.BL)

    @property
    def l(self):
        '''The division point of left/right canonical forms.'''
        return len(self.AL)

    @property
    def ML(self):
        '''A|B list in a string.'''
        ML=self.AL+self.BL
        l=self.l
        hndim=self.hndim
        is_sparse=sps.issparse(ML[0][0])
        if l!=0:
            if self.is_ket:
                if is_sparse:
                    ML[l-1]=[mi.multiply(self.S) for mi in ML[l-1]]
                else:
                    ML[l-1]=ML[l-1]*self.S
            else:
                if is_sparse:
                    ML[l-1]=[mi.multiply(self.S.T) for mi in ML[l-1]]
                else:
                    ML[l-1]=ML[l-1]*self.S[:,newaxis]
        else:
            if self.is_ket:
                if is_sparse:
                    ML[l]=[mi.multiply(self.S.T) for mi in ML[l]]
                else:
                    ML[l]=self.S[:,newaxis]*ML[l]
            else:
                if is_sparse:
                    ML[l]=[mi.multiply(self.S) for mi in ML[l]]
                else:
                    ML[l]=self.S*ML[l]
        return ML

    @property
    def UL(self):
        '''U list in a string.'''
        if self.is_ket:
            ULA=[swapaxes(ai,0,1).reshape([-1,ai.shape[-1]]) for ai in self.AL]
            ULB=[swapaxes(bi,0,1).reshape([bi.shape[0],-1]) for bi in self.BL]
        else:
            ULA=[swapaxes(ai,0,1).reshape([ai.shape[0],-1]) for ai in self.AL]
            ULB=[swapaxes(bi,0,1).reshape([-1,bi.shape[-1]]) for bi in self.BL]
        UL=ULA+ULB
        l=self.l
        if l!=0:
            UL[l-1]=UL[l-1]*self.S
        else:
            UL[l]=self.S[:,newaxis]*UL[l]
        return UL

    @property
    def state(self):
        '''The state representation of this MPS'''
        hndim=self.hndim
        S=identity(1)
        t0=time.time()
        for mi in self.ML[::-1]:
            S=array([mi[j].dot(S) for j in xrange(hndim)])
            S=swapaxes(S,0,1).reshape([S.shape[1],-1])
            #S=swapaxes(S,0,1).reshape([S.shape[1],-1])
        return S.ravel()

    @property
    def hconj(self):
        '''
        Get the hermitian conjugate of this ket or bra.
        '''
        if self.is_ket:
            return BMPS([swapaxes(ai,2,1).conj() for ai in self.AL],[swapaxes(bi,2,1).conj() for bi in self.BL],self.S.conj(),token=self.token)
        else:
            return KMPS([swapaxes(ai,2,1).conj() for ai in self.AL],[swapaxes(bi,2,1).conj() for bi in self.BL],self.S.conj(),token=self.token)

    def check_canonical(self,tol=1e-5):
        '''
        Check the unitaryness of the M matrix in a MPS.
        '''
        res=[]
        hndim=self.hndim
        for i in xrange(self.l):
            mi=self.ML[i]
            res.append(all(abs(sum([mi[j].T.conj().dot(mi[j]) for j in xrange(hndim)],axis=0)-identity(mi.shape[-1]))<tol))
        for i in xrange(self.l,self.nsite):
            mi=self.ML[i]
            res.append(all(abs(sum([mi[j].dot(mi[j].T.conj()) for j in xrange(hndim)],axis=0)-identity(mi.shape[1]))<tol))
        return res


    def query(self,serie):
        '''Query the magnitude of a state.'''
        state=identity(1)
        for si,Mi in zip(serie,self.ML)[::-1]:
            state=Mi[si].dot(state)
        return state[0,0]

    def show(self,offset=(0,0)):
        '''
        Display this MPS instance graphically.

        offset:
            The global displace of the MPS chain.
        '''
        r=0.25
        rs=0.35
        rlc=0.5 if self.is_ket else -0.5
        contourlw=1
        color_A='none'
        color_S=color_A
        color_B=color_A
        edgecolor='#000000'
        baroffset=r if self.is_ket else -r
        barlen=r if self.is_ket else -r
        textoffset=baroffset+barlen+(0.2 if self.is_ket else -0.2)
        ax=gca()
        lc=[]
        for i in xrange(self.l):
            xi,yi=i+offset[0],0+offset[1]
            ax.add_patch(patches.Circle(xy=(xi,yi),radius=r,facecolor=color_A,lw=contourlw,edgecolor=edgecolor,hatch='---'))
            lc.append([(xi,yi+baroffset),(xi,yi+baroffset+barlen)])
            lc.append([(xi+r,yi),(xi+1-(rs if i==self.l-1 else r),yi)])
            text(xi,yi+textoffset,r'$\sigma_%s$'%i,horizontalalignment='center',verticalalignment='center')
        xS,yS=self.l+offset[0],offset[1]
        ax.add_patch(patches.Polygon(xy=[(xS,yS+rs),(xS+rs,yS),(xS,yS-rs),(xS-rs,yS)],lw=contourlw,facecolor=color_S))
        text(xS,yS+textoffset,r'$S$',horizontalalignment='center')
        for i in xrange(self.l+1,self.nsite+1):
            xi,yi=i+offset[0],0+offset[1]
            ax.add_patch(patches.Circle(xy=(xi,yi),radius=r,lw=contourlw,facecolor=color_B,edgecolor=edgecolor,hatch='|||'))
            lc.append([(xi,yi+baroffset),(xi,yi+baroffset+barlen)])
            lc.append([(xi-1+(rs if i==self.l+1 else r),yi),(xi-r,yi)])
            text(xi,yi+textoffset,r'$\sigma_%s$'%(i-1),horizontalalignment='center',verticalalignment='center')
        lc=LineCollection(lc,color=edgecolor,lw=2)
        ax.add_collection(lc)
        ax.autoscale()
        axis('equal')

    def compress(self,tol=1e-5):
        '''
        Compress this state.
        '''
        nsite,l=self.nsite,self.l
        self>>(nsite-l,tol)
        self<<(nsite,tol)
        self>>(l,tol)

class KMPS(MPS):
    '''
    Ket MPS.
    '''
    def __init__(self,AL,BL,S,token=None):
        if token is None: token=id(self)
        super(KMPS,self).__init__(AL,BL,S,is_ket=True,token=token)

class BMPS(MPS):
    '''
    Bra MPS.
    '''
    def __init__(self,AL,BL,S,token=None):
        if token is None: token=id(self)
        super(BMPS,self).__init__(AL,BL,S,is_ket=False,token=token)
