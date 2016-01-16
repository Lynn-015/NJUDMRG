#!/usr/bin/env python
"""The Kane-Mele Model"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt
import Density_of_States as dos

def band_energy(k,t1=0.62,t2=0,t3=-0.38,a=1/np.sqrt(3)):
    '''
    Get the band energy
    t1:nearest hopping
    t2:second nearest hopping
    t3:third nearest hopping
    '''
    a1=np.array([3/2*a,-np.sqrt(3)/2*a])
    a2=np.array([3/2*a,np.sqrt(3)/2*a])
    '''a1,a2 are base vectors of the lattice'''
    b1=np.array([2*np.pi/(3*a),-2*np.pi/(np.sqrt(3)*a)])
    b2=np.array([2*np.pi/(3*a),2*np.pi/(np.sqrt(3)*a)])
    '''b1,b2 are base vectors of the reciprocal lattice'''
    
    delta1_1=2/3*a1-1/3*a2
    delta1_2=-1/3*a1+2/3*a2
    delta1_3=-1/3*a1-1/3*a2
    delta2_1=a1
    delta2_2=a2
    delta2_3=-a1+a2
    delta3_1=2/3*a1+2/3*a2
    delta3_2=-4/3*a1+2/3*a2
    delta3_3=2/3*a1-4/3*a2
    
    T_1=-t1*(np.exp(1j*np.dot(delta1_1,k))\
            +np.exp(1j*np.dot(delta1_2,k))\
            +np.exp(1j*np.dot(delta1_3,k)))
    T_2=-t2*(np.exp(1j*np.dot(delta2_1,k))\
            +np.exp(1j*np.dot(delta2_2,k))\
            +np.exp(1j*np.dot(delta2_3,k)))
    T_3=-t3*(np.exp(1j*np.dot(delta3_1,k))\
            +np.exp(1j*np.dot(delta3_2,k))\
            +np.exp(1j*np.dot(delta3_3,k)))
    
    Ek_plus=T_2+np.conjugate(T_2)+abs(T_1+T_3)
    Ek_minus=T_2+np.conjugate(T_2)-abs(np.conjugate(T_1+T_3))
    
    return {'Hamiltonian':[[T_2+np.conjugate(T_2),T_1+T_3],
            [np.conjugate(T_1+T_3),T_2+np.conjugate(T_2)]],
            'High':Ek_plus ,'Low':Ek_minus}   
               
def band_plot(foot_step=0.1,b=4*np.pi/np.sqrt(3)):
    '''Plot the band energy'''
    k_G_to_K=np.arange(0,b/np.sqrt(3),foot_step)
    k1_G_to_K=k_G_to_K*np.sqrt(3)/2
    k2_G_to_K=k_G_to_K/2

    k_K_to_M=np.arange(b/np.sqrt(3),b*np.sqrt(3)/2,foot_step)
    k1_K_to_M=0*k_K_to_M+2*np.pi/np.sqrt(3)
    k2_K_to_M=2*np.pi/3-(k_K_to_M-4*np.pi/3)

    k_M_to_G=np.arange(b*np.sqrt(3)/2,b*(np.sqrt(3)+1)/2,foot_step)
    k1_M_to_G=2*np.pi/np.sqrt(3)-(k_M_to_G-2*np.pi)
    k2_M_to_G=0*k_M_to_G

    k=np.append(np.append(k_G_to_K,k_K_to_M),k_M_to_G)
    k1=np.append(np.append(k1_G_to_K,k1_K_to_M),k1_M_to_G)
    k2=np.append(np.append(k2_G_to_K,k2_K_to_M),k2_M_to_G)
    
    Ek_high=band_energy([k1,k2])['High']
    Ek_low=band_energy([k1,k2])['Low']
    plt.plot(k,Ek_high)
    plt.plot(k,Ek_low)

    dos.density_of_states_plot(k,Ek_low,N=100)
    dos.density_of_states_plot(k,Ek_high,N=100)
    '''Plot the density of states'''

    plt.plot(4*np.pi/3,band_energy([2*np.pi/np.sqrt(3),2*np.pi/3])['Low'],'o')
    plt.annotate('   K',xy=(4*np.pi/3,band_energy([2*np.pi/np.sqrt(3),2*np.pi/3])['Low']))
    plt.plot(0,band_energy([0,0])['Low'],'o')
    plt.annotate('   G',xy=(0,band_energy([0,0])['Low']))
    plt.plot(2*np.pi,band_energy([2*np.pi/np.sqrt(3),0])['Low'],'o')
    plt.annotate('   M',xy=(2*np.pi,band_energy([2*np.pi/np.sqrt(3),0])['Low']))
    plt.plot(2*np.pi*(1+1/np.sqrt(3)),band_energy([0,0])['Low'],'o')
    plt.annotate('   G',xy=(2*np.pi*(1+1/np.sqrt(3)),band_energy([0,0])['Low']))
    
    print 'Chemical Potential=',np.real(np.sort(Ek_low)[1/2*len(k)])
    return np.real(np.sort(Ek_low)[1/2*len(k)])

def spectral_function(H,mu,eita=0.01):
    '''Get the spectral function'''
    return -2*np.imag(np.linalg.inv((1j*eita+mu)*np.identity(2)-H))

def fermi_surface_plot(mu,foot_step=0.1):
    '''Plot the fermi surface'''
    kx=np.arange(-6,6,foot_step)
    ky=np.arange(-6,6,foot_step)
    x,y=np.meshgrid(kx,ky)
    s=x*y
    for i in range(len(kx)):
        for j in range(len(ky)):
            s[i][j]=np.trace(spectral_function(band_energy([kx[i],ky[j]])['Hamiltonian'],mu))         
    plt.pcolor(x,y,s) 

if __name__ == "__main__":
    mu=band_plot()
    plt.show()
    fermi_surface_plot(mu)
    plt.show()
