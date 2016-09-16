#!/usr/bin/env python
"""
Tight-binding chain
e0:on-site energy
t:hopping t
N:chain length N.
"""

import numpy as np
import matplotlib.pyplot as plt

def band_energy(k,t=1.0,e0=0.2,a=1.0):
    """The function of energy with respect to k."""
    return e0-t*np.exp(1j*k*a)-t*np.exp(-1j*k*a)

def band_plot(N=400,a=1.0):
    """Plot the band in k-space."""
    foot_step=2*np.pi/N
    x=np.arange(0.0,2*np.pi/a,foot_step)
    y=band_energy(x)
    plt.plot(x,y) 

def density_of_state_plot(N=400,a=1.0,eita=0.01):
    """Plot the density_of_state respect to E."""
    foot_step=2*np.pi/N
    k=np.arange(0.0,2*np.pi/a,foot_step)
    Ek=band_energy(k)
    E=np.arange(-3.0,3.0,0.01)
    Ek.shape=(N,1)
    E.shape=(1,600)
    """Reshape E and Ek series with broadcasting method."""
    dirac_function=np.imag(np.true_divide(1/np.pi,np.subtract(E-Ek,1j*eita)))
    D=np.sum(np.true_divide(dirac_function,N),axis=0)
    """Calculate the density of state with lorentzian broadenning method.""" 
    E.shape=(600)
    plt.plot(D,E)

if __name__ == "__main__":
    band_plot()
    density_of_state_plot()
    plt.show()
