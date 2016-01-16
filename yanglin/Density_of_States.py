#!/usr/bin/env python
"""Density of State"""
from __future__ import division
import numpy as np
import matplotlib.pyplot as plt

def density_of_states_plot(k,Ek,N,eita=0.01):
    length=len(Ek)
    E=np.sort(Ek)
    E.shape=(1,length)
    Ek.shape=(length,1)
    dirac_function=np.imag(np.true_divide(1/np.pi,np.subtract(E-Ek,1j*eita)))
    D=np.sum(np.true_divide(dirac_function,N),axis=0)
    E.shape=(length)
    plt.plot(D,E)
