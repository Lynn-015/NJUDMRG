#!/usr/bin/python
'''
Utilities.
'''
from numpy import *

def bcast_dot(A,B):
    '''
    dot product broadcast version.
    '''
    return einsum('...jk,...kl->...jl', A, B)

def formatstr(s,fmt):
    '''
    Get colored string for terminal.

    Parameters
    ------------------
    s:
        The target string.
    fmt:
        The format.

    Return
    --------------
    colored string.
    '''
    ENDC = '\033[0m'
    fmtdict=dict(
    HEADER = '\033[95m',
    OKBLUE = '\033[94m',
    OKGREEN = '\033[92m',
    WARNING = '\033[93m',
    FAIL = '\033[91m',
    BOLD = '\033[1m',
    UNDERLINE = '\033[4m',
    )
    return fmtdict.get(fmt)+str+ENDC
