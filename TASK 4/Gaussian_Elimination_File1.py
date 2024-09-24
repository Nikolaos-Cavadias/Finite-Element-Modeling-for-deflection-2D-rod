# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 16:36:55 2024

@author: Nikolas Cavadias
"""
import numpy as np
def Gaussian_Elimination(A,B,C,f):
    
    n=len(A)
    x=np.zeros(n)
    #forward elimination
    for i in range (1,n,1):
        temp = B[i-1]/A[i-1]
        A[i]=A[i]-temp*C[i-1]
        f[i]=f[i]-temp*f[i-1]
            
    x[n-1]=f[n-1]/A[n-1]
        
    #backward elimination
    for i in range (n-2,-1,-1):
        x[i]=(f[i]-C[i]*x[i+1])/A[i]

    return x
