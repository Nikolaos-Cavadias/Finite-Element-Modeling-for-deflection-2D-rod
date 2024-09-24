# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 16:38:40 2024

@author: Nikolas Cavadias
"""

import numpy as np

def Jacobi_Iteration(A, B, C, f, x0, epst):
    N = len(A)
    x = np.copy(x0)  # Initial guess for the solution
    x_new = np.zeros(N)
    
    iter_count = 0  # Initialize iteration counter
    
    L = np.diag(A) + np.diag(B,-1) + np.diag(C,1)
    
    residual_0 = np.linalg.norm(np.dot(L, x) - f)  # Initial residual norm
    residual_norm=residual_0
    
    while residual_norm > epst*residual_0:
        for i in range (N):
            if i == 0:
                x_new[i] = (f[i] - C[i] * x[i + 1]) / A[i]

            elif i == N-1:
                x_new[i] = (f[i] - B[i - 1] * x[i - 1]) / A[i]

            else:
                x_new[i] = (f[i] - B[i - 1] * x[i - 1] - C[i] * x[i + 1]) / A[i]

        x = x_new.copy()  # Update the solution vector
        residual_norm = np.linalg.norm(np.dot(L, x) - f)  # Calculate the new residual norm
        iter_count += 1  # Increment iteration counter

    residual = np.dot(L, x) - f  # Calculate the final residual vector
    rnk = np.linalg.norm(residual)  # Calculate the norm of the final residual
    
    return x, rnk, iter_count

