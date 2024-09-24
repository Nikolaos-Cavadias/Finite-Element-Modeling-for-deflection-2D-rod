# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 13:25:19 2024

@author: Nikolas Cavadias
"""

"""
N = no of nodes = Ne+1
e(i)= j = end points of each elements 1,Ne

i= node index = 0, 1      
h= element size x(j+1)-x(i)
quad_type= 1 or 2
lower bound =x(j)
Upper bound =x(j+1)

a[nxn] matrix

"""
from sympy import *


def Numerical_Quadrature(e,i,quad_type,h, integrantfunc):
    xx = symbols('x')
    def f(j):
        return integrantfunc.subs(xx, j)
    
    if i==0 :
        φb=0
        φa=1
    elif i==1:
        φb=1
        φa=0
    else:
        print("Error")
   
    def Trapezoidal_function(a, b, φa, φb):
        I=h/2*(f(a)*φa + f(b)*φb)
        return I
        
    def simpsons_function(a, b, φa, φb):
        I=h/6*(f(a)*φa+2*f((a+b)/2)+f(b)*φb)
        return I
    
   
    a=(e-1)*h
    b=e*h


    if quad_type == 2 :
        I = simpsons_function(a, b, φa, φb)

    elif quad_type == 1: 
        I = Trapezoidal_function(a, b, φa, φb)
    
    else:
        print("Invalid. Please enter either '1' or '2' ")
    
    
    return I