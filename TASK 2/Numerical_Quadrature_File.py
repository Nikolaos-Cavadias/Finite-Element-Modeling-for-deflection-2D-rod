# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 13:25:19 2024

@author: Nikolas Cavadias
"""
from sympy import *

xx = symbols('x')

integrantfunc = None  # Define the variable globally

def input_func():
    return sympify(input ("f = "))

def f(j):
    if integrantfunc is None:
        raise ValueError("Function not defined yet")
    return integrantfunc.subs(xx, j)

def Numerical_Quadrature(e,i,quad_type,h):
    
    global integrantfunc  # Access the global variable
    
    if integrantfunc is None:
        integrantfunc = input_func()
    else:
        pass

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