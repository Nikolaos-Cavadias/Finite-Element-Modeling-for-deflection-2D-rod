import numpy as np
import sys
import matplotlib.pyplot as plt
from Numerical_Quadrature_File1 import Numerical_Quadrature
from sympy import *
##
#    Input the number of finite elements and boundary conditions
##


xx = symbols('x')
integrantfunc = sympify(input ("f = "))      # -12*x^2+18*x+6
def f(j):
    return integrantfunc.subs(xx, j)

ε1_iter=[]
ε2_iter=[]
ε3_iter=[]
Norm_ε1=[]
Norm_ε2=[]
Norm_ε3=[]


y0 = np.zeros(1)
y1 = np.zeros(1)
y0[0] = float(input("Input the BC at x=0: "))             #  BC at x=0
y1[0] = float(input("Input the BC at x=1: "))              #  BC at x=1  

quad_type=1

Size = [16, 32, 64, 128, 256, 512] 
for Ne in Size :
    

    ##
    #    Calculate the number of nodes and the element size
    ##  
    N = Ne + 1                  #  number of nodes
    h = 1 / Ne                  #  element size 
    ##
    #    Initialize the global matrix and the right-hand side vector
    ##
    A = np.zeros(N)             #  coefficient matrix main diagonal
    B = np.zeros(N - 1)         #  coefficient matrix sub-diagonal
    C = np.zeros(N - 1)         #  coefficient matrix super-diagonal
    f = np.zeros(N)             #  right-hand side vector
    ##
    #    Assemble the element matrices and right-hand sides
    ##
    for e in range(Ne):                         #  loop over finite elements
        Ae = np.zeros([2,2])                    #  element coefficient matrix
        fe = np.zeros([2])                      #  element right-hand side
        for i in range(2):                      #  loop over the nodes in an element
            ## ----------------------------------------------------------------------------- ##
            ##  HERE IS THE CALL TO THE NUMERICAL QUADRATURE FUNCTION THAT YOU NEED TO WRITE ##
            fe[i]=Numerical_Quadrature(e,i,quad_type,h,integrantfunc);  #  elelemt right-hand side vector contribution  
            ##  ---------------------------------------------------------------------------- ##
            for j in range(2):                  #  loop over nodes in an element
                Ae[i,j]=(-1)**(i+j)/h           #  element matrix contribution
    ##
    #    Impose the boundary conditions
    ##            
        if e==0:                                #  the first element at x=0
           Ae[0,0]=1 
           Ae[0,1]=0 
           fe[0]=y0[0];                         #  y=y0 at x=0;
           fe[1]=fe[1]-Ae[1,0]*y0[0]            #  modify the second equation
           Ae[1,0]=0 
        if e==Ne-1:                             #  the last element at x=1
           Ae[1,1]=1 
           Ae[1,0]=0 
           fe[1]=y1[0]                          #  y=y1 at x=1 
           fe[0]=fe[0]-Ae[0,1]*y1[0]            #  modify the second equation
           Ae[0,1]=0 
    ##
    #    Transfer the elemnt matrix and the element right-hand side to global data structures
    ##      
        A[e]=A[e]+Ae[0,0]                       #  diagonal element Ae[0,0]
        A[e+1]=A[e+1]+Ae[1,1]                   #  diagonal element Ae[1,1]
        B[e]=Ae[1,0]                            #  sub diagonal element Ae[1,0]
        C[e]=Ae[0,1]                            #  super diagonal element Ae[0,1] 
        f[e]=f[e]+fe[0]                         #  rhs element fe[0]
        f[e+1]=f[e+1]+fe[1]                     #  rhs element fe[1]
    ##
    #    Solve the linear system
    ##
    
    AA=A.copy()
    BB=B.copy()
    CC=C.copy()
    ff=f.copy()
    
           #  direct solver
    from Gaussian_Elimination_File1 import Gaussian_Elimination
       ## ----------------------------------------------------------------------------- ##
       ## HERE IS THE CALL TO GUSSIAN ELIMINATION FUNCTION FROM LAB 1                   ##
    xD=Gaussian_Elimination(A,B,C,f)
       ## ----------------------------------------------------------------------------- ##
           #  Jacobi iteration
    from Jacobi_Iteration_File1 import Jacobi_Iteration
    x0=np.zeros(N)
    Tol = [1e-4, 1e-6, 1e-8]
    
    for i, epst in enumerate(Tol):
   ## ----------------------------------------------------------------------------- ##
   ##  HERE IS THE CALL TO THE jacobi method FUNCTION THAT YOU NEED TO WRITE        ##
       
       [xJ,rnk,iter]=Jacobi_Iteration(AA,BB,CC,ff,x0,epst)   
       if i==0:
           ε1_iter.append(iter)
           Norm_ε1.append(np.linalg.norm(xD - xJ))
       elif i==1:
           ε2_iter.append(iter)
           Norm_ε2.append(np.linalg.norm(xD - xJ))               
       elif i==2:
           ε3_iter.append(iter)
           Norm_ε3.append(np.linalg.norm(xD - xJ))
       else:
           print('Invalid')
       
       
       
   ## ----------------------------------------------------------------------------- ##
    print("Method converged in ",iter," iterations")
    print("Norm of the final residual %.6E" % rnk)
              


##
#    Plot the solution   9*sin(3.14*x)
##

plt.plot(Size, ε1_iter, 'g', linewidth=1.5, label='epsilon = 1e-4')
plt.plot(Size, ε2_iter, 'b', linewidth=1.5, label='epsilon = 1e-6')
plt.plot(Size, ε3_iter, 'r', linewidth=1.5, label='epsilon = 1e-8')      #  plot the solution
plt.xlabel("Discrete problem size N")                            #  label x-axis
plt.ylabel("No of iterations (N)")                        #  label y-axis
plt.title('No of iterations - Problem size')       #  figure title
plt.xlim((min(Size)-16,max(Size)+50))                     #  x-axis limits
plt.legend()
plt.ylim((min(ε3_iter)-0.5,max(ε3_iter)+1000))          #  y-axis limits
plt.show()


plt.plot(Size, Norm_ε1, 'g', linewidth=1.5, label='epsilon = 1e-4')
plt.plot(Size, Norm_ε2, 'b', linewidth=1.5, label='epsilon = 1e-6')
plt.plot(Size, Norm_ε3, 'r', linewidth=1.5, label='epsilon = 1e-8')      #  plot the solution
plt.xlabel("Discrete problem size N")                            #  label x-axis
plt.ylabel("Norms of the error  xD - xJ (N)")                        #  label y-axis
plt.title('Norms of the error  xD - xJ (N) - Problem size')       #  figure title
plt.xlim((min(Size)-16,max(Size)+50))                     #  x-axis limits
plt.legend()
plt.ylim((min(Norm_ε3)-0.5,max(Norm_ε1)+2))          #  y-axis limits
plt.show()





