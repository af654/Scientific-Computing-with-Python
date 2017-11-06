import numpy as np
from numpy import linalg as LA
from numpy import dot
import math
import matplotlib.pyplot as plt
from scipy import optimize

# LU decomposition of square systems
def Gaussian_Elimination(A):
   A = np.float32(A)
   n=A.shape[1]
   m=A.shape[0]
   if (m != n):
       print('Matrix is not square!')
       return
   for k in range(0,n-1):
       if A[k,k] == 0:
           return
       for i in range(k+1,n):
           A[i,k]=A[i,k]/A[k,k]
       for j in range(k+1,n):
           for i in range(k+1,n):
               A[i,j]-=A[i,k]*A[k,j]
   s = (n,n)
   L = np.zeros(s)
   U = np.zeros(s)
   for i in range(n):
       L[i,i] = 1
       for k in range(i):
           L[i,k] = A[i,k]

       for j in range(i,n):
           U[i,j] = A[i,j]
   return (L,U)

# Forward substitution for lower triangular system
def Forward_Substitution(A,b):
   m=A.shape[0]
   n=A.shape[1]
   if(m!=n):
       print('Matrix is not square!')
       return
   x = np.zeros(n)
   for j in range(0,n):
       if A[j,j] == 0:
           print('Matrix is singular!')
           return          # matrix is singular
       x[j] = b[j]/A[j,j]
       for i in range(j+1,n):
           b[i] = b[i] - A[i,j]*x[j]
   return x

#Back substitution for upper triangular system
def Back_Substitution(A, b):
  n = A.shape[0]
  x = np.zeros(n)
  for j in range(n - 1, -1, -1):
      if A[j, j] == 0:
          print('Matrix is singular!')
          return  # matrix is singular
      x[j] = b[j] / A[j, j]
      for i in range(0, j):
          b[i] = b[i] - A[i, j] * x[j]
  return x

def plotfigure(y_hat):
  x=np.transpose(np.matrix([0.0,  1.0,  2.0,  3.0,  4.0,  5.0]))
  y=np.transpose(np.matrix([1.0,2.7,5.8,6.6,7.5,9.9]))
  for i in range(6):
    print("polynomial %d"%i)
    A = np.matrix(np.zeros((6,i+1)))
    #print(A)
    for j in range(i+1):
      A[:,j] = np.power(x,j)
    # (A^T)Ax=(A^T)b
    a = np.linalg.inv(np.transpose(A)*A)*np.transpose(A)*y 
    fx = A[:,:i+1] * a
    res =  fx- y
    #print(res)
    ssq = np.sum(np.power(res,2))
    print("Sum Squared error: %e\n"%ssq)
    plt.figure()
    plt.plot(x,y,'.',x,fx)
    plt.show()

def main():
    #normal equation least squares
    A = np.float64(np.matrix([[1,0],[1,1],[1,2],[1,3],[1,4],[1,5]]))
    b = np.float64(np.array([1,2.7,5.8,6.6,7.5,9.9]))
    #compute the condition number of a matrix
    LA.cond(A)
    print A
    A_transpose = A.transpose()
    print A_transpose

    #numpy matrix vector multiplication
    print A_transpose.dot(A)
    print A_transpose.dot(b)

    #A_new = np.float64(np.matrix([[20516,4425,979,225],[4425,979,225,55],[979,225,55,15],[225,55,15,5]]))
    A_new = np.float64(np.matrix([[6,15],[15,55]]))
    #b_new = np.float64(np.array([33.5,113.6,4528,19448,87376,404096]))
    b_new = np.float64(np.array([33.5,113.6]))
    print A_new
    print b_new
    x_hat = Part_a(A_new, b_new)
    print("The answer to the least squares solution using normal equation is\n", x_hat)
    print("You can see that the values are not close to the correct answers")

    print("check that this is the actual answer \n")
    print A.dot(x_hat)

    y_hat = A.dot(x_hat)
    
    r1 = np.float32(b_new - np.dot(A_new,x_hat))
    print("residual is",r1)

    #compute the condition number of a matrix
    condition_number = math.pow((LA.cond(A_new)),2)

    print("conditional number is", condition_number)

    plotfigure();

def Part_a(A, b):
   (L, U) = Gaussian_Elimination(A)
   y = Forward_Substitution(L, b)
   y1 = np.array(y.copy())     # y1 is correct
   x1 = Back_Substitution(U, y1.T)
   return np.float32(x1)

if __name__ == "__main__":
    main()
