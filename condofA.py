import numpy as np
import math
from scipy.linalg import lu

def regularnorm(A):
	A = np.float64(A)
	#print A
	n=A.shape[1]
	m=A.shape[0]
	value = 0
	added_values = 0
	#mxn matrix
	#the array column_sums stores each column sum of the matrix A
	column_sums = np.zeros(n)

	for col in range(A.shape[1]):
    		#print A[:,col]
    		added_values=0
    		for row in range(A.shape[0]):
    			value = A[row,col]
    			added_values += abs(value)
    		column_sums[col]=added_values

	maximum = column_sums[0]
	for value in column_sums:
		if maximum < value:
			maximum = value
	
	print ('The L1 norm of the matrix A is: '), maximum
	return maximum

# Forward substitution for lower triangular system UT
def Forward_Substitution(UT,c):
	UT = np.float64(UT)
   	m=UT.shape[0]
   	n=UT.shape[1]
   	if(m!=n):
       		print('Matrix is not square!')
       		return
   	v = np.zeros(n)
   	for j in range(0,n):
   		if UT[j,j] == 0:
   			print('Matrix is singular!')
   			return          # matrix is singular
   		v[j] = c[j]/UT[j,j]
   		for i in range(j+1,n):
   			c[i] = c[i] - UT[i,j]*v[j]
   	return v

def solveforv(U_transpose):
	UT = np.float64(U_transpose)
	m=UT.shape[0]
	n=UT.shape[1]
	c=np.zeros(n)

	for component in range(0,n):
		c[component] = -1
	print c
	
	v = Forward_Substitution(UT,c)

	#returns the v vector for future usage for system LTy=v
	return v

def solvefory(L,U):
	L_transpose = L.transpose()
	U_transpose = U.transpose()
	v = solveforv(U_transpose)
	print v
	#compute the y in the system LTy=v
	y = np.linalg.solve(L_transpose,v)

	return y

def main():
	A = np.matrix([[92,66,25],[-73,78,24],[-80,37,10]])
	#compute the norm of A
	regularnorm(A);

	#compute the LU decomposition of the linear system Az=y, where A=LU
	P,L,U = lu(A,permute_l=False)

	print L
	print U
	y = solvefory(L,U)
	#y=np.random.rand(3,1)

	print ('The vector y in the linear system Az=y is: '), y
	z = np.linalg.solve(A,y)

	print z

	z = np.float64(np.linalg.norm(z))
	y = np.linalg.norm(y)

	A_inverse = z/np.float64(y)

	print ('The norm of the inverse of matrix A using the approximation method is: '), A_inverse

	true_value = np.linalg.cond(A,1)
	print true_value


if __name__ == "__main__":
    main()
