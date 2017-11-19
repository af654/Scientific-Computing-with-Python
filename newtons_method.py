import numpy as np 
import math

def newton(x):
	#eulers = 2.718281
	regular = math.pow(x,3)-3*math.pow(x,2)+3*x-1
	derivative = 3*math.pow(x,2)-6*x+3
	#use newtons to approximate the x1 for the next iteration
	x1 = np.float(x-(regular/derivative))
	return x1

def rateofconvergence(list):
	length = len(list) - 1
	newlist = []

	#populate the list
	for i in range(length):
		if(list[i]-list[length]) != 0 :
			newlist.append(list[i] - list[length])

	A = np.zeros((len(newlist)-1,2))
	b = np.zeros((len(newlist)-1,1))

	for j in range(len(newlist)-1):
		A[j,0]= np.log(abs(newlist[j]))
		A[j,1]= 1.0
		b[j,0] = np.log(abs(newlist[j+1]))

	#solve for the slope using least squares
	m = np.linalg.solve(np.dot(A.T, A), np.dot(A.T, b))
	return m[0,0]

def main():
	x=2
	xnew = 1.5
	counter = 1
	list = []

	while (np.abs(x-xnew) > 0.0000001):
	  xnew = x 
	  x = newton(x)
	  list.append(x)
	  print "Iteration of" , counter , "gives x value of:" , x , "from the value of:", xnew
	  counter+=1

	convergence = rateofconvergence(list)
	print convergence

if __name__ == "__main__":
    main()


