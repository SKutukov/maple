import  numpy as np
import math
from Norm import norm_inf 
from Norm import mul 
from Norm import spectr_radius
A= np.matrix('8.29381 0.995516 -0.560617 ; 0.995516 6.298198 0.595772 ;-0.560617 0.595772 4.997407')
b=np.matrix(' 0.766522 ; 3.844422 ; 3.844422 ')
print "A : "
print A
print "b : "
print b
#print "1)Solve gaus: "
x_ext= np.linalg.solve(A,b)
#print x_ext
#################################################################
def transformation(A,b):
	num_rows, num_cols=A.shape
	H=np.empty(A.shape)
	g=np.empty([num_rows])
	i=0
	while(i<num_rows):
		j=0
		while(j<num_cols):
			if(i!=j):
				H[i,j]=-A[i,j]/A[i,i]
			else:
				H[i,j]=0
			j=j+1
		i=i+1
	i=0
	while(i<num_rows):
		g[i]=b[i]/A[i,i]
		i=i+1
	return H,np.matrix(g).T
#print "2)transformation to x=Hx+g : "
H,g=transformation(A,b)
x_0=x_ext+np.matrix(' 0.2  0.2  0.2 ').T
eps=0.000001
def sum1(H,x,i):
	sum=0
	num_rows, num_cols=H.shape
	j=0
	while(j<=i):
		sum=sum+H[i,j]*x[j]
		j=j+1
	return sum

def sum2(H,x,i):
	sum=0
	num_rows, num_cols=H.shape
	j=i+1
	while(j<num_rows):
		sum=sum+H[i,j]*x[j]
		j=j+1
	return sum
print "7) up relacs: "
r=spectr_radius(H)
print r

def relacs(H,g,x_0,eps,max_k,q):
	num_rows, num_cols=H.shape
	dx=2*eps
	x=x_0
	k=0
	x1=x
	while(dx>eps):
		i=0
		while (i<num_rows):
			sum_1=sum1(H,x,i)
			sum_2=sum2(H,x1,i)
			x[i]=x[i]+q*(sum_1+sum_2-x[i]+g[i])
			i=i+1
		x1=x	
		dx=norm_inf(x_ext-x)
		k=k+1
	return x,k	

def print_relacs(H,g,x_0,eps,k_max,q):
	x_relacs,k_relacs=relacs(H,g,x_0,eps,k_max,q)
	print "k_relaks",k_relacs
#	print "x_relacs: "	
#	print x_relacs
	print "ERROR relacs:"
	print x_relacs-x_ext
q=2/(1+math.sqrt(1-r*r))
q=q-0.1
print "q: ",q
print_relacs(H,g,x_0,eps,25,q)

#print np.linalg.eig(H)[0]
	
