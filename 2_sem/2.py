import  numpy as np
import math
from Norm import norm_inf 
from Norm import mul 
from Norm import spectr_radius
A= np.matrix('8.29381 0.995516 -0.560617 ; 0.995516 6.298198 0.595772 ;-0.560617 0.595772 4.997407')
b=np.matrix(' 0.766522 ; 3.844422 ; 3.844422 ')

print "1)Solve gaus: "
x_ext= np.linalg.solve(A,b)
print x_ext

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
print "2)transformation to x=Hx+g : "
H,g=transformation(A,b)
print "H: "
print H
print "g: "
print g
print "Norm of H: "
print norm_inf(H)
eps=0.0001
print "3) find k, ||x_*-x^k||<eps"
def find_K(H,g,x_0,eps):
	d=norm_inf(x_0)+norm_inf(g)/(1-norm_inf(H))
	eps_d=eps/d
	return math.ceil(np.log(eps/d)/np.log(norm_inf(H)))
x_0=x_ext+np.matrix(' 0.2  0.2  0.2 ').T
k_teor=find_K(H,g,x_0,eps)	
print "k_teor: ", k_teor
print "4) simple iteration"
def iteration(H,g,x_0,eps):
	dx=2*eps
	x=x_0
	k=0
	h=norm_inf(H)/(1-norm_inf(H))
	while(dx>eps):
		#print(dx)
		x1=x
		x=H*x+g
		dx=h*norm_inf(x1-x)
		k=k+1
	return x,k,dx,x1
x_iter,k_ext,dx,x_k_1= iteration(H,g,x_0,eps)
print x_iter
print "k_ext: ", k_ext
print "ERROR x_ext-x_iter: "
print x_ext-x_iter
print "ERROR aprior: "
def aprior(H,g,x_0,k):
	norm=norm_inf(H)
	d=norm_inf(x_0)+np.linalg.norm(g)/(1-norm)
	return math.pow(norm,k)*(d)
print aprior(H,g,x_0,k_teor)
print "ERROR apostorior: "
print dx
r=spectr_radius(H)
x_lust=x_k_1+(x_iter-x_k_1)/(1-r)
print "lusternick :",x_lust
print "ERROR lust: "
print x_lust-x_ext

print "5) Zeidel"
def transform_H(H):
	num_rows, num_cols=H.shape
	i=0
	H_l=np.empty(H.shape)
	H_r=np.empty(H.shape)
	#H_l
	while(i<num_rows):
		j=i
		while(j<num_cols):
			H_r[i,j]=H[i,j]		
			j=j+1
		i=i+1

	i=0
	while(i<num_rows):
		j=0
		while(j<i):
			H_l[i,j]=H[i,j]
			j=j+1
		i=i+1	
	return H_l,H_r
def zeidel(H,g,x_0,x_ext,eps,max_k):
	H_l,H_r=transform_H(H)	
	num_rows, num_cols=H.shape
	E=np.identity(num_rows)
	M=np.linalg.inv(E-H_l)
	mat=np.empty(H.shape) 
	mat=mul(M,H_r)
	c=M*g
	x,k_ext,dx,x_k_1= iteration(mat,c,x_0,eps)
	return x,k_ext

x_zeid,k_zeid=zeidel(H,g,x_0,x_ext,eps,25)
print "k_zeid: ", k_zeid
print x_zeid
print "ERRROR zeid: "
print (x_zeid-x_ext)

print "7) up relacs: "
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
r=spectr_radius(H)
q=2/(1+math.sqrt(1-r*r))
def relacs(H,g,x_0,eps,max_k,q):
	num_rows, num_cols=H.shape
	dx=2*eps
	x=x_0
	k=0
	while((dx>eps) & (k<max_k)):
		#print(dx)
		dx_it=2*eps
		s=0
		while((dx_it>(eps/1000))&(s<max_k)):
			x1=x
			i=0
			while (i<num_rows):
				sum_1=sum1(H,x,i)
				sum_2=sum2(H,x,i)
				x[i]=x[i]+q*(sum_1+sum_2-x[i]+g[i])
				i=i+1
			s=s+1
			dx_it=norm_inf(x-x1)
		#print s
		dx=norm_inf(x-x_ext)
		#print x1-x
		k=k+1
	return x,k
def print_relacs(H,g,x_0,eps,k_max,q):
	x_relacs,k_relacs=relacs(H,g,x_0,eps,k_max,q)
#	print "k_relaks",k_relacs
#	print "x_relacs: "	
#	print x_relacs
	print "ERROR relacs:"
	print x_relacs-x_ext
print "q: ",1			
print_relacs(H,g,x_0,eps,25,1)
print "q: ",q			
print_relacs(H,g,x_0,eps,25,q)
print "q: ",q+0.1
print_relacs(H,g,x_0,eps,25,q+0.1)
print "q: ",q-0.1
print_relacs(H,g,x_0,eps,25,q-0.1)
	
