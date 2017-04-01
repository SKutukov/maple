import  numpy as np
from Norm import mul 
from Norm import max_fabs

##################### load data #######################
import argparse as ap
# Get the path of matrix
parser = ap.ArgumentParser()
#parser.add_argument("-m", "--matrix", help="Path to matrix file", required="True")
parser.add_argument("-e", "--eps", help="Calculation error ", required="True")
args = vars(parser.parse_args())
eps=float(args["eps"])
########################################################
from sympy import *
from sympy import init_printing
init_printing()
x= Symbol('x')
y= Symbol('y') 
H=sin(x*y)
f=x-0.6
a=0
b=1
n=4
def norm(u,v,a,b):
	X=np.empty([3,1])
	X[0,0]=u.subs(x,a)-v.subs(x,a)
	X[1,0]=(u.subs(x,(a+b)/2)-v.subs(x,(a+b)/2))
	X[2,0]=(u.subs(x,b)-v.subs(x,b))
	#print (X)
	return max_fabs(X)
def fred_iter(H,A,X,n,a,b,f):
	g=np.empty([n,1])
	for i in range(0, n):
		g[i]=f.subs(x,X[i])

	D=np.identity(n)
	for i in range(0, n):
		for j in range(0, n):
			D[i,j]=D[i,j]-A[j]*H.subs(y,X[j]).subs(x,X[i])
	
	#print(D)
	C=np.linalg.solve(D,g)
	#print (C)
	u=f
	for i in range(0, n):
		u=u+A[i]*H.subs(y,X[i])*C[i]	
		i=i+1
	return u

def fredgolm(H,f,n,a,b):
	du=2*eps
	u=f
	k=n
	j=1
	while (du>eps):
		temp=u.copy()
		A=[]
		for i in range(0, k):
			A.append(1/k)
		#print ("A: ")		
		#print (A)
		X=[]
		for i in range(0, k):
			c=(a+b)/(2*k)
			X.append(c*(2*i+1))
		#print ("nodes of kvadrature: ")
		#print (X)
		u= fred_iter(H,A,X,k,a,b,f)
		du=norm(u,temp,a,b)
		print (("i=%d,n=%d,error=%e") % (j,k,du))		
		k=2*k
		j=j+1
		
	return u
from sympy.plotting import plot
u=fredgolm(H,f,n,a,b)
print (u)
Y=[u.subs(x,a).evalf(), u.subs(x,b).evalf()]
plot(u,(x, a, b))

