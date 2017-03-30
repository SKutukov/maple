import  numpy as np
from Norm import mul 

##################### load data #######################
import argparse as ap
# Get the path of matrix
parser = ap.ArgumentParser()
#parser.add_argument("-m", "--matrix", help="Path to matrix file", required="True")
parser.add_argument("-e", "--eps", help="Calculation error ", required="True")
args = vars(parser.parse_args())
# Read file
#f= open(args["matrix"],'r')
#content = f.readlines()
#A= np.matrix(content[0])
#f.close()
eps=float(args["eps"])
########################################################
from sympy import *
from sympy import init_printing
init_printing()
x= Symbol('x')
y= Symbol('y') 
H=sin(x*y)
print (H)
f=x-0.6
a=0
b=1
n=4
A=[]
for i in range(0, n):
	A.append(1/n)
X=[]
for i in range(0, n):
	c=(a+b)/(2*n)
	X.append(c*(2*i+1))
print ("nodes of kvadrature: ")
print (X)

def fred(H,A,X,n,a,b,f):
	g=np.empty([n,1])
	for i in range(0, n):
		g[i]=f.subs(x,X[i]).evalf()

	D=np.identity(n)
	for i in range(0, n):
		for j in range(0, n):
			D[i,j]=D[i,j]-A[j]*H.subs(y,X[j]).subs(x,X[i]).evalf()
	
	print(D)
	C=np.linalg.solve(D,g)
	print (C)
	u=f
	for i in range(0, n):
		u=u+A[i]*H.subs(y,X[i]).evalf()*C[i]	
		i=i+1
	return u
from sympy.plotting import plot
u=fred(H,A,X,n,a,b,f)
print (u)
Y=[u.subs(x,a).evalf(), u.subs(x,b).evalf()]
plot(u,(x, a, b))

