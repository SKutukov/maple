import  numpy as np
import math
import array
from Norm import mul 

##################### load data #######################
import argparse as ap
# Get the path of matrix
parser = ap.ArgumentParser()
parser.add_argument("-m", "--matrix", help="Path to matrix file", required="True")
parser.add_argument("-e", "--eps", help="Calculation error ", required="True")
args = vars(parser.parse_args())
# Read file
f= open(args["matrix"],'r')
content = f.readlines()
A= np.matrix(content[0])
f.close()
eps=float(args["eps"])
########################################################
print "matrix A: "
print A	
print "eig numbers: "
print np.linalg.eig(A)
################### find max in updiagonal elem ##########
def find_max(mat):
	col=1
	row=0
	num_rows, num_cols=A.shape
	ret=np.fabs(mat[0,1])
	i=0
	while(i<num_rows):
		j=i+1
		while(j<num_cols):
			if(ret<np.fabs(mat[i,j])):
				ret=np.fabs(mat[i,j])
				row=i
				col=j	
			j=j+1
		i=i+1
	
	return row,col
################### find max in updiagonal elem ##########
def sign(a):
	return math.copysign(1, a)
def calc_c_s(mat,row,col):
	#prepare
	div=mat[row,row] - mat[col,col] 
	fi=math.atan(2*mat[row,col]*mat[row,col]/div)/2
	c=math.cos(fi)
	s=math.sin(fi)	

	#f=np.fabs(div)
	#d=math.sqrt(div*div + 4 * mat[row,col]*mat[row,col] )
	#c=math.sqrt(( 1 +  f/d)/2)
	#s=sign(mat[row,col]*(mat[row,row]-mat[col,col])) * math.sqrt(( 1 -  f/d)/2)
	#change elem
	res=np.identity(mat.shape[0])
	res[row,row]=c
	res[col,col]=c
	res[row,col]=-s
	res[col,row]=s
	return res
###########################################

def Yakobi(A,eps):
	A_k=A.copy()
	row,col=find_max(A_k)
	print row,col
	if (A[row,col]<eps):
		return A_k
	res=calc_c_s(A,row,col)
	A_k=mul(A_k,res)
	A_k=mul(res.T,A_k)
	
	row,col=find_max(A_k)
	print row,col
	if (A[row,col]<eps):
		return A_k
	res=calc_c_s(A,row,col)
	A_k=mul(A_k,res)
	A_k=mul(res.T,A_k)

	row,col=find_max(A_k)
	print row,col
	if (A[row,col]<eps):
		return A_k
	res=calc_c_s(A,row,col)
	A_k=mul(A_k,res)
	A_k=mul(res.T,A_k)

	
	#A[row,col]=A[col,row]=c*A[row,col]+s*A[col,row]

	return A_k
print Yakobi(A,eps)
def opsiteSpectr(A,lmbd):
	return linalg.eig(A)[0]
def Viland(A,eps):
	return linalg.eig(A)[0][0]
