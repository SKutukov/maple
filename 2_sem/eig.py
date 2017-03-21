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

	fi=math.atan(2.0*mat[row,col]/div)/2.0
	c=math.cos(fi)
	s=math.sin(fi)	

#	f=np.fabs(div)
#	d=math.sqrt(div*div + 4 * mat[row,col]*mat[row,col] )
#	c=math.sqrt(( 1 +  f/d)/2)
#	s=sign(mat[row,col]*(mat[row,row]-mat[col,col])) * math.sqrt(( 1 -  f/d)/2)
	return c,s
##########################################

def Yakobi(A,eps):
	A_k=A.copy()
	i=0
	max_elem=2*eps
	while (max_elem>eps):
		row,col=find_max(A_k)
		max_elem=A_k[row,col]
		print max_elem
		c,s=calc_c_s(A,row,col)
		
		a_i_i=A_k[row,row]
		a_j_j=A_k[col,col]
		a_i_j=A_k[row,col]
		a_j_i=A_k[col,row]
		
		A_i=A_k[:,row]
		A_j=A_k[:,col]
		
		#A_k(*,i)=.......
		#A_k[:,row]=c*A_i+s*A_j
		#A_k(*,j)=.......
		#A_k[:,col]=-s*A_i+c*A_j
		#A_k(i,i)=.......
		#A_k[row,row]=c*c*a_i_i+2*c*s*a_i_j+s*s*a_j_j
		#A_k(j,j)=.......
		#A_k[col,col]=s*s*a_i_i-2*c*s*a_i_j+c*c*a_j_j
		#A_k(i,j)=A_k(j,i)=.........		
		A_k[row,col]=A_k[col,row]=(c*c-s*s)*a_i_j+c*s*(a_j_j-a_i_i)	
		
		print A_k
		i=i+1
		
	return A_k
print Yakobi(A,eps)
def opsiteSpectr(A,lmbd):
	return linalg.eig(A)[0]
def Viland(A,eps):
	return linalg.eig(A)[0][0]
