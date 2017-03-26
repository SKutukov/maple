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
print np.linalg.eig(A)[0]
print "eig vectors: "
print np.linalg.eig(A)[1]
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
			#print i,j	
			j=j+1
		i=i+1
	
	return row,col
################### find max in updiagonal elem ##########
def sign(a):
	return math.copysign(1, a)
def calc_c_s(mat,row,col):
	#prepare
	div=mat[row,row] - mat[col,col] 

#	fi=math.atan(2.0*mat[row,col]/div)/2.0
#	c=math.cos(fi)
#	s=math.sin(fi)	

	f=np.fabs(div)
	d=math.sqrt(div*div + 4 * mat[row,col]*mat[row,col] )
	c=math.sqrt(( 1 +  f/d)/2)
	s=sign(mat[row,col]*(mat[row,row]-mat[col,col])) * math.sqrt(( 1 -  f/d)/2)
	res=np.identity(mat.shape[0])
 	res[row,row]=c
	res[col,col]=c
	res[row,col]=-s
	res[col,row]=s
	return res,c,s
##########################################

def Yakobi(A,eps):
	A_k=A.copy()
	i=0
	row,col=find_max(A_k)
	max_elem=np.fabs(A_k[row,col])
	X=np.identity(A.shape[0])
	while (max_elem>eps and i<4):
		#print row,col
		res,c,s=calc_c_s(A,row,col)
		#print res
		X=mul(X,res)
		#A_k=res.T*A_k*res		
		
		a_i_i=A_k[row,row].copy()
		a_j_j=A_k[col,col].copy()
		a_i_j=A_k[row,col].copy()
		a_j_i=A_k[col,row].copy()
		
		A_i=A_k[:,row].copy()
		A_j=A_k[:,col].copy()
		
		#A_k(*,i)=.......
		A_k[row,:]=(c*A_i+s*A_j).T
		A_k[:,row]=(c*A_i+s*A_j)
		#A_k(*,j)=.......
		A_k[col,:]=(-s*A_i+c*A_j).T		
		A_k[:,col]=(-s*A_i+c*A_j)
		#A_k(i,i)=.......
		A_k[row,row]=c*c*a_i_i+2*c*s*a_i_j+s*s*a_j_j
		#A_k(j,j)=.......
		A_k[col,col]=s*s*a_i_i-2*c*s*a_i_j+c*c*a_j_j
		#A_k(i,j)=A_k(j,i)=.........		
		A_k[row,col]=A_k[col,row]=(c*c-s*s)*a_i_j+c*s*(a_j_j-a_i_i)	
		#A_k[row,col]=A_k[col,row]=0;
		row,col=find_max(A_k)
		max_elem=np.fabs(A_k[row,col])
		#print A_k
		#print X
		i=i+1
	return np.diagonal(A_k),X
print "Yakobi : "
D,X= Yakobi(A,eps)
print D
print X
print X[:,0]/np.linalg.norm(X[:,0])
print X[:,1]/np.linalg.norm(X[:,1])
print X[:,2]/np.linalg.norm(X[:,2])
Y_0=np.linalg.eig(A)[1][0]*0.9+np.linalg.eig(A)[1][1]*0.5+np.linalg.eig(A)[1][2]*0.6
def sc(a,b):
	return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]
def scalar(A,Y,eps):
	Y_k=Y.T
	dl=2*eps
	lm=1
	while(dl>eps):
		Y_k_1=A*Y_k
		tm=lm
		lm=sc(Y_k_1,Y_k)/sc(Y_k,Y_k)
		Y_k=Y_k_1
		dl=np.fabs(lm-tm)			
	return lm,Y_k
def opsiteSpectr(A,lmbd,eps):
	E=np.identity(A.shape[0])
	B=A-lmbd*E
	return scalar(B,Y_0,eps)+lmbd
	#if (lmbd<0):
	#	return np.linalg.eig(A)[0][1]	
	#else:
	#	return np.linalg.eig(A)[0][1]
print "opsite Spectr"
print opsiteSpectr(A,np.linalg.eig(A)[0][0],eps)	
def Viland(A,lm_0,eps):
	W=A-lm_0*np.identity(A.shape[0])
	W=np.linalg.inv(W)
	print np.linalg.eig(W)[1]
	Y_0=np.linalg.eig(W)[1][0]*0.9+np.linalg.eig(W)[1][1]*0.5+np.linalg.eig(W)[1][2]*0.6
	lambd,X=scalar(W,Y_0,eps)			
	return 1/lambd+lm_0,X
lm,X = scalar(A,Y_0,eps)
print "scalar: : "
lm=lm[0,0]
print (X/np.linalg.norm(X))
print lm
print "err: "
print lm-np.linalg.eig(A)[0][0]
print "Viland : "
l,X=Viland(A,lm,eps)
print l
X=(X/np.linalg.norm(X))
print X
print "err"
print l-np.linalg.eig(A)[0][0]
