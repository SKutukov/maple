import  numpy as np
from Norm import norm_inf 
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
	return H,g
print "2)transformation to x=Hx+g : "
H,g=transformation(A,b)
print H,g
print norm_inf(H)


#print ('matrix H: ') 
#print (np.linalg.norm(H))


#print ('x_ext: ')
#print x_ext
#x=np.matrix('0.4 ; 0.4 ');
def iteration(H,x,eps):
	dx=2*eps
	while(dx>eps):
		print(dx)
		x1=x
		x=H*x
		dx=np.linalg.norm(x1-x)
	return x

#print iteration(H,x,0.0000001)

