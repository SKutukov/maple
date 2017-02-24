import numpy as np
A= np.matrix('8.29381 0.995516 -0.560617 ; 0.995516 6.298198 0.595772 ;-0.560617 0.595772 4.997407')
A_exp = np.matrix('8.29381 0.995516 -0.560617 0.766522; 0.995516 6.298198 0.595772 3.844422 ; -0.560617 0.595772 4.997407 3.844422 ')
b=np.matrix(' 0.766522 ; 3.844422 ; 3.844422 ')

def norm_inf(mat):
	max=0
	num_rows, num_cols =mat.shape
	i=0
	while(i<num_rows):
		j=0
		sum=np.sum(np.absolute(mat[i]))
		if (max<sum):
			max=sum	
		i=i+1
	return max
eps=0.00001
def solve(mat):
	i=0
	num_rows, num_cols=A.shape
	#strait run
	while(i<num_rows):
		j=i+1
		a_i_i=mat[i,i]
		if (a_i_i<eps):
			print 'warning: div by small number'
		mat[i]=mat[i]/a_i_i
		while(j<num_cols):
			mat[j]=mat[j]-mat[j,i]*mat[i]				
			j=j+1
			
		i=i+1
	#print(mat)
	#reverce run
	i=num_rows-1;
	while(i>=0):
		j=i-1
		while(j>=0):
			mat[j]=mat[j]-mat[i]*mat[j,i]	
			j=j-1
		i=i-1
#	print(mat)
	return mat[:,num_cols]
def find_max(mat,row,col):
	i=0
	num_rows, num_cols=A.shape
	ret=mat[0,0]
	row,col
	while(i<num_rows):
		j=0
		while(j<num_cols):
			if(ret<mat[i,j]):
				ret=mat[i,j]
				row=i
				col=j						
			j=j+1
		i=i+1
	return ret	
def solve_ch(mat):
	i=0
	num_rows, num_cols=A.shape
	#strait run
	while(i<num_rows):
		j=i+1
		max_row=-1
		max_col=-1
		a_i_i=find_max(mat,max_row,max_col)
		tmp_row=mat[max_row]
		tmp_col=mat[:,max_col]
		mat[max_row]=mat[i]
		mat[:,max_col]=mat[:,i]
		mat[i]=tmp_row
		mat[:,i]=tmp_col
		mat[i]=mat[i]/a_i_i
		while(j<num_cols):
			mat[j]=mat[j]-mat[j,i]*mat[i]				
			j=j+1
			
		i=i+1
	#print(mat)
	#reverce run
	i=num_rows-1;
	while(i>=0):
		j=i-1
		while(j>=0):
			mat[j]=mat[j]-mat[i]*mat[j,i]	
			j=j-1
		i=i-1
#	print(mat)
	return mat[:,num_cols]

def inv_matrix(A):
	return np.linalg.inv(A)

print("matrix A: ") 
print(A)
print("vector b: ")
print(b)
print("solve exact")
print(np.linalg.solve(A,b))
print("solve")
x=solve(A_exp)
print x;
print("solve advance")
x=solve_ch(A_exp)
print x; 
A_inv=inv_matrix(A)
#print(A_inv)
norm = norm_inf(A)
print('norm: ')
print(norm)
norm_inv = norm_inf(A_inv)
print('norm inv: ')
print(norm_inv)
cond =norm*norm_inv;
print('cond:');
print(cond);

