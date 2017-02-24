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

def solve(mat):
	i=0
	num_rows, num_cols=A.shape
	#strait run
	while(i<num_rows):
		j=i+1
		a_i_i=mat[i,i]
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
		print(mat)
	return mat[:,num_cols]

def solve_ch(A,b):
	return np.linalg.solve(A,b)

def inv_matrix(A):
	return np.linalg.inv(A)

print("matrix A: ") 
print(A)
print("vector b: ")
print(b)
print("solve")
x=solve(A_exp)
print x;
print("solve advance")
x=solve_ch(A,b)
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

