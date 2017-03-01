import numpy as np
import array
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
def mul(A,B):
	num_rowsA, num_colsA=A.shape
	num_rowsB, num_colsB=B.shape
	if (num_colsA!=num_rowsB):
		print "study algebra!!!!!"
		return
	C=np.empty([num_rowsA,num_colsB])
	i=0
	while (i<num_rowsA):
		j=0
		while(j<num_colsB):
			k=0
			while(k<num_colsA):
				C[i,j]=C[i,j]+A[i,k]*B[k,j]
				k=k+1
			j=j+1
		i=i+1		
	return C
