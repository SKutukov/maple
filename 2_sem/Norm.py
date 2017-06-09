import numpy as np
def norm_inf(mat):
	max=0
	num_rows, num_cols = mat.shape
	for i in range(0, num_rows):
		sum = np.sum(np.absolute(mat[i]))
		if (max < sum):
			max = sum	
	return max
def mul(A,B):
	num_rowsA, num_colsA=A.shape
	num_rowsB, num_colsB=B.shape
	if (num_colsA!=num_rowsB):
		print ("study algebra!!!!!")
		return
	C=np.empty([num_rowsA, num_colsB])
	for i in range(0, num_rowsA):
		for j in range(0, num_colsB):
			for k in range(0, num_colsA):
				C[i,j] = C[i,j]+A[i,k]*B[k,j]
				
			
	return C

def max_fabs(h):
	max_value=np.fabs(h[0])
	for i in range(0,h.shape[0]):
		if(max_value<np.fabs(h[i])):
			max_value=np.fabs(h[i])
	return max_value
	
def spectr_radius(H):
	return max_fabs(np.linalg.eig(H)[0])

