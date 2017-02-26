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
