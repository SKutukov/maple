import numpy as np
A= np.matrix('-401.52 200.16 ; 1200.96  -601.68')
b1=np.matrix('200 ; -600')
b2=np.matrix('199 ; -601')

print("matrix A: ") 
print(A)
print("vector b1 ")
print(b1)
print("vector b2 ")
print(b2)
print("solve 1")
x1=np.linalg.solve(A,b1)
print x1;
print("solve 2")
x2=np.linalg.solve(A,b2)
print x2;

A_inv=np.linalg.inv(A)
#print(A_inv)
norm = np.linalg.norm(A)
norm_inv = np.linalg.norm(A_inv)
cond =norm*norm_inv;
print('cond:');
print(cond);

norm_b=np.linalg.norm(b1);
#print(norm)
#print(norm_b)
db=1;
dA=0.0001;
eps_teor=(cond/(1-cond*dA/norm))*(db/norm_b+dA/norm)
print('teoretical error:')
print(eps_teor)
print('error:')
print(x1-x2)
print('fact error:')
print(np.linalg.norm(x1-x2)/np.linalg.norm(x1))
