import  numpy as np

H = np.matrix('-0.4 0.3; 0.5 -0.6')
zero = np.matrix('0 ; 0')
print ('matrix H: ') 
print (np.linalg.norm(H))

x_ext= np.linalg.solve(H,zero)
print ('x_ext: ')
print x_ext
x=np.matrix('0.4 ; 0.4 ');
def iteration(H,x,eps):
	dx=2*eps
	while(dx>eps):
		print(dx)
		x1=x
		x=H*x
		dx=np.linalg.norm(x1-x)
	return x

print iteration(H,x,0.0000001)

