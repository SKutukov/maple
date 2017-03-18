import  numpy as np
import math
import array


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

def Yakobi(A,eps):
	return np.linalg.eig(A)

def opsiteSpectr(A,lmbd):
	return linalg.eig(A)[0]
def Viland(A,eps):
	return linalg.eig(A)[0][0]
