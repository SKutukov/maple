import  numpy as np
#from Norm import mul 

##################### load data #######################
import argparse as ap
# Get the path of matrix
parser = ap.ArgumentParser()
#parser.add_argument("-m", "--matrix", help="Path to matrix file", required="True")
parser.add_argument("-e", "--eps", help="Calculation error ", required="True")
args = vars(parser.parse_args())
# Read file
#f= open(args["matrix"],'r')
#content = f.readlines()
#A= np.matrix(content[0])
#f.close()
eps=float(args["eps"])

########################################################
from sympy import *
from sympy import init_printing
init_printing()
x = Symbol('x')
print exp(pi * sqrt(163)).evalf(50)
 
