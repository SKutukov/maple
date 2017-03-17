import  numpy as np
import math

##################### load data #######################
import argparse as ap
# Get the path of matrix
parser = ap.ArgumentParser()
parser.add_argument("-m", "--matrix", help="Path to matrix file", required="True")
parser.add_argument("-e", "--eps", help="Path to matrix file", required="True")
args = vars(parser.parse_args())
# Read file
f= open(args["matrix"],'r')
A= np.matrix(f.read())
f.close()
eps=float(args["eps"])
########################################################

print A
