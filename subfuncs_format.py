import csv
import string
import math
import cmath
import numpy as np	#import package numpy and give it the designation "np"
import os

####################################################
def fmtcols(tmp_arrayp, file):
	out=open(file,"w")
	for x in zip(*tmp_arrayp):
		out.write("{1}\t{0}\n".format(*x))
	return ()
####################################################
print("fmtcols imported successfuly.")