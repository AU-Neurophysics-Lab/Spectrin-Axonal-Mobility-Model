import csv
import string
import math
import cmath
import numpy as np	#import package numpy and give it the designation "np"
import os

####################################################
def check_make_dir(dir):
	d = os.path.dirname(dir)
	print("directory is %s\n" % (d) )
	print("Checking for directory\n")

	if not os.path.exists(d):
		os.makedirs(d)
		print("Directory didn't exist, and was made.")
	return ()

####################################################