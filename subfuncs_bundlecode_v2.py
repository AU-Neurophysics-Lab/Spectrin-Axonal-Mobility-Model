import csv
import string
import math
import cmath
import numpy as np	#import package numpy and give it the designation "np"
import os
import random


#---------------------------------------------------------------#
#  Define generator function for bundle			 	#
#	Input: outputfilename, Obstruction Prob, Bundle Size	#
#---------------------------------------------------------------#

#---------------------------------------------------------------#
#	Version 2 creates a 2-dimensional bundle		#
#       where y-lattice is number of MTs			#
#	x-lattice is un-scaled or scaled			#
#---------------------------------------------------------------#

def generate_bundle_v2_1(bundle_flenam, x_max, y_MT, y_PF):

	bundle_output = open(bundle_flenam + "-2d.bnd", 'w')							# Open Output file for static bundle
	MaxBND = y_MT*y_PF*2
	twoD_Bundle_Distribution_static = [[float(0.0) for x in range(int(x_max))] for z in range(0,int(MaxBND))]	# Generate "array" to store bundle data

	return(twoD_Bundle_Distribution_static)
#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#	

def populate_bundle_ends(two_dim_Bundle, x_max, y_MT, y_PF, x_mean, sig_mt, m):
	x_j = 0
	MaxBND = y_MT*y_PF*2
	Start = y_MT*y_PF + m*y_PF
	End = Start + y_PF
	while x_j < x_max:
		num = random.gauss(x_mean, sig_mt)
		if float(num) > float(0.0) : x_j = x_j + float(num)		
		if x_j < x_max: 
			for r in range(Start,int(End)):
				two_dim_Bundle[r][int(x_j)] = 1
			print(int(x_j))
			print(two_dim_Bundle[0][int(x_j)])
	print("\n")
	return(two_dim_Bundle)

#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#	

def populate_bundle_ends_exp(two_dim_Bundle, x_max, y_MT, y_PF, x_mean, sig_mt, m):
	x_j = 0
	MaxBND = y_MT*y_PF*2
	Start = y_MT*y_PF + m*y_PF
	End = Start + y_PF
	while x_j < x_max:
		num = random.expovariate(x_mean)
		if float(num) > float(0.0) : x_j = x_j + float(num)		
		if x_j < x_max: 
			for r in range(Start,int(End)):
				two_dim_Bundle[r][int(x_j)] = 1
			print(int(x_j))
			print(two_dim_Bundle[0][int(x_j)])
	print("\n")
	return(two_dim_Bundle)

#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#	

def populate_bundle_ends_exp_vb(two_dim_Bundle, x_max, y_MT, y_PF, x_mean, sig_mt, m):
	x_j = 0
	MaxBND = y_MT*y_PF*2
	Start = y_MT*y_PF + m*y_PF
	End = Start + 2*y_PF
	while x_j < x_max:
		num = random.expovariate(x_mean)
		if float(num) > float(0.0) : x_j = x_j + float(num)		
		if x_j < x_max: 
			for r in range(Start,int(End)):
				two_dim_Bundle[r][int(x_j)] = 1
			print(int(x_j))
			print(two_dim_Bundle[0][int(x_j)])
	print("\n")
	return(two_dim_Bundle)

#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#	
def populate_bundle_ends_exp_vc(two_dim_Bundle, x_max, y_MT, y_PF, x_mean, sig_mt, m):
	x_j = 0
	MaxBND = y_MT*y_PF*2
	Start = y_MT*y_PF + m*y_PF
	End = Start + (4-m)*y_PF
	while x_j < x_max:
		num = random.expovariate(x_mean)
		if float(num) > float(0.0) : x_j = x_j + float(num)		
		if x_j < x_max: 
			for r in range(Start,int(End)):
				two_dim_Bundle[r][int(x_j)] = 1
			print(int(x_j))
			print(two_dim_Bundle[0][int(x_j)])
	print("\n")
	return(two_dim_Bundle)

#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#	

# Here is how this works: (1) at each lattice site determine if tau will be put down by 50% chance
#			  (2) If tau is put down throw a random number to determine which PF it's on
#			  (3) At the end of the bundle, if there is tau left in the pool, repeate (1)(2)
#			   *Note don't put down tau if something already exists at the lattice site

def populate_bundle_tau(two_dim_Bundle, x_max, y_MT, y_PF, x_mean, sig_mt, m, Pool):
	while Pool > 0:
		for xp in range (0,x_max):
			flip = np.random.rand()		#First generate random #
			if flip < 0.5:			# OK now take 50% chance of having a defect down
				r_pf = int(np.random.rand()*10) + y_PF*y_MT
				if two_dim_Bundle[int(r_pf)][int(xp)] < 1 and two_dim_Bundle[int(r_pf)][int(xp)] < y_PF*y_MT and Pool > 0: 
					two_dim_Bundle[int(r_pf)][int(xp)] = 1
					Pool = Pool - 1
					print(Pool)
		if Pool > 0: xp = 0
		

	return(two_dim_Bundle)

#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#	

def populate_bundle_ends_unweighted(two_dim_Bundle, x_max, y_MT, y_PF, x_mean, sig_mt, m):
	x_j = 0
	MaxBND = y_MT*y_PF*2
	Start = y_MT*y_PF + m*y_PF
	End = Start + y_PF
	for k in range(0,x_max):
		num = np.random.rand()
		if float(num) < float(1.0/x_mean): 
			for r in range(Start,int(End)):
				two_dim_Bundle[r][int(k)] = 1
			print(int(k))
			print(two_dim_Bundle[0][int(k)])
	print("\n")
	return(two_dim_Bundle)

#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#

#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#	

def populate_bundle_ends_periodic(two_dim_Bundle, x_max, y_MT, y_PF, x_mean, sig_mt, m):
	x_j = 0
	MaxBND = y_MT*y_PF*2
	Start = y_MT*y_PF + m*y_PF
	End = Start + y_PF
	for k in range(0,x_max):
		if (float(k)/float(x_mean)).is_integer(): 
			for r in range(Start,int(End)):
				two_dim_Bundle[r][int(k)] = 1
			print(int(k))
			print(two_dim_Bundle[0][int(k)])
	print("\n")
	return(two_dim_Bundle)

#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#

def unscaled_2Dbundle(two_dim_Bundle, x_max, y_MT, y_PF):

	Bd = [[float(0.0) for x in range(int(x_max))] for z in range(0,y_MT)]

	for i in range(0,int(y_MT)):							# Check each MT in the bundle seperately
		for j in range(0,x_max):						# Check only the first PF on the MT, note the second portion of the array keeps track of obstructions the first half of the array is a place-holder legacy code
			if two_dim_Bundle[i*int(y_PF)+int(y_MT)*int(y_PF)+1][j] > 0:	# If a MT end exists then we should add it to the new array
				Bd[i][j] = 1
		


	return(Bd)


def scaled_2Dbundle(two_dim_Bundle, x_max, y_MT, y_PF, x_scale):

	Bd = [[float(0.0) for x in range(int(x_max))] for z in range(0,y_MT)]

	for i in range(0,int(y_MT)):							# Check each MT in the bundle seperately
		for j in range(0,x_max):						# Check only the first PF on the MT, note the second portion of the array keeps track of obstructions the first half of the array is a place-holder legacy code
			if two_dim_Bundle[i][j] > 0:	# If a MT end exists then we should add it to the new array
				Bd[i][j] = 1
		


	return(Bd)

