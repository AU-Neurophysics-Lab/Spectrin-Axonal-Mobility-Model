
#########################################################
#	These functions perform post-simulation		#
#	analyses					#
#########################################################

import csv
import sys
import string
import math
import cmath
import numpy as np
import os

#########################################################
#	This function coarse-grains the simulated	#
#	motors and outputs the position versus time	#
#	version 2 accounts for motor reversals		#
#	version 3 removes the extra time due to 	#
#	reversal					#
#########################################################
#########################################################
# Version 3 accounts for multiple lattice site steps    #
# per time-step. This will average the coarse graining  #
# differently.                                          #
#########################################################
#########################################################
# Version 4 accounts for multiple lattice site steps    #
# per time-step. This will average the coarse graining  #
# differently.                                          #
#########################################################
#########################################################
# Version 4 accounts outputs number of runs with data   #
#########################################################

def coarse_grain_v4(dir, file, dX, Ncols, y_MT, res):

	motor_in = open(dir + file, "r")		# Open current motor file
	num_lines = 0
	#readline = motor_in.readline()   		# Use if there is a header
	while motor_in.readline():
		num_lines += 1
#	print(num_lines)

    ### Now run the post-simulation coarse graining ###
	tmp_ct = 0.0					# This will keep track of runs
	a = 0.0						# This will keep track of the current coarse-grained time-step
	lucg = 1.0					# This keeps track of the last un-coarse-grained lattice site
	tlucg = 1.0					# This keeps track of the last un-coarse-grained time
	bund_tran_grained = [[float(0.0) for x in range(0,int(num_lines))] for z in range(0, Ncols)]
	Ph_grained = 0.0				# This keeps track of the number of ends within the coarse-grained lattice-site
	PT1 = 0.0					# This keeps track of how long a motor has paused at lattice site
	PT2 = 0.0					# This keeps track of Pausing Memory
	PT3 = 0.0					# This keeps track of Pausing Memory

	motor_in = open(dir + file, "r")		# Open current motor file
	#readline = motor_in.readline()   		# Use if there is a header
	while True:
		check = motor_in.readline()
		if len(check) > 0:
			cols = check.split()		# Split up the columns so they can be used independantly

	#-------- Add any MT-ends to the Ph count, Note, there is an extra trap condition because of double counting obstructions after waiting at an obstruction
			if int(float(cols[11])) > 0.0 and int(float(cols[1])) != PT1 and int(float(cols[1]))-PT1  > 1.0:
				Ph_grained = Ph_grained + float(cols[11])/y_MT
				PT1 = int(float(cols[1])) 		# Mark this location down as already counted
				

	#-------------------------------------------------------------#
			if int(float(cols[0])) == (dX + dX*a):


				bund_tran_grained[0][int(a)] = float(cols[0]) - float(tlucg)	# Add the time in lattice site to the output
				tlucg = float(cols[0])						# Mark this time down as the last time counted

				bund_tran_grained[3][int(a)] = float(Ph_grained)		# Add the Number of MT-ends to the output
				Ph_grained = 0.0						# Reset the Obstruction Probability count

				if int(float(cols[17])) != tmp_ct:
					bund_tran_grained[2][int(a)] = int(float(cols[17])) # Add the distance to the output
					tmp_ct = int(float(cols[17]))

				bund_tran_grained[1][int(a)] = int(float(cols[1])/res)		# Add the distance to the output
				a = a + 1							# Increment the Obstruction
				

				

#				print(a)

		if not check: break



	return (bund_tran_grained)

####################################################

