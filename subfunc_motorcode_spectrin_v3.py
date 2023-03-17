import csv
import string
import math
import cmath
import numpy as np	#import package numpy and give it the designation "np"
import os
import random


def four_param_motorv5_3(lattice, x_max, Pd, Pa, Pr, Lambda, Ph, velocity, STDEV_vel, location, Num_Motors, y_MT, MTs, Latt_conv, Time_conv, P_ON, Mod_P):



#---------------------------#
#   Loop of endosomes	    #

	for i in range(0, Num_Motors):
		curr_out = open(location + "/Trap-model_motor-position_v2_1-#"+str(i)+".txt","w")		# Open current output file


	# Reset all values to starting value
		x1 = 50.0			# Initialize motor at initial x-position. This keeps track of x-position
		x_conversion = 0.0 
		t = 0				# Initial simulation start time
		tmp_mod = int(np.random.normal(Mod_P,1))			# Initialize modulo counter
		if tmp_mod == 0:
			tmp_mod = int(np.random.normal(Mod_P,1))
		
		y1 = 0.0			# This keeps track of the MT the motor is on during the simulation
		Start_prob = np.random.rand()   #Pull a random MT start location probability between (0,1)
		kk = 0
		while kk <= y_MT:
			kk = kk +1
			if Start_prob <	kk/y_MT and Start_prob > (kk-1)/y_MT:
				y1 = kk


		memory1 = 0				#Memory of lattice defects
				# This section will determine if the motor moves N-lattice sites forward during this time-step
			#-------------------------------#
			#  Generate Actual Motor Speed  #
			#-------------------------------#
			# Pull from gaussian generator with mean
		Pforward = np.random.normal(velocity,STDEV_vel)
		if Pforward < 0:                        # We currently do not want backward steps
			Pforward = - Pforward
		

		step_store = 0;				# This marker will keep track of how many steps the motor will take
		PT = 0;					# This marker will keep track of the number of time-steps to pause
		P_MOD = 0;				# This Marker keeps track of whether pauses occur after N-bouts of motion
		while x1 < x_max:			# Run the simulation until the end of the lattice
			t += 1
			

			curr_position = x1 			# Make conversion if desired
		#----------------------------------#
		#Get random numbers
			probability = np.random.rand()   #Pull a single random number detachment
			re_prob = np.random.rand()   #Pull a single random number for re-attachment
			rev_prob = np.random.rand()   #Pull a single random number for MT switching
			dir_rev_prob = np.random.rand()   #Pull a single random number for MT switching

			tunnel_prob = np.random.rand()   #Pull a random number for defect tunneling

			
			curr_velocity = Pforward	        # Motor moves forward an integer number of lattice sites based on its speed

		#----------------------------------#
		#    Reflecting Barrier at x=0 	   #
			if curr_position < 0:
				curr_velocity = abs(curr_velocity)



	#Detachment/reatachment choices
		#Detchament			
			if float(probability)<=float(Pd):
		#Reatachment			
				if float(re_prob)>float(Pa): break
#				if float(rev_prob)>=float(Pr): curr_velocity = -curr_velocity
				y1_old = y1
				jump_prob = np.random.rand()   #Pull a random MT start location probability between (0,1)
				kk = 0
				while kk <= y_MT:
					kk = kk + 1
					if jump_prob <	kk/y_MT and jump_prob > (kk-1)/y_MT:
						y1 = kk

					# Reversal check
						if MTs[int(kk - 1)] == 1: curr_velocity = -curr_velocity
						if MTs[int(kk - 1)] == 0: curr_velocity = curr_velocity

					if kk == y_MT : break


				if y1_old != y1:
					memory1 = 0

# After decisions have been made, now determine if motor moves forward or not

# First Check if the motor needs to pause
			if PT > 0:							# Is the motor currently waiting at the same lattice site
				PT = PT - 1							# Decrement the pause time counter
				
				
			if PT == 0:							# Now check if the pause time counter has run down and the motor can move
				if step_store == 0:						# First, check if the marker storing steps is empty
                                        
					step_store = curr_velocity
				if step_store > 1:						# Second, check if the marker storing steps is greater than a single lattice site
					next_pos = int(curr_position) + int(step_store)			# Define the next lattice position, based on current velocity
					step_store = 0							# Decrement the marker storing steps

					
				if step_store < 1:						# Third, check if the marker storing steps is less than a single lattice site
					next_pos = int(curr_position) 					# Define the next lattice position as the current position
					step_store = step_store + curr_velocity				# Increment the step counter until it is greater than unity


# Check for defect at next position (x1) on ANY MT-track and record it
			dsp = 0.0
			for qq in range(0,y_MT):
				dsp = dsp + lattice[int(qq)][int(next_pos)]

# Check for defect at next position (x1) on the current MT-track (y1), and that the motor is not paused; if not then run motor options	
			if lattice[int(y1-1)][int(next_pos)] < 1 and PT == 0: 
				x1 += int(step_store)
				memory1 = 0;

				if Mod_P > 0:							# Check if the motor pauses every N-bout
					P_MOD = P_MOD + 1;						# Increment the counter 

				if P_MOD == int(tmp_mod):					# Is this a bout where a pause time should be drawn
					PT = int(np.random.exponential(P_ON))					# Draw a new pause time to wait after moving forward
					P_MOD = 0						# Reset the pause counter
					tmp_mod = int(np.random.normal(Mod_P,1))		# Draw a new travel distance
					if tmp_mod == 0:
						tmp_mod = int(np.random.normal(Mod_P,1))


# Check for defect at next position, if yes and motor is blocked, then store in memory and wait a time-step
			if lattice[int(y1-1)][int(next_pos)] == 1: 
				memory1 = memory1 + 1			

			curr_position = x1
			curr_out.write("%s\t%s\t%s\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n" % (t, x1, y1, curr_position, float(Pd), float(probability),float(Pa),float(re_prob),float(Pr),float(rev_prob), float(lattice[int(y1-1)][int(next_pos)]), float(dsp), PT, Pforward, int(step_store), step_store, P_MOD, tmp_mod ) )
			if next_pos >= x_max : break
			if t > 10000 : break

	return ()

#--------------------------------------------------------------#
#--------------------------------------------------------------#
#--------------------------------------------------------------#









