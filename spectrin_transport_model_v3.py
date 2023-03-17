#######################################################
###    This motility code is adapted from a model   ###
###    of SVP motility on an axonal bundle.         ###
###    MW Gramlich, S Balseiro-Gómez, SMA Tabei,    ###
###    M Parkes, S Yogev; Scientific reports 11 (1) ###
###    , 1-11                                       ###
###    This code also includes random bundle        ###
###    orientation considered for in vitro bundles  ###
#######################################################
###     This main routine is used to generate       ###
###     simulated transport of spectrin on an axon  ###
###     bundled lattice.                            ###
#######################################################
###     Written By: M.W. Gramlich                   ###
###     Department of Physics Auburn University     ###
###     Last Date Modified: 03/2023                 ###
#######################################################

#######################################################
###    This routine operates in the following order ###
###    [1] Initiate all subroutines                 ###
###    [2] Initiate variables and constraints       ###
###    [3] Initialize the axonal bundle array       ###
###    [4] Initialize the directory used to output  ###
###     all files. This directory will be named     ###
###     based on the values of the variables.       ###
###    [5] Initiate the spectrin simulations        ###
###    [6] The code then Coarse-Grains the 	    ###
###	   simulated motility to reproduce	    ###
###	   experimental results.		    ###
#######################################################

import csv
import string
import math
import cmath
import numpy as np	#import package numpy and give it the designation "np"
import os
import random

from subfuncs_os import check_make_dir
from subfunc_motorcode_spectrin_v3 import four_param_motorv5_3
print ("four_param_motorv2_1 imported successfully")
from subfuncs_bundlecode_v2 import generate_bundle_v2_1
from subfuncs_bundlecode_v2 import populate_bundle_ends_periodic
from subfuncs_bundlecode_v2 import populate_bundle_ends
from subfuncs_bundlecode_v2 import unscaled_2Dbundle
from subfuncs_bundlecode_v2 import scaled_2Dbundle
from subfuncs_postsimanalyses_v1 import coarse_grain_v4

# Note: Define the version number of the model here
ver_num = "v6_1"


#---------------------#
#  Initial Variables  #


detachment_prob = 0.00016
reatachment_prob = 0.0
reversal_prob = 1.01
motor_stop = 399
N_motors = 500
Latt_conv = 4           # This is the nanometers per lattice site
Time_conv = 30          # This is the msec per time-step
velocity = 5.25         # This number should be given in lattice sites per time-step
STDEV_vel = 1.5       	# This number should be given in lattice sites per time-step
Resolution = 25         # This is the resolution in number of lattice-sites/pixel
P_ON = 310;             # This is a parameter for the number of time-steps a motor will pause
Mod_P = 3;		# This parameter allows for consecutive runs between pausing
#--------------------#





dir = os.getcwd() + "\\2D-motor simulation data\\obstructions-test\\"

#---------------------------------------------------------------------------------------------------#
#	                          Initialize bundle lattice                                         #
#                       Note: This section only matters if MT-ends are being considered             #
#                             explicitly, rather than a mean-field approach.                        #

Bundle_size = 12504
Lambda = 10.0	# Define periodicity of MT-ends if used
mean_Ph = 0.0
N_MT = 5        # Define number of MTs in the bundle
N_PF = 13       # Define the number of protofilaments per microtubule
re_scale = 3	# Minimum re-scale of bundle in lattice units


motor_stop = Bundle_size/re_scale # Define when the motor

# Set how many times the simulations will run for different bundles for the same endosome motility conditions
for i in range(1,2):	
	Bundle = generate_bundle_v2_1("tmp_bundleb", Bundle_size, N_MT, N_PF)

	# Skip generating MT-ends for the spectrin model
	#for j in range(0,N_MT): Bundle = populate_bundle_ends(Bundle, Bundle_size, N_MT, N_PF, 450, 190, j)


#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#

	Bundle_Sum = unscaled_2Dbundle(Bundle, Bundle_size, N_MT, N_PF)
	Bundle_Sum_scale = scaled_2Dbundle(Bundle_Sum, Bundle_size, N_MT, N_PF, re_scale)
	
        # This section considers the possiblity of randomly oriented microtobules.
	#MTd = [0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1]     # Use for anti-parallel microtubules
	MTd = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]

	
#---------------------------------------------------------------------------------------------------#
#---------------------------------------------------------------------------------------------------#


	# Set the directory used to store all output files for later analysis
	new_dir = dir + "\\NMT=" + str(N_MT) + "\\Pd" + str(format(detachment_prob, ".2f")) + "_Pa" + str(format(reatachment_prob, ".2f"))+ "_Pr" + str(format(reversal_prob, ".2f")) + "_lambda" + str(format(Lambda, ".2f")) + "_Ph" + str(format(mean_Ph, ".2f"))+ "_PON" + str(format(P_ON, ".2f")) + "-" + str(i) + "\\"
	check_make_dir(new_dir)

	# Output the bundle geometry for reference

	bundle_out = open(new_dir + "trap-" + ver_num + "-bundle-geometry_nonobstruction.bnd","w")		# Open bundle output file
	bundle_out.write("x\t")
	for y in range(0, int(N_MT*N_PF*2)):
		if y < int(N_MT*N_PF): bundle_out.write("PF\t" )
		if y >= int(N_MT*N_PF): bundle_out.write("Obs\t" )
		if y == int(N_MT*N_PF*2): bundle_out.write("\n")
	bundle_out.write("\n")
			#--------------------##--------------------##--------------------#
			#--------------------##--------------------##--------------------#
	for x in range(0, Bundle_size):
		bundle_out.write("%f\t" % (float(x)) )

		for y in range(0, N_MT*N_PF*2):
			bundle_out.write("%f\t" % (Bundle[int(y)][int(x)]) )
		bundle_out.write("\n")
			#--------------------##--------------------##--------------------#
			#--------------------##--------------------##--------------------#
	sum_bundle_out = open(new_dir + "trap-" + ver_num + "-bundle-geometry_summed.bnd","w")
	for x in range(0, Bundle_size):
		sum_bundle_out.write("%f\t%f\t%f\t%f\t%f\t%f\n" % (x, Bundle_Sum[0][int(x)], Bundle_Sum[1][int(x)], Bundle_Sum[2][int(x)], Bundle_Sum[3][int(x)], Bundle_Sum[4][int(x)]) )

			#--------------------##--------------------##--------------------#
			#--------------------##--------------------##--------------------#
	sum_scale_bundle_out = open(new_dir + "trap-" + ver_num + "-bundle-geometry_summed_scaled.bnd","w")
	for x in range(0, int(Bundle_size/re_scale)):
		sum_scale_bundle_out.write("%f\t%f\t%f\t%f\t%f\t%f\n" % (x, Bundle_Sum_scale[0][int(x)], Bundle_Sum_scale[1][int(x)], Bundle_Sum_scale[2][int(x)], Bundle_Sum_scale[3][int(x)], Bundle_Sum_scale[4][int(x)]) )

#---------------------------------------------------------------------------------------------------#
#		      Run the simulated endosome motility routine             			    #

	r = four_param_motorv5_3(Bundle_Sum_scale, motor_stop, detachment_prob, reatachment_prob, reversal_prob, Lambda, mean_Ph, velocity, STDEV_vel, new_dir, N_motors, N_MT, MTd, Latt_conv, Time_conv, P_ON, Mod_P)
	
#---------------------------------------------------------------------------------------------------#
#		      Coarse Grain the simulated endosome motility             			    #

	# This file outputs the Speed and Runlength of each endosome
	VR = open(new_dir + "Speed_Runlength_v2.txt","w")
	VR.write("#Time\tDistance\n")


	for qqq in [10]:

		# Make the sub-directory for the Coarse-Graining
		CG_dir = new_dir + str(qqq) + "\\"
		check_make_dir(CG_dir)

		MSL = open(new_dir + "binned_PT" + str(qqq) + "_v2.txt","w")
		MSL.write("#Time\tPh\n")



		for i in range(0,N_motors):

			t = coarse_grain_v4(new_dir, "Trap-model_motor-position_v2_1-#" + str(i) + ".txt", qqq, 5, N_MT, Resolution)			# Now run the coarse graining

			#-----------------   Output only the Pause-times for binning and analysis	---------------#
			for lll in range(1,len(t[0])):
				if float(t[2][int(lll)]) > 0.0 and float(t[2][int(lll)]) != float(t[2][int(lll-1)]):
					MSL.write("%f\t%f\t%f\n" % (float(t[0][int(lll)]), float(t[1][int(lll)]), float(t[2][int(lll)]) ) )	
				if float(t[2][int(lll)]) > 0.0 and float(t[2][int(lll)]) == float(t[2][int(lll-1)]):
					MSL.write("%f\t%f\t%f\n" % (float(t[0][int(lll)]), float(t[1][int(lll)]), 0 ) )	


			#-----------------   Output all data for this particular motor simulation	---------------#

			coarsegrained_bundle_out = open(CG_dir + "trap-coarsegrained-bund_" + str(i) + "-v_" + ver_num + "_v1.txt","w")
			chkz = 0				# Counter for end of simulation
			chktmp = 0
			count_tmp = 1
			pause_tmp = 0				# Counter for how long the coarse-grained spectrin has paused
			pause_num = 0				# Counter for how many pauses occured
			pause_last_pos = t[2][int(1)]		# Last pause position
			pause_last_time = int(0)		# Last pause time
			for x in range(1, len(t[0])):
				if t[2][int(x)] != pause_last_pos and chkz == 0:			# Is the current position different than the last position
					pause_num = pause_num + 1					# Increment the Pause Counter
					pause_tmp = pause_tmp + (int(x)-pause_last_time)/(Mod_P-1)	# Add the current pause time to the current average; Note: If multiple runs occur between pauses, then this undercounts the pause time by that amount 
					pause_last_pos = t[2][int(x)]					# Make the current position the new "Last pause position"
					pause_last_time = int(x)					# Make the current time the new "Last pause time"
					
				if t[2][int(x)] != 0:
					chktmp = chktmp + float(t[2][int(x)])
					count_tmp = count_tmp + 1
				if t[0][int(x)] == 0 and chkz == 0 and int(x-1)>0:
					tt = float(float(chktmp)/float(count_tmp))			# Calculate Total Travel Length
					AP = 0
					if float(pause_num) > 0:
						AP = float(float(pause_tmp)/float(pause_num))		# Calculate Average Pause Time
						
					VR.write("%f\t%f\t%f\t%f\t%f\n" % (int(x-1), t[1][int(x-1)], float(t[1][int(x-1)])/int(x-1)*0.33, tt,AP))
					chkz = 1
				coarsegrained_bundle_out.write("%f\t%f\t%f\n" % (t[0][int(x)], t[1][int(x)], t[2][int(x)]) )
VR.close()	

input("Press Enter to Exit.")
