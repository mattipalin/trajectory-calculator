# Descent trajectory calculator

import os.path
import sys
import math

r_Mars 			= 3386000.0				# Radius of Mars (in m)
m_Mars 			= 6.41693e23				# Mass of Mars (in kg)
#coriolis		= 1.0					# The Coriolis effect; -1 = going against the rotation of Mars; 0 = no effect; +1 = landing in the direction of rotation
deltat_init		= 0.01					# Time step in integration
wrtint 			= 20					# Writing interval (number of time steps)
write_to_file		= True					# Write to file (True/False)
silent			= False					# Silent operation (True/False)
gammasearch		= False					# Normal operation (False) or Gamma Search (True) 

# Input phase
os.system('clear')

#------------------ Mission parameters ------------------------------------
m 			= 22.0					# Vehicle mass, in kg			
velocity_init 		= 4000.0				# Initial velocity, m/s
gamma_init 		= 17*3.14159/180.0			# Initial trajectory angle (gamma), in rad
altitude_init 		= 120000.0				# Initial altitude, in km
C_D 			= 0.7					# Drag coefficient
C_L 			= 0.0					# Lift coefficient (for ex. 0)
A 			= 0.7864				# Reference area
beta 			= m/(C_D*A)*3.71			# Ballistic coefficient, related to the weight (=m*g) of the vehicle


# ------------------ Scalar functions -------------------------------------

def g(altitude): 			return 6.674e-11*M_Mars/((R_Mars + altitude)**2)			# Acceleration of gravity
def temperature(altitude): 		return 226.0-4.5*math.log(1.0+altitude/1000)**2				# Temperature on Mars
def viscosity(temperature):		return 1.92e-7*(temperature**(1.8))/(temperature+60)			# Air dynamic viscosity on Mars
def dynpres(altitude,velocity):		return 0.5*rho(altitude)*velocity**2					# Dynamic pressure
def speed_of_sound(altitude):		return 236.4*math.exp(-0.002470*altitude/1000.0)			# Speed of sound on Mars
def density(altitude): 												# Medium density
	if altitude<46000:		return math.exp(-4.113-0.095*altitude/1000)
	else:				return math.exp(-2.320-0.134*altitude/1000)

def update_deltat(gamma): 
	global deltat_init
	try:				return min(deltat_init*max(1.0/math.cos(math.sqrt(gamma**2)+5.0*3.14159/180.0)**6.0,1.0),0.25)			# Calculates new time step (to make the calculation faster
	except:				return deltat_init


# --------------- Vectors ------------------------------

def abs(vector):						return math.sqrt(sum(c**2 for c in vector))					# Calculates the absolute value of a vector

def cross(a, b):						return [a[1]*b[2] - a[2]*b[1],a[2]*b[0] - a[0]*b[2],a[0]*b[1] - a[1]*b[0]]	# Calculates the cross product
def calculate_altitude(xd,yd):	global r_Mars;			return abs([xd,yd])-r_Mars

# ------------------- Derivatives   ------------------------------------------

def dUxdt(xd,yd,Ux,Uy):
	global C_D, m, m_Mars, A
	altitude = calculate_altitude(xd,yd)
	return -C_D/m*0.5*density(altitude)*A*math.sqrt(Ux**2 + Uy**2 )*Ux  - xd*6.674e-11*m_Mars/((xd**2+yd**2)**(3.0/2.0))

def dUydt(xd,yd,Ux,Uy):
	global C_D, m, m_Mars, A
	altitude = calculate_altitude(xd,yd)
	return -C_D/m*0.5*density(altitude)*A*math.sqrt(Ux**2 + Uy**2 )*Uy  - yd*6.674e-11*m_Mars/((xd**2+yd**2)**(3.0/2.0))


# ----------------- Functions for printing to file --------------------------

def writeline(velocity_vector,gamma,altitude,time,flight_range, x, y):				# Defines how each line is written into the result file. 
	global A

	velocity=abs(velocity_vector)
	Ma	= velocity/speed_of_sound(altitude)
	Re	= density(altitude)*math.sqrt(A*4.0/3.14159)*abs(velocity_vector)/viscosity(temperature(altitude))
	Kn	= Ma/Re*math.sqrt(1.35*3.141592/2)

	entry00 = time									# Time (seconds)
	entry01 = velocity								# Velocity (m/s)
	entry02 = gamma*180.0/3.14159							# Gamma (degrees)
	entry03 = altitude								# Altitude (m)
	entry04 = flight_range								# Flight range
	entry05 = Ma									# Mach number
	entry06 = Re									# Reynolds number
	entry07 = Kn									# Knudsen number
	entry08 = density(altitude)							# Density
	entry09	= (r_Mars + 0		)*math.cos(flight_range/r_Mars)			# X-coordinate of Mars surface in Cartesian
	entry10 = (r_Mars + 0		)*math.sin(flight_range/r_Mars)			# Y-coordinate of Mars surface in Cartesian
	entry11 = x									# X-coordinate in Cartesian
	entry12 = y									# Y-coordinate in Cartesian
	return '{0:.4e}      {1:.5e}        {2:.5e}        {3:.5e}      {4:.5e}        {5:.5e}       {6:.5e}       {7:.5e}        {8:.5e}        {9:.5e}           {10:.5e}        {11:.5e}        {12:.5e}'.format(entry00,entry01,entry02,entry03,entry04,entry05,entry06,entry07,entry08,entry09,entry10,entry11,entry12)


# ---------------   Main program ---------------------------------------------

def fly():
	global altitude_init, silent, write_to_file, r_Mars
	# Phase 1: Initial values
	time			= 0
	r 			= [	r_Mars + altitude_init			,	0					]
	U 			= [	-velocity_init*math.sin(gamma_init)	,	velocity_init*math.cos(gamma_init)	]
	altitude		= altitude_init
	rho			= density(altitude_init)
	deltat 			= deltat_init
	theta			= 0.0
	flight_range		= 0.0

	iteration_number	= 0

	# Phase 2: Prepare output file

	if not(silent):
		print "This is the Cartesian Descent Trajectory Calculator (CADETRAC) - created by Matti Palin. 2016."
		print " ------------------------------------ "
		if write_to_file == True: 	print "The results will be written to a file."
		else:				print "The results are not written to a file."

	f = open("Results.txt",'w')

	if write_to_file:
		f.write("This is the result file of TrajectoryCalculation.py \n")
		f.write("---------------------------------------------------------------- \n")
		f.write("Initial values: \n")
		f.write("Vehicle mass:               %5.3f" %m)						; f.write(" kg \n")
		f.write("Initial velocity:           %5.3f" %velocity_init)				; f.write(" m/s \n")
		f.write("Initial flight path angle:  %5.3f" %(float(gamma_init)*180.0/3.14159))		; f.write(" degrees \n")
		f.write("Initial altitude:           %5.3f" %altitude_init)				; f.write(" m \n")
		f.write("Drag coefficient:           %5.3f" %C_D)					; f.write(" \n")
		f.write("Lift coefficient:           %5.3f" %C_L)					; f.write(" \n")
		f.write("Reference area A:           %5.3f" %A)						; f.write(" m^2 \n")
		f.write("Time step:                  %5.3e" %deltat)					; f.write(" s \n")
		f.write("---------------------------------------------------------------- \n")
		f.write("Time             Velocity          Gamma              Altitude          r                 Ma               Re               Kn                   Rho                 SCart_X             SCart_Y              Cart_X          Cart_Y")
		f.write("\n")

	# Phase 3: Calculation loop
	while altitude > 0:
		Ux		= U[0]
		Uy		= U[1]
		theta		= math.atan(r[1]/r[0])
		flight_range	= r_Mars*theta
		gamma		= math.atan(-Ux/Uy) - theta
		altitude	= calculate_altitude(r[0],r[1])

		#print "time		=", time
		#print "Delta t		=", deltat
		#print "altitude         =", altitude
		#print "rho		=", rho
		#print "x		=", r[0]
		#print "y		=", r[1]
		#print "Ux		=", U[0]
		#print "Uy		=", U[1]
		#print "dUxdt            =", dUxdt(r[0],r[1],Ux,Uy)
		#print "dUydt            =", dUydt(r[0],r[1],Ux,Uy)
		#print "gamma		=", gamma*180.0/3.141592
		#print "theta		=", theta*180.0/3.141592
		#raw_input()

		if altitude > altitude_init: return False
		r		= [	r[0] + U[0]                     *deltat		, r[1] + U[1]                   *deltat		]
		U		= [	U[0] + dUxdt(r[0],r[1],Ux,Uy)   *deltat		, U[1] + dUydt(r[0],r[1],Ux,Uy) *deltat		]
		time		+= deltat




		if (write_to_file and iteration_number % wrtint == 0):
			f.write(writeline(U,gamma,altitude,time,flight_range, r[0], r[1]))
			f.write("\n")
		iteration_number += 1
		time		+= deltat
		deltat		= update_deltat(gamma)

	# Resume of results
	if not(silent):
		print "Total steps: ", iteration_number
		print "I recommend changing parameter wrtint to approx. ", int(iteration_number/150.0)
		print "Calculation successful!"
		print "Results have been written into Results.txt. Have a good day!"
		print 'Distance travelled: {0:.2f} km'.format(flight_range/1000)
		print 'Final trajectory angle: {0:.1f} degrees'.format(gamma*180.0/3.14159)
	f.close()
	return True



print fly()