# Descent trajectory calculator

import os.path
import sys
import math
# The needed constants
G = 6.67e-11

R_Earth = 6367444.0
R_Mars = 3386000.0
M_Earth = 5.97212e24
M_Mars = 6.41693e23
deltat = 0.0002
deltah = -0.2

wrtint = 200

# Input phase
os.system('clear')
def user_inputs():
	print "Please provide the vehicle weight."
	W = raw_input() ; W = float(W)
	print "Please give the initial velocity (in m/s)"
	velocity_init = raw_input(); velocity_init = float(velocity_init)
	print "Please provide the initial vehicle flight path in degrees. (0 = horizontal to surface)"
	gamma_init = raw_input()
	gamma_init = float(gamma_init)*3.14159/180.0 # Conversion to radians
	print "Please provide the initial altitude in meters."
	altitude_init = raw_input(); altitude_init = float(altitude_init)
	print "Please give the drag coefficient C_D of the vehicle."
	C_D = raw_input() ; C_D = float(C_D)
	print "Please give the lift coefficient C_L of the vehicle."
	C_L = raw_input() ; C_L = float(C_L)
	print "Please give the vehicle reference area A."
	A = raw_input() ; A = float(A)
	beta = W/(C_D*A) # Ballistic coefficient
	print "In which planet do you wish to perform the calculation? (e= Earth, m = Mars)"
	planet = ""
	while planet != "e" and planet != "m":
		planet = raw_input()
		if planet == "e":
			R = R_Earth
		elif planet == "m":
			R = R_Mars
		else:
			print "Could not understand input. Please write again."
	print "Do you wish to perform the calculation as a function of time or altitude? (t/a)"
	deltamethod = raw_input()
	while deltamethod != "t" and deltamethod != "a":
		print "Could not understand input. Please write again."
		deltamethod = raw_input()
	return W, velocity_init, gamma_init, altitude_init, C_D, C_L, A, beta, planet, R, deltamethod
def defaults():
	W = 232.52
	velocity_init = 5000.0
	gamma_init = 77.5*3.14159/180.0
	altitude_init = 120000.0
	C_D = 0.7
	C_L = 0.0
	A = 0.7864
	planet = "m" ; R = R_Mars
	beta = W/(C_D*A)*3.71
	deltamethod = "t" # t or a
	return W, velocity_init, gamma_init, altitude_init, C_D, C_L, A, beta, planet, R, deltamethod
# Scalar functions

def g(planet,altitude): # Calculates the acceleration of gravitation on a given altitude
	if planet == "e":
		return G*M_Earth/((R_Earth + altitude)**2)
	elif planet == "m":
		return G*M_Mars/((R_Mars + altitude)**2)
	else:
		print "An error occured with the calculation of acceleration of gravitation."

def Temp_Earth(altitude):
	T = 228.0
	if altitude<0.:
		print "Negative altitude! Error ..."
		quit()
	if altitude<10000.:
		T = 288.15 -0.0065*altitude
	elif altitude>10000.:
		if altitude<20000.:
			T = 216.65
		else:
			T = 196.65 + 0.001*altitude
	return T

def Temp_Mars(altitude):
	return 226.0-4.5*math.log(1.0+altitude/1000)**2

def viscosity_Mars(temperature):
	return 1.92e-7*(temperature**(1.8))/(temperature+60)

def rho(planet,altitude): # Calculates the medium density at a given altitude
	if planet == "e":
		if altitude <10.0:
			return 1.22
		elif altitude>10 and altitude<17000:
			return 101325.0*0.0289644/(8.31447*Temp_Earth(altitude))*(1.0 - altitude*0.0065/288.15)**(g(planet,altitude)*0.028644/(8.31447*0.0065))
		else:
			return 0.08
	elif planet == "m":
		#return math.exp(-0.1153*altitude/1000-3.9392)
		if altitude<46000:
			return math.exp(-4.113-0.095*altitude/1000)
		else:
			return math.exp(-2.320-0.134*altitude/1000)
def Q(planet,altitude,velocity): # Calculates the dynamic pressure
	return 0.5*rho(planet,altitude)*velocity**2

def speed_of_sound(altitude):
	return 236.4*math.exp(-0.00247*altitude/1000)

# Derivatives
def dVdh(planet,velocity,gamma,altitude):
	return g(planet,altitude)*((Q(planet,altitude,velocity)/beta) - math.sin(gamma))/(velocity*math.sin(gamma))
def dgammadh(planet,velocity,gamma,altitude):
	return (g(planet,altitude)*(Q(planet,altitude,velocity)/beta)*C_L/C_D + math.cos(gamma)*(-g(planet,altitude) + velocity**2/(R+altitude)))/(velocity**2 * math.sin(gamma))
def dtdh(planet,velocity,gamma,altitude):
	if gamma == 0:
		return 100
	if abs(-1.0/(velocity*math.sin(gamma))) > 1000:
		return (abs(-1.0/(velocity*math.sin(gamma))))/abs(-1.0/(velocity*math.sin(gamma)))*1000
	return -1.0/(velocity*math.sin(gamma))
def drdh(planet,velocity,gamma,altitude):
	return -R*math.cos(gamma)/((R+altitude)*math.sin(gamma))

def dVdt(planet,velocity,gamma,altitude):
	return g(planet,altitude)*(-Q(planet,altitude,velocity)/beta + math.sin(gamma))
def dgammadt(planet,velocity,gamma,altitude):
	global W
	return (-Q(planet,altitude,velocity)*g(planet,altitude)/beta*C_L/C_D + math.cos(gamma)*(g(planet,altitude)-velocity**2/(R+altitude)))/velocity -2*W*7.087919e-5
def dhdt(planet,velocity,gamma,altitude):
	return -velocity*math.sin(gamma)
def drdt(planet,velocity,gamma,altitude):
	return R*velocity*math.cos(gamma)/(R + altitude)

def iterate(deltamethod,planet,velocity,gamma,altitude,time,r):
	if altitude > 120000:
		sys.exit("Overshoot")
	if deltamethod == "a":
		velocity = velocity +     dVdh(planet,velocity,gamma,altitude)*deltah
		gamma    = gamma    + dgammadh(planet,velocity,gamma,altitude)*deltah
		time     = time     +     dtdh(planet,velocity,gamma,altitude)*deltah
		r        = r        +     drdh(planet,velocity,gamma,altitude)*deltah
		return velocity,gamma,time,r
	elif deltamethod == "t":
		velocity = velocity +     dVdt(planet,velocity,gamma,altitude)*deltat
		gamma    = gamma    + dgammadt(planet,velocity,gamma,altitude)*deltat
		altitude = altitude +     dhdt(planet,velocity,gamma,altitude)*deltat
		r        = r        +     drdt(planet,velocity,gamma,altitude)*deltat
		return velocity,gamma,altitude,r

def printvalues(deltamethod,velocity,gamma,altitude,time,r):
	if deltamethod == "a":
		f.write('{0:.4f}      {1:.5e}        {2:.5e}        {3:.5e}       {4:.5e}        {5:.5e}       {6:.5e}       {7:.5e}'.format(altitude,velocity,gamma*180.0/3.14159,time,r,velocity/speed_of_sound(altitude),rho(planet,altitude)*math.sqrt(A*4.0/3.14159)*velocity/viscosity_Mars(Temp_Mars(altitude)),(velocity/speed_of_sound(altitude))/(rho(planet,altitude)*math.sqrt(A*4.0/3.14159)*velocity/viscosity_Mars(Temp_Mars(altitude)))*math.sqrt(1.35*3.14159/2)))
	else:
		f.write('{0:.4f}      {1:.5e}        {2:.5e}        {3:.5e}      {4:.5e}        {5:.5e}       {6:.5e}       {7:.5e}'.format(time,velocity,gamma*180.0/3.14159,altitude,r,velocity/speed_of_sound(altitude),rho(planet,altitude)*math.sqrt(A*4.0/3.14159)*velocity/viscosity_Mars(Temp_Mars(altitude)),(velocity/speed_of_sound(altitude))/(rho(planet,altitude)*math.sqrt(A*4.0/3.14159)*velocity/viscosity_Mars(Temp_Mars(altitude)))*math.sqrt(1.35*3.14159/2)))
	f.write("\n")

# Main program starts
print "This is the Descent Trajectory Calculator - created by Matti Palin. 2015."
t_init = 0.0
r = 0.0
W, velocity_init, gamma_init, altitude_init, C_D, C_L, A, beta, planet, R, deltamethod = defaults()
print "Thank you. The calculation will now commence." 
print " ------------------------------------ "
# Beginning of calculation
time = t_init
velocity = velocity_init
gamma = gamma_init
altitude = altitude_init

#First lines of output file
f = open("Results.txt",'w')
f.write("This is the result file of TrajectoryCalculation.py \n")
f.write("---------------------------------------------------------------- \n")
if planet == "e":
	f.write("Simulation on Earth.\n")
elif planet == "m":
	f.write("Simulation on Mars.\n")
f.write("Initial values: \n")
f.write("Vehicle mass:               %5.3f" %W)
f.write(" kg \n")
f.write("Initial velocity:           %5.3f" %velocity_init)
f.write(" m/s \n")
f.write("Initial flight path angle:  %5.3f" %(float(gamma_init)*180.0/3.14159))
f.write(" degrees \n")
f.write("Initial altitude:           %5.3f" %altitude_init)
f.write(" m \n")
f.write("Drag coefficient:           %5.3f" %C_D)
f.write(" \n")
f.write("Lift coefficient:           %5.3f" %C_L)
f.write(" \n")
f.write("Reference area A:           %5.3f" %A)
f.write(" m^2 \n")
if deltamethod == "a":
	f.write("Altitude step:              %5.3f" %abs(deltah))
	f.write(" m \n")
elif deltamethod == "t":
	f.write("Time step:                  %5.3f" %deltat)
	f.write(" s \n")
f.write("---------------------------------------------------------------- \n")
if deltamethod == "a":
	f.write("Altitude        Velocity          Gamma          Time            r         Mach number      Reynolds number      Knudsen number")
	f.write("\n")
else:
	f.write("Time        Velocity          Gamma          Altitude        r         Mach number      Reynolds number      Knudsen number")
	f.write("\n")
# Calculation loop
calc1 = 0
calc2 = 0
calc3 = 0
calc4 = 0
if deltamethod == "a":
	if velocity_init>200.0 and planet == "e":
		deltah = 0.1*deltah
	deltah = deltah*0.1
	while abs(altitude_init - altitude)<10:
		printvalues(deltamethod,velocity,gamma,altitude,time,r)
		velocity,gamma,time,r = iterate(deltamethod,planet,velocity,gamma,altitude,time,r)
		altitude = altitude + deltah
		
	
	deltah = deltah*5
	while abs(altitude_init - altitude)<0.02*altitude_init:
		printvalues(deltamethod,velocity,gamma,altitude,time,r)
		velocity,gamma,time,r = iterate(deltamethod,planet,velocity,gamma,altitude,time,r)
		altitude = altitude + deltah
		
	
	deltah = deltah*2.0
	while abs(altitude_init - altitude)<0.25*altitude_init:
		printvalues(deltamethod,velocity,gamma,altitude,time,r)
		velocity,gamma,time,r = iterate(deltamethod,planet,velocity,gamma,altitude,time,r)
		altitude = altitude + deltah
		

	deltah = deltah*20.0
	while altitude>0:
		printvalues(deltamethod,velocity,gamma,altitude,time,r)
		velocity,gamma,time,r = iterate(deltamethod,planet,velocity,gamma,altitude,time,r)
		altitude = altitude + deltah
		
	

elif deltamethod == "t":
	if velocity_init>200.0 and planet == "e":
		deltat = 0.1*deltat
	deltat = 0.1*deltat
	while abs(altitude_init - altitude)<10:
		beta = W/(C_D*A)*g(planet,altitude)
		if ((calc1+calc2+calc3+calc4) % wrtint == 0):
			printvalues(deltamethod,velocity,gamma,altitude,time,r)
		velocity,gamma,altitude,r = iterate(deltamethod,planet,velocity,gamma,altitude,time,r)
		time = time + deltat
		calc1 = calc1 +1
	deltat = deltat*5.0
	print "Phase 1: ", calc1
	while abs(altitude_init - altitude)<0.02*altitude_init:
		beta = W/(C_D*A)*g(planet,altitude)
		if ((calc1+calc2+calc3+calc4) % wrtint == 0):
			printvalues(deltamethod,velocity,gamma,altitude,time,r)
		velocity,gamma,altitude,r = iterate(deltamethod,planet,velocity,gamma,altitude,time,r)
		time = time + deltat
		calc2 = calc2 +1
	print "Phase 2: ", calc2
	deltat = deltat*2.0
	while abs(altitude_init - altitude)<0.25*altitude_init:
		beta = W/(C_D*A)*g(planet,altitude)
		if ((calc1+calc2+calc3+calc4) % wrtint == 0):
			printvalues(deltamethod,velocity,gamma,altitude,time,r)
		velocity,gamma,altitude,r = iterate(deltamethod,planet,velocity,gamma,altitude,time,r)
		time = time + deltat
		calc3 = calc3 +1
	print "Phase 3: ", calc3
	deltat = deltat*5.0
	while altitude>0:
		beta = W/(C_D*A)*g(planet,altitude)
		if ((calc1+calc2+calc3+calc4) % wrtint == 0):
			printvalues(deltamethod,velocity,gamma,altitude,time,r)
		velocity,gamma,altitude,r = iterate(deltamethod,planet,velocity,gamma,altitude,time,r)
		time = time + deltat
		calc4 = calc4 +1
	print "Phase 4: ", calc4
	print "Total steps: ", calc1+calc2+calc3+calc4
	print "I recommend changing parameters WrtInt to approx. ", int((calc1+calc2+calc3+calc4)/150.0)
print "Calculation successful!"
print "Results have been written into Results.txt. Have a good day!"
print 'Distance travelled: {0:.2f} km'.format(r/1000)
print 'Final trajectory angle: {0:.1f} degrees'.format(gamma*180.0/3.14159)
