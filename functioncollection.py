import math
import numpy as np
def calculate_tracking_error_angle(q_ref,q_true):
	# Returns just the angle of the tracking error quaternion.
	# In other words, this should work similarly to calculate_tracking_error_quaternion.
	# It is furtermore assumed that the quaternions are unit quaternions.
	return math.acos(min(sum([q_ref[j]*q_true[j] for j in range(4)]),1))

def calculate_tracking_error_quaternion(q_ref, q_true):
	# Calculates the tracking error quaternion, a measure of the attitude/orientation error.
	# Based on 'A Survey of Attitude Error Representationts', Younes & Turner, 2014.
	# q_ref is the quaternion that is assumed to have some degree of error, and q_true is the true quaternion.
	q_ref_len	= math.sqrt(sum([e**2 for e in q_ref]))
	q_2		= [q_ref[0], -q_ref[1], -q_ref[2], -q_ref[3], ]	# The quaternion inverse of q_ref

	l0		= [q_true[0], -q_true[1], -q_true[2], -q_true[3]]
	l1		= [q_true[1],  q_true[0], -q_true[3],  q_true[2]]
	l2		= [q_true[2],  q_true[3],  q_true[0], -q_true[1]]
	l3		= [q_true[3], -q_true[2],  q_true[1],  q_true[0]]
	out0		= sum(i[0] * i[1] for i in zip(l0, q_2))/q_ref_len
	out1		= sum(i[0] * i[1] for i in zip(l1, q_2))/q_ref_len
	out2		= sum(i[0] * i[1] for i in zip(l2, q_2))/q_ref_len
	out3		= sum(i[0] * i[1] for i in zip(l3, q_2))/q_ref_len
	return [out0, out1, out2, out3]

def calculate_tracking_error_PRA(q_ref, q_true):
	# 2018-05-09. For some reason the two methods above don't work.
	# Therefore, a more intuitive method is proposed. We can use the
	# intuitive Principal Rotation Axis method. We just have to compare
	# two DCMs. The order doesn't matter, due to symmetry of the error.
	M1 = Quaternion_To_TransformationMatrix(q_ref)
	M2 = Quaternion_To_TransformationMatrix(q_true)
	M  = np.dot(np.linalg.inv(M1), M2)
	# Then we must extract the angle
	return math.acos(min(0.5*(np.trace(M) - 1.0), 1.0))
	
def generate_cross_product_matrix(v):
	# Returns the cross product matrix, given a 3-vector v
	return np.matrix([   [0.0, -v[2], v[1]]   ,  [v[2] , 0.0 , -v[0]]   ,   [-v[1] , v[0], 0.0]  ])

def normalize(q):
	# Returns the normalized version of q, i.e with norm = 1.0
	norm = math.sqrt(sum([j**2 for j in q]))
	return [i/norm for i in q]

def quaternion_norm(q):
	# I assume that q is a 4-D vector
	return (sum([j**2 for j in q]))**(0.5)

def interpret_PRA_from_matrix(M):
	# I'll assume M in a normal matrix, a rotation matrix
	trace = sum([M[i,i] for i in range(3)])
	theta = math.acos(0.5*(trace - 1.0))
	return [(M[2,1] - M[1,2])/(2.0*math.sin(theta)), (M[0,2] - M[2,0])/(2.0*math.sin(theta)),(M[1,0] - M[0,1])/(2.0*math.sin(theta))]
	

def Get_EulerTransformation_matrix(phi_d, theta_d, psi_d, direction, rotation_case, intrinsic):
	# Takes the Euler angles and returns the rotation matrix associated with those angles
	S1, C1	= math.sin(  phi_d), math.cos(  phi_d)
	S2, C2	= math.sin(theta_d), math.cos(theta_d)
	S3, C3	= math.sin(  psi_d), math.cos(  psi_d)

	if intrinsic.lower() == 'intrinsic':
		print "Get_EulerTransformation_matrix: intrinsic initiated"
		rotation_case = rotation_case[::-1]
		S1, C1	= math.sin(  psi_d), math.cos(  psi_d)
		S2, C2	= math.sin(theta_d), math.cos(theta_d)
		S3, C3	= math.sin(  phi_d), math.cos(  phi_d)
	if rotation_case== "xyz":
		A11, A12, A13 =  C2*C3			,  C2*S3		,  -S2
		A21, A22, A23 =  C3*S1*S2 - C1*S3	, C1*C3 + S1*S2*S3	, C2*S1
		A31, A32, A33 =  C1*C3*S2 + S1*S3	, C1*S2*S3 - C3*S1	, C1*C2

	elif rotation_case== "yxz":
		A11, A12, A13 =  C1*C3 - S1*S2*S3	, C3*S1*S2 + C1*S3	, -C2*S1
		A21, A22, A23 = -C2*S3 			, C2*C3			, S2
		A31, A32, A33 =  C3*S1 + C1*S2*S3	, S1*S3 - C1*C3*S2	, C1*C2
	elif rotation_case== "xzy":
		A11, A12, A13 =  C2*C3			,  S2			, -C2*S3
		A21, A22, A23 =  S1*S3 - C1*C3*S2	, C1*C2			, C3*S1 + C1*S2*S3
		A31, A32, A33 =  C3*S1*S2 + C1*S3	, -C2*S1		, C1*C3 - S1*S2*S3
	elif rotation_case== "yzx":
		A11, A12, A13 =  C1*C2			, C1*C3*S2 + S1*S3	, C1*S2*S3 - C3*S1
		A21, A22, A23 = -S2			, C2*C3			, C2*S3
		A31, A32, A33 =  C2*S1			, C3*S1*S2 - C1*S3	, C1*C3 + S1*S2*S3
	elif rotation_case== "zyx":
		A11, A12, A13 =  C1*C2			, C3*S1 + C1*S2*S3	, S1*S3 - C1*C3*S2
		A21, A22, A23 = -C2*S1			, C1*C3 - S1*S2*S3	,  C3*S1*S2 + C1*S3
		A31, A32, A33 =  S2			, -C2*S3		, C2*C3
	elif rotation_case== "zxy":
		A11, A12, A13 =  C1*C3 + S1*S2*S3	, C2*S1			, C3*S1*S2 - C1*S3
		A21, A22, A23 =  C1*S2*S3 - C3*S1	, C1*C2			, C1*C3*S2 + S1*S3
		A31, A32, A33 =  C2*S3			, -S2			,  C2*C3
	else: print "ERROR: rotation case not understood!"

	if direction == "reverse":	XN, YN, ZN = [A11 , A12 , A13], [A21 , A22 , A23], [A31 , A32 , A33]
	else:				XN, YN, ZN = [A11 , A21 , A31], [A12 , A22 , A32], [A13 , A23 , A33]

        print "Transformation matrix:"
	print XN
	print YN
	print ZN
	return np.matrix([XN, YN, ZN])

def Quaternion_To_TransformationMatrix(q):
	r0	= q[0]
	r1	= q[1]
	r2	= q[2]
	r3	= q[3]
	out0	= [r0**2 + r1**2 - r2**2 - r3**2  , 2.0*(r1*r2 - r0*r3)	  , 2.0*(r1*r3 + r0*r2)]
	out1	= [2.0*(r1*r2 + r0*r3)  ,   r0**2 - r1**2 + r2**2 - r3**2 , 2.0*(r2*r3 - r0*r1)]
	out2	= [2.0*(r1*r3 - r0*r2)  ,  2.0*(r2*r3 + r0*r1)     , r0**2 - r1**2 - r2**2 + r3**2]
	return np.matrix([out0, out1, out2])
	

def convert_euler_to_quaternion(phi_d, theta_d, psi_d, rotation_case):
	# This was revised on Dec 1 2017

	s1, c1	= math.sin(phi_d  /2.0), math.cos(phi_d  /2.0)
	s2, c2	= math.sin(theta_d/2.0), math.cos(theta_d/2.0)
	s3, c3	= math.sin(psi_d  /2.0), math.cos(psi_d  /2.0)

	if rotation_case == "xyz":
		q0	= c3*c2*c1 + s3*s2*s1
		q1	= c3*c2*s1 - s3*s2*c1
		q2	= s3*c2*s1 + c3*s2*c1
		q3	= s3*c2*c1 - c3*s2*s1
	if rotation_case == "yxz":
		q0	= c3*c2*c1 - s3*s2*s1
		q1	= c3*s2*c1 - s3*c2*s1
		q2	= s3*s2*c1 + c3*c2*s1
		q3	= s3*c2*c1 + c3*s2*s1
	if rotation_case == "xzy":
		q0	= c3*c2*c1 - s3*s2*s1
		q1	= c3*c2*s1 + s3*s2*c1
		q2	= s3*c2*c1 + c3*s2*s1
		q3	=-s3*c2*s1 + c3*s2*c1
	if rotation_case == "yzx":
		q0	= c3*c2*c1 + s3*s2*s1
		q1	= s3*c2*c1 - c3*s2*s1
		q2	= c3*c2*s1 - s3*s2*c1
		q3	= s3*c2*s1 + c3*s2*c1
	if rotation_case == "zyx":
		q0	= c3*c2*c1 - s3*s2*s1
		q1	= s3*c2*c1 + c3*s2*s1
		q2	= c3*s2*c1 - s3*c2*s1
		q3	= s3*s2*c1 + c3*c2*s1
	if rotation_case == "zxy":
		q0	= c3*c2*c1 + s3*s2*s1
		q1	= c3*s2*c1 + s3*c2*s1
		q2	= s3*c2*c1 - c3*s2*s1
		q3	=-s3*s2*c1 + c3*c2*s1

	if q0<0.0:
		return [-q0, -q1, -q2, -q3]
	else:
		return [q0, q1, q2, q3]

def convert_quaternion_to_euler(q, rotation_case):
	# This is a revised version, made Dec 4th 2017 - in order to check
	# if the equations are correct.
	p0	= q[0]
	p1	= q[1]
	p2	= q[2]
	p3	= q[3]
	
	if rotation_case == 'xyz':
		phi	= math.atan2( 2.0*(p0*p1 + p2*p3),(p0**2 - p1**2 - p2**2 + p3**2) )
		theta	= math.asin( max( min( 2.0*(p0*p2 - p1*p3 ) ,1.0 ) ,-1.0 ) )
		psi	= math.atan2( 2.0*(p0*p3 + p1*p2),(p0**2 + p1**2 - p2**2 - p3**2) )
	elif rotation_case == 'yxz':
		phi	= math.atan2( 2.0*(p0*p2 - p1*p3),(p0**2 - p1**2 - p2**2 + p3**2) )
		theta	= math.asin( max( min( 2.0*(p0*p1 + p2*p3 ) ,1.0 ) ,-1.0 ) )
		psi	= math.atan2( 2.0*(p0*p3 - p1*p2),(p0**2 - p1**2 + p2**2 - p3**2) )
	elif rotation_case == 'xzy':
		phi	= math.atan2( 2.0*(p0*p1 - p2*p3),(p0**2 - p1**2 + p2**2 - p3**2) )
		theta	= math.asin( max( min( 2.0*(p0*p3 + p1*p2 ) ,1.0 ) ,-1.0 ) )
		psi	= math.atan2( 2.0*(p0*p2 - p1*p3),(p0**2 + p1**2 - p2**2 - p3**2) )
	elif rotation_case == 'yzx':
		phi	= math.atan2( 2.0*(p0*p2 + p1*p3),(p0**2 + p1**2 - p2**2 - p3**2) )
		theta	= math.asin( max( min( 2.0*(p0*p3 - p1*p2 ) ,1.0 ) ,-1.0 ) )
		psi	= math.atan2( 2.0*(p0*p1 + p2*p3),(p0**2 - p1**2 + p2**2 - p3**2) )
	elif rotation_case == 'zyx':
		phi	= math.atan2( 2.0*(p0*p3 - p1*p2),(p0**2 + p1**2 - p2**2 - p3**2) )
		theta	= math.asin( max( min( 2.0*(p0*p2 + p1*p3 ) ,1.0 ) ,-1.0 ) )
		psi	= math.atan2( 2.0*(p0*p1 - p2*p3),(p0**2 - p1**2 - p2**2 + p3**2) )
	elif rotation_case == 'zxy':
		phi	= math.atan2( 2.0*(p0*p3 + p1*p2),(p0**2 - p1**2 + p2**2 - p3**2) )
		theta	= math.asin( max( min( 2.0*(p0*p1 - p2*p3 ) ,1.0 ) ,-1.0 ) )
		psi	= math.atan2( 2.0*(p0*p2 + p1*p3),(p0**2 - p1**2 - p2**2 + p3**2) )
	# Parse error: no rotation case given
	else: 
		print "Error! convert_quaternion_to_euler: No rotation case given."
		return 0
	return [phi, theta, psi]

def generate_rotation_vector_line(omega):
	n = 50		# Number of points
	tot_len = math.sqrt(sum([j**2 for j in omega]))
	if abs(tot_len) < 1.0e-16: return [], [], []
	x = [0.0 + omega[0]/tot_len * j / n for j in range(n)]
	y = [0.0 + omega[1]/tot_len * j / n for j in range(n)]
	z = [0.0 + omega[2]/tot_len * j / n for j in range(n)]
	return x, y, z

def generate_airplane(scale, pixel_distance, coordinate_system):
	x_points	= []
	y_points	= []
	z_points	= []
	
	t	= -scale	# Fuselage
	while t < 0.0: #0.5*scale:
		x_points.append(t); y_points.append(0.0); z_points.append(0.0); 	t += pixel_distance

	t	= 0.0	# Right wing
	while t < 0.5*scale: x_points.append(t); y_points.append(0.0); z_points.append(t); t += pixel_distance

	t	= 0.0	# Left wing
	while t < 0.5*scale: x_points.append(t); y_points.append(0.0); z_points.append(-t);t += pixel_distance
	
	t = 0.5*scale	# Rudder
	while t < scale: x_points.append(t - 0.5*scale);	y_points.append(t - 0.5*scale);	z_points.append(0.0);t += pixel_distance

	if coordinate_system == "FINFLO":	return x_points, y_points, z_points
	if coordinate_system == "particle":	return [-j for j in x_points], z_points , [-j for j in y_points]
	if coordinate_system == "default": 	return [-j for j in x_points], z_points, y_points
	else:	print "Error: Generate_Airplane -- no coordinate system defined!"

def calculate_updated_quaternion_LeapfrogOLD(omega,q_d, dt):
	# This function should have exactly the same output as calculate_updated_quaternion,
	# but here we're only testing the written out form (Equations 5.39-5.41).
	#print "----- Leapfrog starting ---- "
	p	= 0.25*dt*omega[0]	# ! Tarkista kerroin!
	q	= 0.25*dt*omega[1]
	r	= 0.25*dt*omega[2]

	K0	= p**2 + q**2 + r**2 + 1.0
	K1	= p**2 - q**2 - r**2
	K	= K0*(K1-1.0)

	print "p/q/r = ", [round(el,5) for el in [p, q, r]]
	print "K0    = ", K0
	print "K1    = ", K1
	print "K     = ", K

	r0	= [-p**4 + 2.0*p**2 - 4.0*p*q*r + (q**2 + r**2)**2 -1.0 , -2.0*p*(p**2 - q**2 +r**2 -1.0)			, -2.0*q*(K1 + 1.0) - 4.0*p*q	, 2.0*r*K0]
	r1	= [2.0*p*(p**2 + q**2  - r**2 - 1.0)			, - p**4 +2.0*p**2 + 4.0*p*q*r + (q**2 + r**2)**2 - 1.0	, -2.0*r*(K1 - 1.0) + 4.0*p*q	, -2.0*q*K0]
	r2	= [-2.0*q*K0			, -2.0*r*K0		, -(K1+1.0)*K0	,	2.0*p*K0		]
	r3	= [ 2.0*r*(K1-1.0)+ 4.0*p*q	, -2.0*q*(K1-1) + 4.0*p*r	, 2.0*p*(K1 + 1.0)	, -K0*(K1+1.0)	]

	r0	= [e/K for e in r0]
	r1	= [e/K for e in r1]
	r2	= [e/K for e in r2]
	r3	= [e/K for e in r3]


	M	= np.matrix([r0,r1,r2,r3])
	print M
	qout	= [sum([M[i,j]*q_d[i] for i in range(4)]) for j in range(4)]
	#print "---- Leapfrog ending ------"
	if qout[0]<0: 	return [-e for e in qout]
	else:		return qout

def calculate_updated_quaternion_Omelyan0(omega,q_d, dt):
	# Written 2018-08-08 by Matti Palin; replaces calculate_updated_quaternion_Leapfrog.
	p	= 0.25*dt*omega[0]	# ! Tarkista kerroin!
	q	= 0.25*dt*omega[1]
	r	= 0.25*dt*omega[2]
	L	= p**2 + q**2 +r**2

	l0	= [1.-L, -2.*p, -2.*q, -2.*r]
	l1	= [2.*p, 1. -L, -2.*r,  2.*q]
	l2	= [2.*q,  2.*r, 1.- L, -2.*p]
	l3	= [2.*r, -2.*q,  2.*p, 1.- L]

	r0	= [el/(L + 1.) for el in l0]
	r1	= [el/(L + 1.) for el in l1]
	r2	= [el/(L + 1.) for el in l2]
	r3	= [el/(L + 1.) for el in l3]

	M	= np.matrix([r0,r1,r2,r3])
	qout0	= np.dot(M, np.transpose(np.matrix(q_d)))
	qout  = [qout0.item(i) for i in range(4)]

	if qout[0]<0: 	return [-e for e in qout]
	else:		return qout

def calculate_updated_quaternion_Omelyan1(omega,q_d, dt):
	# Written 2018-08-08 by Matti Palin; calculate_updated_quaternion_Omelyan1.
	Lsquared	 = (dt/4.0)**2 * sum([el**2 for el in omega])
	
	l0	= [0       , -omega[0], -omega[1], -omega[2]]
	l1	= [omega[0],         0, -omega[2],  omega[1]]
	l2	= [omega[1],  omega[2],         0, -omega[0]]
	l3	= [omega[2], -omega[1],  omega[0],         0]
	Q	= 0.5*np.matrix([l0, l1, l2, l3])

	inv	= (1.-Lsquared)/(1.+Lsquared)*np.identity(4) + dt/(1.+Lsquared)*Q
	qout0	= np.dot(inv, np.transpose(np.matrix(q_d)))
	qout  = [qout0.item(i) for i in range(4)]

	if qout[0]<0: 	return [-e for e in qout]
	else:		return qout
	

def calculate_quaternion_time_derivative(q_d,omega,dt):
	# Calculates the quaternion time derivative
	p	= omega[0]
	q	= omega[1]
	r	= omega[2]
	pdot0	= 0.5*(0*q_d[0] - p*q_d[1] - q*q_d[2] - r*q_d[3])
	pdot1	= 0.5*(p*q_d[0] + 0*q_d[1] - r*q_d[2] + q*q_d[3])
	pdot2	= 0.5*(q*q_d[0] + r*q_d[1] + 0*q_d[2] - p*q_d[3])
	pdot3	= 0.5*(r*q_d[0] - q*q_d[1] + p*q_d[2] + 0*q_d[3])
	return np.array([pdot0,pdot1,pdot2,pdot3])

def calculate_updated_quaternion_FirstEuler(q_d, omega, dt):
	# This is the 'naive' first-order Euler method; the first thing that comes to your mind
	qdot = calculate_quaternion_time_derivative(q_d,omega,dt)
	return normalize([q_d[i] + dt*qdot[i] for i in range(4)])

def interpret_EulerTransformation_matrix(M, rotation_case):
	# Takes in a 3x3 transformation matrix M and interprets the Euler angles from it
	#print "DET=", np.linalg.det(M)
	#print "test:", [ ( M[2,1]/M[2,2]) , (-M[2,0]) , ( M[1,0]/M[0,0]) ]
	if rotation_case== "xyz": return [ math.atan2( M[2,1],M[2,2]) , math.asin(-M[2,0]) , math.atan2( M[1,0],M[0,0]) ]
	if rotation_case== "yxz": return [ math.atan2(-M[2,0],M[2,2]) , math.asin( M[2,1]) , math.atan2(-M[0,1],M[1,1]) ]
	if rotation_case== "xzy": return [ math.atan2(-M[1,2],M[1,1]) , math.asin( M[1,0]) , math.atan2(-M[2,0],M[0,0]) ]
	if rotation_case== "yzx": return [ math.atan2( M[0,2],M[0,0]) , math.asin(-M[0,1]) , math.atan2( M[2,1],M[1,1]) ]
	if rotation_case== "zyx": return [ math.atan2(-M[0,1],M[0,0]) , math.asin( M[0,2]) , math.atan2(-M[1,2],M[2,2]) ]
	if rotation_case== "zxy": return [ math.atan2( M[1,0],M[1,1]) , math.asin(-M[1,2]) , math.atan2( M[0,2],M[2,2]) ]

def Get_OmegaToEulerDot_matrix(phi, theta, psi, rotation_case):
	#print "Get_OmegaToEulerDot_matrix: interpreting rotation_case=", rotation_case
	if rotation_case == "xyz":
		# Matin versio
		A00, A01, A02	= math.cos(psi)			, math.sin(psi)			, 0.0
		A10, A11, A12	=-math.sin(psi)*math.cos(theta)	, math.cos(psi)*math.cos(theta)	, 0.0
		A20, A21, A22	= math.cos(psi)*math.sin(theta)	, math.sin(psi)*math.sin(theta)	, math.cos(theta)		# June 26 2018: Fixed 1.0 -> math.cos(theta)

		# FINFLO version
		#A00, A01, A02	= math.cos(theta)	, math.sin(phi)*math.sin(theta) , math.cos(phi)*math.sin(theta)
		#A10, A11, A12	= 0.0			, math.cos(phi)*math.cos(theta) ,-math.sin(phi)*math.cos(theta)
		#A20, A21, A22	= 0.0			, math.sin(phi)			, math.cos(phi)

	if rotation_case == "yxz":
		A00, A01, A02	=-math.sin(psi)			, math.cos(psi)			, 0.0
		A10, A11, A12	= math.cos(theta)*math.cos(psi)	, math.cos(theta)*math.sin(psi)	, 0.0
		A20, A21, A22	= math.sin(theta)*math.sin(psi)	,-math.sin(theta)*math.cos(psi)	, math.cos(theta)

	if rotation_case == "xzy":
		A00, A01, A02	= math.cos(psi)			, 0.0				,-math.sin(psi)
		A10, A11, A12	= math.cos(theta)*math.sin(psi)	, 0.0				, math.cos(theta)*math.cos(psi)
		A20, A21, A22	=-math.sin(theta)*math.cos(psi)	, math.cos(theta)		, math.sin(theta)*math.sin(psi)

	if rotation_case == "yzx":
		A00, A01, A02	=	0.0			, 		math.cos(psi)	, math.sin(psi)
		A10, A11, A12	=	0.0			, -math.cos(theta)*math.sin(psi), math.cos(theta)*math.cos(psi)
		A20, A21, A22	=	math.cos(theta)	,  math.sin(theta)*math.cos(psi), math.sin(theta)*math.sin(psi)
		
	if rotation_case == "zyx":
		A00, A01, A02	=	0.0			,-math.sin(psi)			, math.cos(psi)
		A10, A11, A12	=	0.0			, math.cos(theta)*math.cos(psi)	, math.cos(theta)*math.sin(psi)
		A20, A21, A22	=	math.cos(theta)	, math.sin(theta)*math.sin(psi)	,-math.sin(theta)*math.cos(psi)

	if rotation_case == "zxy":
		A00, A01, A02	= math.sin(psi)			, 	0.0			, math.cos(psi)
		A10, A11, A12	= math.cos(theta)*math.cos(psi)	,	0.0			,-math.cos(theta)*math.sin(psi)
		A20, A21, A22	= math.sin(theta)*math.sin(psi) , 	math.cos(theta)		, math.sin(theta)*math.cos(psi)

	#print "---------- EulerDot matrix:--------"
	#print [round(i, 4) for i in [A00, A01, A02]]
	#print [round(i, 4) for i in [A10, A11, A12]]
	#print [round(i, 4) for i in [A20, A21, A22]]
	#print " ----------------------------------"

	l1  = [ A00/math.cos(theta) , A01/math.cos(theta) , A02/math.cos(theta) ]
	l2  = [ A10/math.cos(theta) , A11/math.cos(theta) , A12/math.cos(theta) ]
	l3  = [ A20/math.cos(theta) , A21/math.cos(theta) , A22/math.cos(theta) ]
	return np.matrix([l1, l2, l3])

def intrinsic_switch(phi_d, theta_d, psi_d, rotation_case):
	# Switches the rotation case and Euler angles so that they correspond
	# to the intrinsic version of the rotation case
	# How to use: phi, theta, psi, rotation_case = intrinsic_switch(phi, theta, psi, rotation_case)
	return psi_d, theta_d, phi_d, rotation_case[::-1]

def LC(a,b,c):		# The Levi-Civita symbol
	i = a+1; j = b+1 ; k = c+1
	if i==1 and j == 2 and k==3	: return 1
	if i==3 and j == 1 and k==2	: return 1
	if i==2 and j == 3 and k==1	: return 1
	if k==1 and j == 2 and i==3	: return -1
	if k==3 and j == 1 and i==2	: return -1
	if k==2 and j == 3 and i==1	: return -1
	if (i-j)*(j-k)*(k-i) == 0	: return 0

def Kronecker_delta(i,j):
	if i-j == 0: return 1
	else: return 0

def generate_PRA_matrix(theta,n):
	# theta is angle in radians
	# and n is a unit vector (3-array)

	nx	= n[0]
	ny	= n[1]
	nz	= n[2]
	c	= math.cos(theta)
	s	= math.sin(theta)

	r0	= [	c + nx*nx*(1.0-c)	, nx*ny*(1.0-c) - nz*s		, ny*s	+nx*nz*(1.0-c)	]
	r1	= [  nz*s + nx*ny*(1.0-c)	, c + ny*ny*(1.0-c)		,-nx*s	+ny*nz*(1.0-c)	]
	r2	= [ -ny*s + nx*nz*(1.0-c)	, nx*s + ny*nz*(1.0-c)		, c + nz*nz*(1.0-c)	]
	return np.matrix([r0,r1,r2])	

def get_EulerDot(EulerAngles_d, omega_x, rotation_case):
	phi_d, theta_d, psi_d = EulerAngles_d.item(0), EulerAngles_d.item(1), EulerAngles_d.item(2)
	OmegaDotMatrix = Get_OmegaToEulerDot_matrix(phi_d, theta_d, psi_d, rotation_case)
	return np.dot(OmegaDotMatrix, np.transpose(np.matrix(omega_x)))

def Run_Static_RungeKutta(function,initial_value,h):
	# Runs a time integration using the 4th order Runge-Kutta method.
	# In this case, it is assumed that the function does not depend on
	# time, only the dependent variable.
	# h is the time step, should be equal to delta t etc.
	k1	= function(initial_value)
	k2	= function(np.add(initial_value , np.multiply(0.5*h,np.transpose(k1))))
	k3	= function(np.add(initial_value , np.multiply(0.5*h,np.transpose(k2))))
	k4	= function(np.add(initial_value , np.multiply(1.0*h,np.transpose(k3))))
	# We have to figure out the dimension of the k's in order to format the output
	try:
		return initial_value[1] + h/6.0 * (k1  + 2.0 * k2 +  2.0 * k3 + k4)
	except ValueError:
		return np.add(initial_value , np.multiply(h/6.0,np.add(np.add(np.transpose(k1) , np.multiply(2.0,np.transpose(k2))) , np.add(np.multiply(2.0,np.transpose(k3)),np.transpose(k4)))))




def Run_TwoVariable_RungeKutta(function,initial_value1, initial_value2,h):
	# Runs a time integration using the 4th order Runge-Kutta method.
	# In this case, 'function' takes two input values: initial_value1, initial_value2
	# On the other hand, function return exactly one scalar (y).
	k1	= h*function(initial_value1        , initial_value2         )
	k2	= h*function(initial_value1 + 0.5*h, initial_value2 + 0.5*k1)
	k3	= h*function(initial_value1 + 0.5*h, initial_value2 + 0.5*k2)
	k4	= h*function(initial_value1 + 1.0*h, initial_value2 + 1.0*k3)
	return initial_value2 + 1.0/6.0*(k1 + 2.0*k2 + 2.0*k3 + k4)
	

def Run_Static_ButcherKutta(function,initial_value,h):
	# Runs a time integration using Butcher's Fifth-order Runge-Kutta
	# In this case, it is assumed that the function does not depend on
	# time, only the dependent variable (usually denoted 'y').
	# h is the time step, should be equal to delta t etc.
	k1	= function(initial_value)
	k2	= function(np.add(initial_value, np.multiply(1.0/4.0*h,k1)))
	k3	= function(np.add(initial_value, np.add(  np.multiply( 1.0/8.0 * h, k1) , np.multiply(1.0/8.0*h, k2) )))
	k4	= function(np.add(initial_value, np.add(  np.multiply(-1.0/2.0 * h, k2) , np.multiply(        h, k3) )))
	k5	= function(np.add(initial_value, np.add(  np.multiply(3.0/16.0 * h, k1) , np.multiply(3.0/16.0*h,k4) )))
	k6	= function(
		   np.add( initial_value,
			np.add( np.multiply(-3.0/7.0*h, k1)
			      , np.add(  np.multiply(  2.0/7.0*h, k2)
			      , np.add(  np.multiply( 12.0/7.0*h, k3)
			      , np.add(  np.multiply(-12.0/7.0*h, k4)
			      , np.multiply(8.0/7.0*h, k5) ))))))
	term0	= np.add( np.multiply(32.0/90.0*h, k5)	,  np.multiply( 7.0/90.0*h, k6))
	term1	= np.add( np.multiply(32.0/90.0*h, k3)  ,  np.multiply(12.0/90.0*h, k4))
	term2   = np.add( term1                         , term0	)
	term3	= np.add( np.multiply(7.0/90.0*h, k1)   , term2 )
	return np.add(initial_value, np.transpose( term3)  )

def iterative_matrix_exp(M):
	# Calculates the exponential of 3x3 matrix M by taking terms of the Taylor polynomial.
	# M is a numpy matrix.
	n_terms = 35
	I	= np.matrix([[1.0,0.0,0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]])
	out	= I
	power	= I
	for i in range(n_terms):
		power	= np.dot(M, power)
		sign	= (1.0)/(math.factorial(i+1))
		out = np.add(  out ,  np.multiply(sign, power) )
	return out

def calculate_OmegaExpAlt(M, deltat):
	omega	=[ M.item(5), -M.item(2), M.item(1)]
	L	= math.sqrt(sum([e**2 for e in omega]))
	E	= [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
	for j in range(3):
		for i in range(3):
			E[i][j] = Kronecker_delta(i,j)*math.cos(L*deltat) + (1.0-math.cos(L*deltat))*omega[i]*omega[j]/L**2 - math.sin(L*deltat)/L*sum([LC(i,j,k)*omega[k] for k in range(3)])
	return np.matrix(E)

def calculate_OmegaExpVecAlt(omega, deltat):
	L	= math.sqrt(sum([e**2 for e in omega]))
	E	= [[0.0,0.0,0.0],[0.0,0.0,0.0],[0.0,0.0,0.0]]
	for j in range(3):
		for i in range(3):
			E[i][j] = Kronecker_delta(i,j)*math.cos(L*deltat) + (1.0-math.cos(L*deltat))*omega[i]*omega[j]/L**2 - math.sin(L*deltat)/L*sum([LC(i,j,k)*omega[k] for k in range(3)])
	return np.matrix(E)

def quaternion_update_Rabault(q_old,omega,dt, rotation_case):
	import time
	# Direct time stepping of orientation quaternions.
	# Based on 'A brief summary about quaternion use fordescribing rotations in 3D space' by J. Rabault, 2017

	# 2018-07-11: Observation: This works! But I only used a constant omega, so I don't know
	# if it works for time-dependent omega.
	l	= math.sqrt(sum([e**2 for e in omega]))
	theta	= dt*l
	omega0	= omega[0]
	omega1	= omega[1]
	omega2  = omega[2]

	startR00 = time.time()
	r0	= math.cos(theta/2.0)
	r1	= math.sin(theta/2.0)*omega0/l
	r2	= math.sin(theta/2.0)*omega1/l
	r3	= math.sin(theta/2.0)*omega2/l

	q0	= q_old[0]
	q1	= q_old[1]
	q2	= q_old[2]
	q3	= q_old[3]

	# The output is q_old multiplied by q from the right. And by 'multiplying' I mean quaternion multiplication.
	out0	= q0*r0 - q1*r1 - q2*r2 - q3*r3
	out1	= q1*r0 + q0*r1 - q3*r2 + q2*r3
	out2	= q2*r0 + q3*r1 + q0*r2 - q1*r3
	out3	= q3*r0 - q2*r1 + q1*r2 + q0*r3
	return normalize(np.array([out0,out1,out2,out3]))

def euler_coordinate_transform(x_input, y_input, z_input, rotation_case, phi_d, theta_d, psi_d, direction, intrinsic):
	M = Get_EulerTransformation_matrix(phi_d, theta_d, psi_d, direction, rotation_case, intrinsic)
	if isinstance(x_input, list):
		assert len(x_input) == len(y_input)
		assert len(y_input) == len(z_input)

		x_out = [0.0 for k in x_input]
		y_out = [0.0 for k in y_input]
		z_out = [0.0 for k in z_input]
	
		for i in range(len(x_input)):
			xi = np.matrix([[x_input[i]], [y_input[i]], [z_input[i]]])
			xii= np.dot(M, xi)
			x_out[i] = xii[0]
			y_out[i] = xii[1]
			z_out[i] = xii[2]

	else: 	# If the input is three scalars
		xii= np.dot(M, np.matrix([[x_input], [y_input], [z_input]]))
		x_out = [xii[0]]
		y_out = [xii[1]]
		z_out = [xii[2]]
	return x_out, y_out, z_out

def euler_coordinate_transform_vector(vector, rotation_case, phi_d, theta_d, psi_d, direction, intrinsic):
	M = Get_EulerTransformation_matrix(phi_d, theta_d, psi_d, direction, rotation_case, intrinsic)
	xii= np.dot(M, np.matrix([[vector[0]], [vector[1]], [vector[2]]]))
	x_out = xii.item(0)
	y_out = xii.item(1)
	z_out = xii.item(2)
	return [x_out, y_out, z_out]

def euler_filter(Euler_angles):
	# An attempt to make all (output) Euler angles unique.
	# When put through this filter, the numbers may change
	# but the orientation they represent shouldn't change (*praying*)
	phi_d 	= Euler_angles[0]
	theta_d = Euler_angles[1]
	psi_d 	= Euler_angles[2]
	out0	= phi_d
	out1	= theta_d
	out2	= psi_d
	#if abs(out0) > math.pi: out0 = 2.0*math.pi - out0
	#if abs(out0) > math.pi: out2 = 2.0*math.pi - out2
	if abs(out2) > math.pi: out2 = 2.0*math.pi - out2
	if max([abs(math.pi - (theta_d % (2.0*math.pi))),2.0*abs(phi_d), 2.0*abs(psi_d)])<math.pi/2:
		print "Filter applied!"
		out0 = (phi_d % (2.0*math.pi)) - math.pi
		out1 = math.pi - (theta_d % (2.0*math.pi))
		out2 = (psi_d % (2.0*math.pi)) - math.pi

	return [(out0 % (2.0*math.pi)), (out1%(2.0*math.pi)), (out2%(2.0*math.pi))]

def euler_filter2(angles):
	# I believe this is the one that I should use ...
	if abs(angles[1])>math.pi/2.0:
		#print "Filter applied!"
		return [(angles[0])%(2.0*math.pi) -math.pi , (-angles[1])%(2.0*math.pi) -math.pi ,(angles[2] )%(2.0*math.pi) -math.pi,]
	else: return [(angles[i]+math.pi)%(2.0*math.pi) - math.pi for i in range(3)]

def euler_filter_alt(Euler_angles):
	# An attempt to make all (output) Euler angles unique.
	# When put through this filter, the numbers may change
	# but the orientation they represent shouldn't change (*praying*)
	# This is a much more simple attempt, limiting the
	return [((i-math.pi)%(2.0*math.pi))-math.pi for i in Euler_angles]

def generate_aux_curve(phi_d, theta_d, psi_d, rotation_case, coordinate_system):
	dphi 	= math.radians(2.0)
	dtheta	= math.radians(2.0)
	dpsi	= math.radians(2.0)

	if phi_d < 0:	dphi = -dphi
	if theta_d < 0:	dtheta= -dtheta
	if psi_d < 0:	dpsi = -dpsi

	# This will draw red auxiliary lines on top of the "sphere" in order to further understand the Euler angles.
	if rotation_case == "xyz":
		A11, A12, A13 	= 1.0, 			0.0,            0.0
		A21, A22, A23	= 0.0,        math.cos(dphi), -math.sin(dphi)
		A31, A32, A33	= 0.0,        math.sin(dphi),  math.cos(dphi)

		B11, B12, B13 	= math.cos(dtheta),    	 0.0,  math.sin(dtheta)
		B21, B22, B23	= 0.0		 ,       1.0,  	           0.0
		B31, B32, B33	= -math.sin(dtheta),      0.0,  math.cos(dtheta)

		C11, C12, C13 	= math.cos(dpsi)	,-math.sin(dpsi),  0.0
		C21, C22, C23	= math.sin(dpsi)	, math.cos(dpsi),  0.0
		C31, C32, C33	= 0.0		,           0.0,  1.0

	if rotation_case == "xzy":
		A11, A12, A13 	= 1.0, 			0.0,            0.0
		A21, A22, A23	= 0.0,        math.cos(dphi), -math.sin(dphi)
		A31, A32, A33	= 0.0,        math.sin(dphi),  math.cos(dphi)

		B11, B12, B13 	= math.cos(dtheta)	,-math.sin(dtheta),  0.0
		B21, B22, B23	= math.sin(dtheta)	, math.cos(dtheta),  0.0
		B31, B32, B33	= 0.0		,           0.0,            1.0

		C11, C12, C13 	= math.cos(dpsi),    	 0.0,  math.sin(dpsi)
		C21, C22, C23	= 0.0		 ,       1.0,  	           0.0
		C31, C32, C33	= -math.sin(dpsi),        0.0,  math.cos(dpsi)

	if rotation_case == "yxz":
		A11, A12, A13 	= math.cos(dphi),    	 0.0,  math.sin(dphi)
		A21, A22, A23	= 0.0		 ,       1.0,  	           0.0
		A31, A32, A33	= -math.sin(dphi),        0.0,  math.cos(dphi)

		B11, B12, B13 	= 1.0, 			0.0,            0.0
		B21, B22, B23	= 0.0,        math.cos(dtheta), -math.sin(dtheta)
		B31, B32, B33	= 0.0,        math.sin(dtheta),  math.cos(dtheta)

		C11, C12, C13 	= math.cos(dpsi)	,-math.sin(dpsi),  0.0
		C21, C22, C23	= math.sin(dpsi)	, math.cos(dpsi),  0.0
		C31, C32, C33	= 0.0		,           0.0,  	   1.0

	if rotation_case == "yzx":
		A11, A12, A13 	= math.cos(dphi),    	 0.0,  math.sin(dphi)
		A21, A22, A23	= 0.0		 ,       1.0,  	           0.0
		A31, A32, A33	= -math.sin(dphi),        0.0,  math.cos(dphi)

		B11, B12, B13 	= math.cos(dtheta)	,-math.sin(dtheta),  0.0
		B21, B22, B23	= math.sin(dtheta)	, math.cos(dtheta),  0.0
		B31, B32, B33	= 0.0		,                   0.0,    1.0

		C11, C12, C13 	= 1.0, 			0.0,            0.0
		C21, C22, C23	= 0.0,        math.cos(dpsi), -math.sin(dpsi)
		C31, C32, C33	= 0.0,        math.sin(dpsi),  math.cos(dpsi)

	if rotation_case == "zyx":
		A11, A12, A13 	= math.cos(dphi)	,-math.sin(dphi),  0.0
		A21, A22, A23	= math.sin(dphi)	, math.cos(dphi),  0.0
		A31, A32, A33	= 0.0		,           0.0,           1.0

		B11, B12, B13 	= math.cos(dtheta),    	 0.0,  math.sin(dtheta)
		B21, B22, B23	= 0.0		 ,       1.0,  	           0.0
		B31, B32, B33	= -math.sin(dtheta),      0.0,  math.cos(dtheta)

		C11, C12, C13 	= 1.0, 			0.0,            0.0
		C21, C22, C23	= 0.0,        math.cos(dpsi), -math.sin(dpsi)
		C31, C32, C33	= 0.0,        math.sin(dpsi),  math.cos(dpsi)


	if rotation_case == "zxy":
		A11, A12, A13 	= math.cos(dphi)	,-math.sin(dphi),  0.0
		A21, A22, A23	= math.sin(dphi)	, math.cos(dphi),  0.0
		A31, A32, A33	= 0.0		,           0.0,           1.0

		B11, B12, B13 	= 1.0, 			0.0,            0.0
		B21, B22, B23	= 0.0,        math.cos(dtheta), -math.sin(dtheta)
		B31, B32, B33	= 0.0,        math.sin(dtheta),  math.cos(dtheta)

		C11, C12, C13 	= math.cos(dpsi),    	 0.0,  math.sin(dpsi)
		C21, C22, C23	= 0.0		 ,       1.0,  	           0.0
		C31, C32, C33	= -math.sin(dpsi),        0.0,  math.cos(dpsi)

	# Start drawing the dots - or making a list of their coordinates
	if coordinate_system == "default" or coordinate_system == "particle" : starting_vector = [1.0,0.0, 0.0]
	if coordinate_system == "FINFLO": starting_vector = [-1.0,0.0, 0.0]
	aux_x_history,   aux_y_history,   aux_z_history		= [], [], []

	# First we have phi_d
	phi_d_d 	= 0.0
	theta_d_d	= 0.0
	psi_d_d		= 0.0
	current_vector = [i for i in starting_vector]
	while abs(phi_d_d) < abs(phi_d):
		current_vector[0] = A11*current_vector[0] + A12*current_vector[1] + A13*current_vector[2]
		current_vector[1] = A21*current_vector[0] + A22*current_vector[1] + A23*current_vector[2]
		current_vector[2] = A31*current_vector[0] + A32*current_vector[1] + A33*current_vector[2]

		aux_x_history.append(current_vector[0])
		aux_y_history.append(current_vector[1])
		aux_z_history.append(current_vector[2])
		phi_d_d = phi_d_d + dphi

	# Then theta_d
	while abs(theta_d_d) < abs(theta_d):
		current_vector[0] = B11*current_vector[0] + B12*current_vector[1] + B13*current_vector[2]
		current_vector[1] = B21*current_vector[0] + B22*current_vector[1] + B23*current_vector[2]
		current_vector[2] = B31*current_vector[0] + B32*current_vector[1] + B33*current_vector[2]

		aux_x_history.append(current_vector[0])
		aux_y_history.append(current_vector[1])
		aux_z_history.append(current_vector[2])
		theta_d_d = theta_d_d + dtheta

	# Then psi_d
	while abs(psi_d_d) < abs(psi_d):
		current_vector[0] = C11*current_vector[0] + C12*current_vector[1] + C13*current_vector[2]
		current_vector[1] = C21*current_vector[0] + C22*current_vector[1] + C23*current_vector[2]
		current_vector[2] = C31*current_vector[0] + C32*current_vector[1] + C33*current_vector[2]

		aux_x_history.append(current_vector[0])
		aux_y_history.append(current_vector[1])
		aux_z_history.append(current_vector[2])
		psi_d_d = psi_d_d + dpsi

	return aux_x_history, aux_y_history, aux_z_history



def WindComponent(position, velocity, time):
    """
    Calculate wind force component for Mars flight trajectory simulation.
    
    Parameters:
    position (np.array): 3D position vector [x, y, z] in km
    velocity (np.array): 3D velocity vector [vx, vy, vz] in km/s
    time (float): Current simulation time in seconds
    
    Returns:
    np.array: 3D wind force vector [Fx, Fy, Fz] in Newtons
    """
    # Mars atmospheric density model (simplified)
    # Exponential decay with altitude
    altitude = np.linalg.norm(position) - 3389.5  # Mars radius is ~3389.5 km
    density = 0.020 * np.exp(-altitude / 11.1)  # kg/m^3
    
    # Wind speed varies with altitude and time
    # Adding some random oscillations for variety
    wind_base = np.array([
        30 * np.sin(time / 1000 + altitude / 50),  # x-component
        25 * np.cos(time / 800 + altitude / 40),   # y-component
        15 * np.sin(time / 1200)                   # z-component
    ])
    
    # Dust storm effect (random occurrence)
    dust_storm = np.random.random() < 0.01  # 1% chance of dust storm
    if dust_storm:
        wind_base *= 3.5
    
    # Calculate relative wind velocity
    relative_wind = wind_base - velocity
    
    # Simple drag equation: F = 1/2 * rho * v^2 * Cd * A
    drag_coefficient = 0.47  # Approximate sphere drag coefficient
    reference_area = 10.0    # m^2
    
    wind_force = (0.5 * density * np.linalg.norm(relative_wind) * 
                 drag_coefficient * reference_area * 
                 relative_wind)
    
    # Add turbulence
    turbulence = np.random.normal(0, 0.1, 3) * np.linalg.norm(wind_force)
    wind_force += turbulence
    
    return wind_force