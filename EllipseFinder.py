def countNumberOfPointInsideEllipse(xAxis,yAxis):
	global x_coordinates
	global y_coordinates
	global number_of_points
	global x_center
	global y_center

	successes = 0
	for i in range(number_of_points):
		conditionNumber = ( (x_coordinates[i]-x_center)/xAxis )**2 + ( (y_coordinates[i]/1000.0-y_center)/yAxis )**2
		if conditionNumber<=1.0: successes += 1
	return successes
	


def main():
	global x_coordinates
	global y_coordinates
	global number_of_points
	global x_center
	global y_center
	filename = ''

	# First we read in the file
	print("Make sure that x-coordinates are in kilometres and ")
	print("y-coordinates are in metres.")
	print("Please enter the data file name:")
	filename = input()
	f = open(filename, "r")

	for line in f:
		line_split = line.lstrip().split('\t')
		x_coordinates.append(float(line_split[0]))
		y_coordinates.append(float(line_split[1]))

	# Now we need to know the aspect ratio of the ellipse
	print("Please enter the aspect ratio of the ellipse (x/y):")
	aspect_ratio = float(input())

	# Figure out the center of the ellipse
	number_of_points =len(x_coordinates)
	x_center = sum(x_coordinates)/number_of_points
	y_center = sum(y_coordinates)/(number_of_points*1000.0)

	print("x center [km] = ", x_center)
	print("y center [km] = ", y_center)

	# Then we can start guessing the y-direction radius of the ellipse!
	y_radius = 0.5
	while True:
		ratio = countNumberOfPointInsideEllipse(y_radius*aspect_ratio,y_radius)/number_of_points
		if ratio>0.95:break
		y_radius = y_radius*1.025

	print("x radius [km] = ", y_radius*aspect_ratio)
	print("y radius [km] = ", y_radius)
	
x_coordinates = []
y_coordinates = []
number_of_points = 0
x_center = 0.0
y_center = 0.0

main()
