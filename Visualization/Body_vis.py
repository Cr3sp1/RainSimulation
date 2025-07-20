from vpython import *


def plot_geometric_bodies(filename):

    # Read the data from the file and visualize

    with open(filename, 'r') as file:

        for line in file:

            parts = line.strip().split(',')

            body_type = parts[0]

            params = list(map(float, parts[1:]))



            if body_type == 'S':  # Sphere

                center = vector(params[0], params[1], params[2])

                radius = params[3]

                sphere(pos=center, radius=radius, color=color.blue)



            elif body_type == 'P':  # Parallelepiped
                center = vector(params[0], params[1], params[2])
                axis1 = vector(params[3], params[4], params[5])
                axis2 = vector(params[6], params[7], params[8])
                axis3 = vector(params[9], params[10], params[11])

                # Define the box
                box(pos=center, size=vector(mag(axis1), mag(axis2), mag(axis3)), 
                    axis=axis1, up=axis2, texture='checkers.png')




            elif body_type == 'C':  # Capsule
                point1 = vector(params[0], params[1], params[2])

                point2 = vector(params[3], params[4], params[5])

                radius = params[6]

                axis = point2 - point1

                cylinder(pos=point1, axis=axis, radius=radius, texture='stripes.png')

                sphere(pos=point1, radius=radius, texture='stripes.png')

                sphere(pos=point2, radius=radius, texture='stripes.png')

# plot_geometric_bodies('../data/Temp/Vit.dat')
plot_geometric_bodies('../data/Run/Status/Run0.750000.dat')


# Adjust camera position and orientation
walk_pos = vector(0.446092, -2, 1.2)
run_pos = vector(0.7, -2, 1.2)
scene.camera.pos = run_pos
scene.camera.axis = vector(0, 1, 0)
scene.camera.up = vector(0, 0, 1)

# Change background color
scene.background = color.white

# Add a distant light
distant_light(direction=vector(0, -1, 0.3), color=vector(1,1,1))
distant_light(direction=vector(1, 0, 1), color=vector(0.7,0.7,0.7))

# Set the resolution of the scene
scene.width = 800  # Set the width in pixels
scene.height = 600  # Set the height in pixels
