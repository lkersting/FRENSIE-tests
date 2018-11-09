import math as m

d = 100.0
minor_radius = 0.1
pi = m.acos(-1.0)

# First angle (0.15)
angle = 0.15
mu = m.cos(angle)
z_position = d*mu
major_radius = d*m.sin(angle)

cmd1 = "Create Torus Major " + str(major_radius) + " Minor 0.1"
cmd2 = "Vol 1 Move Z " + str(z_position)
cubit.cmd(cmd1)
cubit.cmd(cmd2)

# Intermediate angles (11 steps of pi/12)
for i in range(1,12):
    angle = pi/12*i
    mu = m.cos(angle)
    z_position = d*mu
    major_radius = d*m.sin(angle)
    cmd1 = "Create Torus Major " + str(major_radius) + " Minor 0.1"
    cmd2 = "Vol " + str(i+1) + " Move Z " + str(z_position)
    cubit.cmd(cmd1)
    cubit.cmd(cmd2)

# Finale angle (pi-0.15)
angle=pi-0.15
mu = m.cos(angle)
z_position = d*mu
major_radius = d*m.sin(angle)

cmd1 = "Create Torus Major " + str(major_radius) + " Minor 0.1"
cmd2 = "Vol 13 Move Z " + str(z_position)
cubit.cmd(cmd1)
cubit.cmd(cmd2)

# Create the estimators
cubit.cmd("Group \"estimator_1.cell.tl.flux.p\" add vol 1 2 3 4 5 6 7 8 9 10 11 12 13")
cubit.cmd("Group \"estimator_2.cell.tl.flux.p\" add vol 1 2 3 4 5 6 7 8 9 10 11 12 13")

# Create the cylinder where collisions will be forced
cubit.cmd("Create Cylinder Height 0.001 Radius 0.1")
cubit.cmd("Group \"mat_1_rho_-0.001\" add vol 14")

# Create the termination cell
cubit.cmd("Create Brick X 300 Y 300 Z 300")
cubit.cmd("Create Brick X 301 Y 301 Z 301")
cubit.cmd("Subtract Vol 15 from Vol 16")
cubit.cmd("Group \"termination.cell\" add vol 17")


