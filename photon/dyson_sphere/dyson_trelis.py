import math as m

d = 100.0
minor_radius = 0.1

pi = m.acos(-1.0)
angle_start = 0.15
angle_end = pi - angle_start
angle_step = (angle_end-angle_start)/800

for i in range(0,801):
    angle = angle_start + i*angle_step;
    mu = m.cos(angle)
    z_position = d*mu
    major_radius = d*m.sin(angle)
    print z_position
    cmd1 = "Create Torus Major "
    cmd1 += str(major_radius)
    cmd1 += " Minor 0.1"
    print cmd1
    cmd2 = "Vol "
    cmd2 += str(i+1)
    cmd2 += " Move Z "
    cmd2 += str(z_position)
    print cmd2
    cubit.cmd(cmd1)
    cubit.cmd(cmd2)
    cubit.cmd("Group \"estimator_1.cell.tl.flux.p\" add vol "+str(i+1))

cubit.cmd("Create Cylinder Height 0.001 Radius 1.0")
cubit.cmd("Group \"mat_1_rho_-0.001\" add vol 802")

cubit.cmd("Create Brick X 300 Y 300 Z 300")
cubit.cmd("Create Brick X 301 Y 301 Z 301")
cubit.cmd("Subtract Vol 803 from Vol 804")
cubit.cmd("Group \"termination.cell\" add vol 805")
