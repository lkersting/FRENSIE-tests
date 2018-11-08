import math as m

d = 100.0
minor_radius = 0.1

pi = m.acos(-1.0)
angle_start = 0.15
angle_end = pi - angle_start
steps = 800
angle_step = (angle_end-angle_start)/steps

vol_counter = 1

# i = 0   ==> phi == 0.15
# i = 179 ==> phi ~ pi/4
# i = 253 ==> phi ~ pi/3
# i = 326 ==> phi ~ 5pi/12
# i = 400 ==> phi ~ pi/2
# i = 474 ==> phi ~ 7pi/12
# i = 547 ==> phi ~ 2pi/3
# i = 621 ==> phi ~ 3pi/4
# i = 800 ==> phi == pi - 0.15
for i in range(0,steps+1):
    if i == 0 or i == 179 or i == 253 or i== 326 or i == 400 or i == 474 or i == 547 or i == 621 or i == steps:
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
        cmd2 += str(vol_counter)
        cmd2 += " Move Z "
        cmd2 += str(z_position)
        print cmd2
        cubit.cmd(cmd1)
        cubit.cmd(cmd2)
        cubit.cmd("Group \"estimator_1.cell.tl.flux.p\" add vol "+str(vol_counter))
        vol_counter += 1

cubit.cmd("Create Cylinder Height 0.001 Radius 0.1")
cubit.cmd("Group \"mat_1_rho_-0.001\" add vol "+ str(vol_counter))

cubit.cmd("Create Brick X 300 Y 300 Z 300")
cubit.cmd("Create Brick X 301 Y 301 Z 301")
cubit.cmd("Subtract Vol "+str(vol_counter+1)+" from Vol "+str(vol_counter+2))
cubit.cmd("Group \"termination.cell\" add vol "+str(vol_counter+3))
