import math as m

d = 100.0
minor_radius = 0.1

pi = m.acos(-1.0)
angle_start = 0.15
angle_end = pi - angle_start
steps = 800
angle_step = (angle_end-angle_start)/steps

cell_lines = []
surface_lines = []

for i in range(0,steps+1):
    angle = angle_start + i*angle_step;
    mu = m.cos(angle)
    z_position = d*mu
    major_radius = d*m.sin(angle)
    
    surface_lines.append(str(i+1) + " tz 0 0 " + str(z_position) + " " + str(major_radius) + " 0.1 0.1")
    cell_lines.append(str(i+1) + "  0        -"+str(i+1)+"  imp:p=1")

surface_lines.append(str(steps+2)+" pz -0.0005")
surface_lines.append(str(steps+3)+" pz  0.0005")
surface_lines.append(str(steps+4)+" cz  0.1")
surface_lines.append(str(steps+5)+" so  10000.0")

cell_lines.append(str(steps+2)+"  1  -0.001 "+str(steps+2)+" -"+str(steps+3)+" -"+str(steps+4)+"  imp:p=1")
cell_lines.append(str(steps+3)+"  0  -"+str(steps+5)+" (-"+str(steps+2)+":"+str(steps+3)+":"+str(steps+4)+") ")

for i in range(0,steps+1):
    if len(cell_lines[-1]) > 70:
        cell_lines[-1] += " &"
        cell_lines.append("  ")
    
    cell_lines[-1] += str(i+1) + " "

cell_lines[-1] += "  imp:p=1"
cell_lines.append(str(steps+4)+"  0         "+str(steps+5)+"  imp:p=0")

print "dyson sphere"

for i in range(0,len(cell_lines)):
    print cell_lines[i]

print ""

for i in range(0,len(surface_lines)):
    print surface_lines[i]

print ""

data_lines = []
data_lines.append("mode p")
data_lines.append("nps   1e8")
data_lines.append("sdef  pos = 0 0 -0.1  dir = 1  vec = 0 0 1  erg = 0.1")
data_lines.append("c")
data_lines.append("m1    1000.12p 1.0")
data_lines.append("c")
data_lines.append("fcl:p ")

for i in range(0,steps+1):
    if len(data_lines[-1]) > 70:
        data_lines[-1] += " &"
        data_lines.append("  ")
    
    data_lines[-1] += "0 "

data_lines[-1] += "-1 0 0"
                  
data_lines.append("f04:p ")

for i in range(0,steps+1):
    if len(data_lines[-1]) > 70:
        data_lines[-1] += " &"
        data_lines.append("  ")
        
    data_lines[-1] += str(i+1) + " "

data_lines.append("e0 1e-3 998i 0.1")
data_lines.append("ft4 INC")
data_lines.append("fu4 0 1 10")
data_lines.append("phys:p 100 1 0 0 0")
data_lines.append("prdmp  j  1e7  1   1")

for i in range(0,len(data_lines)):
    print data_lines[i]


