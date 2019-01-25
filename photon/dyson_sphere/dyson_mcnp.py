import math as m

d = 100.0
minor_radius = 0.1
pi = m.acos(-1.0)

cell_lines = []
surface_lines = []

# First angle (0.15)
angle = 0.15
mu = m.cos(angle)
z_position = d*mu
major_radius = d*m.sin(angle)

surface_lines.append("1 tz 0 0 " + str(z_position) + " " + str(major_radius) + " 0.1 0.1")
cell_lines.append("1  0         -1  imp:p=1")

# Intermediate angles (11 teps of pi/12)
for i in range(1,12):
    angle = pi/12*i
    mu = m.cos(angle)
    z_position = d*mu
    major_radius = d*m.sin(angle)
    
    surface_lines.append(str(i+1) + " tz 0 0 " + str(z_position) + " " + str(major_radius) + " 0.1 0.1")
    cell_lines.append(str(i+1) + "  0        -"+str(i+1)+"  imp:p=1")

# Finale angle (pi-0.15)
angle=pi-0.15
mu = m.cos(angle)
z_position = d*mu
major_radius = d*m.sin(angle)

surface_lines.append("13 tz 0 0 " + str(z_position) + " " + str(major_radius) + " 0.1 0.1")
cell_lines.append("13  0        -13  imp:p=1")

# Create the cylinder where collisions will be forced
surface_lines.append("14 pz -0.0005")
surface_lines.append("15 pz  0.0005")
surface_lines.append("16 cz  0.1")
cell_lines.append("14  1  -0.001  14 -15 -16  imp:p=1")

# Create the complement cell and the termination cell
surface_lines.append("17 so  10000.0")

cell_lines.append("15  0         -17 (-14:15:16) 1 2 3 4 5 6 7 8 9 10 11 12 13  imp:p=1")
cell_lines.append("16  0         17  imp:p=0")

# Print the cell card
print "dyson sphere"

for i in range(0,len(cell_lines)):
    print cell_lines[i]

print ""

# Print the surface card
for i in range(0,len(surface_lines)):
    print surface_lines[i]

print ""

# Create the data card lines
data_lines = []
data_lines.append("mode p")
data_lines.append("nps   1e9")
data_lines.append("sdef  pos = 0 0 -0.1  dir = 1  vec = 0 0 1  erg = 0.1")
data_lines.append("c")
data_lines.append("m1    1000.12p 1.0")
data_lines.append("c")
data_lines.append("fcl:p 0 0 0 0 0 0 0 0 0 0 0 0 0 -1 0 0")

data_lines.append("f04:p 1 2 3 4 5 6 7 8 9 10 11 12 13")
data_lines.append("f14:p 1 2 3 4 5 6 7 8 9 10 11 12 13")
data_lines.append("e04 1e-3 998i 0.1")
data_lines.append("e14 1e-3 9980i 0.1")
data_lines.append("ft4 INC")
data_lines.append("fu4 0 1")
data_lines.append("phys:p 100 1 0 0 0")
data_lines.append("prdmp  j  1e7  1   1")

# Print the data card
for i in range(0,len(data_lines)):
    print data_lines[i]


