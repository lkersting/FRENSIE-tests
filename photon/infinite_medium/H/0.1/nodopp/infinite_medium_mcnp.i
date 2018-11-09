10.0 MeV photons in room temp water
100   1    -1.0       -10         imp:n=1
101   1    -1.0       -11 10      imp:n=1
102   1    -1.0       -12 11      imp:n=1
103   1    -1.0       -13 12      imp:n=1
104   1    -1.0       -14 13      imp:n=1
105   1    -1.0       -15 14      imp:n=1
999   0               +15         imp:n=0

10   so    10.0
11   so    20.0
12   so    30.0
13   so    40.0
14   so    50.0
15   so    10000.0

mode p
nps    1e8
sdef   pos = 0 0 0   erg = 0.1
c
m1     1000.12p   1.0
c
f01:p  10 11 12 13 14
f02:p  10 11 12 13 14
e0     1e-3 1.5e-3 198i 0.1
c
phys:p j  1  0 0  1
prdmp  j  1e7  1   1
