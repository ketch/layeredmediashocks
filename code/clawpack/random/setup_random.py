#Generate randomized interface locations and write to a file
#The layer half-widths are normally distributed with mean delta/2

from numpy import random, savetxt

nlayer = 300
delta = 1.
dh = delta/2.

scale=dh/1.

deviations=random.normal(0.,scale,(nlayer))

savetxt("randint.data",deviations)
