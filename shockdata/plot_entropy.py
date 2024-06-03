
from pylab import *

figure(6)

clf()
t,entropy = loadtxt('_output/fort.69',unpack=True)
plot(t,entropy)
title("Entropy")
xlabel("time")
ylabel("entropy")
