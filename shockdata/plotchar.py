
from pylab import *
from pyclaw.data import Data

#xmin = 50.
#tmin = 30.
#xmin = 100.
#tmin = 90.
xmin = 0.
tmin = 0.
setprob_data = Data('setprob.data')
KA = setprob_data.K_A
KB = setprob_data.K_B
rhoA = setprob_data.rho_A
rhoB = setprob_data.rho_B
c = 2.*sqrt(1./((1./KA + 1./KB)*(rhoA+rhoB)))
print "c = ",c

chardata = loadtxt('_output/fort.31')
tchar = chardata[:,0]
xchar = chardata[:,1:]
figure(5,figsize=(10,12))
clf()
subplot(211)
for j in range(xchar.shape[1]):
    plot(xchar[:,j], tchar,'b')
xlabel("x")
ylabel("t")
title("Characteristics")

xtmaxdata = loadtxt('_output/fort.68')
plot(xtmaxdata[:,1],xtmaxdata[:,0],'r')
#plot(xtmaxdata[:,2],xtmaxdata[:,0],'m')
#xlim([100,180])
#ylim([110,200])
xlim(xmin=xmin)
ylim(ymin=tmin)

if 0:
    for x in range(100,201):
        plot([x,x],[0,200],'g')


subplot(212)
for j in range(xchar.shape[1]):
    plot(xchar[:,j]-tchar*c, tchar,'b')
plot(xtmaxdata[:,1]-xtmaxdata[:,0]*c,xtmaxdata[:,0],'r')
title("Shifted characteristics")
xlabel("x-ct, c = %6.4f" % c, fontsize=15)
ylabel("t", fontsize=15)
ylim(ymin=tmin)
xlim(xmax=40)

figure(6)

clf()
t,entropy = loadtxt('_output/fort.69',unpack=True)
plot(t,entropy)
title("Entropy")
xlabel("time")
ylabel("entropy")
xlim(xmin=tmin)
