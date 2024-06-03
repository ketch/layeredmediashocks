#Test entropy evolution after time reversal
import os
import setrun
#from setrun import setrun
from entropy import entropy
import pylab as pl
from pyclaw.runclaw import runclaw
import numpy as np

reload(setrun)
#grids = [450,900,1800,3600]
#Ks = np.array([1.1,1.3,1.5,1.7,1.9,2.0,2.1,2.3,2.5,2.7,3.9])
Ks = np.array([1.0,1.5,2.0,2.5,3.0,3.5,4.0,4.5,5.0,5.5])
#Ks = np.array([3.5])
ntrp = np.zeros(np.shape(Ks))

for i,K in enumerate(Ks):
  rundata=setrun.setrun()
  rundata.probdata.K_B=K
  rundata.probdata.rho_B=K
  rundata.probdata.K_A=1.
  rundata.probdata.rho_A=1.
  print K

  rundata.probdata.a1=0.1
  rundata.probdata.trtime = 250.
  rundata.clawdata.tfinal = 500.

  #Rescale time and IC:
  #cmean = 2*np.sqrt(K/((K+1.)*(rho+1.)))
  #Kmean = 2*K/(K+1.)
  #rhomean = 0.5*(rho+1.)
  #rundata.probdata.a1 = rundata.probdata.a1/np.sqrt(rhomean)
  #rundata.clawdata.tfinal = 500/(cmean/0.8)
  rundata.write()
  runclaw(xclawcmd='xsclaw',outdir='./_output')

  ent=entropy(54,6)
  ntrp[i]=ent[-1]
  pl.plot(Ks,np.log10(abs(1.-ntrp)))
  pl.draw()

pl.plot(Ks,np.log10(abs(1.-ntrp)))
pl.draw()
