import numpy as np
import matplotlib
matplotlib.rcParams['font.size'] = 15

import pylab as pl

def plotnrg(nrgs):
    t=np.arange(80,510,10)
    print t.shape
    print nrgs[0].shape
    for nrg in nrgs:
        nrg=nrg/nrg[0]
        pl.plot(t,nrg,linewidth=3)
        pl.hold(True)

    pl.xlabel('t',fontsize=20)
    pl.ylabel('$\eta$',fontsize=30)
    pl.legend(['Z=4','Z=3','Z=2','Z=1'],loc='best')
    pl.axis([80,500,0.5,1.1])
    pl.hold(False)
