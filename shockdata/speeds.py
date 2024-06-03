
from numpy import *

def speeds(KA,rhoA,KB,rhoB,sigmaL):

    Kharm = 2./(1./KA + 1./KB)
    rhobar = 0.5*(rhoA + rhoB)
    cA = sqrt(KA/rhoA)
    cB = sqrt(KB/rhoB)
    ZA = rhoA*cA
    ZB = rhoB*cB

    # homogenized linearized speeds:
    cl = sqrt(Kharm*(1.+sigmaL) / rhobar)
    cr = sqrt(Kharm*(1.+0.) / rhobar)

    epslA = log(sigmaL+1)/KA
    epslB = log(sigmaL+1)/KB

    # homogenized shock speed:
    seff = sqrt(2.*sigmaL/(epslA+epslB) / rhobar)

    # speed of fastest shock:
    cmax = 2./(1./cA + 1./cB)

    print "cA = %s, ZA = %s" %(cA,ZA)
    print "cB = %s, ZB = %s" %(cB,ZB)
    print "cl = %s, cr = %s" %(cl, cr)
    print "effective shock speed:  seff = %s" % seff
    print "maximum   shock speed:  cmax = %s" % cmax
    
