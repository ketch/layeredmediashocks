from numpy import sqrt, log

def cshock(sig,rhoB=4,KB=4,rhoA=1,KA=1):
    """Compute the homogenized shock speed.
       Assumes right state is zero.
    """

    rhobar = 0.5*(rhoA+rhoB)
    #sigmapA= KA*(sig+1)
    #sigmapB= KB*(sig+1)

    epsA=log(sig+1.)/KA
    epsB=log(sig+1.)/KB
    sigmapA=sig/epsA
    sigmapB=sig/epsB
    sigmapbar = 2./(1./sigmapA + 1./sigmapB)

    cshockeff = sqrt(sigmapbar/rhobar)

    return cshockeff

def makeplot(Ks,rhos,sigma_max=5):
    """
    This is used to make Figure \ref{fig:cshock} (Figure 3.1) of the paper, as follows:

        >>> makeplot((1.,2.,3.,4.),(1.,2.,3.,4.),2.)
    """
    import numpy as np
    import pylab as pl

    assert len(Ks)==len(rhos), 'input arrays must be of same length'

    linestyles=['-k','--k','-.k',':k','-b']

    sigs=np.linspace(0.01,sigma_max,400)

    ii=0

    for KB,rhoB in zip(Ks,rhos):
        cbars=cshock(sigs,rhoB,KB)
        cA = 1.#sqrt(KA/rhoA) -- it is assumed that KA=rhoA=1.
        cB = sqrt(KB/rhoB)
        cmeanlin = 2./(1./cA + 1./cB)

        pl.plot(sigs,cbars/cmeanlin,linestyles[ii],linewidth=3)
        pl.hold(True)
        ii=ii+1


    #legstr=['$(rho_B,K_B)=('+str(rhoB)+','+str(KB)+')$' for KB,rhoB in zip(Ks,rhos)]
    legstr=['$Z_B='+str(int(rhoB))+'$' for rhoB in rhos]
    pl.legend(legstr,loc='best')
    pl.plot([min(sigs),max(sigs)],[1.,1.],'--k')
    pl.xlabel('$\sigma_l$',fontsize=30); pl.ylabel('$s_{eff}$',fontsize=30)
    pl.draw()
    pl.hold(False)
