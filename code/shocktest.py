"""
Scripts to test propagation of shocks in periodic layered media.
The initial condition used is a discontinuity that satisfies a
"homogenized" Rankine-Hugoniot condition.
"""
def runtests(mx=7200,shock_amplitude=[0.1,0.2,0.4,0.8,1.2],rhos=[2.0,4.0,9.0],Ks=[1.0,1.0,1.0]):
    """
    Run the test for a given grid size and a sequence of shock
    amplitudes and material densities and bulk moduli.
    Note that the script doesn't run all combinations from (rhos,Ks) but just
    all sequential pairs.
    """
    import setrun
    from pyclaw.runclaw import runclaw
    import numpy as np

    assert len(rhos)==len(Ks)

    reload(setrun)

    for rhoB,KB in zip(rhos,Ks):
        for a2 in shock_amplitude:
            rundata=setrun.setrun('sharpclaw')

            #Material parameters
            rundata.probdata.K_B=KB
            rundata.probdata.rho_B=rhoB
            rundata.probdata.K_A=1.
            rundata.probdata.rho_A=1.

            #No boundary motion and don't switch to periodic BCs
            rundata.probdata.t1 =1000.0
            rundata.probdata.tw1=1000.0
            rundata.probdata.a1=0.0
            
            #The initial condition is a constant-stress region
            #and an unstressed region, separated by a front.
            rundata.probdata.ic=5

            #a2 is the magnitude of the stress behind the front
            rundata.probdata.a2=a2

            #We ought to rerun these with a larger domain and longer times
            #to see more clearly the entropy loss.
            rundata.probdata.trtime = 1000.
            rundata.clawdata.tfinal = 100.
            rundata.clawdata.xupper=150.0

            #Reflecting at left, outflow at right
            rundata.clawdata.mthbc_xlower = 3
            rundata.clawdata.mthbc_xupper = 1
            rundata.clawdata.mx=mx

            rundata.clawdata.nout = 100

            rundata.write()
            runclaw(xclawcmd='xsclaw',outdir='./shock/'+str(mx)+'_'+str(a2)+'_'+str(rhoB)+'_'+str(KB))

def entropy_change_convergence(grids=[3600,7200,14400,28800],shock_amplitude=1.2,rho=4.0,K=4.0):
    """
    Make a convergence plot for the entropy change versus resolution.
    """
    from entropy import entropy_oneframe
    import cbar
    import matplotlib.pyplot as pl
    import numpy as np

    a2=shock_amplitude
    ent=[]; rate=[]; relentloss=[]

    for i,mx in enumerate(grids):
        ent0=entropy_oneframe(50,
                outdir='./shock/'+str(mx)+'_'+str(a2)+'_'+str(rho)+'_'+str(K),)
        ent.append(entropy_oneframe(100,
                outdir='./shock/'+str(mx)+'_'+str(a2)+'_'+str(rho)+'_'+str(K),)/ent0)
        relentloss.append(abs(1.-ent[-1]))
        if i>0: rate.append(np.log2(relentloss[i-1]/relentloss[i]))

    pl.clf()
    pl.loglog(grids,relentloss)
    print rate
    return grids, ent


def shock_entropy_plot(mx=7200,shock_amplitude=[0.1,0.2,0.4,0.8],rho=2.0,K=2.0):
    """
    Plot entropy evolution for a range of initial shock amplitudes.
    This is used to generate Figures 3.3 and 3.4.
    """
    from entropy import entropy
    import cbar
    import matplotlib.pyplot as pl

    pl.clf()
    linestyles=['-k','--k','-.k',':k']
    speeds=[cbar.cshock(amp,rho,K) for amp in shock_amplitude]
    speedstring=[str(speed)[:4] for speed in speeds]
    for i,a2 in enumerate(shock_amplitude):
        entropy(0,100,outdir='./shock/'+str(mx)+'_'+str(a2)+'_'+str(rho)+'_'+str(K),
                                style=linestyles[i])
        pl.hold(True)

    speedstring=['$S_{eff}='+str(round(speed,2))+'$' for speed in speeds]
    print speedstring
    pl.legend(speedstring,loc='best')
    pl.xlabel('Time')
    pl.ylabel('Relative total entropy')
    pl.hold(False)


def shock_location_plot(mx=7200,shock_amplitude=[0.4,0.8],rho=2.0,K=2.0):
    """
    Plot predicted and observed location of shock versus time.
    Note that the method for finding the shock location is a bit ad hoc and may fail.
    """
    from pyclaw.plotters.data import ClawPlotData
    import matplotlib.pyplot as pl
    import numpy as np
    import cbar

    pl.clf()

    for a2 in shock_amplitude:
        outdir='./shock/'+str(mx)+'_'+str(a2)+'_'+str(rho)+'_'+str(K)

        plotdata = ClawPlotData()
        plotdata.outdir=outdir
        aux = np.loadtxt(outdir+'/fort.aux')

        t=[]
        xshock=[]
        xpreds=[]

        cshock=cbar.cshock(a2,rho,K)

        for i in range(101):
            #Read in the solution
            dat = plotdata.getframe(i)
            t.append(dat.t)

            xpred = 30. + dat.t*cshock #Predicted location
            xpreds.append(xpred)
            ipred = int(np.round(xpred/dat.d)[0])
            print xpred,ipred

            sig = dat.q[:,0]
            dsig = np.diff(sig)
            ishock = np.argmax(-sig[(ipred-1000):(ipred+1000)])+ipred-1000
            print ishock
            xshock.append(0.5*(dat.c_center[0][ishock]+dat.c_center[0][ishock+1]))

        pl.plot(t,xshock,t,xpreds)



def compute_entropy(mx=7200,a2s=[0.1,0.2,0.4,0.8,1.2],rhos=[2.0,4.0,9.0],Ks=[1.0,1.0,1.0]):
    """
    So far I'm not using this function.
    """
    import pylab as pl
    from entropy import entropy
    import numpy as np
    from cbar import cshock

    assert len(rhos)==len(Ks)

    ntrp = []  #Storage for final entropy ratios
    cshockhom = []  #storage for effective shock speeds

    for rhoB,KB in zip(rhos,Ks):
        for a2 in a2s:
            entropy_time_sequence,t = entropy(0,10,outdir='./shock/'+str(mx)+'_'+str(a2)+'_'+str(rhoB)+'_'+str(KB),doplot=False)
            ntrp.append(entropy_time_sequence[-1])

            cA = 1.#np.sqrt(KA/rhoA)
            cB = np.sqrt(KB/rhoB)
            cmeanlin = 2./(1./cA + 1./cB)

            cshockhom.append(cshock(a2/2.,rhoB,KB)/cmeanlin)
            pl.clf()
            pl.plot(cshockhom,ntrp,'o')
            pl.draw()
            print np.log10(abs(1.-np.array(ntrp)))
#            except:
#                print 'no files found for K='+str(KB)+', rho='+str(rhoB)+' and a2='+str(a2)

    ntrp=np.array(ntrp)
    cshockhom=np.array(cshockhom)
    pl.clf()
    pl.plot(cshockhom,ntrp,'o')
    pl.xlabel('Homogenized shock speed')
    pl.ylabel('Final Entropy/Initial Entropy')
