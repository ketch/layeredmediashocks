def runtests(mx=7200,a2s=[0.1,0.2,0.4,0.8,1.2,1.6,2.0],rhos=[1.0,2.0,3.0,4.0,6.0,9.0],Ks=[1.0,2.0,3.0,4.0,6.0,9.0]):
    """
    Test time-reversed entropy change for a front.
    This is the test used for figure 3.2 of the paper.

    a2s is an array of initial condition amplitudes
    rhos,Ks are arrays of densities and bulk moduli

    Note that the script doesn't run all combinations from (rhos,Ks) but just
    all sequential pairs.
    """
    import setrun
    from pyclaw.runclaw import runclaw
    import numpy as np

    assert len(rhos)==len(Ks)

    reload(setrun)

    for rhoB,KB in zip(rhos,Ks):
        for a2 in a2s:
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
            rundata.probdata.ic=4

            #a2 is the magnitude of the stress behind the front
            rundata.probdata.a2=a2

            #We ought to rerun these with a larger domain and longer times
            #to see more clearly the entropy loss.
            rundata.probdata.trtime = 100.
            rundata.clawdata.tfinal = 200.
            rundata.clawdata.xupper=300.0

            #Reflecting at left, outflow at right
            rundata.clawdata.mthbc_xlower = 3
            rundata.clawdata.mthbc_xupper = 1
            rundata.clawdata.mx=mx

            rundata.write()
            runclaw(xclawcmd='xsclaw',outdir='./threshold/'+str(mx)+'_'+str(a2)+'_'+str(rhoB)+'_'+str(KB))


def compute_tr_entropy(mx=7200,a2s=[0.1,0.2,0.4,0.8,1.2,1.6,2.0,3.0,4.0,6.0,9.0],rhos=[1.0,2.0,3.0,4.0,6.0,9.0],Ks=[1.0,2.0,3.0,4.0,6.0,9.0]):
    """
    Compute entropy change for time-reversed front evolution problem (the function above).
    This is used to generate Figure 3.2 of the paper.
    """
    import pylab as pl
    from entropy import entropy_oneframe
    import numpy as np
    from cbar import cshock

    ntrp = []  #Storage for final entropy ratios
    cshockhom = []  #storage for effective shock speeds
    allrhos=[]
    allKs=[]
    alla2s=[]

    pl.clf()
    for rhoB in rhos:
        for KB in Ks:
            for a2 in a2s:
                sigmaL = a2/2.
                try:
                    if cshock(sigmaL,rhoB,KB)>2.5: raise Exception('This will go past the end of the domain')

                    ent0 = entropy_oneframe(1,outdir='./threshold/'+str(mx)+'_'+str(a2)+'_'+str(rhoB)+'_'+str(KB))
                    ent1 = entropy_oneframe(9,outdir='./threshold/'+str(mx)+'_'+str(a2)+'_'+str(rhoB)+'_'+str(KB))
                    ntrp.append(ent1/ent0)

                    cA = 1.#np.sqrt(KA/rhoA)
                    cB = np.sqrt(KB/rhoB)
                    cmeanlin = 2./(1./cA + 1./cB)

                    cshockhom.append(cshock(sigmaL,rhoB,KB)/cmeanlin)

                    # An alternative hypothesis:
                    #Kharm = 2./(1. + 1./KB)
                    #rhobar = 0.5*(1. + rhoB)
                    #cl = np.sqrt(Kharm*(1.+sigmaL) / rhobar)
                    #cr = np.sqrt(Kharm*(1.+0.) / rhobar)
                    #cshockhom.append(0.5*(cl+cr)/cmeanlin)

                    colorstring = str(1.-np.min([1,1.5*np.log(np.sqrt(KB*rhoB))/np.log(27.)]))
                    pl.plot(cshockhom[-1],ntrp[-1],'o',color=colorstring)
                    pl.draw()
                    pl.hold(True)
                    allrhos.append(rhoB)
                    allKs.append(KB)
                    alla2s.append(a2)
                except:
                    print 'no files found for K='+str(KB)+', rho='+str(rhoB)+' and a2='+str(a2)

    ntrp=np.array(ntrp)
    cshockhom=np.array(cshockhom)
    pl.plot([1.0,1.0],[0.7,1.1],'--k')
    pl.xlabel('Relative effective shock speed')
    pl.ylabel('Final Entropy/Initial Entropy')
    pl.hold(False)
    return ntrp,cshockhom,allrhos,allKs,alla2s
