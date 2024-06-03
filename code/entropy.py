def entropy_oneframe(frame,outdir='./_output',sigfun='exp'):
    """
    Compute the total entropy for a single output file.
    """
    import numpy as np
    from pyclaw.plotters.data import ClawPlotData
    #Get the material parameters
    aux = np.loadtxt(outdir+'/fort.aux')
    rho = aux[:,0]; K   = aux[:,1]

    plotdata = ClawPlotData()
    plotdata.outdir=outdir

    #Read in the solution
    dat = plotdata.getframe(frame)

    #Compute the entropy
    u = dat.q[1,:]
    nrg=rho*u**2/2.
    if sigfun=='exp':
      # For exponential relation sigma(eps) = e^(K*eps) - 1:
      # YOU MUST use the version of out1nel.f that OUTPUTS SIGMA!
      sigma = dat.q[0,:]
      sigint = (sigma-np.log(sigma+1.))/K 

    elif sigfun=='quad':
      # quadratic: sigma(eps) = K*eps + 0.5*eps**2:
      # YOU MUST use the version of out1nel.f that OUTPUTS EPSILON!
      eps = dat.q[0,:]
      sigint = K*eps/2. + eps**3/6.

    else:
      raise Exception('unrecognized stress-strain relation in entropy_oneframe()')

    return np.sum(nrg+sigint)*dat.grid.d[0]

def entropy(nstart,nend,outdir='./_output',sigfun='exp',doplot=True,linewidth=3,style='-k',firstlast=False):
    r"""Return total entropy values for a simulation, starting with frame
        nstart and ending with frame nend.  Optionally also produces a 
        plot.  The output is normalized so that the value at frame nstart is 1.
        In general, the entropy is given by
        $$\eta = 0.5 * \rho * u^2 + \int_0^epsilon sigma(s) ds.$$
    """
    import numpy as np
    import pylab as pl
    from pyclaw.plotters.data import ClawPlotData

    #Get the material parameters
    aux = np.loadtxt(outdir+'/fort.aux')
    rho = aux[:,0]; K   = aux[:,1]

    plotdata = ClawPlotData()
    plotdata.outdir=outdir
    totalentropy=np.zeros([nend-nstart+1])

    t=[]

    for i in range(nstart,nend+1):
        #Read in the solution
        dat = plotdata.getframe(i)
        t.append(dat.t)

        #Compute the entropy
        u = dat.q[:,1]
        nrg=rho*u**2/2.
        if sigfun=='exp':
          # For exponential relation sigma(eps) = e^(K*eps) - 1:
          # MUST modify out1nel.f to OUTPUT SIGMA!
          sigma = dat.q[:,0]
          sigint = (sigma-np.log(sigma+1.))/K 

        elif sigfun=='quad':
          # quadratic: sigma(eps) = K*eps + 0.5*eps**2:
          # MUST modify out1nel.f to OUTPUT EPS!
          eps = dat.q[:,0]
          sigint = K*eps/2. + eps**3/6.

        else:
          print 'unrecognized stress-strain relation'
          raise

        totalentropy[i-nstart]=np.sum(nrg+sigint)*dat.grid.d[0]
        print totalentropy[i-nstart]

    #Normalize:
    totalentropy=totalentropy/totalentropy[0]

    if doplot: 
        #pl.clf(); 
        pl.plot(t,totalentropy,style,linewidth=linewidth); 
        pl.draw()

    return totalentropy,t

def entropy_limiter_tests(grids=[1800,3600,7200,14400],lims=[1,2,3,4,5],outdir='./_output'):
    r"""Run a sequence of entropy tests.  This is used to produce
        Figure \ref{fig:entropy_lims} (Figure 2.5)."""
    import setrun
    from entropy import entropy
    import pylab as pl
    from pyclaw.runclaw import runclaw
    import numpy as np

    reload(setrun)
    ntrp = np.zeros([len(lims),len(grids)])

    for i,lim in enumerate(lims):
        if lim==5: 
            claw_pkg='sharpclaw'
            xclawcmd='xsclaw'
        else: 
            claw_pkg='classic'
            xclawcmd='xclaw'

        #Set up default run parameters
        rundata=setrun.setrun(claw_pkg)
        rundata.clawdata.mthlim = [lim,lim]

        #Set material parameters
        rundata.probdata.K_B=4.0
        rundata.probdata.rho_B=4.0
        rundata.probdata.K_A=1.0
        rundata.probdata.rho_A=1.0

        rundata.clawdata.tfinal=500.
        rundata.clawdata.xupper=150.
        rundata.clawdata.nout=10
        rundata.probdata.trtime=250.
        rundata.probdata.a2=0.
        rundata.probdata.ic=1
        rundata.probdata.a1=0.1
        rundata.probdata.t1=10.
        rundata.probdata.tw1=10.

        for j,mx in enumerate(grids):
            print 'Grid: ',mx,',  Limiter: ',lim

            rundata.clawdata.mx = mx
            rundata.write()
            runclaw(xclawcmd,outdir=outdir)

            ent,t=entropy(1,rundata.clawdata.nout-1,outdir=outdir)
            ntrp[i,j]=ent[-1]
            print ntrp


    #Now plot the results

    linestyles=['-k','--k','-.k',':k','-b']
    limnames=[lim_name(ilim) for ilim in lims]
    nlayers = 150. #Number of layers in domain

    pl.hold(False)
    for i,lim in enumerate(lims):
        pl.loglog(np.array(grids)/nlayers,abs(ntrp[i,:]-1.),linestyles[i],linewidth=3)
        pl.hold(True)

    pl.legend(limnames)

    for i,lim in enumerate(lims):
        for j,mx in enumerate(grids):
            if ntrp[i,j]>1: pl.plot(mx/150.,abs(ntrp[i,j]-1.),'+k',markersize=15)
            else: pl.plot(mx/150.,abs(ntrp[i,j]-1.),'ok',markersize=7)

    #pl.axis([4,150,1.e-4,10])
    pl.xlabel('Grid cells per layer',fontsize=20)
    pl.ylabel('$\Delta \eta/\eta_0$',fontsize=20)
    pl.title('Entropy Convergence',fontsize=20)
    pl.hold(False)
    pl.draw()

def lim_name(ilim):
    if ilim==1: return 'Minmod'
    elif ilim==2: return 'Superbee'
    elif ilim==3: return 'van Leer'
    elif ilim==4: return 'MC'
    elif ilim==5: return 'WENO5'


def entropy_Z_tests(Zs=[1.,2.,4.],mx=7200):
    """Run a sequence of entropy tests varying Z_B with $K_B=rho_B=Z_B$.
       This is used to generate figure \ref{fig:entropy_1} (Figure 2.4)."""
    import setrun
    import matplotlib.pyplot as pl
    from pyclaw.runclaw import runclaw
    import numpy as np

    reload(setrun)
    ntrp = np.zeros([len(Zs)])
    styles=['--k','.k','-k']

    for i,Z in enumerate(Zs):
        print Z
        rundata=setrun.setrun()
        rundata.probdata.K_B=Z
        rundata.probdata.rho_B=Z
        rundata.probdata.K_A=1.
        rundata.probdata.rho_A=1.

        rundata.probdata.a1 = 0.1
        rundata.probdata.ic = 1
        rundata.probdata.t1 = 10
        rundata.probdata.tw1 = 10
        rundata.probdata.trtime = 100000
        rundata.clawdata.tfinal = 500.
        rundata.clawdata.nout = 100
        rundata.clawdata.mx = mx
        rundata.clawdata.xupper = 150.
        rundata.write()
        runclaw(outdir='./_output')
            
        ent,t=entropy(rundata.clawdata.nout/10,rundata.clawdata.nout,doplot=True,linewidth=3,style=styles[i])
        pl.hold(True)

    pl.ylabel('$\eta(t)/\eta(0)$',fontsize=30)
    pl.xlabel('$time$',fontsize=30)
    pl.title('Entropy evolution',fontsize=25)
    pl.legend(['$Z_B =$'+str(Z) for Z in Zs],loc=3)
    pl.draw()
    pl.hold(False)

def entropy_K_and_rho_tests():
    """Run a sequence of entropy tests varying K and rho.
       The results of these tests didn't go into the paper but are quite interesting.
    """
       
    import setrun
    import pylab as pl
    from pyclaw.runclaw import runclaw
    import numpy as np

    reload(setrun)
    Ks=np.array([1,1.5,2,3,4,6,8,16,32,64])
    rhos=np.array([1,1.5,2,3,4,6,8,16,32,64])
    ntrp = np.zeros([len(Ks),len(rhos)])

    for i,K in enumerate(Ks):
        for j,rho in enumerate(rhos):
            rundata=setrun.setrun()
            rundata.probdata.K_B=K
            rundata.probdata.rho_B=rho
            rundata.probdata.K_A=1.
            rundata.probdata.rho_A=1.
            print K,rho

            #Rescale time and IC:
            cmean = 2*np.sqrt(K/((K+1.)*(rho+1.)))
            Kmean = 2*K/(K+1.)
            rhomean = 0.5*(rho+1.)
            rundata.probdata.a1 = rundata.probdata.a1/np.sqrt(rhomean)
            rundata.clawdata.tfinal = 500/(cmean/0.8)
            rundata.write()
            runclaw(outdir='./_output')
            
            ent,t=entropy(rundata.clawdata.nout)
            ntrp[i,j]=ent[-1]
            pl.pcolor(np.log2(Ks),np.log2(rhos),np.log10(abs(1.-ntrp)),cmap='gray')
            pl.colorbar(); pl.axis('image')
            pl.colorbar()
            pl.draw()

    pl.pcolor(np.log2(Ks),np.log2(rhos),np.log10(abs(1.-ntrp)),cmap='gray')
    pl.colorbar(); pl.axis('image')
    pl.xlabel('$\log_2(K_B)$',fontsize=20)
    pl.ylabel('$\log_2(p_B)$',fontsize=20)
    pl.title('$\log_{10}|\eta-\eta_0|$',fontsize=25)
    pl.draw()
