import pylab as pl

def test_tr(grids=[1800,3600,7200],Z_B=4.0,claw_pkg='sharpclaw',frame=3,outdir='./output'):
    r"""Run a sequence of time-reversal tests.
        This is used to generate the values in Table \ref{tbl:trconv} (Table 2.1)
        and the figures in Figures \ref{fig:stego_tr} and \ref{fig:stego_ntr}
        (Figures 2.2 and 2.3).
        To make the figures, use also the function plot_tr() below."""
    import setrun
    from pyclaw.runclaw import runclaw
    reload(setrun) #This way if we change anything it gets updated
    errs = []

    if   claw_pkg=='sharpclaw': xclawcmd='xsclaw'
    elif claw_pkg=='classic' : xclawcmd='xclaw'

    for i,mx in enumerate(grids):
        rundata=setrun.setrun(claw_pkg)

        #Set material parameters
        rundata.probdata.K_B=Z_B
        rundata.probdata.rho_B=Z_B
        rundata.probdata.K_A=1.0
        rundata.probdata.rho_A=1.0

        rundata.probdata.trtime=600.
        rundata.probdata.t1= 10.
        rundata.probdata.a1= 0.1
        rundata.probdata.tw1=10.

        rundata.clawdata.mx=mx
        rundata.clawdata.tfinal=1200.
        rundata.clawdata.nout=60
        rundata.write()
        runclaw(xclawcmd=xclawcmd,outdir=outdir)

        sig0,sig1,u0,u1,xc = get_timereversal_frames(frame)

        discrepancy = max(abs(u1-u0))
        errs.append(discrepancy)

    return grids,errs

def get_timereversal_frames(frame=3,nframes=60,outdir='./output'):
    """Get corresponding early and late time frames from a 
       time-reversal test.
    """
    from pyclaw.plotters.data import ClawPlotData

    plotdata = ClawPlotData()
    plotdata.outdir=outdir

    #"Exact" (early-time) solution
    dat0 = plotdata.getframe(frame)
    sig0 = dat0.q[0,:]
    u0   = dat0.q[1,:]

    #Time-reversed (late-time) solution
    dat1 = plotdata.getframe(nframes-frame)
    sig1 = dat1.q[0,:]
    u1   = dat1.q[1,:]

    xc = dat0.c_center
    print xc[0].shape

    return sig0,sig1,u0,u1,xc[0]

def plot_tr(frame=3,nframes=60):
    r"""Plot corresponding early- and late-time frames from a time-reversal
        test.  This is used to generate the plots in Figures 
        \ref{fig:stego_tr} and \ref{fig:stego_ntr}.
    """

    sig0,sig1,u0,u1,xc = get_timereversal_frames(frame,nframes)

    pl.plot(xc,sig0,'-k',linewidth=4)
    pl.hold(True)
    N=len(xc); skip=3
    sub=range(0,N,skip)
    pl.plot(xc[sub],sig1[sub],'bs',linewidth=4)
    pl.hold(False)
    pl.axis([30,60,-0.1,0.7])
    pl.xlabel('x',fontsize=30)
    pl.ylabel('$\sigma$',fontsize=30)
    pl.title('Stress at time 1140')
