def compmyframe(frame,m,outdirs):
  from pyclaw.plotters.data import ClawPlotData
  import pylab as pl
  import numpy as np

  plotdata = ClawPlotData()
  styles=['o','-k','x']

  for i,dir in enumerate(outdirs):
    plotdata.outdir=dir

    dat0 = plotdata.getframe(frame)
    u = np.squeeze(dat0.q[:,m])

    xc = np.squeeze(dat0.c_center)

    pl.plot(xc,u,styles[i],linewidth=2)
    pl.hold(True)

  pl.xlabel('x',fontsize=20)
  pl.ylabel('$\sigma$',fontsize=20)
  pl.hold(False)

  return None


