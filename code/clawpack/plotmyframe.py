def plotmyframe(frame,m):
  from pyclaw.plotters.data import ClawPlotData
  import pylab as pl
  import numpy as np

  plotdata = ClawPlotData()
  plotdata.outdir='./_output'

  dat0 = plotdata.getframe(frame)
  u = np.squeeze(dat0.q[:,m])

  xc = np.squeeze(dat0.c_center)

  print np.shape(xc),np.shape(u)
  pl.plot(xc,u,'-k',linewidth=3)
  pl.xlabel('x',fontsize=20)
  pl.ylabel('$\sigma$',fontsize=20)

  return None


