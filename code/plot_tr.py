from pyclaw.plotters.data import ClawPlotData
from pylab import hold, clf, plot

clf()
plotdata = ClawPlotData()
plotdata.outdir = './_output'

# Figure:
plotfigure = plotdata.new_plotfigure(name='Solution', figno=1)
plotfigure.kwargs = {'figsize':[5,3]}

# Axes:
plotaxes = plotfigure.new_plotaxes(name='Stress')
plotaxes.xlim = [30,55]

plotitem = plotaxes.new_plotitem(name='SharpClaw 3600', plot_type='1d')
plotitem.plot_var = 0       # q[2] is the stress
plotitem.plotstyle = '-'
plotitem.color = 'k'   # could use 'r' or 'red' or '[1,0,0]'
plotitem.kwargs = {'linewidth':3,'markersize':5}

#plotitem = plotaxes.new_plotitem(name='ClawPack 3600', plot_type='1d')
#plotitem.outdir = '/users/ketch/research/claw42/fwave2/tr28800'
#plotitem.plot_var = 2       # q[2] is the stress
#plotitem.plotstyle = '-'
#plotitem.color = 'k'
#plotitem.kwargs = {'linewidth':3,'markersize':5}

#plotitem = plotaxes.new_plotitem(name='ClawPack 28800', plot_type='1d')
#plotitem.outdir = '/users/ketch/research/claw42/fwave2/'
#plotitem.plot_var = 0       # q[2] is the stress
#plotitem.plotstyle = '-'
#plotitem.color = 'k'
#plotitem.kwargs = {'linewidth':3}


plotdata.plotframe(2)
hold(False)
#plotitem.plotstyle = 's'
#plotitem.color = 'b'   
plotitem.plotstyle = 'o'
plotitem.color = 'r' 
plotdata.plotframe(48)

q1=plotitem.getframe(2).q
q2=plotitem.getframe(48).q
xc=plotitem.getframe(2).grids[0].x.center
clf()
plot(xc,abs(q1[:,0]-q2[:,0]),'-k')
