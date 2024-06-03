
""" 
Set up the plot figures, axes, and items to be done for each frame.

This module is imported by the plotting routines and then the
function setplot is called to set the plot parameters.
    
""" 

#--------------------------
def setplot(plotdata):
#--------------------------
    
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of pyclaw.plotters.data.ClawPlotData.
    Output: a modified version of plotdata.
    
    """ 


    plotdata.clearfigures()  # clear any old figures,axes,items data

#    plotdata.afterframe=zoom_to_shock

    # Figure for q[0]
    plotfigure = plotdata.new_plotfigure(name='Stress', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'#[0.,150.]
    plotaxes.ylimits = 'auto'#[-22.1,12.4]
    plotaxes.title = 'Stress'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 0
    plotitem.plotstyle = '-o'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':2,'markersize':5}
    


    # Figure for q[1]
    plotfigure = plotdata.new_plotfigure(name='Velocity', figno=2)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = 'auto'
    plotaxes.ylimits = 'auto'#[-1.1,0.5]
    plotaxes.title = 'Velocity'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d_plot')
    plotitem.plot_var = 1
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.show = True       # show on plot?
    plotitem.kwargs = {'linewidth':3,'markersize':5}
    

    # Parameters used only when creating html and/or latex hardcopy
    # e.g., via pyclaw.plotters.frametools.printframes:

    plotdata.printfigs = True                # print figures
    plotdata.print_format = 'png'            # file format
    plotdata.print_framenos = 'all'          # list of frames to print
    plotdata.print_fignos = 'all'            # list of figures to print
    plotdata.html = True                     # create html files of plots?
    plotdata.html_homelink = '../README.html'
    plotdata.latex = True                    # create latex file of plots?
    plotdata.latex_figsperline = 2           # layout of plots
    plotdata.latex_framesperline = 1         # layout of plots
    plotdata.latex_makepdf = False           # also run pdflatex?

    return plotdata

def zoom_to_shock(current_data):
    import matplotlib.pyplot as pl
    pl.xlim([current_data.t*1.192+28,current_data.t*1.192+32])

def dif(current_data):
    import numpy as np
    d=np.empty(np.size(current_data.q[0,:]))
    d[1:] = np.diff(current_data.q[0,:])/current_data.grid.d[0]
    return d

def w1(current_data):
    """First Riemann invariant"""
    import numpy as np
    sig = current_data.q[0,:]
    u   = current_data.q[1,:]
    rho=1.5
    K=4./3.
    return rho*(u-2./np.sqrt(rho*K) * np.sqrt(sig+1))

def w2(current_data):
    """Second Riemann invariant"""
    import numpy as np
    sig = current_data.q[0,:]
    u   = current_data.q[1,:]
    rho=1.5
    K=4./3.
    return rho*(u+2./np.sqrt(rho*K) * np.sqrt(sig+1))
