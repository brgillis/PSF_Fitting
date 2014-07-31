import numpy as np
import pylab

def draw_whisker_shear(g1, g2, x, y, size=None, title=None, filename=None, magnify=60, adjust_power=1):
    """Draw shear whisker plots from an array of g1, g2 values.  Write to file `outfile`.
    
       Requires: g1 [<float>]
                 g2 [<float>]
                 x [<float>]
                 y [<float>] (arrays of shears and positions)
       Optional: size [<float>] (a size measurement for each object)
                 title <string> (title for the plot)
                 filename <string> (filename to output to. If not specified, will display in window)
                 magnify <float> (magnification factor to apply to size of bars, for manual tuning)
                 adjust_power <float> (for fine-tuning the relative sizes of bars. For instance, 0.5
                               will plot based on the square-root of the shear. The length of the
                               median shear will not be affected by adjustments of this, so if the
                               median length is fine, it is not necessary to simultaneously adjust
                               this value and magnify in order to tune the display sizes.)
                               
       Returns: (nothing)
       
       Side-effects: Overwrites filename if specified, otherwise opens a window to display plot
       
       Except: Will raise an exception if:
                   filename is passed and cannot be written to
                   
    """
    
    g = np.power(np.add(np.power(g1,2),np.power(g2,2)),0.5*adjust_power)
    mag = magnify * np.power(np.median(g),-adjust_power)
    
    theta = 0.5*np.arctan2(g2, g1)
    gx = np.multiply(g,np.cos(theta))
    gy = np.multiply(g,np.sin(theta))
    pylab.figure()
    unused_res=pylab.quiver(x, y, mag*gx, mag*gy, size, scale=16, headwidth=0, pivot='middle',
                     units='width')
    pylab.colorbar()
    unused_str = str(np.median(g))
    if title is not None:
        if adjust_power==0:
            titlestr = ""
        else:
            titlestr = ", median shear = {0:.3}".format(np.power(np.median(g),1./adjust_power))
        pylab.title(title+titlestr)
    pylab.xlabel('X position [pixels]')
    pylab.ylabel('Y position [pixels]')
    if filename is None:
        pylab.show()
    else:
        pylab.savefig(filename)