import numpy as np
from scipy import interpolate, signal, stats
from astropy.convolution import Gaussian1DKernel, Box1DKernel, Gaussian2DKernel
from scipy.signal import convolve as scipy_convolve
from astropy.convolution import convolve as astropy_convolve

def pmetrics(x, y, z, MWL=-0.2, MHW=1.28, zrange=0.3):
    """
    Find points on beach profile
    Assumes x starts somewhere offshore and is positive landward

    Input:
       x - distance on cross-shore profile from arbitrary point (array; m)
       y - distance alongshore (float, m) - used to customize algorithm for particular stretch of beach
       z - elevation on cross-shore profile at x locations(array same size as x, datum as MWL and MHW...normally NAVD88; m)
       MWL - elevation of mean water level (m)
       MHW - elevation of mean high water level (m)
    """
    # dict to hold result
    m = {}
    m['phi'] = np.NaN
    m['xphi'] = np.NaN
    m['dhi'] = np.NaN
    m['xdhi'] = np.NaN
    m['xMHW']=np.NaN
    m['xMWL']=np.NaN
    m['dtoe']=np.NaN
    m['xdtoe']=np.NaN

    # hand-crafted max. beach widths specific to Sandwich
    bwmx = 100.
    if(y>600. and y<=1000.):
        bwmx = 50.
    if(y>1000):
        bwmx = 20.

    # determine spacing of profile points (assume they are uniform)
    try:
        dx = np.median(np.diff(x))
    except:
        print("Warning: problem calculating dx")
    if(np.isnan(dx)):
        print("Warning: dx = ",dx," changing to dx = 1")
        dx = 1.

    print("dx = ",dx)


    # smoothed profile
    #  size of Gauss kernal (meters)
    gkn = int(3/dx)
    #  amount to remove at ends
    gkno2 = int(round(gkn/2.))
    gauss_kernel = Gaussian1DKernel(gkn)
    zf = convolve(np.array(z), gauss_kernel)

    # derivatives of unsmoothed profile
    zip = np.array((0))
    dz =  np.diff( np.append(zip, z) )
    ddz = np.diff( np.append(zip, dz) )

    # derivative of smoothed profile
    dzf = np.diff( np.append(zip, zf))

    # eliminate big jumps
    dz[np.abs(ddz)>20.]=np.NaN
    ddz[np.abs(ddz)>20.]=np.NaN
    dzf[np.abs(dzf)>20.]=np.NaN

    # find highest point on profile: phi
    # max of z
    iphi = np.nanargmax(z)
    m['phi'] = z[iphi]
    m['xphi'] = x[iphi]

    if(m['phi']>MWL):
        # we can then call dhi = phi
        m['dhi'] = z[iphi]
        m['xdhi'] = x[iphi]

        # find MWL by fitting line to nearby points
        idx = np.argwhere(np.logical_and(np.abs(zf-MWL)<zrange, x<175))
        if len(idx) >= 3:
            slope, intercept, r_value, p_value, stderr = \
               stats.linregress( np.squeeze(x[idx]),np.squeeze(z[idx]) )
            #print("slope, intercept, r:",slope, intercept, r_value)
            xMWL = (MWL-intercept)/slope
            if np.any(x<=xMWL):
                iMWL = np.argwhere(x<=xMWL)[-1]
                m['xMWL']=xMWL

        # find MHW the same way
        if(m['phi']>MHW):
            idx = np.argwhere(np.logical_and(np.abs(zf-MHW)<zrange, x<175))
            if len(idx) >= 3:
                slope, intercept, r_value, p_value, stderr = \
                   stats.linregress( np.squeeze(x[idx]),np.squeeze(z[idx]) )
                #print("slope, intercept, r:",slope, intercept, r_value)
                xMHW = (MHW-intercept)/slope
                if np.any(x<=xMHW):
                    iMHW = np.asscalar(np.argwhere(x<=xMHW)[-1])
                    m['xMHW']=xMHW

                    # dune toe = max ddz between MHW and 30 m
                    idtoe = iMHW+np.nanargmax(ddz[iMHW:iMHW+int(30/dx)])-2
                    m['dtoe'] = z[idtoe]
                    m['xdtoe'] = x[idtoe]

    m['bslope']=np.nan
    if(~np.isnan(m['xMHW']) and ~np.isnan(m['xMWL'])):
        # beach slope between MWL and MHW
        m['bslope'] = (MHW-MWL)/(xMHW-xMWL)

    return m
