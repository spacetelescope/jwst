from . import poly
import numpy as np 
import os
from astropy.io import fits
from astropy.table import Table
from scipy.interpolate import interp1d

class interp1d_picklable(object):
    """ class wrapper for piecewise linear function
    """

    def __init__(self, xi, yi, **kwargs):
        self.xi = xi
        self.yi = yi
        self.args = kwargs
        self.f = interp1d(xi, yi, **kwargs)

    def __call__(self, xnew):
        return self.f(xnew)

    def __getstate__(self):
        return self.xi, self.yi, self.args

    def __setstate__(self, state):
        self.f = interp1d(state[0], state[1], **state[2])

class Config(object):
    """Class to read and hold GRISM configuration info"""
    def __init__(self,GRISM_CONF,DIRFILTER=None):
        """Return a Config object
        
        Parameters
        ----------
        GRISM_CONF : str
            The full path and name to a grism configuration file
            
        DIRFILTER : str
            The name of the direct filter so that filter wedge offsets can be included.
            Should match the filter used when a direct image was used in the same visit as 
            the grism observations
            
        Returns
        -------
        C : Config class object
        """

        self.GRISM_CONF = open(GRISM_CONF).readlines()
        self.GRISM_CONF_PATH = os.path.dirname(GRISM_CONF)
        self.GRISM_CONF_FILE = os.path.basename(GRISM_CONF)

        self.orders = self._get_orders()
        self._DISPX_data = {}
        self._DISPY_data = {}
        self._DISPL_data = {}

        self._DISPX_polyname = {}
        self._DISPY_polyname = {}
        self._DISPL_polyname = {}

        self.SENS = {}
        self.SENS_data = {}

        # Wavelength range of the grism
        self.WRANGE = {}

        # Extent of FOV in detector pixel
        self.XRANGE = {}
        self.YRANGE = {}

        if DIRFILTER!=None:
            # We get the wedge offset values for this direct filter
            r = self._get_value("WEDGE_%s" % (DIRFILTER),type=float)
            self.wx = r[0]
            self.wy = r[1]
        else:
            self.wx = 0.
            self.wy = 0.

        # Get physical size of detector
        self.NAXIS = self._get_value("NAXIS",type=int)

        # Load the name of a POM file
        self.POM = None
        try:
            self.POM = os.path.join(self.GRISM_CONF_PATH,self._get_value("POM"))
        except:
            pass

        for order in self.orders:    
            self._DISPX_data[order] = self._get_parameters("DISPX",order)
            self._DISPY_data[order] = self._get_parameters("DISPY",order)
            self._DISPL_data[order] = self._get_parameters("DISPL",order)
            self.SENS[order] = self._get_sensitivity(order)
            
            self._DISPX_polyname[order] = np.shape(self._DISPX_data[order])
            self._DISPY_polyname[order] = np.shape(self._DISPY_data[order])
            self._DISPL_polyname[order] = np.shape(self._DISPL_data[order])

            self.SENS_data[order] = self._get_sensitivity(order)

            vg = self.SENS_data[order][1]>np.max(self.SENS_data[order][1])*1e-3
            wmin = np.min(self.SENS_data[order][0][vg])
            wmax = np.max(self.SENS_data[order][0][vg])
            self.WRANGE[order] = [wmin,wmax]

            self.SENS[order] = interp1d_picklable(self.SENS_data[order][0],self.SENS_data[order][1],bounds_error=False,fill_value=0.)
            

            self.XRANGE[order] = self._get_value("XRANGE_%s" % (order),type=float)
            self.YRANGE[order] = self._get_value("YRANGE_%s" % (order),type=float)

    @staticmethod
    def _rotate_coords(dx, dy, theta=0, origin=[0,0]):
        """Rotate cartesian coordinates CW about an origin
        
        Parameters
        ----------
        dx, dy : float or `~numpy.ndarray`
            x and y coordinages
                    
        theta : float
            CW rotation angle, in radians
            
        origin : [float,float]
            Origin about which to rotate
        
        Returns
        -------
        dxr, dyr : float or `~numpy.ndarray`
            Rotated versions of `dx` and `dy`
            
        """
        _mat = np.array([[np.cos(theta), -np.sin(theta)],
                     [np.sin(theta), np.cos(theta)]])

        rot = np.dot(np.array([dx-origin[0], dy-origin[1]]).T, _mat)
        dxr = rot[:,0]+origin[0]
        dyr = rot[:,1]+origin[1]
        return dxr, dyr

    def DISPL(self,order,x0,y0,t):
        """Returns the wavelength corresponding to a value t for an object at posittion x0,y0
        
        Parameters
        ----------
        x0, y0 : float or `~numpy.ndarray`
            x and y coordinates in the direct image
                    
        t : float
            Value of the t variable (usually 0<t<1)
        
        Returns
        -------
        wav : float or `~numpy.ndarray`
            wavelength value
            
        """
        return poly.POLY[self._DISPL_polyname[order]](self._DISPL_data[order],x0,y0,t)

    def DDISPL(self,order,x0,y0,t):
        """Returns the first derivate of the wavelength (wrt to t) for a value t for an object at posittion x0,y0
        
        Parameters
        ----------
        x0, y0 : float or `~numpy.ndarray`
            x and y coordinates in the direct image
                    
        t : float
            Value of the t variable (usually 0<t<1)
        
        Returns
        -------
        dwav : float or `~numpy.ndarray`
            First derivative of the wavelength with respect to 't', as a function of 't'
        """
        return poly.DPOLY[self._DISPL_polyname[order]](self._DISPL_data[order],x0,y0,t)

    def DISPXY(self, order, x0, y0, t, theta=0):
        """Return both `x` and `y` coordinates of a rotated trace
        
        Parameters
        ----------
        order : str
            Order string
            
        x0, y0 : float
            Reference position (i.e., in direct image)

        t : float or `~numpy.ndarray`
            Parameter where to evaluate the trace
            
        theta : float
            CW rotation angle, in radians
        
        Returns
        -------
        dxr, dyr : float or `~np.ndarray`
            Rotated trace coordinates as a function of `t`

        """
        dx = -self.wx + poly.POLY[self._DISPX_polyname[order]](self._DISPX_data[order],x0,y0,t)
        dy = -self.wy + poly.POLY[self._DISPY_polyname[order]](self._DISPY_data[order],x0,y0,t)

        if theta != 0:
            dxr, dyr = self._rotate_coords(dx, dy, theta=theta, origin=[0,0])
            return dxr, dyr
        else:
            return dx, dy
            
    
    def INVDISPXY(self, order, x0, y0, dx=None, dy=None, theta=0, t0=np.linspace(0,1,10)):
        """Return independent variable `t` along rotated trace
        
        Parameters
        ----------
        order : str
            Order string
            
        x0, y0 : float
            Reference position (i.e., in direct image)

        dx : float, `~numpy.ndarray` or None
            `x` coordinate in *rotated* trace where to evaluate the trace
            independent variable `t`.
            
        dy : float, `~numpy.ndarray` or None
            Same as `dx` but evaluate along 'y' axis.
        
        t0 : `~np.ndarray`
            Independent variable location where to evaluate the rotated trace.
            For low-order trace shapes, this can be coarsely sampled as 
            in the default.
            
        Returns
        -------
        tr : float or `~np.ndarray`
            Independent variable `t` evaluated on the rotated trace at
            `dx` or `dy`.

        .. note::
        
        Order of execution is first check if `dx` supplied.  If not, then
        check `dy`.  And if both are None, then return None (do nothing).
        
        """
        if dx is not None:
            xr, yr = self.DISPXY(order, x0, y0, t0, theta=theta)
            so = np.argsort(xr)
            f = interp1d_picklable(xr[so], t0[so])
            tr = f(dx)
            return tr
        
        if dy is not None:
            xr, yr = self.DISPXY(order, x0, y0, t0, theta=theta)
            so = np.argsort(yr)
            f = interp1d_picklable(yr[so], t0[so])
            tr = f(dy)
            return tr
                 
        return None
        
    def DISPX(self,order,x0,y0,t):
        """Returns the x offset x'-x = DISPL(x0,y0,t) where x0,y0 is the 
        position on the detector, x'-x is the difference between direct and grism image x-coordinates and 0<t<1
        
        Parameters
        ----------
        order : str
            Order string
            
        x0, y0 : float
            Reference position (i.e., in direct image)

        t : float or `~numpy.ndarray`
            Parameter where to evaluate the trace
        
        Returns
        -------
        dx : float or `~np.ndarray`
            Trace x-coordinates as a function of `t`

        """
        dx = -self.wx + poly.POLY[self._DISPX_polyname[order]](self._DISPX_data[order],x0,y0,t)
            
        return  dx

    def DDISPX(self,order,x0,y0,t):
        """Returns the first derivative of the x offset (x'-x) wrt to t, where x0,y0 is the 
        position on the detector, x'-x is the difference between direct and grism image x-coordinates and 0<t<1
        
        Parameters
        ----------
        order : str
            Order string
            
        x0, y0 : float
            Reference position (i.e., in direct image)

        t : float or `~numpy.ndarray`
            Parameter where to evaluate the trace
        
        Returns
        -------
        dxdt : float or `~np.ndarray`
            First derivative of the trace x-coordinates with respect to 't', as a function of `t`

        """
        return  poly.DPOLY[self._DISPX_polyname[order]](self._DISPX_data[order],x0,y0,t)

    def DISPY(self,order,x0,y0,t):
        """Returns the x offset (y'-y) wrt to t, where x0,y0 is the 
        position on the detector, y'-y is the difference between direct and grism image y-coordinates and 0<t<1
        
        Parameters
        ----------
        order : str
            Order string
            
        x0, y0 : float
            Reference position (i.e., in direct image)

        t : float or `~numpy.ndarray`
            Parameter where to evaluate the trace
        
        Returns
        -------
        dydt : float or `~np.ndarray`
            First derivative of the trace y-coordinates with respect to 't', as a function of `t`

        """
        return  -self.wy + poly.POLY[self._DISPY_polyname[order]](self._DISPY_data[order],x0,y0,t)

    def DDISPY(self,order,x0,y0,t):
        """Returns the first derivative of the y offset (y'-y) wrt to t, where x0,y0 is the 
        position on the detector, y'-y is the difference between direct and grism image x-coordinates and 0<t<1
        
        Parameters
        ----------
        order : str
            Order string
            
        x0, y0 : float
            Reference position (i.e., in direct image)

        t : float or `~numpy.ndarray`
            Parameter where to evaluate the trace
        
        Returns
        -------
        dydt : float or `~np.ndarray`
            First derivative of the trace y-coordinates with respect to 't', as a function of `t`

        """
        return poly.DPOLY[self._DISPY_polyname[order]](self._DISPY_data[order],x0,y0,t)

    def INVDISPL(self,order,x0,y0,l):
        """Returns the value of 't' that corresponds to a given wavelength for a source at position x0,y0
        
        Parameters
        ----------
        order : str
            Order string
            
        x0, y0 : float
            Reference position (i.e., in direct image)

        l : float or `~numpy.ndarray`
            Wavelength
        
        Returns
        -------
        t : float or `~np.ndarray`
            `t` value

        """
        return poly.INVPOLY[self._DISPL_polyname[order]](self._DISPL_data[order],x0,y0,l)

    def INVDISPX(self,order,x0,y0,dx,t0=np.linspace(0,1,10)):
        """Returns the value of 't' that corresponds to a given x-offset for a source at position x0,y0
        
        Parameters
        ----------
        order : str
            Order string
            
        x0, y0 : float
            Reference position (i.e., in direct image)

        dx: float or `~numpy.ndarray`
            x-offset between source and a given pixel
        t0 : `~np.ndarray`
            Independent variable location where to evaluate the trace.
            For low-order trace shapes, this can be coarsely sampled as 
            in the default.
        Returns
        -------
        t : float or `~np.ndarray`
            `t` value

        """
        if self._DISPX_polyname[order] in poly.INVPOLY.keys():
            return poly.INVPOLY[self._DISPX_polyname[order]](self._DISPX_data[order],x0,y0,dx+self.wx)
        else:
            xr = self.DISPX(order, x0, y0, t0)
            so = np.argsort(xr)
            f = interp1d_picklable(xr[so], t0[so])
            tr = f(dx)
            return tr

        

    def INVDISPY(self,order,x0,y0,dy,t0=np.linspace(0,1,10)):
        """Returns the value of 't' that corresponds to a given y-offset for a source at position x0,y0
        
        Parameters
        ----------
        order : str
            Order string
            
        x0, y0 : float
            Reference position (i.e., in direct image)

        dy: float or `~numpy.ndarray`
            y-offset between source and a given pixel
        t0 : `~np.ndarray`
            Independent variable location where to evaluate the trace.
            For low-order trace shapes, this can be coarsely sampled as 
            in the default.
        Returns
        -------
        t : float or `~np.ndarray`
            `t` value

        """        
        if self._DISPY_polyname[order] in poly.INVPOLY.keys():
            return poly.INVPOLY[self._DISPY_polyname[order]](self._DISPY_data[order],x0,y0,dy+self.wy)
        else:      
            xr, yr = self.DISPXY(order, x0, y0, t0)
            so = np.argsort(yr)
            f = interp1d_picklable(yr[so], t0[so])
            tr = f(dy)
            return tr

    def _get_orders(self):
        """Returns all the know orders in Config
        
        Parameters
        ----------
        
        Returns
        -------
        orders: `array`
            List of orders

        """        

        orders = []

        # Get orders 
        for l in self.GRISM_CONF:
            k = "BEAM_"
            if l[0:len(k)]==k:
                ws = l.split()
                order = ws[0].split("_")[-1]
                orders.append(order)
        return orders

    def _get_parameters(self,name,order,str_fmt="%s_%s_"):
        """Return the 2D polynomial array stored in the config file"""
        str = str_fmt % (name,order)
        # Find out how many we have to store
        n = 0
        m = 0
        for l in self.GRISM_CONF:
            if l[0]=="#": continue
            ws = l.split()
            if len(ws)>0 and str in ws[0]:
                i = ws[0].split(str)[-1]
                n = n + 1
                m = len(ws)-1

        arr = np.zeros((n,m))

        for l in self.GRISM_CONF:
            ws = l.split()
            if len(ws)>0 and str in ws[0]:
                i = int(ws[0].split(str)[-1])
                if len(ws)-1 !=m:
                    print("ERROR: Wrong format for ",GRISM_CONF,name,order)
                    return None
                vals = [float(ww) for ww in ws[1:]]
                arr[i,0:m] = vals

        return arr            

    def _get_value(self,str,type=None):
        """Helper function to simply return the value for a simple keyword parameters
        in the config file."""
        
        for l in self.GRISM_CONF:
            ws = l.split()
            if len(ws)>0 and ws[0]==str:
                if len(ws)==2:
                    if type==None:
                        return ws[1]
                    elif type==float:
                        return float(ws[1])
                    elif type==int:
                        return int(ws[1])
                else:
                    if type==None:
                        return ws[1:]
                    elif type==float:
                        return [float(x) for x in ws[1:]]
                    elif type==int:
                        return [int(x) for x in ws[1:]]
        return None


    def _get_sensitivity(self,order):
        """Helper function that looks for the name of the sensitivity file, reads it and
        stores the content in a simple list [WAVELENGTH, SENSITIVITY]."""
        fname = os.path.join(self.GRISM_CONF_PATH,self._get_value("SENSITIVITY_%s" % (order)))
        fin = fits.open(fname)
        wavs = fin[1].data.field("WAVELENGTH")[:]
        sens = fin[1].data.field("SENSITIVITY")[:]
        fin.close()        
        # Fix for cases where sensitivity is not zero on edges
        sens[0:2] = 0.
        sens[-2:] = 0.
                
        return [wavs,sens]    
