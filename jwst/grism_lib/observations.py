import os
import numpy as np
from astropy.io import fits
from scipy import sparse
from scipy.interpolate import interp1d
from jwst.grism_lib import grismconf
from jwst.grism_lib.disperse import dispersed_pixel

import logging

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


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


def helper(vars):
    # parse the input list of vars; ID is a dummy number here
    x0s, y0s, flux, order, config, ID, extrapolate_SED, xoffset, yoffset = vars

    # use the list of vars to compute dispersed attributes for the pixel
    p = dispersed_pixel(x0s, y0s, flux, order, config, ID,
                        extrapolate_SED=extrapolate_SED, xoffset=xoffset, yoffset=yoffset)

    # unpack the results
    xs, ys, areas, lams, counts, ID = p
    IDs = [ID] * len(xs)

    # return results as numpy array
    pp = np.array([xs, ys, areas, lams, counts, IDs])
    return pp


class observation():
    # This class defines an actual observation. It is tied to a single image and a single config file
    
    def __init__(self, direct_images, segmentation_data, config, mod="A", order="+1",
                 plot=0, max_split=100, SED_file=None, extrapolate_SED=False, max_cpu=1,
                 ID=0, SBE_save=None, boundaries=[], renormalize=True):

        """
        direct_images: List of file name(s) containing direct imaging data
        segmentation_data: an array of the size of the direct images, containing 0 and 1's, 0 being pixels to ignore
        config: The path and name of a GRISMCONF NIRCAM configuration file
        mod: Module, A or B
        order: The name of the spectral order to simulate, +1 or +2 for NIRCAM
        max_split: Number of chunks to compute instead of trying everything at once.
        SED_file: Name of HDF5 file containing datasets matching the ID in the segmentation file
                  and each consisting of a [[lambda],[flux]] array.
        SBE_save: If set to a path, HDF5 containing simulated stamps for all obsjects will be saved.
        boundaries: a tuple containing the coordinates of the FOV within the larger seed image.
        """

        # This loads all the info from the configuration file for this grism mode;
        # things like the name of the POM mask file, the dispersion coeffs, etc.
        self.Cfg = grismconf.Config(config)
        self.ID = ID
        self.IDs = []
        self.dir_image_names = direct_images
        self.seg = segmentation_data
        self.dims = np.shape(self.seg)
        self.order = order
        self.SED_file = SED_file   # should always be NONE for baseline pipeline (use flat SED)
        self.SBE_save = SBE_save   # may want to have a new product in which simulated spectra are saved (NOT HDF5 format!)
        self.max_cpu = max_cpu
        self.cache = False
        self.POM_mask = None
        self.POM_mask01 = None
        self.renormalize = renormalize

        # Current NIRCAM config files DO have POM defined, so this branch is used
        if self.Cfg.POM is not None:
            #print("Using POM mask",self.Cfg.POM)
            log.info("Using POM mask", self.Cfg.POM)
            with fits.open(self.Cfg.POM) as fin:
                self.POM_mask = fin[1].data
                self.POM_mask01 = fin[1].data * 1. # A version of the mask that only contains 0-1 values. Pixels with labels/values>10000 are set to 1.
                self.POM_mask01[self.POM_mask01>=10000] = 1.
                self.Pxstart = int(fin[1].header["NOMXSTRT"])
                self.Pxend = int(fin[1].header["NOMXEND"])
                self.Pystart = int(fin[1].header["NOMYSTRT"])
                self.Pyend = int(fin[1].header["NOMYEND"])

                if len(fin)>2:
                    #print("Loading extra optical element transmission curves from POM file")
                    log.info("Loading extra optical element transmission curves from POM file")
                    self.POM_transmission = {}
                    for i in range(2,len(fin)):
                        w = fin[i].data["Wavelength"]
                        t = fin[i].data["Transmission"]
                        indx = int(fin[i].header["POMINDX"])
                        self.POM_transmission[indx] = interp1d_picklable(w,t,bounds_error=False,fill_value=0.)

        if len(boundaries)!=4:           
            xpad = (np.shape(segmentation_data)[1]-self.Cfg.NAXIS[0])//2
            ypad = (np.shape(segmentation_data)[0]-self.Cfg.NAXIS[1])//2
            self.xstart = 0 + xpad
            self.xend = xpad + self.Cfg.NAXIS[0]-1
            self.ystart = 0 + ypad
            self.yend = ypad + self.Cfg.NAXIS[1]-1
            #print("No boundaries passed. Assuming symmetrical padding of {} {} pixels and a final size of {} {} .".format(xpad,ypad,self.xend+1-self.xstart,self.yend+1-self.ystart))
            log.info("No boundaries passed. Assuming symmetrical padding of {} {} pixels and a final size of {} {} .".format(xpad,ypad,self.xend+1-self.xstart,self.yend+1-self.ystart))
        else:
            self.xstart,self.xend,self.ystart,self.yend = boundaries

        
        self.extrapolate_SED = extrapolate_SED # Allow for SED extrapolation
        if self.extrapolate_SED:
            #print("Warning: SED Extrapolation turned on.")
            log.info("Warning: SED Extrapolation turned on.")

        self.apply_POM()
        self.create_pixel_list()   # this will create pixel lists for ALL sources in segmentation map
        
        self.p_l = []
        self.p_a = []


    def apply_POM(self):
        """Account for the finite size of the POM and remove pixels in segmentation files which should not
        be dispersed.
        If a POM mask array was loaded then we make sure it is modified to have the same size as the input seg map
        and with the detector FOV starting at the same pixel locations as the input seg map."""

        if self.POM_mask is None:
            x0 = int(self.xstart+self.Cfg.XRANGE[self.Cfg.orders[0]][0] + 0.5)
            x1 = int(self.xend+self.Cfg.XRANGE[self.Cfg.orders[0]][1] + 0.5)
            y0 = int(self.ystart+self.Cfg.YRANGE[self.Cfg.orders[0]][0] + 0.5)
            y1 = int(self.yend+self.Cfg.YRANGE[self.Cfg.orders[0]][1] + 0.5)

            if x0<0: x0 = 0
            if y0<0: y0 = 0
            ymax,xmax = np.shape(self.seg)
            if x0+1>xmax: x0 = xmax-1
            if y0+1>ymax: y0 = ymax-1

            from astropy.io import fits
            print("POM footprint applied: [{}:{},{}:{}]".format(x0,x1,y0,y1))
            print("Pixels outside of this region of the input data ([{},{}]) will not be dispersed.".format(xmax,ymax))
            self.seg[:,:x0+1] = 0
            self.seg[:,x1:] = 0
            self.seg[:y0+1,:] = 0
            self.seg[y1:,:] = 0
        else:
            self.Pxstart - self.xstart
            self.Pystart - self.ystart

            ys1,xs1 = np.shape(self.POM_mask)
            ys2,xs2 = np.shape(self.seg)

            #print("xstart ettc:",self.xstart,self.xend,self.ystart,self.yend)
            #print("POM start etc:",self.Pxstart ,self.Pxend,self.Pystart ,self.Pyend)
            #print(self.Pystart - self.ystart,self.yend-self.ystart,self.Pxstart - self.xstart+1,self.xend-self.xstart+1)
            if (xs1>xs2) or (ys1>ys2):
                print("Warning: The input seg map is smaller ({},{}) than the POM mask ({},{}) for this mode.".format(xs2,ys2,xs1,ys1))
                raise Exception("Invalid seg map size.")

            # Crate a new POM mask that is the size of the input seg and data arrays
            # Compute the offset between where the detector is in the original POM mask and the input sef file
            xoff = self.xstart - self.Pxstart
            yoff = self.ystart - self.Pystart

            if (xoff<0) or (yoff<0):
                print("Warning: the input seg map and POM mask are not compatible. The seg map needs a larger border around the detector FOV than the POM mask.")
                print("seg map detector FOV starts at {} {}".format(self.xstart,self.ystart))
                print("POM mask detector FOV starts at {} {}".format(self.Pxstart,self.Pystart))
                raise Exception("Invalid seg map size.")
 
            POM = np.zeros(np.shape(self.seg))
            POM[yoff:yoff+ys1,xoff:xoff+xs1] = self.POM_mask
            self.POM_mask = POM*1

            # POM = np.zeros(np.shape(self.seg))
            # POM[yoff:yoff+ys1,xoff:xoff+xs1] = self.POM_mask
            POM[POM>1] = 1
            self.POM_mask01 = POM*1

            self.Pxstart  = self.xstart
            self.Pxend  = self.xend
            self.Pystart  = self.ystart
            self.Pyend  = self.ystart

            # We apply the POM mask to the seg file, but only as a mask, romving pixels getting no lights. We keep the POM mask with its orginal values otherwise 
            # to be able to account for partial transmission later on.
    

            from astropy.io import fits
            #print("POM size:",np.shape(POM))
            #print("seg size:",np.shape(self.seg))
            #fits.writeto("seg_org.fits",self.seg,overwrite=True)
            self.seg = self.seg * self.POM_mask01 
            #fits.writeto("POM_msk.fits",self.POM_mask,overwrite=True)
            #fits.writeto("POM.fits",self.POM_mask01 ,overwrite=True)
            #fits.writeto("seg_msk.fits",self.seg,overwrite=True)

        #sys.exit(1)


    def create_pixel_list(self):
        # Create a list of pixels to dispersed, grouped per object ID

        # In the examples this is called with ID=0, so processing is
        # done for all sources in the segmentation map ("seg" attribute)
        if self.ID==0:
            # This creates a huge list of all x,y pixel indices that
            # have a non-zero value in the seg map, sorted by those
            # indices belonging to a particular source ID
            self.xs = []
            self.ys = []
            all_IDs = np.array(list(set(np.ravel(self.seg))))
            all_IDs = all_IDs[all_IDs>0]
            #print("We have ",len(all_IDs),"Objects")
            log.info("We have ",len(all_IDs),"Objects")
            for ID in all_IDs:
                ys,xs = np.nonzero(self.seg==ID)

                if (len(xs)>0) & (len(ys)>0):
                    self.xs.append(xs)
                    self.ys.append(ys)
                    self.IDs = all_IDs
        else:
            vg = self.seg==self.ID
            ys,xs = np.nonzero(vg)            
           
            if (len(xs)>0) & (len(ys)>0):    
                self.xs.append(xs)
                self.ys.append(ys)
                self.IDs = [self.ID]

        self.fs = {}
        for dir_image_name in self.dir_image_names:
            #print("dir image:",dir_image_name)
            log.info("dir image:",dir_image_name)
            if self.SED_file==None:
                # Pipeline will use SED_file=None, so this code will try to read
                # photometry keyword values from the direct images.
                # Unfortunately it's trying to read the old HST-style keywords
                # PHOTPLAM (pivot wavelength) and PHOTFLAM, which DO NOT EXIST
                # for JWST products. We have PHOTMJSR and PHOTUJA2, so this will
                # need to be modified in some way to comply.
                try:
                    pivlam = fits.getval(dir_image_name,'PHOTPLAM') / 10000. # in Angsrrom and we want Micron now
                except:
                    print("WARNING: unable to find PHOTPLAM keyword in {}".format(dir_image_name))
                    sys.exit()

                try:
                    photflam = fits.getval(dir_image_name,'photflam')
                except:
                    print("WARNING: unable to find PHOTFLAM keyword in {}".format(dir_image_name))
                    sys.exit()
                #print("Loaded",dir_image_name, "wavelength:",pivlam,"micron")
                log.info("Loaded",dir_image_name, "wavelength:",pivlam,"micron")
            try:
                dimage = fits.open(dir_image_name)[1].data
            except:
                dimage = fits.open(dir_image_name)[0].data

            # If we do not use an SED file then we use photometry to get fluxes
            # Otherwise, we assume that objects are normalized to 1.
            if self.SED_file==None:
                # This next line uses the PHOTPLAM wavelength value, but we don't have that
                # for JWST images. Will need to figure out wavelength some other way, e.g.
                # by simply parsing from the filter name.
                self.fs[pivlam] = []
                dnew = dimage
                if self.POM_mask01 is not None:
                    dnew = dimage * self.POM_mask01  # Apply POM transmission mask to the data pixels
                for i in range(len(self.IDs)):
                    # This multiplies direct image pixels for each source, which apparently are
                    # still in units of countrate, by PHOTFLAM, to convert them to fluxes.
                    # So "fs" is the list of fluxes for all source pixels, at a given wavelength "pivlam"
                    self.fs[pivlam].append(dnew[self.ys[i],self.xs[i]] * photflam)
            else:
                # Need to normalize the object stamps              
                for ID in self.IDs:
                    vg = self.seg==ID
                    dnew = d
                    if self.POM_mask01 is not None:
                        #print("Applying POM transmission to data")
                        dnew = d * self.POM_mask01 # Apply POM transmission mask to the data pixels. This is a single grey correction for the whole object.
                    if self.renormalize is True:
                        sum_seg = np.sum(d[vg]) # But normalize by the whole flux
                        if sum_seg!=0.:
                            d[vg] = d[vg]/sum_seg
                    else:
                        print("not renormlazing sources to unity")

                self.fs["SED"] = []
                for i in range(len(self.IDs)):
                    self.fs["SED"].append(dnew[self.ys[i],self.xs[i]])
    

    def disperse_all(self,cache=False):

        if cache:
            print("Object caching ON")
            self.cache = True
            self.cached_object = {}

        self.simulated_image = np.zeros(self.dims,np.float)

        # This tqdm stuff just creates a progress bar while
        # the remaining code is running
        #for i in tqdm.tqdm(range(len(self.IDs))):
        for i in range(len(self.IDs)):
            if self.cache:
                self.cached_object[i] = {}
                self.cached_object[i]['x'] = []
                self.cached_object[i]['y'] = []
                self.cached_object[i]['f'] = []
                self.cached_object[i]['w'] = []
                self.cached_object[i]['minx'] = []
                self.cached_object[i]['maxx'] = []
                self.cached_object[i]['miny'] = []
                self.cached_object[i]['maxy'] = []

            # Call disperse for this object "i"
            this_object = self.disperse_chunk(i)

            if self.SBE_save != None:
                # If SBE_save is enabled, we create an HDF5 file containing the stamp of this simulated object
                # order is in self.order
                # We just save the x,y,f,w arrays as well as info about minx,maxx,miny,maxy
    
                # We trim the stamp to avoid padding area
                this_SBE_object =  this_object[self.ystart:self.yend+1,self.xstart:self.xend+1]
                
                yss,xss = np.nonzero(this_SBE_object>0)
                
                if len(xss)<1:
                    continue 

                minx = np.min(xss)
                maxx = np.max(xss)
                miny = np.min(yss)
                maxy = np.max(yss)

                this_SBE_object = this_SBE_object[miny:maxy+1,minx:maxx+1]

                if os.path.isfile(self.SBE_save):
                    mode = "a"
                else:
                    mode = "w"

                with h5py.File(self.SBE_save,mode) as fhdf5:
                    dset = fhdf5.create_dataset("%d_%s" % (self.IDs[i],self.order),data=this_SBE_object,dtype='f',compression="gzip",compression_opts=9)
                    dset.attrs[u'minx'] = minx
                    dset.attrs[u'maxx'] = maxx
                    dset.attrs[u'miny'] = miny
                    dset.attrs[u'maxy'] = maxy
                    dset.attrs[u'units'] = 'e-/s'


    def disperse_background_1D(self,background):
        """Method to create a simple disperse background, obtained by dispersing a full row or column.
        We assume no field dependence in the cross dispersion direction and create a full 2D image by tiling a single dispersed row or column"""

        # Create a fake object, line in middle of detector
        Cfg = self.Cfg
        naxis = [self.xend-self.xstart+1 ,self.yend-self.ystart+1] 
        xpos,ypos = naxis[0]//2,naxis[1]//2

        # Find out if this an x-direction or y-direction dispersion
        dydx = np.array(Cfg.DISPXY(self.order,1000,1000,1))-np.array(Cfg.DISPXY(self.order,1000,1000,0))
        if np.abs(dydx[0])>np.abs(dydx[1]):
            print("disperse_background_1D: x-direction")
            direction = "x"
            xs = np.arange(self.Cfg.XRANGE[self.order][0]+0,self.Cfg.XRANGE[self.order][1]+naxis[0])
            ys = np.zeros(np.shape(xs))+ypos
        else:
            print("disperse_background_1D: y-direction")
            direction = "y"
            ys = np.arange(self.Cfg.YRANGE[self.order][0]+0,self.Cfg.YRANGE[self.order][1]+naxis[0])
            xs = np.zeros(np.shape(ys))+xpos

        print(xpos,ypos)
        
        
        lam = background[0]
        fnu = background[1]

        fnu = fnu/4.25e10 # MJy/arcsec^2
        fnu = fnu*1e6 # Jy/arcsec^2
        fnu = fnu * (0.065**2) # Jy/pixel

        fnu = fnu*1e-23
        c = 299792458.* 1e10 # A
        wa = lam*10000
        flam = fnu/(wa**2/c)

        f = [lam,flam]
                
        pars = []        
        for i in range(len(xs)):
            ID = 1
            xs0 = [xs[i],xs[i]+1,xs[i]+1,xs[i]]
            ys0 = [ys[i],ys[i],ys[i]+1,ys[i]+1]
            pars.append([xs0,ys0,f,self.order,Cfg,ID,False,self.xstart,self.ystart])


        #from multiprocessing import Pool
        import time
        time1 = time.time()
        #mypool = Pool(self.max_cpu) # Create pool
        #all_res = mypool.imap_unordered(helper,pars) # Stuff the pool
        all_res = []
        for i in range(len(pars)):
            all_res.append(helper(pars[i]))
        #mypool.close()

        bck = np.zeros(naxis,np.float)
        for i,pp in enumerate(all_res, 1): 
            if np.shape(pp.transpose())==(1,6):
                continue
            x,y,w,f = pp[0],pp[1],pp[3],pp[4]


            vg = (x>=0) & (x<naxis[0]) & (y>=0) & (y<naxis[1]) 

            x = x[vg]
            y = y[vg]
            f = f[vg]
            w = w[vg]
            
            if len(x)<1:
                continue

            

            minx = int(min(x))
            maxx = int(max(x))
            miny = int(min(y))
            maxy = int(max(y))
            a = sparse.coo_matrix((f, (y-miny, x-minx)), shape=(maxy-miny+1,maxx-minx+1)).toarray()
            bck[miny:maxy+1,minx:maxx+1] = bck[miny:maxy+1,minx:maxx+1] + a

        if direction=="x":
            bck = np.sum(bck,axis=0)
            bck = np.tile(bck,[naxis[1],1])
        else:
            bck = np.sum(bck,axis=1)
            bck = np.tile(bck,[naxis[0],1]).transpose()

        return bck


    def disperse_chunk(self,c):
        """Method that handles the dispersion. To be called after create_pixel_list()"""
        from multiprocessing import Pool
        import time

        print(" dispersing object %d" % c)
        log.info(" dispersing object %d" % c)

        # We won't use SED file in pipeline, so this gets skipped
        if self.SED_file!=None:
            # We use an input spectrum file
            import h5py
            with h5py.File(self.SED_file,'r') as h5f:
                pars = []
                ID = int(self.seg[self.ys[c][0],self.xs[c][0]])

                tmp = h5f["%s" % (ID)][:]
                for i in range(len(self.xs[c])):
                    
                    lams = tmp[0]
                    fffs = tmp[1]*self.fs["SED"][c][i]

                    # trim input spectrum 
                    try:
                        ok = (lams>self.Cfg.WRANGE["+1"][0]) & (lams<self.Cfg.WRANGE["+1"][1])
                    except:
                        ok = (lams>self.Cfg.WRANGE["A"][0]) & (lams<self.Cfg.WRANGE["A"][1])
                    lams = lams[ok]
                    fffs = fffs[ok]

                    if self.POM_mask is not None:
                        POM_value = self.POM_mask[self.ys[c][i],self.xs[c][i]]
                    else:
                        POM_value = 1.
                    if POM_value>=10000:
                        #print("Applying additional transmission of :",self.xs[c][i],self.ys[c][i],POM_value)
                        trans = self.POM_transmission[POM_value](lams)
                        fffs = fffs * trans
                    else:
                        fffs = fffs * POM_value

                    f = [lams,fffs]
                    xs0 = [self.xs[c][i],self.xs[c][i]+1,self.xs[c][i]+1,self.xs[c][i]]
                    ys0 = [self.ys[c][i],self.ys[c][i],self.ys[c][i]+1,self.ys[c][i]+1]
                    pars.append([xs0,ys0,f,self.order,self.Cfg,ID,self.extrapolate_SED,self.xstart,self.ystart])
        else:
            # This is the branch we'll use in pipeline
            # No spectrum passed
            pars = []  # initialize params for this object
            # Loop over all pixels in list for object "c"
            for i in range(len(self.xs[c])):
                # Here "i" and "ID" are just indexes into the pixel list for the object
                # being processed, as opposed to the ID number of the object itself
                ID = i
                xs0 = [self.xs[c][i],self.xs[c][i]+1,self.xs[c][i]+1,self.xs[c][i]]
                ys0 = [self.ys[c][i],self.ys[c][i],self.ys[c][i]+1,self.ys[c][i]+1]
                lams = np.array(list(self.fs.keys()))  # these are the wavelengths previously stored in flux list
                flxs = np.array([self.fs[l][c][i] for l in self.fs.keys()])  # these are the direct image pixel fluxes
                ok = flxs!=0  # We avoid any pixel containing pure 0's
                if len(flxs[ok])==0: continue
                flxs = flxs[ok]
                lams = lams[ok]
                ok = np.argsort(lams)
                flxs = flxs[ok]  # list of good fluxes for this pixel
                lams = lams[ok]  # list of wavelengths for this pixel

                if self.POM_mask is not None:
                    POM_value = self.POM_mask[self.ys[c][i],self.xs[c][i]]
                else:
                    POM_value = 1.
                if POM_value>1:
                    #print("Applying additional transmission of :",self.xs[c][i],self.ys[c][i],POM_value)
                    trans = self.POM_transmission[POM_value](lams)
                    flxs = flxs * trans
                else:
                    flxs = flxs * POM_value

                f = [lams,flxs]
                pars.append([xs0,ys0,f,self.order,self.Cfg,ID,self.extrapolate_SED,self.xstart,self.ystart])
                # now have full pars list for all pixels for this object



        time1 = time.time()
        #mypool = Pool(self.max_cpu) # Create pool
        #all_res = mypool.imap_unordered(helper,pars) # Stuff the pool
        all_res = []
        for i in range(len(pars)):
            all_res.append(helper(pars[i]))
        #mypool.close() # No more work

        this_object = np.zeros(self.dims,np.float)

        for i,pp in enumerate(all_res, 1): 

            if np.shape(pp.transpose())==(1,6):
                continue

            x,y,w,f = pp[0],pp[1],pp[3],pp[4]

            vg = (x>=0) & (x<self.dims[1]) & (y>=0) & (y<self.dims[0]) 

            x = x[vg]
            y = y[vg]
            f = f[vg]
            w = w[vg]
            
            if len(x)<1:
                continue

            

            minx = int(min(x))
            maxx = int(max(x))
            miny = int(min(y))
            maxy = int(max(y))
            a = sparse.coo_matrix((f, (y-miny, x-minx)), shape=(maxy-miny+1,maxx-minx+1)).toarray()
            self.simulated_image[miny:maxy+1,minx:maxx+1] = self.simulated_image[miny:maxy+1,minx:maxx+1] + a
            this_object[miny:maxy+1,minx:maxx+1] = this_object[miny:maxy+1,minx:maxx+1] + a

            if self.cache:
                self.cached_object[c]['x'].append(x)
                self.cached_object[c]['y'].append(y)
                self.cached_object[c]['f'].append(f)
                self.cached_object[c]['w'].append(w)
                self.cached_object[c]['minx'].append(minx)
                self.cached_object[c]['maxx'].append(maxx)
                self.cached_object[c]['miny'].append(miny)
                self.cached_object[c]['maxy'].append(maxy)



        time2 = time.time()

        return this_object

    def disperse_all_from_cache(self,trans=None):
        if not self.cache:
            print("No cached object stored.")
            return

        self.simulated_image = np.zeros(self.dims,np.float)

        #for i in tqdm.tqdm(range(len(self.IDs)),desc="Dispersing from cache"):
        for i in range(len(self.IDs)):
            this_object = self.disperse_chunk_from_cache(i,trans=trans)


    def disperse_chunk_from_cache(self,c,trans=None):
        """Method that handles the dispersion. To be called after create_pixel_list()"""
        
        import time

        if not self.cache:
            print("No cached object stored.")
            return

        time1 = time.time()
        
        this_object = np.zeros(self.dims,np.float)

        if trans!=None:
                print("Applying a transmission function...")
        for i in range(len(self.cached_object[c]['x'])): 
            x = self.cached_object[c]['x'][i]
            y = self.cached_object[c]['y'][i]
            f = self.cached_object[c]['f'][i]*1.
            w = self.cached_object[c]['w'][i]

            if trans!=None:
                f *= trans(w)

            minx = self.cached_object[c]['minx'][i]
            maxx = self.cached_object[c]['maxx'][i]
            miny = self.cached_object[c]['miny'][i]
            maxy = self.cached_object[c]['maxy'][i]
   
            a = sparse.coo_matrix((f, (y-miny, x-minx)), shape=(maxy-miny+1,maxx-minx+1)).toarray()
            self.simulated_image[miny:maxy+1,minx:maxx+1] = self.simulated_image[miny:maxy+1,minx:maxx+1] + a
            this_object[miny:maxy+1,minx:maxx+1] = this_object[miny:maxy+1,minx:maxx+1] + a

        time2 = time.time()

        print(time2-time1,"s.")
        return this_object

    def show(self):
        import matplotlib.pyplot as plt
        plt.ion()

        xx = self.p_x - min(self.p_x)
        yy = self.p_y - min(self.p_y)
        a = sparse.coo_matrix((self.p_f, (yy, xx)), shape=(max(yy)+1, max(xx)+1)).toarray()

        im = plt.imshow(a)
        im.set_clim(0,1)

        plt.draw()
        raw_input("...")

