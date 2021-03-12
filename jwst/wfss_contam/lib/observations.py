import time
import numpy as np
from multiprocessing import Pool
from astropy.io import fits
from scipy import sparse
from scipy.interpolate import interp1d
from .disperse import dispersed_pixel

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
    x0s, y0s, flux, order, wmin, wmax, seg_wcs, grism_wcs, ID, extrapolate_SED, xoffset, yoffset = vars

    # use the list of vars to compute dispersed attributes for the pixel
    p = dispersed_pixel(x0s, y0s, flux, order, wmin, wmax, seg_wcs, grism_wcs, ID,
                        extrapolate_SED=extrapolate_SED, xoffset=xoffset, yoffset=yoffset)

    # unpack the results
    xs, ys, areas, lams, counts, ID = p
    IDs = [ID] * len(xs)

    # return results as numpy array
    pp = np.array([xs, ys, areas, lams, counts, IDs])
    return pp


class observation():
    # This class defines an actual observation. It is tied to a single grism image.

    def __init__(self, direct_images, segmap_model, grism_wcs, wr_ref, filter, config,
                 mod="A", order=1, max_split=100, SED_file=None, extrapolate_SED=False, max_cpu=8,
                 ID=0, SBE_save=None, boundaries=[], renormalize=True):

        """
        direct_images: List of file name(s) containing direct imaging data
        segmentation_data: a data array of the same size as the direct image(s), containing 0's in
            pixels with no source signal and an integer source number for pixels with signal
        config: Full path name of a NIRCam GRISMCONF configuration file
        mod: NIRCam Module: A or B
        order: The spectral order to simulate, +1 or +2 for NIRCAM
        max_split: Number of chunks to compute instead of trying everything at once.
        SED_file: Name of HDF5 file containing datasets matching the ID in the segmentation file
                  and each consisting of a [[lambda],[flux]] array.
        extrapolate_SED:
        max_cpu:
        ID:
        SBE_save: If set to a path, HDF5 containing simulated stamps for all obsjects will be saved.
        boundaries: a tuple containing the coordinates of the FOV within the larger seed image.
        renormalize:
        """

        # This loads all the info from the configuration file for this grism mode;
        # things like the name of the POM mask file, the dispersion coeffs, etc.
        self.seg_wcs = segmap_model.meta.wcs
        self.grism_wcs = grism_wcs
        self.ID = ID
        self.IDs = []
        self.dir_image_names = direct_images
        self.seg = segmap_model.data
        #self.dims = np.shape(self.seg)
        self.dims = [2048, 2048]
        self.filter = filter
        self.order = order
        self.SED_file = SED_file   # should always be NONE for baseline pipeline (use flat SED)
        self.SBE_save = SBE_save   # may want to have a new product in which simulated spectra are saved (NOT HDF5 format!)
        self.max_cpu = max_cpu
        self.cache = False
        self.POM_mask = None
        self.POM_mask01 = None
        self.renormalize = renormalize

        # Get the wavelength range for the grism image
        wavelength_range = wr_ref.get_wfss_wavelength_range(self.filter, [self.order])
        self.wmin = wavelength_range[self.order][0]
        self.wmax = wavelength_range[self.order][1]
        print("wmin, wmax:", self.wmin, self.wmax)

        if len(boundaries) != 4:
            #self.NAXIS = [2048, 2048]
            #xpad = (np.shape(self.seg)[1] - self.NAXIS[0])//2
            #ypad = (np.shape(self.seg)[0] - self.NAXIS[1])//2
            #self.xstart = 0 + xpad
            #self.xend = xpad + self.NAXIS[0] - 1
            #self.ystart = 0 + ypad
            #self.yend = ypad + self.NAXIS[1] - 1
            #print(f"No boundaries passed. Assuming symmetrical padding of {xpad} {ypad} pixels")
            #print(f"and a final size of {self.xend+1-self.xstart} {self.yend+1-self.ystart}.")
            self.xstart = 0
            self.xend = self.xstart + self.dims[0] - 1
            self.ystart = 0
            self.yend = self.ystart + self.dims[1] - 1
            print("No boundaries passed.")
            print(f"Using final size of {self.xend+1-self.xstart} {self.yend+1-self.ystart}.")
        else:
            self.xstart, self.xend, self.ystart, self.yend = boundaries

        # Allow for SED extrapolation
        self.extrapolate_SED = extrapolate_SED
        if self.extrapolate_SED:
            print("Warning: SED Extrapolation turned on.")

        # Apply the POM mask
        #self.apply_POM()

        # Create pixel lists for ALL sources labeled in segmentation map
        self.create_pixel_list()

    def create_pixel_list(self):
        # Create a list of pixels to dispersed, grouped per object ID

        # In the examples this is called with ID=0, so processing is
        # done for all sources in the segmentation map ("seg" attribute)

        # Process all sources defined in seg map
        if self.ID == 0:
            # This creates a huge list of all x,y pixel indices that
            # have a non-zero value in the seg map, sorted by those
            # indices belonging to a particular source ID
            self.xs = []
            self.ys = []
            all_IDs = np.array(list(set(np.ravel(self.seg))))
            all_IDs = all_IDs[all_IDs > 0]
            print(f"We have {len(all_IDs)} objects")
            for ID in all_IDs:
                ys, xs = np.nonzero(self.seg == ID)
                if (len(xs) > 0) & (len(ys) > 0):
                    self.xs.append(xs)
                    self.ys.append(ys)
                    self.IDs = all_IDs

        # Process only the given source ID
        else:
            ys, xs = np.nonzero(self.seg == self.ID)
            if (len(xs) > 0) & (len(ys) > 0):
                self.xs.append(xs)
                self.ys.append(ys)
                self.IDs = [self.ID]

        self.fluxes = {}
        for dir_image_name in self.dir_image_names:
            print(f"dir image: {dir_image_name}")
            if self.SED_file is None:

                # Pipeline will use SED_file=None, so this code will try to read
                # photometry keyword values from the direct images.
                # Unfortunately it's trying to read the old HST-style keywords
                # PHOTPLAM (pivot wavelength) and PHOTFLAM, which DO NOT EXIST
                # for JWST products. We have PHOTMJSR and PHOTUJA2, so this will
                # need to be modified in some way to comply.
                #try:
                #    # convert pivlam from Angstroms to microns
                #    pivlam = fits.getval(dir_image_name, 'PHOTPLAM') / 10000.
                #except KeyError:
                #    print("ERROR: unable to find PHOTPLAM keyword in {}".format(dir_image_name))
                #    return
                # temporary hack to compute an approximate pivlam from filter name
                pivlam = float(self.filter[1:-1]) / 100

                #try:
                #    photflam = fits.getval(dir_image_name, 'photflam')
                #except KeyError:
                #    print("ERROR: unable to find PHOTFLAM keyword in {}".format(dir_image_name))
                #    return
                # JWST data is already flux calibrated, so just set photflam=1
                photflam = 1.0

                print(f"Loaded {dir_image_name} wavelength: {pivlam} micron")
            try:
                dimage = fits.open(dir_image_name)[1].data
            except IndexError:
                dimage = fits.open(dir_image_name)[0].data

            # If we do not use an SED file then we use photometry to get fluxes
            # Otherwise, we assume that objects are normalized to 1.
            if self.SED_file is None:

                # This next line uses the PHOTPLAM wavelength value, but we don't have that
                # for JWST images. Will need to figure out wavelength some other way, e.g.
                # by simply parsing from the filter name.
                self.fluxes[pivlam] = []
                dnew = dimage
                if self.POM_mask01 is not None:
                    dnew = dimage * self.POM_mask01  # Apply POM transmission mask to the data pixels
                for i in range(len(self.IDs)):
                    # This multiplies direct image pixels for each source, which apparently are
                    # still in units of countrate, by PHOTFLAM, to convert them to fluxes.
                    # So "fs" is the list of fluxes for all source pixels, at a given wavelength "pivlam"
                    # Note that JWST data are already flux calibrated, so we've set photflam=1, and
                    # hence the multiplication does nothing.
                    self.fluxes[pivlam].append(dnew[self.ys[i], self.xs[i]] * photflam)

            # Use an SED file
            else:
                # Need to normalize the object stamps
                for ID in self.IDs:
                    vg = self.seg == ID
                    dnew = dimage
                    if self.POM_mask01 is not None:
                        print("Applying POM transmission to data")

                        # Apply POM transmission mask to the data pixels.
                        # This is a single grey correction for the whole object.
                        dnew = dimage * self.POM_mask01

                    if self.renormalize is True:
                        sum_seg = np.sum(dimage[vg])  # But normalize by the whole flux
                        if sum_seg != 0:
                            dimage[vg] /= sum_seg
                    else:
                        print("not renormlazing sources to unity")

                self.fluxes["SED"] = []
                for i in range(len(self.IDs)):
                    self.fluxes["SED"].append(dnew[self.ys[i], self.xs[i]])

    def disperse_all(self, cache=False):

        if cache:
            print("Object caching ON")
            self.cache = True
            self.cached_object = {}

        self.simulated_image = np.zeros(self.dims, np.float)

        for i in range(len(self.IDs)):
            #if i != 70:
            #    continue

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

            if self.SBE_save is not None:
                # If SBE_save is enabled, we create an HDF5 file containing the stamp of this simulated object
                # order is in self.order
                # We just save the x,y,f,w arrays as well as info about minx,maxx,miny,maxy

                # We trim the stamp to avoid padding area
                this_SBE_object = this_object[self.ystart:self.yend + 1, self.xstart:self.xend + 1]

                yss, xss = np.nonzero(this_SBE_object > 0)

                if len(xss) < 1:
                    continue

                minx = np.min(xss)
                maxx = np.max(xss)
                miny = np.min(yss)
                maxy = np.max(yss)

                this_SBE_object = this_SBE_object[miny:maxy + 1, minx:maxx + 1]

                #if os.path.isfile(self.SBE_save):
                #    mode = "a"
                #else:
                #    mode = "w"

                #with h5py.File(self.SBE_save,mode) as fhdf5:
                #    dset = fhdf5.create_dataset("%d_%s" % (self.IDs[i],self.order),
                #                                data=this_SBE_object, dtype='f',
                #                                compression="gzip", compression_opts=9)
                #    dset.attrs[u'minx'] = minx
                #    dset.attrs[u'maxx'] = maxx
                #    dset.attrs[u'miny'] = miny
                #    dset.attrs[u'maxy'] = maxy
                #    dset.attrs[u'units'] = 'e-/s'

    def disperse_background_1D(self, background):
        """
        Method to create a simple disperse background, obtained by dispersing a full row or column.
        We assume no field dependence in the cross-dispersion direction and create a full 2D image
        by tiling a single dispersed row or column
        """

        # Create a fake object, line in middle of detector
        Cfg = self.Cfg
        naxis = [self.xend - self.xstart + 1, self.yend - self.ystart + 1]
        xpos, ypos = naxis[0] // 2, naxis[1] // 2

        # Find out if this an x-direction or y-direction dispersion
        dydx = np.array(Cfg.DISPXY(self.order, 1000, 1000, 1)) - np.array(Cfg.DISPXY(self.order, 1000, 1000, 0))
        if np.abs(dydx[0]) > np.abs(dydx[1]):
            print("disperse_background_1D: x-direction")
            direction = "x"
            xs = np.arange(self.Cfg.XRANGE[self.order][0] + 0, self.Cfg.XRANGE[self.order][1] + naxis[0])
            ys = np.zeros(np.shape(xs)) + ypos
        else:
            print("disperse_background_1D: y-direction")
            direction = "y"
            ys = np.arange(self.Cfg.YRANGE[self.order][0] + 0, self.Cfg.YRANGE[self.order][1] + naxis[0])
            xs = np.zeros(np.shape(ys)) + xpos

        print(xpos, ypos)

        lam = background[0]
        fnu = background[1]

        fnu = fnu / 4.25e10     # MJy/arcsec^2
        fnu = fnu * 1e6         # Jy/arcsec^2
        fnu = fnu * (0.065**2)  # Jy/pixel

        fnu = fnu * 1e-23
        c = 299792458. * 1e10  # A
        wa = lam * 10000
        flam = fnu / (wa**2 / c)

        f = [lam, flam]

        pars = []
        for i in range(len(xs)):
            ID = 1
            xs0 = [xs[i], xs[i]+1, xs[i]+1, xs[i]]
            ys0 = [ys[i], ys[i], ys[i]+1, ys[i]+1]
            pars.append([xs0, ys0, f, self.order, Cfg, ID, False, self.xstart, self.ystart])

        # Compute the dispersed results
        if self.max_cpu > 1:
            mypool = Pool(self.max_cpu)  # Create pool
            all_res = mypool.imap_unordered(helper, pars)  # Stuff the pool
            mypool.close()
        else:
            all_res = []
            for i in range(len(pars)):
                all_res.append(helper(pars[i]))

        bck = np.zeros(naxis, np.float)
        for i, pp in enumerate(all_res, 1):
            if np.shape(pp.transpose()) == (1, 6):
                continue
            x, y, w, f = pp[0], pp[1], pp[3], pp[4]

            vg = (x >= 0) & (x < naxis[0]) & (y >= 0) & (y < naxis[1])

            x = x[vg]
            y = y[vg]
            f = f[vg]
            w = w[vg]

            if len(x) < 1:
                continue

            minx = int(min(x))
            maxx = int(max(x))
            miny = int(min(y))
            maxy = int(max(y))
            a = sparse.coo_matrix((f, (y-miny, x-minx)), shape=(maxy-miny+1, maxx-minx+1)).toarray()
            bck[miny:maxy+1, minx:maxx+1] = bck[miny:maxy+1, minx:maxx+1] + a

        if direction == "x":
            bck = np.sum(bck, axis=0)
            bck = np.tile(bck, [naxis[1], 1])
        else:
            bck = np.sum(bck, axis=1)
            bck = np.tile(bck, [naxis[0], 1]).transpose()

        return bck

    def disperse_chunk(self, c):
        """Method that handles the dispersion. To be called after create_pixel_list()"""

        print(f" dispersing object {c+1}")

        # No spectrum passed
        pars = []  # initialize params for this object

        # Loop over all pixels in list for object "c"
        print("len(xs):", len(self.xs[c]))
        for i in range(len(self.xs[c])):

            # Here "i" and "ID" are just indexes into the pixel list for the object
            # being processed, as opposed to the ID number of the object itself
            ID = i

            # xs0, ys0 form a 2x2 grid of pixel indexes surrounding the direct image
            # pixel index
            xs0 = [self.xs[c][i], self.xs[c][i]+1, self.xs[c][i]+1, self.xs[c][i]]
            ys0 = [self.ys[c][i], self.ys[c][i], self.ys[c][i]+1, self.ys[c][i]+1]

            # "lams" is the list of wavelengths previously stored in flux list
            # and correspond to the central wavelengths of the filters used in
            # the input direct image(s). For the simple case of 1 combined direct image,
            # this contains a single value (e.g. 4.44 for F444W).
            lams = np.array(list(self.fluxes.keys()))
            #print("lams:", lams)

            # "flxs" is the list of pixel values ("fluxes") from the direct image(s).
            # For the simple case of 1 combined direct image, this contains a
            # a single value (just like "lams").
            flxs = np.array([self.fluxes[lam][c][i] for lam in self.fluxes.keys()])
            #print("flxs:", flxs)

            # Only include data for pixels with non-zero fluxes
            ok = flxs != 0
            if len(flxs[ok]) == 0:
                continue
            flxs = flxs[ok]
            lams = lams[ok]
            ok = np.argsort(lams)
            flxs = flxs[ok]  # list of good fluxes for this pixel
            lams = lams[ok]  # list of wavelengths for this pixel

            # Apply POM mask correction to the fluxes
            if self.POM_mask is not None:
                POM_value = self.POM_mask[self.ys[c][i], self.xs[c][i]]
            else:
                POM_value = 1.
            #print("POM_value:", POM_value)
            if POM_value > 1:
                print("Applying additional transmission of:", self.xs[c][i], self.ys[c][i], POM_value)
                trans = self.POM_transmission[POM_value](lams)
                flxs = flxs * trans
            else:
                flxs = flxs * POM_value

            f = [lams, flxs]
            pars.append([xs0, ys0, f, self.order, self.wmin, self.wmax, self.seg_wcs, self.grism_wcs,
                         ID, self.extrapolate_SED, self.xstart, self.ystart])
            # now have full pars list for all pixels for this object

        time1 = time.time()
        if self.max_cpu > 1:
            mypool = Pool(self.max_cpu)  # Create pool
            all_res = mypool.imap_unordered(helper, pars)  # Stuff the pool
            mypool.close()  # No more work
        else:
            all_res = []
            for i in range(len(pars)):
                all_res.append(helper(pars[i]))

        this_object = np.zeros(self.dims, np.float)

        for i, pp in enumerate(all_res, 1):
            if np.shape(pp.transpose()) == (1, 6):
                continue

            x, y, w, f = pp[0], pp[1], pp[3], pp[4]

            vg = (x >= 0) & (x < self.dims[1]) & (y >= 0) & (y < self.dims[0])

            x = x[vg]
            y = y[vg]
            f = f[vg]
            w = w[vg]

            if len(x) < 1:
                continue

            minx = int(min(x))
            maxx = int(max(x))
            miny = int(min(y))
            maxy = int(max(y))
            a = sparse.coo_matrix((f, (y-miny, x-minx)), shape=(maxy-miny+1, maxx-minx+1)).toarray()
            self.simulated_image[miny:maxy+1, minx:maxx+1] = self.simulated_image[miny:maxy+1, minx:maxx+1] + a
            this_object[miny:maxy+1, minx:maxx+1] = this_object[miny:maxy+1, minx:maxx+1] + a

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
        print(time2-time1, "sec")

        return this_object

    def disperse_all_from_cache(self, trans=None):
        if not self.cache:
            print("No cached object stored.")
            return

        self.simulated_image = np.zeros(self.dims, np.float)

        for i in range(len(self.IDs)):
            this_object = self.disperse_chunk_from_cache(i, trans=trans)

        return this_object

    def disperse_chunk_from_cache(self, c, trans=None):
        """Method that handles the dispersion. To be called after create_pixel_list()"""

        if not self.cache:
            print("No cached object stored.")
            return

        time1 = time.time()

        this_object = np.zeros(self.dims, np.float)

        if trans is not None:
            print("Applying a transmission function...")

        for i in range(len(self.cached_object[c]['x'])):
            x = self.cached_object[c]['x'][i]
            y = self.cached_object[c]['y'][i]
            f = self.cached_object[c]['f'][i]*1.
            w = self.cached_object[c]['w'][i]

            if trans is not None:
                f *= trans(w)

            minx = self.cached_object[c]['minx'][i]
            maxx = self.cached_object[c]['maxx'][i]
            miny = self.cached_object[c]['miny'][i]
            maxy = self.cached_object[c]['maxy'][i]

            a = sparse.coo_matrix((f, (y-miny, x-minx)), shape=(maxy-miny+1, maxx-minx+1)).toarray()
            self.simulated_image[miny:maxy+1, minx:maxx+1] = self.simulated_image[miny:maxy+1, minx:maxx+1] + a
            this_object[miny:maxy+1, minx:maxx+1] = this_object[miny:maxy+1, minx:maxx+1] + a

        time2 = time.time()

        print(time2-time1, "sec")
        return this_object
