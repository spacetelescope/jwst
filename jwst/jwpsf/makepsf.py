#! /usr/bin/env python

import sys
import os.path
import time
import math
import numpy as np
import numpy.fft
from astropy.io import fits as pyfits

try:
    from scipy import ndimage
except ImportError:
    try:
        import ndimage
    except:
        raise ImportError


DOUBLE_PREC = 64
SINGLE_PREC = 32

ARCSECtoDEGREES = 1. / 3600.
RADIANStoDEGREES = 180. / math.pi
MICRONStoMETERS = 1.e-6
MICRONStoNANOMETERS = 1000.

def main(args):

    # Note that main() uses only one wavelength, not a finite bandpass.

    if len(args) < 4 or len(args) > 8:
        print("syntax:  makepsf.py pupil phase output " \
              "[wavelength oversample output_size verbose]")
        print("pupil and phase are input FITS files,")
        print("while output is the output FITS file")
        print("defaults:  wavelength=1. (microns), oversample=4, " \
              "output_size=512, verbose=0")
        sys.exit()

    instrument = args[0]
    pupil_file = args[1]
    phase_file = args[2]
    output = args[3]
    # defaults
    wavelength = 1.
    oversample = 4
    output_size = 512
    verbose = False

    if len(args) > 4:
        wavelength = float(args[4])
    if len(args) > 4:
        oversample = int(args[5])
    if len(args) > 5:
        output_size = int(args[6])
    if len(args) > 6:
        verbose = bool(int(args[7]))

    # wavelength in microns, weight.
    filter = ((wavelength,), (1.,))

    psf = MakePsf(pupil_file=pupil_file, phase_file=phase_file,
              output=output, oversample=oversample, filter=filter,
              output_size=output_size, verbose=verbose)


def read_filter(filter_file):
    """Read wavelengths and weights from filter file.

    The function value is a pair of lists in the format expected by
    MakePsf as its 'filter' argument.
    """

    fd = open(filter_file, "r")
    lines = fd.readlines()
    fd.close()

    wavelengths = []
    weights = []
    for line in lines:
        line = line.strip()
        words = line.split()
        wavelengths.append(float(words[0]))
        weights.append(float(words[1]))

    return (wavelengths, weights)


class MakePsf:
    """Create a PSF

    psf = makepsf.MakePsf (instrument=None, pupil_file=None, phase_file=None,
                  output=None,
                  diameter=None, oversample=4, type=numpy.float64,
                  filter=((1.,),(1.,)),
                  output_size=512, pixel_size=0.034, verbose=False)

    arguments:
    instrument    The name of the instrument for which the PSF is constructed
    pupil_file    FITS file containing pupil map image; the values should
                    be non-negative, with zero indicating no throughput.
    phase_file    FITS file containing phase map image; the values are
                    assumed to be in microns.
    output        output FITS file.  if output=None then you should use
                    the writeto() method to save the results.
    diameter      pupil diameter in meters.  if this is specified it
                    overrides the value gotten from the phase_file header
                    (keyword CD1_1 or CDELT1 multiplied by the image width).
    oversample    oversampling factor.  default is twice Nyquist sampling.
    type          data type (float32 or float64) to use for computation;
                    the output will be float32.
    filter        a filter bandpass over which to integrate, specified
                    as a pair (tuple) of arrays, wavelengths and weights;
                    the wavelengths are in microns, and the weights are
                    relative and will be normalized so their sum is one.
    output_size   size of output image (centered on PSF).
    pixel_size    size of detector pixel in arcseconds; size in output image
                    will be 'oversample' times smaller.
                    pixel_size=None is OK if only one wavelength and weight
                    were specified in filter, and in that case the output
                    image will not be resampled.
                    The default value is the pixel size for NIRCam for
                    wavelengths less than 2.5 micron.
    verbose       print info?

    attributes:
    integrated_psf     an array of shape (output_size,output_size) and
                         type float32 containing the computed PSF

    methods:
    writeto (output)   writes the computed PSF to the specified
                         output FITS file
    info()             prints values of some of the attributes
    """

    def __init__(self, instrument=None, pupil_file=None, phase_file=None, output=None,
                  diameter=None, oversample=4, type=np.float64,
                  filter=((1.,), (1.,)),
                  output_size=512, pixel_size=0.034, verbose=False):

        self.instrument = instrument
        self.pupil_file = pupil_file
        self.phase_file = phase_file
        self.D = diameter
        self.oversample = oversample
        if type == np.float32:
            self.type = SINGLE_PREC
        else:
            self.type = DOUBLE_PREC
        self.filter = self.normalizeWeights(filter)
        self.output_size = output_size
        self.pixel_size = pixel_size
        self.verbose = verbose

        # template header, from phase_file or pupil_file
        self.header = None
        # amplitude and phase, read from pupil_file, phase_file
        self.pupil = None
        self.phase = None
        # complex array computed from self.pupil and self.phase, and
        # the same size as those arrays (generally smaller than npix)
        self.opd = None
        # the Fourier transform will be done on a complex array of
        # shape (self.npix, self.npix), larger by the factor
        # self.oversample than self.pupil or self.phase
        self.npix = None

        if pixel_size is None and len(filter[0]) > 1:
            raise ValueError("pixel_size must be specified " \
                "if filter includes more than one wavelength.")

        # Width of an output image pixel in degrees.  If pixel_size
        # was specified, this will be the resampled size; if not, this
        # will be the scale resulting from just taking the Fourier
        # transform of the OPD (expanded pupil and phase) image.  The
        # actual value will be assigned later.
        self.cdelt = None

        # Initialize the integrated PSF to an array of zeros.
        self.integrated_psf = self.initPsf()

        # Read pupil and/or phase from FITS files into memory.
        self.readAmplPhase()

        # For each wavelength, compute the PSF, resample, multiply by
        # the weight, and add it to the result (self.integrated_psf).
        self.createIntegratedPsf()

        self.output_written = False
        if output is None:
            print("Warning:  computed PSF has not been saved;")
            print("  use <psf>.writeto(output)")
        else:
            self.writeto(output)
            if self.output_written and self.verbose:
                print("Computed PSF has been written to", output)

    def info(self):
        """Print some of the attributes."""

        print("pupil file =", self.pupil_file)
        print("phase file =", self.phase_file)
        print("wavelengths and weights =")
        for i in range(len(self.filter[0])):
            print("  %10.5f  %6.4f" % (self.filter[0][i], self.filter[1][i]))
        print("pupil diameter (meters) =", self.D)
        if self.oversample == 2:
            print("oversampling factor = 2 (Nyquist sampling)")
        else:
            r = float(self.oversample) / 2.
            print("oversampling factor = %d (%g * Nyquist sampling)" % \
                        (self.oversample, r))
        if self.type == SINGLE_PREC:
            print("computations will use single precision")
        else:
            print("computations will use double precision")
        print("size of output image =", self.output_size)
        if self.cdelt is not None:
            print("output pixel size (arcsec) =", self.cdelt / ARCSECtoDEGREES)
        if self.output_written:
            print("The computed PSF has been written to the output file.")
        else:
            print("The output file has not been written yet.")

    def normalizeWeights(self, filter):
        """Scale the weights so their sum is 1."""

        (wavelengths, weights) = filter
        weights = np.array(weights, dtype=np.float64)
        sum = weights.sum()
        weights /= sum

        return (wavelengths, weights)

    def initPsf(self, type=np.float32):
        """Create an array for the integrated PSF."""

        return np.zeros((self.output_size, self.output_size), dtype=type)

    def readAmplPhase(self):
        """Read the amplitude and/or phase from FITS files"""

        if self.pupil_file is None and self.phase_file is None:
            raise ValueError("no input specified")

        if self.phase_file is None:
            self.phase = 0.
        else:
            (self.phase, self.header) = self.getData(self.phase_file)
            shape = self.phase.shape
            if len(shape) != 2 or shape[0] != shape[1]:
                raise ValueError("phase image must be 2-D and square")

        if self.pupil_file is None:
            self.pupil = 1.
        else:
            (self.pupil, header) = self.getData(self.pupil_file)
            shape = self.pupil.shape
            if len(shape) != 2 or shape[0] != shape[1]:
                raise ValueError("pupil image must be 2-D and square")
            # normalize the pupil map
            normalize = math.sqrt(self.pupil.sum())
            self.pupil /= normalize
            if self.phase_file is None:
                self.header = header
            elif self.phase.shape != self.pupil.shape:
                raise ValueError("sizes of pupil image and phase image are not the same")

        if self.header["naxis1"] != self.header["naxis2"]:
            raise ValueError("pupil and phase must be square images")

        self.npix = shape[0] * self.oversample

        # If the pupil size was not specified, get it from image header.
        if self.D is None:
            self.getDiameter()

        assert self.npix // 2 * 2 == self.npix

    def getData(self, filename):
        """Return the data and header from a FITS image.

        This converts the data to double precision if that data type
        was specified.
        """

        fd = pyfits.open(filename)
        if fd[0].data is None:
            hdu = fd[1]
        else:
            hdu = fd[0]
        header = hdu.header
        if self.type == DOUBLE_PREC:
            data = hdu.data.astype(np.float64)
        else:
            data = hdu.data.astype(np.float32)
        fd.close()

        return (data, header)

    def getDiameter(self):
        """Get size of input pixel, compute image width in meters."""

        hdr = self.header
        if "cd1_1" in hdr:
            self.D = abs(hdr["cd1_1"]) * hdr["naxis1"]
        elif "cdelt1" in hdr:
            self.D = abs(hdr["cdelt1"]) * hdr["naxis1"]
        else:
            print("Warning:  no coordinate information found in input header;")
            print("  pupil width assumed to be 6.5 meters")
            self.D = 6.5

    def createIntegratedPsf(self):
        """Compute a PSF for each wavelength, and add them up."""

        (wavelengths, weights) = self.filter
        for i in range(len(wavelengths)):

            wavelength = wavelengths[i]
            weight = weights[i]
            self.convertToOpd(wavelength)      # creates self.opd
            opd = self.embedOpd()
            zf = numpy.fft.fft2(opd)
            del opd
            # Compute the amplitude squared.
            # (psf is not really the point spread function yet)
            psf = np.conjugate(zf)
            # psf will now be the point spread function, but still complex
            np.multiply(psf, zf, psf)
            del zf
            # normalize the PSF, and convert to single precision
            psf = psf.real / psf.size
            psf = psf.astype(np.float32)

            self.center(psf)

            # This describes the image scale if no resampling is done.
            cdelt_before_resampling = (wavelength * MICRONStoMETERS) / \
                    (self.D * self.oversample) * RADIANStoDEGREES
            if self.pixel_size is None:
                # we won't resample the output image
                self.cdelt = cdelt_before_resampling
                # Extract a subset.
                if self.output_size < self.npix:
                    o_npix = self.output_size
                    n0 = (self.npix - o_npix) // 2
                    self.integrated_psf += \
                        (psf[n0:n0 + o_npix, n0:n0 + o_npix] * weight)
                else:
                    self.integrated_psf += (psf * weight)
            else:
                # we'll resample to this image scale
                self.cdelt = self.pixel_size / self.oversample * ARCSECtoDEGREES
                # These three parameters are only used by mapPsf and for
                # normalizing the weight after resampling.
                self.rescale = self.cdelt / cdelt_before_resampling
                self.input_center = (self.npix + 1) // 2
                self.output_center = (self.output_size + 1) // 2
                sub_psf = np.zeros((self.output_size, self.output_size),
                                   dtype=np.float32)
                # Do the resampling, writing the output to sub_psf.
                ndimage.geometric_transform(psf, self.mapPsf,
                        output_shape=(self.output_size, self.output_size),
                        output=sub_psf, prefilter=True)
                weight = weight * self.rescale**2
                self.integrated_psf += (sub_psf * weight)
                del sub_psf

            if self.verbose:
                print("PSF for wavelength %g has been computed" % wavelength)

    def mapPsf(self, locn):
        """Compute the point in psf corresponding to a point in sub_psf.

        This function is used by ndimage.geometric_transform.
        """

        x = (locn[1] - self.output_center) * self.rescale + self.input_center
        y = (locn[0] - self.output_center) * self.rescale + self.input_center

        return (y, x)


    def convertToOpd(self, wavelength):
        """Create optical path difference array

        The input wavelength and the phase image are in microns.
        """

        if wavelength is None:
            scale = 1.
        else:
            scale = 2. * math.pi / wavelength
        self.opd = self.pupil * np.exp(1.j * self.phase * scale)

    def embedOpd(self):
        """Embed the OPD in a larger array.

        Copy into the middle of the full array.  It's OK to not be
        centered on [0,0] because we're just going to use the amplitude
        of the Fourier transform of this array.
        """

        shape = self.opd.shape
        if self.type == DOUBLE_PREC:
            type = np.complex128
        else:
            type = np.complex64
        opd = np.zeros(shape=(self.npix, self.npix), dtype=type)

        n0 = (self.npix - shape[1]) // 2
        opd[n0:n0 + shape[0], n0:n0 + shape[1]] = self.opd

        return opd

    def writeto(self, output):
        """Write the integrated PSF to an output FITS file."""

        hdu = pyfits.PrimaryHDU(data=self.integrated_psf)
        (year, month, day, hour, minute, second, weekday, DOY, DST) = \
                    time.gmtime()
        hdu.header.update("DATE", "%4d-%02d-%02dT%02d:%02d:%02d" %
                           (year, month, day, hour, minute, second))
        hdu.header.update("FILENAME", os.path.basename(output),
                           comment="Name of this file")
        hdu.header.update("INSTRUME", self.instrument, "Instrument name")

        # Copy some specific keywords from the input header.
        ihdr = self.header
        if "BUNIT" in ihdr:
            hdu.header.update("BUNIT", ihdr.get("BUNIT"))
        if "ERR_BUDG" in ihdr:
            hdu.header.update("ERR_BUDG", ihdr.get("ERR_BUDG"),
                               comment="Optical error budget version number")
        if "SI_FP" in ihdr:
            hdu.header.update("SI_FP", ihdr.get("SI_FP"),
                               comment="Focal plane for OPD calculation")
        if "OPD_WFE" in ihdr:
            hdu.header.update("OPD_WFE", ihdr.get("OPD_WFE"),
                               comment="OPD wavefront error (nm)")
        if "W" in ihdr:
            hdu.header.update("W", ihdr.get("W"),
                               comment="Flat width of hex segment (m)")
        if "GAP" in ihdr:
            hdu.header.update("GAP", ihdr.get("GAP"),
                               comment="Gap width between hex segments (m)")
        if "EDGE" in ihdr:
            hdu.header.update("EDGE", ihdr.get("EDGE"),
                               comment="Edge roll off (m)")
        if "SW" in ihdr:
            hdu.header.update("SW", ihdr.get("SW"),
                               comment="Obscuring strut width (m)")
        if "HTS" in ihdr:
            hdu.header.update("HTS", ihdr.get("HTS"),
                               comment="Height of segment isogrid")
        if "HT2" in ihdr:
            hdu.header.update("HT2", ihdr.get("HT2"),
                               comment="Height of secondary isogrid")
        if "HT3" in ihdr:
            hdu.header.update("HT3", ihdr.get("HT3"),
                               comment="Height of tertiary isogrid")
        if "FL" in ihdr:
            hdu.header.update("FL", ihdr.get("FL"),
                               comment="Focal length (m)")

        # Add some keywords.
        if self.phase_file is not None:
            hdu.header.update("PHASE", os.path.basename(self.phase_file),
                               "Name of phase image file")
        if self.pupil_file is not None:
            hdu.header.update("PUPIL", os.path.basename(self.pupil_file),
                               "Name of pupil image file")
        hdu.header.update("OVERSAMP", self.oversample, "Oversampling factor")
        hdu.header.update("CALCTYPE", self.type,
                           "32 = single precision, 64 = double precision")
        hdu.header.update("DIAMETER", self.D, "pupil diameter (meters)")
        hdu.header.update("ORIG_NX", self.header["naxis1"],
                           "NAXIS1 in input image")
        hdu.header.update("ORIG_NY", self.header["naxis2"],
                           "NAXIS2 in input image")

        self.putCoordInfo(hdu)

        (wavelengths, weights) = self.filter
        if len(wavelengths) >= 99:
            root_wln = "WAV"
            root_wgt = "WGT"
        else:
            root_wln = "WAVELN"
            root_wgt = "WEIGHT"
        for i in range(len(wavelengths)):
            keyword = "%s%d" % (root_wln, i + 1)
            hdu.header.update(keyword, wavelengths[i],
                    "wavelength in microns")
            keyword = "%s%d" % (root_wgt, i + 1)
            hdu.header.update(keyword, weights[i], "weight")

        ofd = pyfits.HDUList(hdu)
        try:
            ofd.writeto(output)
        except IOError as message:
            print("ERROR:  Output file has NOT been written; " \
                  "use <psf>.writeto(output)")
            print(message)
            return
        self.output_written = True

    def putCoordInfo(self, hdu):
        """Assign coordinate parameters in output header."""

        hdu.header.update("CDELT1", self.cdelt, "pixel size, degrees")
        hdu.header.update("CDELT2", self.cdelt, "pixel size, degrees")
        shape = hdu.data.shape
        crpix1 = (shape[0] + 1) // 2 + 1
        crpix2 = (shape[1] + 1) // 2 + 1
        hdu.header.update("CRPIX1", float(crpix1), "reference pixel")
        hdu.header.update("CRPIX2", float(crpix2), "reference pixel")

    def center(self, x):
        """Shift the contents of x by half its width.

        x is assumed to be rank 2, and the length of each axis is assumed
        to be even.  x will be modified in-place.
        """

        shape = x.shape
        nx = shape[1]
        ny = shape[0]
        hnx = nx // 2
        hny = ny // 2

        temp = x[0:hny, 0:hnx].copy()
        x[0:hny, 0:hnx] = x[hny:ny, hnx:nx].copy()
        x[hny:ny, hnx:nx] = temp

        temp = x[0:hny, hnx:nx].copy()
        x[0:hny, hnx:nx] = x[hny:ny, 0:hnx].copy()
        x[hny:ny, 0:hnx] = temp

if __name__ == "__main__":

    main(sys.argv[1:])
