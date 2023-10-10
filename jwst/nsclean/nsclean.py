import logging
import numpy as np

from astropy.stats import sigma_clipped_stats

from jwst import datamodels
from stdatamodels.jwst.datamodels import dqflags
from jwst.nsclean.util import make_lowpass_filter

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


class NSClean:
    """
    NSClean is the base class for removing residual correlated
    read noise from JWST NIRSpec images.  It is intended for use
    on Level 2a pipeline products, i.e. IRS$^2$ corrected slope
    images. All processing is done in detector coordinates, with
    arrays transposed and flipped so that the IRS$^2$ "zipper"
    appears along the bottom. For most users, NSClean's `clean`
    method automatically transposes and flips the data as
    necessary.
    """
    
    # Class variables. These are the same for all instances.
    ny     = 2048     # Number of lines in detector space
    nx     = 2048     # Number of columns in detector space
    sigrej = 3.0      # Standard deviation threshold for flagging
                      #   statistical outliers.
        
    # I hope to automate this one...
    weights_kernel_sigma = 32 # Used to assign weights for MASK mode fitting

    def __init__(self, detector, mask, fc=1/(2*2048/16), kw=1/(2*2048/16)/4,
                     buffer_sigma=1.5):
        """
        JWST NIRSpec background modeling and subtraction -AKA "clean" (NSClean)

        Parameters: detector, string
                      A string selected from {'NRS1','NRS2'}
                    mask, boolean array
                      A 2048x2048 pixel boolean array. The background model
                      is fitted to pixels =True. Pixels =False are ignored.
                    fc, real number
                      Critical frequency. This is the 1/2 power point for NSClean's
                      low pass filter. The units are pixels$^{-1}$
                    kw, real number
                      "Kill width". The low pass filter's cutoff is 1/2 of a cosine.
                      The filter function goes from 1x to 0x gain over this bandwidth.
                      The units of kw are pixels$^-1$
                    buffer_sigma, real number
                      Standard deviation of the buffer's Gaussian smooth. This is an
                      optional 1-dimensional smooth in the spectral dispersion direction
                      only. If selected, buffing is done after modeling the background
                      but before subtracting it.
        
        """
        self.detector = detector
        self.M = mask
        self.fc = fc
        self.kw = kw
        self.buffer_sigma = buffer_sigma
        
        # Other handy definitions
        self.rfftfreq = np.fft.rfftfreq(self.nx) # FFT frequencies
                    
        # Transpose and flip mask to detector coordinates with the IRS$^2$
        # zipper running along the bottom as displayed in ds9.
        if self.detector=='NRS1':
            # NRS1 requires transpose only
            self.M = self.M.transpose()
        else:
            # NRS2 requires transpose and flip
            self.M = self.M.transpose()[::-1]

        # Extract other necessary information from model_parameters and build the
        # low pass filter. This alse sets the number of Fourier vectors that
        # we need to project out.
        self.apodizer = np.array(make_lowpass_filter(self.fc, self.kw, self.nx))
        self.nvec = np.sum(self.apodizer > 0) # Project out this many frequencies

        # Only MASK mode uses a weighted fit. Compute the weights here. The aim is
        # to weight by the reciprocal of the local background sample density along
        # each line. Roughly approximate the local density, P, using the reciprocal of
        # convolution by a Gaussian kernel for now. For now, hard code the kernel. We
        # will optimize this later.
        W = np.zeros((self.ny, self.nx), dtype=np.float32) # Build the kernel here
        _x = np.arange(self.nx)
        _mu = self.nx//2 + 1
        _sigma = self.weights_kernel_sigma # Hard code for now
        W[self.ny//2+1] = np.exp(-(_x - _mu)**2 / _sigma**2/2) / _sigma / np.sqrt(2*np.pi)
        FW = np.fft.rfft2(np.fft.ifftshift(W))
        with np.errstate(divide='ignore'):
            self.P = 1 / np.fft.irfft2(np.fft.rfft2(np.array(self.M,
                                       dtype=np.float32)) * FW, (self.ny,self.nx))
        self.P = np.where(self.M==True, self.P, 0.) # Illuminated areas carry no weight

        # Build a 1-dimensional Gaussian kernel for "buffing". Buffing is in the
        # dispersion direction only. In detector coordinates, this is axis zero. Even though
        # the kernel is 1-dimensional, we must still use a 2-dimensional array to 
        # represent it. I tried broadcasting a vector, but that made a kernel 2048 
        # columns wide (in detector space).
        _y = np.arange(self.ny) # Row indices
        _mu = self.ny//2 + 1 # Center of kernel
        _sigma = self.buffer_sigma # Standard deviation of kernel
        _gkern = np.exp(-((_y-_mu) / _sigma)**2 / 2) / _sigma / np.sqrt(2*np.pi) # Centered kernel as a vector
        gkern = np.zeros((self.ny, self.nx), dtype=np.float32) # 2D kernel template
        gkern[:, _mu] = _gkern # Copy in the kernel. Normalization is already correct.
        gkern = np.fft.ifftshift(gkern) # Shift for Numpy
        self.fgkern = np.array(np.fft.rfft2(gkern), dtype=np.complex64) # FFT for fast convolution
        
        
    def fit(self, D):
        """
        Fit a background model to the supplied frame of data.
        
            Parameters: D::CuArray{Float32,2}
                          A NIRSpec image. The image must be in the detector-space
                          orientation with the IRS$^2$ zipper running along the bottom as
                          displayed in SAOImage DS9.
               Returns: B::CuArray{Float32,2}
                          The fitted background model.
        Notes:
        1) Fitting is done line by line because the matrices get very big if one
           tries to project out Fourier vectors from the entire 2K x 2K image area.
        """
        model = np.zeros((self.ny, self.nx), dtype=np.float32) # Build the model here
        for y in np.arange(self.ny)[4:-4]:
            
            # Get data and weights for this line
            d = D[y][self.M[y]] # Data
            p = np.diag(self.P[y][self.M[y]]) # Weights
            
            # Fill statistical outliers with line median. We know that the rolling
            # median cleaning technique worked reasonably well, so this is fast and not crazy.
            _mu = np.median(d) # Robust estimate of mean
            _sigma = 1.4826 * np.median(np.abs(d - _mu)) # Robust estimate of standard deviation
            d = np.where(np.logical_and(_mu - self.sigrej * _sigma <= d,
                                        d <= _mu + self.sigrej * _sigma), d, _mu) # Fill outiers
            
            # Build the Fourier basis matrix for this line
            m = np.arange(self.nx)[self.M[y]].reshape((-1, 1)) # Must be a column vector to broadcast
            k = np.arange(self.nvec).reshape((1,-1)) # Must be a row vector to broadcast. We can optimize
                                                     # the code later by putting this into object instantiation
                                                     # since it is the same for every line. For now, leave it
                                                     # here since the cost is negligible and it may aid 
                                                     # comprehension.
            B = np.array(np.exp(2 * np.pi * 1J * m * k / self.nx) / m.shape[0], dtype=np.complex64) # Build the basis matrix

            # Compute the Moore-Penrose inverse of A = P*B.
            #     $A^+ = (A^H A)^{-1} A^H$
            A = np.matmul(p, B)
            AH = np.conjugate(A.transpose()) # Hermitian transpose of A
            pinv_PB = np.matmul(np.linalg.inv(np.matmul(AH, A)), AH)
            
            # Solve for the Fourier transform of this line's background samples.
            # The way that we have done it, this multiplies the input data by the 
            # number of samples used for the fit.
            rfft = np.zeros(self.nx//2 + 1, dtype=np.complex64)
            rfft[:k.shape[1]] = np.matmul(np.matmul(pinv_PB, p), d)
                        
            # Numpy requires that the forward transform multiply
            # the data by n. Correct normalization.
            rfft *= self.nx / m.shape[0]
            
            # Apodize if necessary
            if self.kw > 0:
                rfft[:self.nvec] *= self.apodizer[:self.nvec]

            # Invert the FFT to build the background model for this line
            model[y] = np.fft.irfft(rfft, self.nx)
        
        # Done!
        return(model)
    

    def clean(self, D, buff=True):
        """
        "Clean" NIRspec images by fitting and subtracting the
        instrumental background. This is intended to improve the
        residual correlated noise (vertical banding) that is
        sometimes seen.  Because the banding is not seen by the
        reference pixels, normal IRS^2 processing does not
        remove it.
        
        This is an ad-hoc correction. We model the background by
        fitting it in Fourier space. There is an option
        to "improve" the result by "buffing" in the spectral
        dispersion direction. "Buffing" is light smoothing using
        a 1-dimensional Gaussian.
        
        Parameters: D::array_like
                      The input data. This should be the normal end result of
                      Stage 1 processing.
           Returns: D::array_like
                      The data, but with less striping and the background subtracted.
                    buff, bool
                      "Buff" the fitted spectrum by applying a slight Gaussian blur
                      in the spectral dispersion direction.
        """
        
        # Transform the data to detector space with the IRS^2
        # zipper running along the bottom.
        if self.detector == 'NRS2':
            D = D.transpose()[::-1] # Transpose and flip for NRS2
        else:
            D = D.transpose()       # No flip required for NRS1
            
        # Fit, optionally buff, and subtract the background model
        B = self.fit(D) # Background model
        if buff is True:
            B = np.fft.irfft2(np.fft.rfft2(B) * self.fgkern, s=B.shape) # Buff
        D -= B # Subtract background model
        
        # Transform back to DMS space
        if self.detector=='NRS2':
            D = D[::-1].transpose() # Go back to DMS space
        else:
            D = D.transpose()       # No flip required for NRS1
            
        # Done
        return(D)


def create_mask(image, dq, n_sigma):
    """Create the pixel mask needed for setting which pixels to use
    for measuring 1/f noise
    """
    # Initialize mask with everything set to False (do not use)
    Mask = np.full(np.shape(dq), False)

    # Set NON-SCIENCE pixels to True (use)
    dark_pix = np.where(dq & dqflags.pixel['NON_SCIENCE'])
    Mask[dark_pix] = True

    # Set pixels affected by failed-open shutters to False (don't use)
    open_pix = np.where(dq & dqflags.pixel['MSA_FAILED_OPEN'])
    Mask[open_pix] = False

    # Reset NaN pixels and flag as don't use
    nan_pix = np.where(np.isnan(image))
    image[nan_pix] = 0
    Mask[nan_pix] = False

    # Don't use fixed-slit area pixels; uses hardwired indexes
    Mask[922:1116, :] = False

    # Use left/right reference pixel columns (first and last 4)
    Mask[:, 0:5] = True
    Mask[:, -5:] = True  # keep one extra column on the right (always empty)

    # Get stats to flag outliers
    _, median, sigma = sigma_clipped_stats(image, mask=~Mask, mask_value=0, sigma=5.0)
    outliers = np.where(image > (median + n_sigma * sigma))
    Mask[outliers] = False

    # Return the mask and the record of which pixels were NaN in the input;
    # it'll be needed later
    return Mask, nan_pix


def do_correction(input_model, n_sigma):

    """Execute the NSClean 1/f noise correction

    Parameters
    ----------
    input_model : data model object
        science data to be corrected

    n_sigma : float
        n-sigma rejection level for finding outliers

    Returns
    -------
    output_model : jwst.datamodel.JwstDataModel
        corrected data
    """

    detector = input_model.meta.instrument.detector.upper()
    exp_type = input_model.meta.exposure.type
    log.info(f'Input exposure type is {exp_type}, detector={detector}')
    output_model = input_model.copy()

    image = np.float32(input_model.data)
    #dq = input_model.dq

    # Hardwired override for now, until we get updated DQ values from input
    cal_model = datamodels.open("/Users/bushouse/JWST/softw/NIRSpec_NSClean/devel/jw01270003001_02101_00001_nrs2_cal.fits")
    dq = cal_model.dq
    cal_model.close

    # Create the pixel mask that'll be used to indicate which pixels
    # to include in the 1/f noise measurements
    Mask, nan_pix = create_mask(image, dq, n_sigma)

    # Instantiate an NSClean object
    cleaner = NSClean(detector, Mask)

    # Clean it
    log.info(f"cleaning image {input_model.meta.filename}")
    try:
        cleaned_image = cleaner.clean(image, buff=True)
    except np.linalg.LinAlgError:
        log.warning("Error cleaning image; step will be skipped")
        output_model.meta.cal_step.nsclean = 'SKIPPED'
        return output_model

    # Restore NaN's from original image
    cleaned_image[nan_pix] = np.nan
    output_model.data = cleaned_image

    # Set completion status
    output_model.meta.cal_step.nsclean = 'COMPLETE'

    return output_model
