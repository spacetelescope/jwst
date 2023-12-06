import numpy as np


class NSClean:
    """
    NSClean is the base class for removing residual correlated
    read noise from JWST NIRSpec images.  It is intended for use
    on Level 2a pipeline products, i.e. IRS2 corrected slope
    images. All processing is done in detector coordinates, with
    arrays transposed and flipped so that the IRS2 "zipper"
    appears along the bottom. For most users, NSClean's `clean`
    method automatically transposes and flips the data as
    necessary.
    """
    
    def __init__(self, detector, mask, fc=1/(2*2048/16), kill_width=1/(2*2048/16)/4,
                 buffer_sigma=1.5, sigrej=3.0, weights_kernel_sigma=32):
        """
        JWST NIRSpec background modeling and subtraction -AKA "clean" (NSClean)

        Parameters
        ----------
        detector : str
            A string selected from {'NRS1','NRS2'}

        mask : bool array
            The background model is fitted to pixels set to True.
            Pixels set to False are ignored.

        fc : float
            Critical frequency. This is the 1/2 power point for NSClean's
            low pass filter. The units are 1/pixels.

        kill_width : float
            "Kill width". The low pass filter's cutoff is 1/2 of a cosine.
            The filter function goes from 1x to 0x gain over this bandwidth.
            The units are 1/pixels.

        buffer_sigma : float
            Standard deviation of the buffer's Gaussian smooth. This is an
            optional 1-dimensional smooth in the spectral dispersion direction
            only. If selected, buffing is done after modeling the background
            but before subtracting it.

        sigrej : float
            Standard deviation threshold for flagging statistical outliers,
            to exclude them from use in fitting the background model.

        weights_kernel_sigma : int
            Used to assign weights for MASK mode fitting.
        """
        self.detector = detector
        self.mask = mask
        self.fc = fc
        self.kill_width = kill_width
        self.buffer_sigma = buffer_sigma
        self.sigrej = sigrej
        self.weights_kernel_sigma = weights_kernel_sigma
        
        # Transpose and flip mask to detector coordinates with the IRS2
        # zipper running along the bottom as displayed in ds9.
        if self.detector == 'NRS1':
            # NRS1 requires transpose only
            self.mask = self.mask.transpose()
        else:
            # NRS2 requires transpose and flip
            self.mask = self.mask.transpose()[::-1]

        self.ny = self.mask.shape[0]
        self.nx = self.mask.shape[1]

        # FFT frequencies
        self.rfftfreq = np.fft.rfftfreq(self.ny)
                    
        # Extract other necessary information from model_parameters and build the
        # low pass filter. This alse sets the number of Fourier vectors that
        # we need to project out.
        self.apodizer = np.array(make_lowpass_filter(self.fc, self.kill_width, self.nx))
        self.nvec = np.sum(self.apodizer > 0)  # Project out this many frequencies

        # Only MASK mode uses a weighted fit. Compute the weights here. The aim is
        # to weight by the reciprocal of the local background sample density along
        # each line. Roughly approximate the local density, P, using the reciprocal of
        # convolution by a Gaussian kernel for now. For now, hard code the kernel. We
        # will optimize this later.
        W = np.zeros((self.ny, self.nx), dtype=np.float32)  # Build the kernel here
        _x = np.arange(self.nx)
        _mu = self.nx//2 + 1
        _sigma = self.weights_kernel_sigma
        W[self.ny//2+1] = np.exp(-(_x - _mu)**2 / _sigma**2/2) / _sigma / np.sqrt(2*np.pi)
        FW = np.fft.rfft2(np.fft.ifftshift(W))
        with np.errstate(divide='ignore'):
            self.P = 1 / np.fft.irfft2(np.fft.rfft2(np.array(self.mask, dtype=np.float32)) * FW, (self.ny,self.nx))
        self.P = np.where(self.mask, self.P, 0.) # Illuminated areas carry no weight

        # Build a 1-dimensional Gaussian kernel for "buffing". Buffing is in the
        # dispersion direction only. In detector coordinates, this is axis zero. Even though
        # the kernel is 1-dimensional, we must still use a 2-dimensional array to 
        # represent it. I tried broadcasting a vector, but that made a kernel 2048 
        # columns wide (in detector space).
        _y = np.arange(self.ny)
        _mu = self.nx//2 + 1
        _sigma = self.buffer_sigma  # Standard deviation of kernel
        _gkern = np.exp(-((_y-_mu) / _sigma)**2 / 2) / _sigma / np.sqrt(2*np.pi)  # Centered kernel as a vector
        gkern = np.zeros((self.ny, self.nx), dtype=np.float32)  # 2D kernel template
        gkern[:, _mu] = _gkern  # Copy in the kernel. Normalization is already correct.
        gkern = np.fft.ifftshift(gkern)  # Shift for Numpy
        self.fgkern = np.array(np.fft.rfft2(gkern), dtype=np.complex64)  # FFT for fast convolution
        
        
    def fit(self, data):
        """
        Fit a background model to the supplied frame of data.
        
        Parameters
        ----------
        data : float array
            A NIRSpec image. The image must be in the detector-space
            orientation with the IRS2 zipper running along the bottom as
            displayed in SAOImage DS9.

        Returns
        -------
        Bkg : float array
            The fitted background model.

        Notes
        -----
        Fitting is done line by line because the matrices get very big if one
        tries to project out Fourier vectors from the entire 2K x 2K image area.
        """
        model = np.zeros((self.ny, self.nx), dtype=np.float32)  # Build the model here
        for y in np.arange(self.ny)[4:-4]:
            
            # Get data and weights for this line
            d = data[y][self.mask[y]]  # unmasked (useable) data
            p = np.diag(self.P[y][self.mask[y]])  # Weights

            # If none of the pixels in this line is useable (all masked out),
            # skip and move on to the next line.
            if len(d) == 0:
                continue
            
            # Fill statistical outliers with line median. We know that the rolling
            # median cleaning technique worked reasonably well, so this is a fast
            # justifiable approximation.
            _mu = np.median(d)  # Robust estimate of mean
            _sigma = 1.4826 * np.median(np.abs(d - _mu))  # Robust estimate of standard deviation

            # Fill outliers
            d = np.where(np.logical_and(_mu - self.sigrej * _sigma <= d,
                                        d <= _mu + self.sigrej * _sigma), d, _mu)
            
            # Build the Fourier basis matrix for this line
            m = np.arange(self.nx)[self.mask[y]].reshape((-1, 1))  # Must be a column vector to broadcast
            k = np.arange(self.nvec).reshape((1,-1))  # Must be a row vector to broadcast. We can optimize
                                                      # the code later by putting this into object instantiation
                                                      # since it is the same for every line. For now, leave it
                                                      # here since the cost is negligible and it may aid 
                                                      # comprehension.

            # Build the basis matrix
            B = np.array(np.exp(2 * np.pi * 1J * m * k / self.nx) / m.shape[0], dtype=np.complex64)

            # Compute the Moore-Penrose inverse of A = P*B.
            #     $A^+ = (A^H A)^{-1} A^H$
            A = np.matmul(p, B)
            AH = np.conjugate(A.transpose())  # Hermitian transpose of A
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
            if self.kill_width > 0:
                rfft[:self.nvec] *= self.apodizer[:self.nvec]

            # Invert the FFT to build the background model for this line
            model[y] = np.fft.irfft(rfft, self.nx)
        
        # Done!
        return model
    

    def clean(self, data, buff=True):
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
        
        Parameters
        ----------
        data : array_like
            The input data. This should be the normal end result of Stage 1 processing.

        buff : bool
            "Buff" the fitted spectrum by applying a slight Gaussian blur
            in the spectral dispersion direction.

        Returns
        -------
        data : array_like
            The data, but with less striping and the background subtracted.
        """
        
        # Transform the data to detector space with the IRS2 zipper running along the bottom.
        if self.detector == 'NRS2':
            # Transpose and flip for NRS2
            data = data.transpose()[::-1]
        else:
            # Transpose (no flip) for NRS1
            data = data.transpose()
            
        # Fit the background model
        Bkg = self.fit(data)  # Background model

        # Buff, if requested
        if buff:
            Bkg = np.fft.irfft2(np.fft.rfft2(Bkg) * self.fgkern, s=Bkg.shape)

        # Subtract the background model from the data
        data -= Bkg
        
        # Transform back to DMS space
        if self.detector=='NRS2':
            data = data[::-1].transpose()
        else:
            data = data.transpose()
            
        # Done
        return data 


def make_lowpass_filter(f_half_power, w_cutoff, n, d=1.0):
    """
    Make a lowpass Fourier filter
    
    Parameters
    ----------
    f_half_power : float
        Half power frequency

    w_cutoff : float
        Width of cosine cutoff. The response transitions
        from 1x to 0x over this range of frequencies

    n : int
        Number of samples in the timeseries

    d : float
       Sample spacing (inverse of the sampling rate). Defaults to 1.

    Returns
    -------
    filt : array
       Filter array
    """
    
    # Make frequencies vector
    freq = np.fft.rfftfreq(n, d=d)
    
    # Build a cosine wave that is approriately shifted
    cos = (1 + np.cos(np.pi * (freq - f_half_power) / w_cutoff + np.pi / 2)) / 2
    
    # Construct low-pass filter with cosine rolloff
    filt = np.where(freq <= f_half_power - w_cutoff / 2, 1, cos)
    filt = np.where(freq <= f_half_power + w_cutoff / 2, filt, 0)
    
    # Done
    return filt


def med_abs_deviation(d, median=True):
    """
    Median absolute deviation
        
    Computes the median and the median absolute deviation (MAD). For normally
    distributed data, multiply the MAD by 1.4826 to approximate standard deviation. 
    If median=True (the default), this returns a tuple containing (median, MAD).
    Otherwise, only the MAD is returned.
    
    Parameters
    ----------
    d : ndarray
        The input data

    median : bool
        Return both the median and MAD

    Returns
    -------
    m : float
        median

    mad : float
        median absolute deviation
    """
    d = d[np.isfinite(d)]  # Exclude NaNs
    m = np.median(d)
    mad = np.median(np.abs(d - m))
    if median is True:
        return(m, mad)
    else:
        return(mad)


class NSCleanSubarray:
    """
    NSCleanSubarray is the base class for removing residual correlated
    read noise from generic JWST near-IR Subarray images.  It is
    intended for use on Level 2a pipeline products, i.e. slope images.
    """
    
    # Class variables. These are the same for all instances.
    nloh = np.int32(12)       # New line overhead in pixels
    tpix = np.float32(10.e-6) # Pixel dwell time in seconds
    sigrej = np.float32(4.0)  # Standard deviation threshold for flagging
                              #   statistical outliers.
        
    def __init__(self, data, mask, fc=(1061, 1211, 49943, 49957),
                 exclude_outliers=True, weights_kernel_sigma=None):
        """
        Background modeling and subtraction for generic JWST near-IR subarrays.

        Parameters
        ----------
        data : float array
            The 2D input image data array to be operated on.

        mask : bool array
            The background model is fitted to pixels set to True.
            Pixels set to False are ignored.

        fc : tuple
            Apodizing filter definition. These parameters are tunable. They
            happen to work well for NIRSpec BOTS exposures.
              1) Unity gain for f < fc[0]
              2) Cosine roll-off from fc[0] to fc[1]
              3) Zero gain from fc[1] to fc[2]
              4) Cosine roll-on from fc[2] to fc[3]

        exclude_outliers : bool
            Exclude statistical outliers and their nearest neighbors
            from the background pixels mask.

        weights_kernel_sigma : float
            Standard deviation of the 1-dimensional Gaussian kernel
            that is used to approximate background sample density. This
            is ad-hoc. See the NSClean journal article for more information. The
            default for subarrays results in nearly equal weighting of all background
            samples.

        Notes:
        1) NSCleanSubarray works in detector coordinates. Both the data and mask
           need to be transposed and flipped so that slow-scan runs from bottom
           to top as displayed in SAOImage DS9. The fast scan direction is 
           required to run from left to right.
        """
        # Definitions
        self.data = np.array(data, dtype=np.float32)
        self.mask = np.array(mask, dtype=np.bool_)
        self.ny = np.int32(data.shape[0])  # Number of pixels in slow scan direction
        self.nx = np.int32(data.shape[1])  # Number of pixels in fast scan direction
        self.fc = np.array(fc, dtype=np.float32)
        self.n = np.int32(self.ny * (self.nx + self.nloh))  # Number of ticks in clocking pattern
        self.rfftfreq = np.array(np.fft.rfftfreq(self.n, self.tpix),
                                 dtype=np.float32)  # Fourier frequencies in clocking pattern
        
        # We will weight by the inverse of the local sample density in time. We compute the local
        # sample density by convolution using a Gaussian. Define the standard deviation of the
        # Gaussian here. This is ad-hoc.
        if weights_kernel_sigma is None:
            self.weights_kernel_sigma = 1 / ((fc[0]+fc[1]) / 2) / self.tpix / 2 / 4
        else:
            self.weights_kernel_sigma = weights_kernel_sigma
        
        # The mask potentially contains NaNs. Exclude them.
        self.mask[np.isnan(self.data)] = False
        
        # The mask potentially contains statistical outliers.
        # Optionally exclude them.
        if exclude_outliers is True:
            m, s = med_abs_deviation(self.data[self.mask])  # Compute median and median absolute deviation
            s *= 1.4826  # Convert MAD to std
            vmin = m - self.sigrej*s  # Minimum value to keep
            vmax = m + self.sigrej*s  # Maximum value to keep

            # Temporarily change NaNs to inf, so that they don't cause problems.
            self.data[np.isnan(self.data)] = np.inf

            # Flag statistical outliers
            bdpx = np.array(np.where(np.logical_or(self.data < vmin, self.data > vmax),
                                     1, 0), dtype=np.float32)
            self.data[np.isinf(self.data)] = np.nan  # Restore NaNs

            bdpx[np.logical_not(self.mask)] = 0  # We don't need to worry about non-background pixels
            # Also flag 4 nearest neighbors
            bdpx = bdpx +\
                        np.roll(bdpx, (+1,0), axis=(0,1)) +\
                        np.roll(bdpx, (-1,0), axis=(0,1)) +\
                        np.roll(bdpx, (0,+1), axis=(0,1)) +\
                        np.roll(bdpx, (0,-1), axis=(0,1))
            # bdpx now contains the pixels to exclude from the background pixels
            # mask. Exclude them.
            self.mask[bdpx != 0] = False
            self.data -= m  # STUB - Median subtract
        
        # Build the apodizing filter. This has unity gain at low frequency to capture 1/f. It 
        # also has unity gain at Nyquist to capture alternating column noise.
        # At mid frequencies, the gain is zero.  Unity gain at low frequencies.
        self.apodizer = np.zeros(len(self.rfftfreq), dtype=np.float32)
        self.apodizer[self.rfftfreq < self.fc[0]] = 1.0
        # Cosine roll-off
        here = np.logical_and(self.fc[0] <= self.rfftfreq, self.rfftfreq < self.fc[1])
        self.apodizer[here] = 0.5*(np.cos((np.pi/(self.fc[1]-self.fc[0])) *
                                          (self.rfftfreq[here] - self.fc[0]))+1)
        # Cosine roll-on
        here = np.logical_and(self.fc[2] <= self.rfftfreq, self.rfftfreq < self.fc[3])
        self.apodizer[here] = 0.5*(np.cos((np.pi/(self.fc[3]-self.fc[2])) *
                                          (self.rfftfreq[here] - self.fc[2]) + np.pi)+1)
        # Unity gain between f[3] and end
        self.apodizer[self.rfftfreq >= self.fc[3]] = 1.0

        
    def fit(self, return_fit=False, weight_fit=False):
        """
        Fit a background model to the data.
        
        Parameters
        ----------
        return_fit : bool
            Return the Fourier transform.

        weight_fit : bool
            Use weighted least squares as described in the NSClean paper.
            Turn off by default. For subarrays it is TBD if this is necessary.

        Returns
        -------
        rfft : numpy array
            The computed Fourier transform.
        """
        
        # To build the incomplete Fourier matrix, we require the index of each
        # clock tick of each valid pixel in the background samples. For consistency with
        # numpy's notation, we call this 'm' and require it to be a column vector.
        _x = np.arange(self.nx).reshape((1,-1))
        _y = np.arange(self.ny).reshape((-1,1))
        m = (_y*(self.nx+self.nloh) + _x)[self.mask].reshape((-1,1))
        
        # Define which Fourier vectors to fit. For consistency with numpy, call this k.
        k = np.arange(len(self.rfftfreq))[self.apodizer>0.].reshape((1,-1))
        
        # Build the incomplete Fourier matrix
        B = np.array(np.exp(2*np.pi*1J*m*k/self.n)/m.shape[0], dtype=np.complex64)
        
        # Weighted NSClean fitting
        if weight_fit:
            
            # Build the weight matrix. Weight by the reciprocal of the local background
            # sample density in time. Roughly approximate the local density, P, using
            # the reciprocal of convolution by a Gaussian kernel.
            _x = np.arange(self.n, dtype=np.float32) # x-values for building a 1-D Gaussian
            _mu = np.float32(self.n//2+1) # Center point of Gaussian
            _sigma = np.float32(self.weights_kernel_sigma) # Standard deviation of Gaussian
            W = np.exp(-(_x-_mu)**2/_sigma**2/2)/_sigma/np.sqrt(2*np.pi) # Build centered Gaussian
            FW = np.fft.rfft(np.fft.ifftshift(W)) # Forward FFT
            _M = np.hstack((self.mask, np.zeros((self.ny,self.nloh), dtype=np.bool_))).flatten() # Add new line overhead to mask
            with np.errstate(divide='ignore'):
                P = 1/np.fft.irfft(np.fft.rfft(np.array(_M, dtype=np.float32)) * FW, self.n) # Compute weights
            P = P[_M] # Keep only background samples
        
            # NSClean's weighting requires the Moore-Penrose invers of A = P*B.
            #     $A^+ = (A^H A)^{-1} A^H$
            A = P.reshape((-1,1)) * B # P is diagonal. Hadamard product is most RAM efficient
            AH = np.conjugate(A.transpose()) # Hermitian transpose of A
            pinv_PB = np.matmul(np.linalg.inv(np.matmul(AH, A)), AH)
            
        else:
            # Unweighted fit
            pinvB = np.linalg.pinv(B)
        
        # Solve for the (approximate) Fourier transform of the background samples.
        rfft = np.zeros(len(self.rfftfreq), dtype=np.complex64)
        if weight_fit is True:
            rfft[self.apodizer>0.] = np.matmul(pinv_PB * P.reshape((1,-1)), self.data[self.mask])
        else:
            rfft[self.apodizer>0.] = np.matmul(pinvB, self.data[self.mask])
        
        # Numpy requires that the forward transform multiply
        # the data by n. Correct normalization.
        rfft *= self.n / m.shape[0]
        
        # Invert the apodized Fourier transform to build the background model for this integration
        self.model = np.fft.irfft(rfft*self.apodizer, self.n).reshape((self.ny,-1))[:,:self.nx]
        
        # Done
        if return_fit:
            return(rfft)
        
        
    def clean(self, weight_fit=True):
        """
        Clean the data
        
        Parameters
        ----------
        weight_fit : bool
            Use weighted least squares as described in the NSClean paper.
            Otherwise, it is a simple unweighted fit.

        Returns
        -------
        data : 2D float array
            The cleaned data array.
        """ 
        self.fit(weight_fit=weight_fit)  # Fit the background model
        self.data -= self.model  # Overwrite data with cleaned data
        return(self.data)
