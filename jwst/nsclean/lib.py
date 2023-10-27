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
    
    def __init__(self, detector, mask, fc=1/(2*2048/16), kw=1/(2*2048/16)/4,
                 buffer_sigma=1.5, sigrej=3.0, weights_kernel_sigma=32):
        """
        JWST NIRSpec background modeling and subtraction -AKA "clean" (NSClean)

        Parameters
        ----------
        detector : string
            A string selected from {'NRS1','NRS2'}

        mask : boolean array
            The background model is fitted to pixels set to True.
            Pixels set to False are ignored.

        fc : float
            Critical frequency. This is the 1/2 power point for NSClean's
            low pass filter. The units are 1/pixels.

        kw : float
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

        weights_kernel_sigma : integer
            Used to assign weights for MASK mode fitting.
        """
        self.detector = detector
        self.M = mask
        self.fc = fc
        self.kw = kw
        self.buffer_sigma = buffer_sigma
        self.sigrej = sigrej
        self.weights_kernel_sigma = weights_kernel_sigma
        
        # Transpose and flip mask to detector coordinates with the IRS2
        # zipper running along the bottom as displayed in ds9.
        if self.detector == 'NRS1':
            # NRS1 requires transpose only
            self.M = self.M.transpose()
        else:
            # NRS2 requires transpose and flip
            self.M = self.M.transpose()[::-1]

        self.ny = self.M.shape[0]
        self.nx = self.M.shape[1]

        # FFT frequencies
        self.rfftfreq = np.fft.rfftfreq(self.ny)
                    
        # Extract other necessary information from model_parameters and build the
        # low pass filter. This alse sets the number of Fourier vectors that
        # we need to project out.
        self.apodizer = np.array(make_lowpass_filter(self.fc, self.kw, self.nx))
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
            #self.P = 1 / np.fft.irfft2(np.fft.rfft2(np.array(self.M, dtype=np.float32)) * FW, (self.ny,self.nx))
            temp = np.fft.rfft2(np.array(self.M, dtype=np.float32))
            self.P = 1 / np.fft.irfft2(temp * FW, (self.ny, self.nx))
        self.P = np.where(self.M==True, self.P, 0.) # Illuminated areas carry no weight

        # Build a 1-dimensional Gaussian kernel for "buffing". Buffing is in the
        # dispersion direction only. In detector coordinates, this is axis zero. Even though
        # the kernel is 1-dimensional, we must still use a 2-dimensional array to 
        # represent it. I tried broadcasting a vector, but that made a kernel 2048 
        # columns wide (in detector space).
        _y = np.arange(self.ny)
        #_mu = self.ny//2 + 1
        _mu = self.nx//2 + 1
        _sigma = self.buffer_sigma  # Standard deviation of kernel
        _gkern = np.exp(-((_y-_mu) / _sigma)**2 / 2) / _sigma / np.sqrt(2*np.pi)  # Centered kernel as a vector
        gkern = np.zeros((self.ny, self.nx), dtype=np.float32)  # 2D kernel template
        gkern[:, _mu] = _gkern  # Copy in the kernel. Normalization is already correct.
        gkern = np.fft.ifftshift(gkern)  # Shift for Numpy
        self.fgkern = np.array(np.fft.rfft2(gkern), dtype=np.complex64)  # FFT for fast convolution
        
        
    def fit(self, Data):
        """
        Fit a background model to the supplied frame of data.
        
        Parameters
        ----------
        Data : float array
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
            d = Data[y][self.M[y]]  # unmasked (useable) data
            p = np.diag(self.P[y][self.M[y]])  # Weights

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
            m = np.arange(self.nx)[self.M[y]].reshape((-1, 1))  # Must be a column vector to broadcast
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
            if self.kw > 0:
                rfft[:self.nvec] *= self.apodizer[:self.nvec]

            # Invert the FFT to build the background model for this line
            model[y] = np.fft.irfft(rfft, self.nx)
        
        # Done!
        return model
    

    def clean(self, Data, buff=True):
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
        Data : array_like
            The input data. This should be the normal end result of Stage 1 processing.

        buff : boolean
            "Buff" the fitted spectrum by applying a slight Gaussian blur
            in the spectral dispersion direction.

        Returns
        -------
        Data : array_like
            The data, but with less striping and the background subtracted.
        """
        
        # Transform the data to detector space with the IRS2 zipper running along the bottom.
        if self.detector == 'NRS2':
            # Transpose and flip for NRS2
            Data = Data.transpose()[::-1]
        else:
            # Transpose (no flip) for NRS1
            Data = Data.transpose()
            
        # Fit the background model
        Bkg = self.fit(Data)  # Background model

        # Buff, if requested
        if buff:
            Bkg = np.fft.irfft2(np.fft.rfft2(Bkg) * self.fgkern, s=Bkg.shape)

        # Subtract the background model from the data
        Data -= Bkg
        
        # Transform back to DMS space
        if self.detector=='NRS2':
            Data = Data[::-1].transpose()
        else:
            Data = Data.transpose()
            
        # Done
        return Data 


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

    n : integer
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
