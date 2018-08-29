#! /usr/bin/env python

import sys
import os
import os.path
import string
import types
import numpy as np
import math
from astropy.io import fits as pyfits

from tkinter import *
import tkinter.simpledialog as tkSimpleDialog
import tkinter.filedialog as tkFileDialog
import tkinter.messagebox as tkMessageBox

from . import makepsf

try:
    from matplotlib import pylab
    plot_enabled = True
except:
    plot_enabled = False

ANGSTROMStoMICRONS = 0.0001
NANOMETERStoMICRONS = 0.001
METERStoMICRONS = 1.e6

# for specifying the data type for calculations in makepsf.MakePsf
DOUBLE_PREC = 64
SINGLE_PREC = 32


if os.path.dirname(os.path.abspath(__file__)).split('/')[-1] == 'lib':
    USEWORKINGDIR = "../data"
else:
    USEWORKINGDIR = str(os.path.dirname(__file__)) + "/data"

class CreatePSF:

    def __init__(self, root):

        self.root = root

        self.spectrum_file = None       # name of file containing spectrum
        self.spectrum = None            # tuple of wavelength and flux
        self.filter_file = None         # name of file containing filter
        self.filter_name = None         # filter name (e.g. F090W)
        self.filter = None              # tuple of wavelength and throughput

        # These two are used for getting the spectrum from a Kurucz model,
        # specified by spectral type.
        self.spectral_type = None       # spectral type, if relevant
        self.kurucz_directory = None    # directory for Kurucz data

        self.psf = None                 # makepsf.MakePsf object
        self.output = None              # output file, if specified

        # this dictionary may have these entries:
        #   "wl_info":  (wl_min, wl_max, nbands)
        #      wl_min is the minimum wavelength (microns) in the bandpass
        #      wl_max is the maximum wavelength (microns) in the bandpass
        #      nbands is the number of slices into which the bandpass
        #          will be divided
        #   "bandpass":  (wavelengths, weights)
        #      wavelengths is an array of wavelengths (microns)
        #      weights is an array of weights (spectrum * throughput)
        self.integration_parameters = {
            "bandpass": ((1.,), (1.,))}             # default

        # default values
        self.psf_parameters = {
            "diameter": 6.5,
            "oversample": 4,
            "type": DOUBLE_PREC,
            "output_size": 512,
            "pixel_size": 0.034,
            "verbose": True}

        self.instrument_pixels = {
            "NIRCam": 0.034,
            "NIRSpec": 0.10,
            "MIRI": 0.11}

        self.frame = Frame(root)
        self.frame.grid()

        self.get_instrument = Button(self.frame, text="specify instrument",
                                     command=self.getInstrument)
        self.get_instrument.grid(row=0, column=0, sticky=W)

        self.get_spectrum = Button(self.frame, text="specify spectrum",
                                    command=self.getSpectrum, state=DISABLED)
        self.get_spectrum.grid(row=1, column=0, sticky=W)

        self.get_filter = Button(self.frame, text="specify filter",
                                  command=self.getFilter, state=DISABLED)
        self.get_filter.grid(row=2, column=0, sticky=W)

        # Option to set the wavelength range.
        self.integ_param = Button(self.frame,
                                   text="set wavelength range",
                                   command=self.setIntegParam, state=DISABLED)
        self.integ_param.grid(row=3, column=0, sticky=W)

        # Option to set the parameters for computing the PSF.
        self.set_psf_param = Button(self.frame,
                                     text="set PSF parameters",
                                     command=self.setPsfParam, state=DISABLED)
        self.set_psf_param.grid(row=4, column=0, sticky=W)

        # include a dummy label to separate sections of the frame
        for_spacing = Label(self.frame, text="")
        for_spacing.grid(row=5, column=0, sticky=W)

        # The QUIT button.
        self.button = Button(self.frame, text="QUIT", fg="red",
                              command=self.exitQuit)
        self.button.grid(row=6, column=0, sticky=W)

        # Status line.
        self.status_label = Label(self.frame,
                                   text="First choose science instrument", fg="red")
        self.status_label.grid(row=7, column=0, columnspan=3, sticky=W)

        self.plot_spectrum = Button(self.frame, text="plot spectrum",
                                     command=self.plotSpectrum)
        self.plot_spectrum.grid(row=0, column=1, sticky=W)

        self.plot_filter = Button(self.frame, text="plot filter",
                                   command=self.plotFilter)
        self.plot_filter.grid(row=1, column=1, sticky=W)

        self.plot_bandpass = Button(self.frame, text="plot bandpass",
                                     command=self.plotBandpass)
        self.plot_bandpass.grid(row=2, column=1, sticky=W)

        self.print_info = Button(self.frame, text="print info",
                                  command=self.printInfo)
        self.print_info.grid(row=3, column=1, sticky=W)

        # Create the PSF.
        self.create_psf = Button(self.frame, text="create the PSF",
                                  command=self.createJwstPsf)
        self.create_psf.grid(row=0, column=2, sticky=W)

        # Save computed PSF.
        self.save_psf = Button(self.frame, text="write PSF to file",
                                command=self.writePSF)
        self.save_psf.grid(row=1, column=2, sticky=W)

        # Display the PSF.
        self.show_psf = Button(self.frame, text="display the PSF",
                                command=self.displayPSF)
        self.show_psf.grid(row=2, column=2, sticky=W)

        # Display the PSF profiles.
        self.show_psf = Button(self.frame, text="display profiles",
                                command=self.displayprof)
        self.show_psf.grid(row=3, column=2, sticky=W)

    def exitQuit(self):
        if self.psf is not None and not self.psf.output_written:
            save_it = tkMessageBox.askyesno("PSF not saved",
                                             "The PSF has not been saved to a file.\n" \
                                             "Do you want to save it?", default="yes")
            if save_it:
                self.writePSF()
        self.frame.quit()

    def showStatus(self, message):
        """Display message on the status line."""
        self.status_label.config(text=message)
        self.status_label.update_idletasks()

    def getInstrument(self):
        '''Choose science instrument. This modifies choice of OPD files
        and filter selection by pointing to different data folders'''
        self.instop = Toplevel(self.frame)
        self.insframe = Frame(self.instop)
        self.insframe.place(x=0, y=5)
        self.insframe.grid()

        self.insname = StringVar()
        self.insname.set(None)

        Label(self.insframe, text="Choose a science instrument").grid(row=0, column=0, sticky=W)

        Radiobutton(self.insframe, text="NIRCam", variable=self.insname,
                    value="NIRCam", command=self.ichosen).grid(row=1, column=0, sticky=W)

        Radiobutton(self.insframe, text="NIRSpec", variable=self.insname,
                    value="NIRSpec", command=self.ichosen).grid(row=2, column=0, sticky=W)

        Radiobutton(self.insframe, text="MIRI", variable=self.insname,
                    value="MIRI", command=self.ichosen).grid(row=3, column=0, sticky=W)

        Button(self.insframe, text="OK", fg='red', command=self.accept).grid(row=4, column=0)

    def ichosen(self):
        self.instrument = self.insname.get()
        print("You chose", self.instrument)

    def accept(self):
        self.instrument = self.insname.get()
        if self.instrument in ["NIRCam", "NIRSpec", "MIRI"]:
            print("Value ", self.instrument, " accepted")
            self.get_spectrum.config(state=NORMAL)
            self.get_filter.config(state=NORMAL)  # Activate choose filter button
            self.set_psf_param.config(state=NORMAL)
            self.integ_param.config(state=NORMAL)
            self.set_psf_param.config(state=NORMAL) # and PSF parameters
            self.psf_parameters["pixel_size"] = self.instrument_pixels[self.instrument]
            self.instop.destroy()
            self.showStatus("Now choose other parameters")
            if plot_enabled:
                pylab.close('all') # Clear plotting area
        else:
            print("No instrument chosen")





    def getSpectrum(self):
        """Get spectrum name, and read into memory."""

        # Select how the spectrum is to be specified.
        self.spectop = Toplevel(self.frame)
        self.specframe = Frame(self.spectop)
        self.specframe.grid()
        self.specframe.place(x=0, y=0)

        Button(self.specframe, text="specify a spectral type",
                command=self.specFromSpectralType).grid(row=0, column=0, sticky=W)
        Button(self.specframe, text="specify a file name",
                command=self.specFromFileName).grid(row=1, column=0, sticky=W)
        # Button (self.specframe, text="analytic spectrum",
        #     command=self.analyticSpectrum).grid (row=2, column=0, sticky=W)
        Button(self.specframe, text="CANCEL", fg="red",
                command=self.getSpecCancel).grid(row=3, column=0, sticky=W)

        self.spectop.grab_set()
        self.root.wait_window(self.spectop)

        if self.spectrum:
            if self.filter:
                self.plotBandpass()
            else:
                self.plotSpectrum()

        if self.spectrum and self.filter is None:
            self.showStatus("specify filter")
        else:
            self.showStatus("")

    def specFromSpectralType(self):
        """Get spectrum specified by spectral type."""

        # User specifies spectral type and Kurucz directory.
        if self.setSpectralType() is None:
            self.spectop.destroy()
            self.showStatus("spectrum not specified")
            return

        # The keys are spectral types; the values are tuples of the
        # table name and column name.
        Kurucz_filename = {
            "O3V": ("kp00_50000.fits", "g50"),
            "O5V": ("kp00_45000.fits", "g50"),
            "O6V": ("kp00_40000.fits", "g45"),
            "O8V": ("kp00_35000.fits", "g40"),
            "O5I": ("kp00_40000.fits", "g45"),
            "O6I": ("kp00_40000.fits", "g45"),
            "O8I": ("kp00_34000.fits", "g40"),
            "B0V": ("kp00_30000.fits", "g40"),
            "B3V": ("kp00_19000.fits", "g40"),
            "B5V": ("kp00_15000.fits", "g40"),
            "B8V": ("kp00_12000.fits", "g40"),
            "B0III": ("kp00_29000.fits", "g35"),
            "B5III": ("kp00_15000.fits", "g35"),
            "B0I": ("kp00_26000.fits", "g30"),
            "B5I": ("kp00_14000.fits", "g25"),
            "A0V": ("kp00_9500.fits", "g40"),
            "A5V": ("kp00_8250.fits", "g45"),
            "A0I": ("kp00_9750.fits", "g20"),
            "A5I": ("kp00_8500.fits", "g20"),
            "F0V": ("kp00_7250.fits", "g45"),
            "F5V": ("kp00_6500.fits", "g45"),
            "F0I": ("kp00_7750.fits", "g20"),
            "F5I": ("kp00_7000.fits", "g15"),
            "G0V": ("kp00_6000.fits", "g45"),
            "G5V": ("kp00_5750.fits", "g45"),
            "G0III": ("kp00_5750.fits", "g30"),
            "G5III": ("kp00_5250.fits", "g25"),
            "G0I": ("kp00_5500.fits", "g15"),
            "G5I": ("kp00_4750.fits", "g10"),
            "K0V": ("kp00_5250.fits", "g45"),
            "K5V": ("kp00_4250.fits", "g45"),
            "K0III": ("kp00_4750.fits", "g20"),
            "K5III": ("kp00_4000.fits", "g15"),
            "K0I": ("kp00_4500.fits", "g10"),
            "K5I": ("kp00_3750.fits", "g05"),
            "M0V": ("kp00_3750.fits", "g45"),
            "M2V": ("kp00_3500.fits", "g45"),
            "M5V": ("kp00_3500.fits", "g50"),
            "M0III": ("kp00_3750.fits", "g15"),
            "M0I": ("kp00_3750.fits", "g00"),
            "M2I": ("kp00_3500.fits", "g00")}

        startype = self.spectral_type.upper()
        (fname, gcol) = Kurucz_filename[startype]
        fullname = os.path.join(USEWORKINGDIR, self.kurucz_directory, fname[0:4], fname)
        self.spectrum_file = fullname
        fd = pyfits.open(fullname)
        hdu = fd[1]
        data = hdu.data
        fd.close()
        factor = self.wavelengthUnits(hdu, "WAVELENGTH")
        wave = data.field('WAVELENGTH') * factor
        flux = data.field(gcol)
        self.spectrum = (wave, flux)
        self.spectral_type = startype

        status = self.computeBandpass()

        self.spectop.destroy()

    def setSpectralType(self):
        """User specifies spectral type.

        This uses getSpectralType.  If the user cancels, this function
        will return None; otherwise, the function value will be the
        spectral type.  Both the spectral type and directory for
        Kurucz data (a default value) will be assigned to attributes.
        """

        global kurucz_dir

        kurucz_dir = "k93models"
        self.kurucz_directory = kurucz_dir

        #x = GetKuruczDir (self.frame,
        #    title="specify directory with Kurucz data")
        #if x.kurucz_directory is None:
        #    return None
        #self.kurucz_directory = x.kurucz_directory

        self.getSpectralType()

        return self.spectral_type

    def getSpectralType(self):

        self.star_top = Toplevel(self.frame)
        self.starframe = Frame(self.star_top)
        self.starframe.grid()
        # self.starframe.place (x=0, y=0)

        self.star = StringVar()
        self.star.set("G0V")

        Radiobutton(self.starframe, text="O3V", variable=self.star,
                    value="O3V", command=self.chosen).grid(row=0, column=0, sticky=W)
        Radiobutton(self.starframe, text="O5V", variable=self.star,
                    value="O5V", command=self.chosen).grid(row=0, column=1, sticky=W)
        Radiobutton(self.starframe, text="O6V", variable=self.star,
                    value="O6V", command=self.chosen).grid(row=0, column=2, sticky=W)
        Radiobutton(self.starframe, text="O8V", variable=self.star,
                    value="O8V", command=self.chosen).grid(row=0, column=3, sticky=W)
        Radiobutton(self.starframe, text="O5I", variable=self.star,
                    value="O5I", command=self.chosen).grid(row=0, column=4, sticky=W)
        Radiobutton(self.starframe, text="O6I", variable=self.star,
                    value="O6I", command=self.chosen).grid(row=0, column=5, sticky=W)
        Radiobutton(self.starframe, text="O8I", variable=self.star,
                    value="O8I", command=self.chosen).grid(row=0, column=6, sticky=W)

        Radiobutton(self.starframe, text="B0V", variable=self.star,
                    value="B0V", command=self.chosen).grid(row=1, column=0, sticky=W)
        Radiobutton(self.starframe, text="B3V", variable=self.star,
                    value="B3V", command=self.chosen).grid(row=1, column=1, sticky=W)
        Radiobutton(self.starframe, text="B5V", variable=self.star,
                    value="B5V", command=self.chosen).grid(row=1, column=2, sticky=W)
        Radiobutton(self.starframe, text="B8V", variable=self.star,
                    value="B8V", command=self.chosen).grid(row=1, column=3, sticky=W)
        Radiobutton(self.starframe, text="B0III", variable=self.star,
                    value="B0III", command=self.chosen).grid(row=1, column=4, sticky=W)
        Radiobutton(self.starframe, text="B5III", variable=self.star,
                    value="B5III", command=self.chosen).grid(row=1, column=5, sticky=W)
        Radiobutton(self.starframe, text="B0I", variable=self.star,
                    value="B0I", command=self.chosen).grid(row=1, column=6, sticky=W)
        Radiobutton(self.starframe, text="B5I", variable=self.star,
                    value="B5I", command=self.chosen).grid(row=1, column=7, sticky=W)

        Radiobutton(self.starframe, text="A0V", variable=self.star,
                    value="A0V", command=self.chosen).grid(row=2, column=0, sticky=W)
        Radiobutton(self.starframe, text="A5V", variable=self.star,
                    value="A5V", command=self.chosen).grid(row=2, column=1, sticky=W)
        Radiobutton(self.starframe, text="A0I", variable=self.star,
                    value="A0I", command=self.chosen).grid(row=2, column=2, sticky=W)
        Radiobutton(self.starframe, text="A5I", variable=self.star,
                    value="A5I", command=self.chosen).grid(row=2, column=3, sticky=W)

        Radiobutton(self.starframe, text="F0V", variable=self.star,
                    value="F0V", command=self.chosen).grid(row=3, column=0, sticky=W)
        Radiobutton(self.starframe, text="F5V", variable=self.star,
                    value="F5V", command=self.chosen).grid(row=3, column=1, sticky=W)
        Radiobutton(self.starframe, text="F0I", variable=self.star,
                    value="F0I", command=self.chosen).grid(row=3, column=2, sticky=W)
        Radiobutton(self.starframe, text="F5I", variable=self.star,
                    value="F5I", command=self.chosen).grid(row=3, column=3, sticky=W)

        Radiobutton(self.starframe, text="G0V", variable=self.star,
                    value="G0V", command=self.chosen).grid(row=4, column=0, sticky=W)
        Radiobutton(self.starframe, text="G5V", variable=self.star,
                    value="G5V", command=self.chosen).grid(row=4, column=1, sticky=W)
        Radiobutton(self.starframe, text="G0III", variable=self.star,
                    value="G0III", command=self.chosen).grid(row=4, column=2, sticky=W)
        Radiobutton(self.starframe, text="G5III", variable=self.star,
                    value="G5III", command=self.chosen).grid(row=4, column=3, sticky=W)
        Radiobutton(self.starframe, text="G0I", variable=self.star,
                    value="G0I", command=self.chosen).grid(row=4, column=4, sticky=W)
        Radiobutton(self.starframe, text="G5I", variable=self.star,
                    value="G5I", command=self.chosen).grid(row=4, column=5, sticky=W)

        Radiobutton(self.starframe, text="K0V", variable=self.star,
                    value="K0V", command=self.chosen).grid(row=5, column=0, sticky=W)
        Radiobutton(self.starframe, text="K5V", variable=self.star,
                    value="K5V", command=self.chosen).grid(row=5, column=1, sticky=W)
        Radiobutton(self.starframe, text="K0III", variable=self.star,
                    value="K0III", command=self.chosen).grid(row=5, column=2, sticky=W)
        Radiobutton(self.starframe, text="K5III", variable=self.star,
                    value="K5III", command=self.chosen).grid(row=5, column=3, sticky=W)
        Radiobutton(self.starframe, text="K0I", variable=self.star,
                    value="K0I", command=self.chosen).grid(row=5, column=4, sticky=W)
        Radiobutton(self.starframe, text="K5I", variable=self.star,
                    value="K5I", command=self.chosen).grid(row=5, column=5, sticky=W)

        Radiobutton(self.starframe, text="M0V", variable=self.star,
                    value="M0V", command=self.chosen).grid(row=6, column=0, sticky=W)
        Radiobutton(self.starframe, text="M2V", variable=self.star,
                    value="M2V", command=self.chosen).grid(row=6, column=1, sticky=W)
        Radiobutton(self.starframe, text="M5V", variable=self.star,
                    value="M5V", command=self.chosen).grid(row=6, column=2, sticky=W)
        Radiobutton(self.starframe, text="M0III", variable=self.star,
                    value="M0III", command=self.chosen).grid(row=6, column=3, sticky=W)
        Radiobutton(self.starframe, text="M0I", variable=self.star,
                    value="M0I", command=self.chosen).grid(row=6, column=4, sticky=W)
        Radiobutton(self.starframe, text="M2I", variable=self.star,
                    value="M2I", command=self.chosen).grid(row=6, column=5, sticky=W)

        self.button = Button(self.starframe, text="OK", fg="red",
                           command=self.spTypeSelected)
        self.button.grid(row=7, column=0)

        self.star_top.grab_set()
        self.root.wait_window(self.star_top)

    def chosen(self):
        print("You wish to use spectrum", self.star.get(), "?")

    def spTypeSelected(self):
        self.spectral_type = self.star.get()
        self.star_top.destroy()

    def specFromFileName(self):
        """Get spectrum from a user-specified file.

        Note that this assumes the first two columns are wavelength and flux.
        """

        self.spectrum_file = tkFileDialog.askopenfilename(
            title="select file containing spectrum",
            parent=self.frame)
        if not self.spectrum_file:
            return

        fd = pyfits.open(self.spectrum_file)
        hdu = fd[1]
        data = hdu.data
        fd.close()
        factor = self.wavelengthUnits(hdu, "WAVELENGTH")
        try:
            wavelength = data.field("WAVELENGTH") * factor
            flux = data.field(1)
        except (IndexError, AttributeError):
            tkMessageBox.showwarning("unknown table format",
                                      "table must contain at least two columns, wavelength and flux")
            return

        self.spectrum = (wavelength, flux)

        status = self.computeBandpass()

        self.spectop.destroy()

    def analyticSpectrum(self):
        """Create a spectrum from user-specified parameters."""

        print("analytic spectrum (not implemented yet)")

        status = self.computeBandpass()

        self.spectop.destroy()

    def getSpecCancel(self):
        self.spectop.destroy()

    def getFilter(self):

        # Select how the filter is to be specified.
        self.filtertop = Toplevel(self.frame)
        self.filterframe = Frame(self.filtertop)
        self.filterframe.grid()
        self.filterframe.place(x=0, y=0)

        Button(self.filterframe, text="specify a filter name",
                command=self.filterFromFilter).grid(row=0, column=0, sticky=W)
        Button(self.filterframe, text="specify a file name",
                command=self.filterFromFileName).grid(row=1, column=0, sticky=W)
        # Button (self.filterframe, text="filter parameters",
        #     command=self.analyticFilter).grid (row=2, column=0, sticky=W)
        Button(self.filterframe, text="CANCEL", fg="red",
                command=self.getFilterCancel).grid(row=3, column=0, sticky=W)

        self.filtertop.grab_set()
        self.root.wait_window(self.filtertop)

        if self.filter:
            if self.spectrum:
                self.plotBandpass()
            else:
                self.plotFilter()

        if self.filter and self.spectrum is None:
            self.showStatus("specify spectrum")
        else:
            self.showStatus("")

    def filterFromFilter(self):
        """Get filter from a user-specified filter name.

        The name of the filter throughput file will be constructed
        from the filter name.
        """

        self.filtname_top = Toplevel(self.frame)
        self.filtnameframe = Frame(self.filtname_top)
        self.filtnameframe.grid()
        # self.filtnameframe.place (x=0, y=0)

        NIRCam_filters = ["F070W", "F090W", "F110W", "F140M", "F150W", "F163M",
                          "F187N", "F200W", "F210M", "F241N", "F250M", "F256N",
                          "F270W", "F300M", "F335M", "F357W", "F360M", "F390M",
                          "F405N", "F430M", "F444W", "F460M", "F469N", "F480M"]

        MIRI_filters = ["F560W", "F770W", "F1000W", "F1130W", "F1280W", "F1500W",
                        "F1800W", "F2100W", "F2550W"]

        NIRSpec_filters = ["F070LP", "F100LP", "F170LP", "F290LP", "F115W", "F140W"]

        self.filtname_var = StringVar()
        self.filtname_var.set("F999W")

        if self.instrument == "NIRCam":
            for f in range(len(NIRCam_filters)):
                fname = NIRCam_filters[f]
                ncol = f / 6
                nrow = f % 6
                Radiobutton(self.filtnameframe, text=fname,
                            variable=self.filtname_var, value=fname,
                            command=self.filtchosen).grid(row=nrow, column=ncol, sticky=W)


        if self.instrument == "MIRI":
            for f in range(len(MIRI_filters)):
                fname = MIRI_filters[f]
                ncol = f / 6
                nrow = f % 6
                Radiobutton(self.filtnameframe, text=fname,
                            variable=self.filtname_var, value=fname,
                            command=self.filtchosen).grid(row=nrow, column=ncol, sticky=W)

        if self.instrument == "NIRSpec":
            for f in range(len(NIRSpec_filters)):
                fname = NIRSpec_filters[f]
                ncol = f / 6
                nrow = f % 6
                Radiobutton(self.filtnameframe, text=fname,
                            variable=self.filtname_var, value=fname,
                            command=self.filtchosen).grid(row=nrow, column=ncol, sticky=W)



        self.button = Button(self.filtnameframe, text="OK", fg="red",
                           command=self.filtNameSelected)
        self.button.grid(row=7, column=0)

        self.filtname_top.grab_set()
        self.root.wait_window(self.filtname_top)

    def filtchosen(self):
        print("You wish to use filter", self.filtname_var.get(), "?")

    def filtNameSelected(self):

        self.filter_name = self.filtname_var.get()
        fullname = os.path.join(USEWORKINGDIR, self.instrument, "filters", self.filter_name + "_thru.fits")
        self.filter_file = fullname
        print("Filter file", fullname)

        # read self.filter_file into memory (self.filter, etc)
        self.readFilter()

        status = self.computeBandpass()
        self.filtname_top.destroy()
        self.filtertop.destroy()


    def filterFromFileName(self):
        """Get filter from a user-specified file.

        Note that this assumes the throughput is the second column
        (rather than assuming a name for throughput).
        """

        self.filter_file = tkFileDialog.askopenfilename(
            title="select file containing filter",
            parent=self.frame)
        if not self.filter_file:
            return

        self.readFilter()

        status = self.computeBandpass()
        self.filtertop.destroy()

    def readFilter(self):

        fd = pyfits.open(self.filter_file)
        hdu = fd[1]
        data = hdu.data
        fd.close()
        factor = self.wavelengthUnits(hdu, "WAVELENGTH")
        try:
            wavelength = data.field("WAVELENGTH") * factor
            throughput = data.field(1)
        except (IndexError, AttributeError):
            tkMessageBox.showwarning("unknown table format",
                                      "table must contain at least two columns, wavelength and throughput")
            return

        self.filter = (wavelength, throughput)

        # Find range of throughput band
        FRACTIONofPEAK = 0.01
        peak = throughput.max()
        nt = len(throughput)
        first = 0
        last = nt - 1
        for i in range(nt):
            if throughput[i] > FRACTIONofPEAK * peak:
                first = i
                break
        for i in range(nt - 1, 0, -1):
            if throughput[i] > FRACTIONofPEAK * peak:
                last = i
                break

        # Override the min and max wavelengths, but leave the number
        # of bands unchanged.
        if "wl_info" in self.integration_parameters:
            (wl_min, wl_max, nbands) = self.integration_parameters["wl_info"]
        else:
            nbands = 5
        wl_min = wavelength[first]
        wl_max = wavelength[last]
        self.integration_parameters["wl_info"] = (wl_min, wl_max, nbands)

    def analyticFilter(self):

        print("filter from parameters (not implemented yet)")

        status = self.computeBandpass()

        self.filtertop.destroy()

    def getFilterCancel(self):
        self.filtertop.destroy()

    def plotSpectrum(self, hold=False):

        if not plot_enabled:
            self.showStatus("plotting is not enabled")
            return

        if self.spectrum is None:
            tkMessageBox.showinfo("no spectrum",
                                   "you must specify a spectrum before you can plot it")
            return

        (wavelength, flux) = self.spectrum
        # normalize the flux to 1 for plotting
        maxval = flux.max()
        flux /= maxval
        pylab.figure(1)
        pylab.plot(wavelength, flux, hold=hold)
        pylab.xlim(0.0, 5.0) # Plot range from 0 to 5 microns
        if self.filter is not None and \
           "wl_info" in self.integration_parameters:
            self.setXlim()
        pylab.ylim(0., 1.05)
        if self.spectral_type is not None:
            pylab.title(self.spectral_type)
        pylab.xlabel('Wavelength (microns)')
        pylab.ylabel(r'$ergs.cm^{-2}.sec^{-1}.A^{-1}$')
        pylab.grid(True)

    def plotFilter(self, hold=False):

        if not plot_enabled:
            self.showStatus("plotting is not enabled")
            return

        if self.filter is None:
            tkMessageBox.showinfo("no filter",
                                   "you must specify a filter before you can plot it")
            return

        (wavelength, throughput) = self.filter
        pylab.figure(1)
        pylab.plot(wavelength, throughput, hold=hold)
        if "wl_info" in self.integration_parameters:
            self.setXlim()
        if not hold:
            maxval = throughput.max()
            pylab.ylim(0., 1.05 * maxval)
        pylab.xlabel('Wavelength (microns)')
        pylab.ylabel('Throughput')
        pylab.grid(True)

    def plotBandpass(self):

        if not plot_enabled:
            self.showStatus("plotting is not enabled")
            return

        if self.computeBandpass() > 0:
            return

        if self.spectrum is not None:
            self.plotSpectrum()
        if self.filter is not None:
            self.plotFilter(hold=True)
        (wavelength, weight) = self.integration_parameters["bandpass"]
        pylab.figure(1)
        pylab.plot(wavelength, weight, hold=True)
        pylab.plot(wavelength, weight, "k+", hold=True)
        if "wl_info" in self.integration_parameters:
            self.setXlim()
        pylab.xlabel('Wavelength (microns)')
        pylab.ylabel('Weight (normalized)')
        pylab.title(self.spectral_type + ' and ' + self.filter_name)
        pylab.grid(True)

    def setXlim(self):
        """Set plot limits in X, based on integration_parameters."""

        if "wl_info" in self.integration_parameters:
            (wl_min, wl_max, nbands) = self.integration_parameters["wl_info"]
            if wl_max > wl_min:
                xmin = wl_min - (wl_max - wl_min) / 2.
                xmin = max(xmin, 0.)
                xmax = wl_max + (wl_max - wl_min) / 2.
            else:
                xmin = 0.95 * wl_max
                xmax = 1.05 * wl_max
            pylab.xlim(xmin, xmax)

    def setIntegParam(self):
        """User specifies wavelength range and spacing.

        This uses GetWlMinMax.
        """

        global wl_min, wl_max, nbands

        if "wl_info" in self.integration_parameters:
            (wl_min, wl_max, nbands) = self.integration_parameters["wl_info"]
        else:
            # arbitrary default values
            (wl_min, wl_max, nbands) = (0.9, 1.1, 20)

        x = GetWlMinMax(self.frame, title="give wavelength info (microns)")
        if x.wl_max is None:
            return

        wl_min = x.wl_min
        wl_max = x.wl_max
        nbands = x.nbands

        self.integration_parameters["wl_info"] = (wl_min, wl_max, nbands)

        if self.computeBandpass() > 0:
            return

        self.showStatus("wavelength parameters set")

    def setPsfParam(self):
        """User specifies parameters for makepsf.MakePsf.

        This uses GetPsfParameters.
        """

        global diameter, oversample, type, output_size, pixel_size, verbose

        diameter = self.psf_parameters["diameter"]
        oversample = self.psf_parameters["oversample"]
        type = self.psf_parameters["type"]
        output_size = self.psf_parameters["output_size"]
        pixel_size = self.psf_parameters["pixel_size"]
        verbose = self.psf_parameters["verbose"]

        x = GetPsfParameters(self.frame, title="give PSF parameters")
        if x.diameter is None:
            return

        self.psf_parameters["diameter"] = x.diameter
        self.psf_parameters["oversample"] = x.oversample
        self.psf_parameters["type"] = x.type
        self.psf_parameters["output_size"] = x.output_size
        self.psf_parameters["pixel_size"] = x.pixel_size
        self.psf_parameters["verbose"] = x.verbose

        self.showStatus("PSF parameters set")

    def printInfo(self):

        print("Version", __version__)

        if "wl_info" in self.integration_parameters:
            (wl_min, wl_max, nbands) = self.integration_parameters["wl_info"]
            delta_wl = (wl_max - wl_min) / nbands
            print("bandpass:  " \
                  "%d bands, from %.4f to %.4f, step %.4f (microns)" % \
                  (nbands, wl_min, wl_max, delta_wl))

        if "bandpass" in self.integration_parameters:
            (wavelength, weight) = self.integration_parameters["bandpass"]
            print("bandpass wavelengths (microns):")
            print(wavelength)
            print("bandpass weights:")
            print(weight)
        if self.spectrum_file is not None:
            print("spectrum was read from", self.spectrum_file)
        if self.spectral_type is not None:
            print("spectrum is for spectral type", self.spectral_type)
        if self.filter_file is not None:
            print("filter was read from", self.filter_file)

        print("PSF parameters:")
        print("  pupil diameter =", self.psf_parameters["diameter"], "(meters)")
        oversample = self.psf_parameters["oversample"]
        if oversample == 2:
            print("  oversampling factor = 2 (Nyquist sampling)")
        else:
            r = float(oversample) / 2.
            print("  oversampling factor = %d (%g * Nyquist sampling)" % \
                  (oversample, r))
        if self.psf_parameters["type"] == SINGLE_PREC:
            print("  computations will use single precision")
        else:
            print("  computations will use double precision")
        print("  size of output image =", self.psf_parameters["output_size"])
        print("  output pixel size =", \
              self.psf_parameters["pixel_size"] / oversample, "(arcsec)")
        print("  verbose =", self.psf_parameters["verbose"])
        if self.psf is not None:
            print("pupil file =", self.psf.pupil_file)
            print("phase file =", self.psf.phase_file)
            if self.output is None:
                print("no output file has been specified")
            else:
                print("output file =", self.output)

        if self.psf is None:
            self.showStatus("PSF has not been created yet")
        elif not self.psf.output_written:
            self.showStatus("Note:  PSF has not been saved")
        else:
            self.showStatus("")

    def createJwstPsf(self):
        """Call the MakePsf routine in makepsf."""

        diameter = self.psf_parameters["diameter"]
        oversample = self.psf_parameters["oversample"]
        if self.psf_parameters["type"] == SINGLE_PREC:
            type = np.float32
        else:
            type = np.float64
        output_size = self.psf_parameters["output_size"]
        pixel_size = self.psf_parameters["pixel_size"]
        verbose = self.psf_parameters["verbose"]

        if self.computeBandpass() > 0:
            return

        self.showStatus("creating the PSF ...")
        self.psf = None

        pupil_file = tkFileDialog.askopenfilename(
            title="select pupil image",
            initialdir=USEWORKINGDIR,
            initialfile="pupil.fits",
            parent=self.frame)
        if not pupil_file:
            self.showStatus("PSF not created")
            return
        phase_file = tkFileDialog.askopenfilename(
            title="Select phase map image",
            initialdir=USEWORKINGDIR + "/" + self.instrument + "/OPD/",
            parent=self.frame)
        if not phase_file:
            self.showStatus("PSF not created")
            return

        # All output will go into a 'Results' directory created in the directory the program was started in.
        if (os.listdir('.')).count('Results') == 0:
            os.mkdir('./Results')
        self.output = tkFileDialog.asksaveasfilename(initialdir="./Results",
                                                      title="Name of file for output PSF", parent=self.frame)

        if not self.output:
            self.output = None
        else:
            # if the user did not specify extension .fits or .fit, append .fits
            if not (self.output.endswith(".fits") or
                    self.output.endswith(".fit")):
                self.output = self.output + ".fits"
        if self.output is not None and os.access(self.output, os.F_OK):
            try:
                os.remove(self.output)
            except:
                tkMessageBox.showwarning("can't delete file",
                                          "can't delete %s" % self.output)
                self.showStatus("PSF not created")
                return

        try:
            self.psf = makepsf.MakePsf(self.instrument, pupil_file, phase_file,
                                        output=self.output,
                                        diameter=diameter, oversample=oversample, type=type,
                                        filter=self.integration_parameters["bandpass"],
                                        output_size=output_size, pixel_size=pixel_size,
                                        verbose=verbose)
        except Exception as message:
            tkMessageBox.showwarning("exception", message)
            self.showStatus("PSF not created")
            return

        if self.psf.output_written:
            self.displayPSF()
            self.showStatus("PSF has been created and saved")
        else:
            self.displayPSF(data=self.psf.integrated_psf)
            self.showStatus("PSF has been created but not saved")

    def writePSF(self):
        if self.psf is None:
            tkMessageBox.showerror("no PSF",
                                    "you must compute a PSF before you can save it")
            return
        self.output = tkFileDialog.asksaveasfilename(title="output file",
                                                      parent=self.frame)
        if not self.output:
            self.showStatus("PSF not saved")
            return
        if os.access(self.output, os.F_OK):
            try:
                os.remove(self.output)
            except:
                tkMessageBox.showwarning("can't delete file",
                                          "can't delete %s" % os.path.basename(self.output))
                self.showStatus("PSF not saved")
                return
        self.psf.writeto(self.output)
        if self.psf.output_written:
            self.showStatus("PSF has been saved")
        else:
            self.showStatus("PSF not saved")
            tkMessageBox.showwarning("PSF not saved",
                                      "Couldn't write PSF to %s" % os.path.basename(self.output))

    def displayPSF(self, data=None):
        """Read the PSF from the file and display it."""

        if not plot_enabled:
            self.showStatus("plotting is not enabled")
            return

        if self.psf is None:
            self.showStatus("create the PSF")
            tkMessageBox.showwarning("PSF hasn't been created",
                                      "PSF must be created before you can display it.")
            return

        if data is None:
            if self.psf.integrated_psf is None:
                fd = pyfits.open(self.output)
                image = fd[0].data
                fd.close()
            else:
                image = self.psf.integrated_psf
        else:
            image = data

        # take the log of the image to show a wider range
        maxval = image.max()
        image = np.log10(image + maxval / 1000.0) - np.log10(1.001 * maxval)
        pylab.figure(1)
        pylab.clf()
        try:         # Originl matplotlib aspect option "preserve"
            pylab.imshow(image, aspect="preserve", interpolation="nearest",
                          origin="lower", hold=False)
        except:      # New option is "equal"
            pylab.imshow(image, aspect="equal", interpolation="nearest",
                          origin="lower", hold=False)

        if self.output is not None:
            pylab.title(self.output)
        pylab.colorbar()
        pylab.figtext(0.8, 0.05, "log base 10")


    def displayprof(self, data=None):
        """Read the PSF from the file and display radial profile etc."""

        if not plot_enabled:
            self.showStatus("plotting is not enabled")
            return

        if self.psf is None:
            self.showStatus("create the PSF")
            tkMessageBox.showwarning("PSF hasn't been created",
                                      "PSF must be created before you can display it.")
            return

        output_size = self.psf_parameters["output_size"]
        pixel_size = self.psf_parameters["pixel_size"]

        if data is None:
            if self.psf.integrated_psf is None:
                fd = pyfits.open(self.output)
                spot = fd[0].data
                headlist = fd[(0)].header
                output_size = headlist['NAXIS1']
                pixel_size = headlist['CDELT1'] * 3600.0
                fd.close()
            else:
                spot = self.psf.integrated_psf

        else:
            spot = data

        print('Pixel size', pixel_size, ' arcsec')
        self.showStatus("computing profiles ...")
        rscale = pixel_size / 0.01   # Bin sum into 0.01 arsec steps
        diag = output_size / math.sqrt(2.0) + 1
        rmax = int(diag * rscale)
        xcen = output_size / 2 - 1
        ycen = output_size / 2 - 1

        radial = np.zeros((rmax), np.float32)  # Make up zeroed array of proper length
        iradial = np.zeros((rmax), np.float32) # Also one for intensity
        count = np.zeros((rmax))

        for x in range(output_size):
            for y in range(output_size):
                r = math.hypot(x - xcen, y - ycen)
                ir = int(r * rscale)
                count[ir] = count[ir] + 1
                radial[ir] = radial[ir] + spot[y, x]  # Build up radial sum
                iradial[ir] = ((count[ir] - 1) * iradial[ir] + spot[y, x]) / count[ir]

        radial = radial / np.sum(radial)
        iradial = iradial / np.sum(iradial)

        # Generate encircled energy plot
        accum = np.zeros((rmax), np.float32)
        accum[0] = radial[0]
        for r in range(1, rmax):
            accum[r] = accum[r - 1] + radial[r]

        xaxis = np.zeros((rmax), np.float32)
        for x in range(rmax):
            xaxis[x] = 0.01 * x  # binsize 0.01 arcsec

        pylab.figure(2)
        pylab.clf()
        pylab.title('Radial Intensity')
        pylab.grid(True)
        pylab.semilogy(xaxis[0:100], iradial[0:100] / iradial[0], lw=2)
        pylab.xlabel('Arcseconds')


        pylab.figure(3)
        pylab.clf()
        pylab.title('Encircled Energy')
        pylab.plot(xaxis[0:100], accum[0:100], lw=2)
        pylab.grid(True)
        pylab.xlabel('Arcseconds')


        pylab.figure(4)
        pylab.clf()
        pylab.title('X and Y Profiles')
        half = output_size / 2
        extent = int(2.0 / pixel_size) # 2 arcsec
        if extent > half - 1:
            extent = half - 1
        x = range(half - extent, half + extent)
        arcsec = np.zeros((output_size), np.float32)
        for ix in x:
            arcsec[ix] = pixel_size * (ix - half)

        pylab.grid(True)
        pylab.semilogy(arcsec[x], spot[half, x] / spot[half, half], lw=2)
        pylab.semilogy(arcsec[x], spot[x, half] / spot[half, half], lw=2)
        pylab.legend(('X Profile', 'Y profile'))
        pylab.xlabel('Arcseconds')

        self.showStatus("profiles completed")


    def computeBandpass(self):
        """Compute the wavelengths and weights for the bandpass.

        The function value should be 0; it will be 1 if the range of
        wavelengths misses either the spectrum or filter.
        """

        if self.spectrum is None:
            # a flat, very wide spectrum
            if self.filter is None:
                # NOTE:  wavelength is assumed to be in microns
                spec_wl = [1.e-9, 1.e9]
                photons = [1., 1.]
            else:
                spec_wl = self.filter[0]
                photons = np.ones(len(spec_wl), dtype=np.float64)
        else:
            spec_wl = self.spectrum[0]
            photons = self.spectrum[0] * self.spectrum[1]

        if self.filter is None:
            if self.spectrum is None:
                # NOTE:  wavelength is assumed to be in microns
                filt_wl = [1.e-9, 1.e9]
                throughput = [1., 1.]
            else:
                filt_wl = self.spectrum[0]
                throughput = np.ones(len(filt_wl), dtype=np.float64)
        else:
            filt_wl = self.filter[0]
            throughput = self.filter[1]

        if "wl_info" in self.integration_parameters:
            (wl_min, wl_max, nbands) = self.integration_parameters["wl_info"]
        else:
            (wl_min, wl_max, nbands) = (0.9, 1.1, 1)            # microns

        bwidth = (wl_max - wl_min) / nbands
        if bwidth <= 0.:
            wavelength = [(wl_min + wl_max) / 2.]
            weight = [1.]
            self.integration_parameters["bandpass"] = (wavelength, weight)
            self.showStatus("bandpass has been computed")
            return 0

        # Make this a little larger than necessary (we need nbands elements).
        temp = np.arange(wl_min, wl_max + 2. * bwidth, bwidth, dtype=np.float64)
        llow = temp[0:nbands]
        lhigh = llow + bwidth
        wavelength = (llow + lhigh) / 2.

        weight = np.zeros(nbands, dtype=np.float64)
        for nb in range(nbands):
            if lhigh[nb] < spec_wl[0] or lhigh[nb] < filt_wl[0] or \
               llow[nb] > spec_wl[-1] or llow[nb] > filt_wl[-1]:
                weight[nb] = 0.
            else:
                weight[nb] = self.sumband(spec_wl, photons,
                                           filt_wl, throughput, llow[nb], lhigh[nb])

        sum = weight.sum()
        if sum == 0.:
            tkMessageBox.showwarning("wavelengths out of range",
                                      "range of wavelengths misses either spectrum or filter")
            self.showStatus("Note:  set wavelength range")
            return 1
        else:
            weight /= weight.sum()
        self.integration_parameters["bandpass"] = (wavelength, weight)

        self.showStatus("bandpass has been computed")
        return 0

    def sumband(self, spec_wl, photons, filt_wl, throughput,
                 lambda1, lambda2):
        """Integrate the product of photon spectrum and throughput

        arguments:
        spec_wl            wavelength array for spectrum
        photons            array of spectrum fluxes, expressed in photons
        filt_wl            wavelength array for filter
        throughput         array of throughputs for filter
        lambda1, lambda2   range of wavelengths for current band
        """

        lenph = len(photons)
        rn = range(lenph)

        r1 = 0
        r2 = lenph - 1
        for r in rn:
            if spec_wl[r] <= lambda1:
                r1 = r
            if spec_wl[r] >= lambda2:
                r2 = r
                break

        lenwl = len(filt_wl)
        rt = range(lenwl)

        s1 = 0
        s2 = lenwl - 1
        for s in rt:
            if filt_wl[s] <= lambda1:
                s1 = s
            if filt_wl[s] >= lambda2:
                s2 = s
                break

        if r2 - r1 >= s2 - s1:      # Source spectrum has finer grid.
            lfine = spec_wl     # Wavelength
            fine = photons      # Value
            lf1 = r1            # Index
            lf2 = r2
            lcoarse = filt_wl
            coarse = throughput
            lc1 = s1
            lc2 = s2
        else:                   # Filter throughput has finer grid
            lfine = filt_wl
            fine = throughput
            lf1 = s1
            lf2 = s2
            lcoarse = spec_wl
            coarse = photons
            lc1 = r1
            lc2 = r2

        # End corrections
        # Values at lambda1
        fine1 = (fine[lf1] * (lfine[lf1 + 1] - lambda1) + fine[lf1 + 1] * (lambda1 - lfine[lf1])) / \
              (lfine[lf1 + 1] - lfine[lf1])
        coarse1 = (coarse[lc1] * (lcoarse[lc1 + 1] - lambda1) + coarse[lc1 + 1] * (lambda1 - lcoarse[lc1])) / \
                (lcoarse[lc1 + 1] - lcoarse[lc1])
        prod1 = fine1 * coarse1

        # Values at lfine[lf1]
        coarself1 = (coarse[lc1] * (lcoarse[lc1 + 1] - lfine[lf1]) + coarse[lc1 + 1] * (lfine[lf1] - lcoarse[lc1])) / \
                  (lcoarse[lc1 + 1] - lcoarse[lc1])
        prodlf1 = fine[lf1] * coarself1
        e1 = prodlf1 * (lfine[lf1 + 1] - lambda1) - prod1 * (lambda1 - lfine[lf1])


        # Values at lambda2
        fine2 = (fine[lf2 - 1] * (lfine[lf2] - lambda2) + fine[lf2] * (lambda2 - lfine[lf2 - 1])) / \
              (lfine[lf2] - lfine[lf2 - 1])
        coarse2 = (coarse[lc2 - 1] * (lcoarse[lc2] - lambda2) + coarse[lc2] * (lambda2 - lcoarse[lc2 - 1])) / \
                (lcoarse[lc2] - lcoarse[lc2 - 1])
        prod2 = fine2 * coarse2

        # Values at lfine[lf2]
        coarself2 = (coarse[lc2 - 1] * (lcoarse[lc2] - lfine[lf2]) + coarse[lc2] * (lfine[lf2] - lcoarse[lc2 - 1])) / \
                  (lcoarse[lc2] - lcoarse[lc2 - 1])
        prodlf2 = fine[lf2] * coarself2


        e2 = prodlf2 * (lambda2 - lfine[lf2 - 1]) - prod2 * (lfine[lf2] - lambda2)

        # xxx Note:  Need to fix this!
        # f and c here can index outside the bounds of the arrays.
        if lf2 - lf1 >= 2:   # normal trapezoidal integration

            sumprod = 0.0  # Trapezoidal integration
            for f in range(lf1 + 1, lf2):
                for c in range(lc1 + 1, lc2 + 2):
                    if lcoarse[c] > lfine[f]: # Do interpolation between m and m-1
                        cint = ((lcoarse[c] - lfine[f]) * coarse[c - 1] + \
                                (lfine[f] - lcoarse[c - 1]) * coarse[c]) / \
                             (lcoarse[c] - lcoarse[c - 1])
                        sumprod = sumprod + cint * fine[f] * \
                                (lfine[f + 1] - lfine[f - 1])
                        break # First qualifying value of c only for each value of f


            integral = 0.5 * (sumprod + e1 + e2)

        else:  #Only two points across integral
            integral = 0.5 * (prod1 + prod2) * (lambda2 - lambda1)
        return integral

    def wavelengthUnits(self, hdu, column):
        """Interpret the units string in the table header.

        The function value will be the multiplicative factor needed
        to convert the wavelengths to microns.  If no units are
        specified (or can't be interpreted) for the wavelength column,
        the units will be assumed to be Angstroms.
        """

        coldefs = hdu.get_coldefs()
        if isinstance(column, int):
            column_units = coldefs.units[column]
        else:
            column = column.lower()
            found = False
            for i in range(len(coldefs.names)):
                column_name = coldefs.names[i].lower()
                if column_name == column:
                    column_units = coldefs.units[i]
                    found = True
                    break
            if not found:
                print("warning:  can't find %s column" % column)
                return ANGSTROMStoMICRONS

        if column_units is None:
            units = "angstrom"          # default
        else:
            units = column_units.lower()
        if units == "a" or units == "angstrom" or units == "angstroms":
            factor = ANGSTROMStoMICRONS
        elif units == "nm" or units == "nanometer" or units == "nanometers":
            factor = NANOMETERStoMICRONS
        elif units == "micron" or units == "microns":
            factor = 1.
        elif units == "m" or units == "meter" or units == "meters":
            factor = METERStoMICRONS
        else:
            print(" wavelength units '%s' not given; " \
                  "Angstroms assumed" % column_units)
            factor = ANGSTROMStoMICRONS

        return factor

class GetKuruczDir (tkSimpleDialog.Dialog):

    def body(self, master):

        global kurucz_dir
        self.kurucz_directory = None

        Label(master, text="Kurucz directory:").grid(row=0)

        self.e1 = Entry(master)

        self.e1.grid(row=0, column=1)

        e1 = StringVar()
        e1.set(kurucz_dir)
        self.e1["textvariable"] = e1

        return self.e1

    def apply(self):

        self.kurucz_directory = self.e1.get()

class GetWlMinMax (tkSimpleDialog.Dialog):

    def body(self, master):

        global wl_min, wl_max, nbands
        self.wl_min = None
        self.wl_max = None
        self.nbands = None

        Label(master, text="lower wavelength:").grid(row=0)
        Label(master, text="upper wavelength:").grid(row=1)
        Label(master, text="number of bands:").grid(row=2)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)

        e1 = StringVar()
        e1.set(str(wl_min))
        self.e1["textvariable"] = e1

        e2 = StringVar()
        e2.set(str(wl_max))
        self.e2["textvariable"] = e2

        e3 = StringVar()
        e3.set(str(nbands))
        self.e3["textvariable"] = e3

        return self.e1

    def apply(self):

        self.wl_min = string.atof(self.e1.get())
        self.wl_max = string.atof(self.e2.get())
        if self.wl_min <= 0. or self.wl_max <= 0.:
            tkMessageBox.showwarning("bad wavelength",
                                      "wavelengths must be greater than zero")
            self.wl_min = 1.            # microns
            self.wl_max = 1.

        self.nbands = string.atoi(self.e3.get())
        if self.nbands < 1:
            tkMessageBox.showwarning("bad nbands",
                                      "number of bands must be greater than zero")
            self.nbands = 1

class GetPsfParameters (tkSimpleDialog.Dialog):

    def body(self, master):

        global diameter, oversample, type, output_size, pixel_size, verbose
        self.diameter = None
        self.oversample = None
        self.type = None
        self.output_size = None
        self.pixel_size = None
        self.verbose = None

        Label(master, text="pupil diameter (meters):").grid(row=0)
        Label(master, text="oversampling factor:").grid(row=1)
        Label(master, text="single or double precision:").grid(row=2)
        Label(master, text="output image size:").grid(row=3)
        Label(master, text="pixel size (arcsec):").grid(row=4)
        Label(master, text="print more info?:").grid(row=5)

        self.e1 = Entry(master)
        self.e2 = Entry(master)
        self.e3 = Entry(master)
        self.e4 = Entry(master)
        self.e5 = Entry(master)
        self.e6 = Entry(master)

        self.e1.grid(row=0, column=1)
        self.e2.grid(row=1, column=1)
        self.e3.grid(row=2, column=1)
        self.e4.grid(row=3, column=1)
        self.e5.grid(row=4, column=1)
        self.e6.grid(row=5, column=1)

        e1 = StringVar()
        e1.set(str(diameter))
        self.e1["textvariable"] = e1

        e2 = StringVar()
        e2.set(str(oversample))
        self.e2["textvariable"] = e2

        e3 = StringVar()
        if type == SINGLE_PREC:
            e3.set("single")
        else:
            e3.set("double")
        self.e3["textvariable"] = e3

        e4 = StringVar()
        e4.set(str(output_size))
        self.e4["textvariable"] = e4

        e5 = StringVar()
        e5.set(str(pixel_size))
        self.e5["textvariable"] = e5

        e6 = StringVar()
        e6.set(str(verbose))
        self.e6["textvariable"] = e6

        return self.e1

    def apply(self):

        self.diameter = float(self.e1.get())
        if self.diameter <= 0.:
            tkMessageBox.showwarning("invalid diameter",
                                      "the pupil diameter must be greater than zero")
            self.diameter = 6.5

        self.oversample = int(self.e2.get())
        if self.oversample < 1:
            self.oversample = 1

        type = self.e3.get().lower()
        if type == "double" or type == "d" or type == "64":
            self.type = DOUBLE_PREC
        elif type == "single" or type == "s" or type == "32":
            self.type = SINGLE_PREC
        else:
            tkMessageBox.showwarning("unknown data type",
                                      "data type should be 'single' or 'double'")
            self.type = DOUBLE_PREC

        self.output_size = int(self.e4.get())
        if self.output_size < 4:
            self.output_size = 4

        self.pixel_size = float(self.e5.get())

        v = self.e6.get()
        if v == "True" or v == "true" or v == "yes" or v == "1":
            self.verbose = True
        elif v == "False" or v == "false" or v == "no" or v == "0":
            self.verbose = False
        else:
            tkMessageBox.showwarning("boolean parameter required",
                                      "'%s' not understood as boolean" % v)
            self.verbose = False

def jwpsf():

    root = Tk()
    root.title('JWPSF')
    psf = CreatePSF(root)
    root.mainloop()

if __name__ == "__main__":

    jwpsf()
