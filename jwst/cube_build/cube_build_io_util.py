""" Read in reference files for the cube_build setp
"""
from .. import datamodels
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

def read_cubepars(par_filename,
                  instrument,
                  weighting,
                  all_channel,
                  all_subchannel,
                  all_grating,
                  all_filter,
                  instrument_info):
    """ Read in cube parameter reference file

    Based on the instrument and channel/subchannels (MIRI) or
    grating/filter(NIRSPEC), read in the appropriate columns in the
    cube parameter reference file and fill in the cooresponding dicitionary in
    instrument_info

    Parameters
    ----------
    par_filename : str
       cube parameter reference filename
    instrument : str
        Either MIRI or NIRSPEC
    weighting : str
        Type of weighting, msm, emem or miripsf
    all_channel : list
        all the channels contained in input data
    all_subchannel : list
        all subchannels contained in input data
    all_grating: list
        all the gratings contained in the input data
    all_filter: list
        all the filters contained in the input data
    instrument_info : dictionary
        holds the defaults spatial scales, spectral scales, roi size,
        weighting parameters, and min and max wavelengths for each
        for each band

    Returns
    -------
    The dictionary, instrument_info, is filled in for each band covered
    by the input data

    """
    if instrument == 'MIRI':
        with datamodels.MiriIFUCubeParsModel(par_filename) as ptab:
            number_bands = len(all_channel)
            # pull out the channels and subchannels that cover the cube
            for i in range(number_bands):
                this_channel = all_channel[i]
                # compare_channel = 'CH'+this_channel
                this_sub = all_subchannel[i]
                # find the table entries for this combination
                for tabdata in ptab.ifucubepars_table:
                    table_channel = tabdata['channel']
                    table_band = tabdata['band'].lower()
                    table_spaxelsize = tabdata['SPAXELSIZE']
                    table_spectralstep = tabdata['SPECTRALSTEP']
                    table_wavemin = tabdata['WAVEMIN']
                    table_wavemax = tabdata['WAVEMAX']
                    # match on this_channel and this_sub
                    if(this_channel == table_channel and this_sub == table_band):
                        instrument_info.SetSpatialSize(table_spaxelsize, this_channel, this_sub)
                        instrument_info.SetSpectralStep(table_spectralstep, this_channel, this_sub)
                        instrument_info.SetWaveMin(table_wavemin, this_channel, this_sub)
                        instrument_info.SetWaveMax(table_wavemax, this_channel, this_sub)
                #  modified shepard method 1/r weighting
                if weighting == 'msm':
                    for tabdata in ptab.ifucubepars_msm_table:
                        table_channel = tabdata['channel']
                        table_band = tabdata['band'].lower()
                        table_sroi = tabdata['ROISPATIAL']
                        table_wroi = tabdata['ROISPECTRAL']
                        table_power = tabdata['POWER']
                        table_softrad = tabdata['SOFTRAD']
                    # match on this_channel and this_sub
                        if(this_channel == table_channel and this_sub == table_band):
                            instrument_info.SetMSM(this_channel, this_sub,
                                                   table_sroi, table_wroi,
                                                   table_power, table_softrad)

                #  modified shepard method e^-r weighting
                elif weighting == 'emsm':
                    for tabdata in ptab.ifucubepars_emsm_table:
                        table_channel = tabdata['channel']
                        table_band = tabdata['band'].lower()
                        table_sroi = tabdata['ROISPATIAL']
                        table_wroi = tabdata['ROISPECTRAL']
                        table_scalerad = tabdata['SCALERAD']
                    # match on this_channel and this_sub
                        if(this_channel == table_channel and this_sub == table_band):
                            instrument_info.SetEMSM(this_channel, this_sub,
                                                    table_sroi, table_wroi,
                                                    table_scalerad)

            #  read in wavelength table for modified shepard method 1/r weighting
            if weighting == 'msm':
                for tabdata in ptab.ifucubepars_multichannel_msm_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    table_sroi = tabdata['ROISPATIAL']
                    table_wroi = tabdata['ROISPECTRAL']
                    table_power = tabdata['POWER']
                    table_softrad = tabdata['SOFTRAD']
                    instrument_info.SetMultiChannelTable(table_wave, table_sroi,
                                                         table_wroi, table_power,
                                                         table_softrad)
            #  read in wavelength table for modified shepard method 1/r weighting
            elif weighting == 'emsm' or weighting == 'miripsf':
                for tabdata in ptab.ifucubepars_multichannel_emsm_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    table_sroi = tabdata['ROISPATIAL']
                    table_wroi = tabdata['ROISPECTRAL']
                    table_scalerad = tabdata['SCALERAD']
                    instrument_info.SetMultiChannelEMSMTable(table_wave, table_sroi,
                                                             table_wroi, table_scalerad)

    # Read in NIRSPEC Values
    elif instrument == 'NIRSPEC':
        with datamodels.NirspecIFUCubeParsModel(par_filename) as ptab:
            number_gratings = len(all_grating)

            for i in range(number_gratings):
                this_gwa = all_grating[i]
                this_filter = all_filter[i]
                for tabdata in ptab.ifucubepars_table:
                    table_grating = tabdata['DISPERSER'].lower()
                    table_filter = tabdata['FILTER'].lower()
                    table_spaxelsize = tabdata['SPAXELSIZE']
                    table_spectralstep = tabdata['SPECTRALSTEP']
                    table_wavemin = tabdata['WAVEMIN']
                    table_wavemax = tabdata['WAVEMAX']

                    if(this_gwa == table_grating and this_filter == table_filter):
                        instrument_info.SetSpatialSize(table_spaxelsize, this_gwa, this_filter)
                        instrument_info.SetSpectralStep(table_spectralstep, this_gwa, this_filter)
                        instrument_info.SetWaveMin(table_wavemin, this_gwa, this_filter)
                        instrument_info.SetWaveMax(table_wavemax, this_gwa, this_filter)

                #  modified shepard method 1/r weighting
                if weighting == 'msm':
                    for tabdata in ptab.ifucubepars_msm_table:
                        table_grating = tabdata['DISPERSER'].lower()
                        table_filter = tabdata['FILTER'].lower()
                        table_sroi = tabdata['ROISPATIAL']
                        table_wroi = tabdata['ROISPECTRAL']
                        table_power = tabdata['POWER']
                        table_softrad = tabdata['SOFTRAD']

                        if(this_gwa == table_grating and this_filter == table_filter):
                            instrument_info.SetMSM(this_gwa, this_filter,
                                                   table_sroi, table_wroi,
                                                   table_power, table_softrad)
                #  modified shepard method e^-r weighting
                elif weighting == 'emsm':
                    for tabdata in ptab.ifucubepars_emsm_table:
                        table_grating = tabdata['DISPERSER'].lower()
                        table_filter = tabdata['FILTER'].lower()
                        table_sroi = tabdata['ROISPATIAL']
                        table_wroi = tabdata['ROISPECTRAL']
                        table_scalerad = tabdata['SCALERAD']

                        if(this_gwa == table_grating and this_filter == table_filter):
                            instrument_info.SetEMSM(this_gwa, this_filter,
                                                    table_sroi, table_wroi,
                                                    table_scalerad)

            # read in wavelength tables
            if weighting == 'msm':
                for tabdata in ptab.ifucubepars_prism_msm_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    table_sroi = tabdata['ROISPATIAL']
                    table_wroi = tabdata['ROISPECTRAL']
                    table_power = tabdata['POWER']
                    table_softrad = tabdata['SOFTRAD']
                    instrument_info.SetPrismTable(table_wave, table_sroi,
                                                  table_wroi, table_power,
                                                  table_softrad)

                for tabdata in ptab.ifucubepars_med_msm_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    table_sroi = tabdata['ROISPATIAL']
                    table_wroi = tabdata['ROISPECTRAL']
                    table_power = tabdata['POWER']
                    table_softrad = tabdata['SOFTRAD']
                    instrument_info.SetMedTable(table_wave, table_sroi,
                                                table_wroi, table_power,
                                                table_softrad)

                for tabdata in ptab.ifucubepars_high_msm_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    table_sroi = tabdata['ROISPATIAL']
                    table_wroi = tabdata['ROISPECTRAL']
                    table_power = tabdata['POWER']
                    table_softrad = tabdata['SOFTRAD']
                    instrument_info.SetHighTable(table_wave, table_sroi,
                                                 table_wroi, table_power,
                                                 table_softrad)

            elif weighting == 'emsm':
                for tabdata in ptab.ifucubepars_prism_emsm_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    table_sroi = tabdata['ROISPATIAL']
                    table_wroi = tabdata['ROISPECTRAL']
                    table_scalerad = tabdata['SCALERAD']
                    instrument_info.SetPrismEMSMTable(table_wave, table_sroi,
                                                      table_wroi, table_scalerad)

                for tabdata in ptab.ifucubepars_med_emsm_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    table_sroi = tabdata['ROISPATIAL']
                    table_wroi = tabdata['ROISPECTRAL']
                    table_scalerad = tabdata['SCALERAD']
                    instrument_info.SetMedEMSMTable(table_wave, table_sroi,
                                                    table_wroi, table_scalerad)

                for tabdata in ptab.ifucubepars_high_emsm_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    table_sroi = tabdata['ROISPATIAL']
                    table_wroi = tabdata['ROISPECTRAL']
                    table_scalerad = tabdata['SCALERAD']
                    instrument_info.SetHighEMSMTable(table_wave, table_sroi,
                                                     table_wroi, table_scalerad)
# _____________________________________________________________________________


def read_resolution_file(resol_filename,
                         all_channel,
                         all_subchannel,
                         instrument_info):
    """ Read in the MIRI Resolution reference file

    If this is MIRI data and the weighting is miripsf then read in the
    MIRI resolition file. Read weighting parameters based on which channel
    and subchannel we are working with and fill in the appropriated values
    into instrument_info.

    Parameters
    ----------
    resol_filename : str
        MIRI resolution reference table
    all_channel : list
        all the channels contained in input data
    all_subchannel : list
        all subchannels contained in input data
    instrument_info : dictionary
        holds the  MIRI psf weighting parameters

    Returns
    -------
    Resolution parameters are filled in in the instrumeent_info
    dectionary

    """

    ptab = datamodels.MiriResolutionModel(resol_filename)
    table_alpha_cutoff = ptab.psf_fwhm_alpha_table['A_CUTOFF']
    table_alpha_a_short = ptab.psf_fwhm_alpha_table['A_A_SHORT']
    table_alpha_b_short = ptab.psf_fwhm_alpha_table['A_B_SHORT']
    table_alpha_a_long = ptab.psf_fwhm_alpha_table['A_A_LONG']
    table_alpha_b_long = ptab.psf_fwhm_alpha_table['A_B_LONG']

    table_beta_cutoff = ptab.psf_fwhm_beta_table['B_CUTOFF']
    table_beta_a_short = ptab.psf_fwhm_beta_table['B_A_SHORT']
    table_beta_b_short = ptab.psf_fwhm_beta_table['B_B_SHORT']
    table_beta_a_long = ptab.psf_fwhm_beta_table['B_A_LONG']
    table_beta_b_long = ptab.psf_fwhm_beta_table['B_B_LONG']

    instrument_info.Set_psf_alpha_parameters(table_alpha_cutoff,
                                             table_alpha_a_short,
                                             table_alpha_b_short,
                                             table_alpha_a_long,
                                             table_alpha_b_long)

    instrument_info.Set_psf_beta_parameters(table_beta_cutoff,
                                            table_beta_a_short,
                                            table_beta_b_short,
                                            table_beta_a_long,
                                            table_beta_b_long)

    number_bands = len(all_channel)
    # pull out the channels and subcahnnels that the cube is make from
    for i in range(number_bands):
        this_channel = all_channel[i]
        this_sub = all_subchannel[i]
        compare_band = this_channel + this_sub
        for tabdata in ptab.resolving_power_table:
            table_sub_band = tabdata['SUB_BAND'].lower()
            table_wave_center = tabdata['R_CENTRE']
            table_res_a_low = tabdata['R_A_LOW']
            table_res_b_low = tabdata['R_B_LOW']
            table_res_c_low = tabdata['R_C_LOW']
            table_res_a_high = tabdata['R_A_HIGH']
            table_res_b_high = tabdata['R_B_HIGH']
            table_res_c_high = tabdata['R_C_HIGH']
            table_res_a_ave = tabdata['R_A_AVG']
            table_res_b_ave = tabdata['R_B_AVG']
            table_res_c_ave = tabdata['R_C_AVG']
            # match on this_channel and this_sub
            if compare_band == table_sub_band:
                instrument_info.Set_RP_Wave_Cutoff(table_wave_center,
                                                   this_channel, this_sub)

                instrument_info.Set_RP_low(table_res_a_low,
                                           table_res_b_low,
                                           table_res_c_low,
                                           this_channel, this_sub)

                instrument_info.Set_RP_high(table_res_a_high,
                                            table_res_b_high,
                                            table_res_c_high,
                                            this_channel, this_sub)

                instrument_info.Set_RP_ave(table_res_a_ave,
                                           table_res_b_ave,
                                           table_res_c_ave,
                                           this_channel, this_sub)
