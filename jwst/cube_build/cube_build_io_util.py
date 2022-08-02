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
    cube parameter reference file and fill in the corresponding dictionary in
    instrument_info

    Parameters
    ----------
    par_filename : str
       cube parameter reference filename
    instrument : str
        Either MIRI or NIRSPEC
    weighting : str
        Type of weighting, msm, emem or drizzle
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
                    if this_channel == table_channel and this_sub == table_band:
                        instrument_info.SetSpatialSize(table_spaxelsize, this_channel, this_sub)
                        instrument_info.SetSpectralStep(table_spectralstep, this_channel, this_sub)
                        instrument_info.SetWaveMin(table_wavemin, this_channel, this_sub)
                        instrument_info.SetWaveMax(table_wavemax, this_channel, this_sub)
                #  modified Shepard method 1/r weighting
                if weighting == 'msm':
                    for tabdata in ptab.ifucubepars_msm_table:
                        table_channel = tabdata['channel']
                        table_band = tabdata['band'].lower()
                        table_sroi = tabdata['ROISPATIAL']
                        table_wroi = tabdata['ROISPECTRAL']
                        table_power = tabdata['POWER']
                        table_softrad = tabdata['SOFTRAD']
                    # match on this_channel and this_sub
                        if this_channel == table_channel and this_sub == table_band:
                            instrument_info.SetMSM(this_channel, this_sub,
                                                   table_sroi, table_wroi,
                                                   table_power, table_softrad)

                #  modified Shepard method e^-r weighting
                elif weighting == 'emsm':
                    for tabdata in ptab.ifucubepars_emsm_table:
                        table_channel = tabdata['channel']
                        table_band = tabdata['band'].lower()
                        table_sroi = tabdata['ROISPATIAL']
                        table_wroi = tabdata['ROISPECTRAL']
                        table_scalerad = tabdata['SCALERAD']
                    # match on this_channel and this_sub
                        if this_channel == table_channel and this_sub == table_band:
                            instrument_info.SetEMSM(this_channel, this_sub,
                                                    table_sroi, table_wroi,
                                                    table_scalerad)

            #  read in wavelength table for multi-band data
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
            #  read in wavelength table for modified Shepard method 1/r weighting
            elif weighting == 'emsm':
                for tabdata in ptab.ifucubepars_multichannel_emsm_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    table_sroi = tabdata['ROISPATIAL']
                    table_wroi = tabdata['ROISPECTRAL']
                    table_scalerad = tabdata['SCALERAD']
                    instrument_info.SetMultiChannelEMSMTable(table_wave, table_sroi,
                                                             table_wroi, table_scalerad)
            elif weighting == 'drizzle':
                for tabdata in ptab.ifucubepars_multichannel_driz_wavetable:
                    table_wave = tabdata['WAVELENGTH']
                    instrument_info.SetMultiChannelDrizTable(table_wave)

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

                    if this_gwa == table_grating and this_filter == table_filter:
                        instrument_info.SetSpatialSize(table_spaxelsize, this_gwa, this_filter)
                        instrument_info.SetSpectralStep(table_spectralstep, this_gwa, this_filter)
                        instrument_info.SetWaveMin(table_wavemin, this_gwa, this_filter)
                        instrument_info.SetWaveMax(table_wavemax, this_gwa, this_filter)

                #  modified Shepard method 1/r weighting
                if weighting == 'msm':
                    for tabdata in ptab.ifucubepars_msm_table:
                        table_grating = tabdata['DISPERSER'].lower()
                        table_filter = tabdata['FILTER'].lower()
                        table_sroi = tabdata['ROISPATIAL']
                        table_wroi = tabdata['ROISPECTRAL']
                        table_power = tabdata['POWER']
                        table_softrad = tabdata['SOFTRAD']

                        if this_gwa == table_grating and this_filter == table_filter:
                            instrument_info.SetMSM(this_gwa, this_filter,
                                                   table_sroi, table_wroi,
                                                   table_power, table_softrad)
                #  modified Shepard method e^-r weighting
                elif weighting == 'emsm':
                    for tabdata in ptab.ifucubepars_emsm_table:
                        table_grating = tabdata['DISPERSER'].lower()
                        table_filter = tabdata['FILTER'].lower()
                        table_sroi = tabdata['ROISPATIAL']
                        table_wroi = tabdata['ROISPECTRAL']
                        table_scalerad = tabdata['SCALERAD']

                        if this_gwa == table_grating and this_filter == table_filter:
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
