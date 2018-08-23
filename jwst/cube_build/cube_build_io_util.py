# Routines used for building cubes
from .. import datamodels
import logging
log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)

#********************************************************************************
# HELPER ROUTINES for CubeData class defined in cube_build.py
# these methods relate to I/O type procedures.
# read_offset_file
# read_cubepars
# read_resolution_file
#********************************************************************************
# Read in dither offset file
# For testing this is useful but possibily this might be useful during flight if the
# images need an additional offset applied to them
def read_offset_file(offset_list):

    ra_offset = []
    dec_offset = []
    f = open(offset_list, 'r')
    i = 0
    for line in f:
        offset_str = line.split()
        offset = [float(xy) for xy in offset_str]
        ra_off = offset[0]
        dec_off = offset[1]

        ra_offset.append(ra_off)
        dec_offset.append(dec_off)
        i = i + 1
    f.close()
    return ra_offset, dec_offset


#********************************************************************************
def read_cubepars(par_filename,
                  instrument,
                  all_channel,
                  all_subchannel,
                  all_grating,
                  all_filter,
                  instrument_info):
#********************************************************************************
    """
    Short Summary
    -------------
    Based on the instrument and channel/subchannels (MIRI) or grating/filter(NIRSPEC)
    that covers the full range of the data, read in the appropriate columns in the
    cube parameter reference file and fill in the cooresponding dicitionary in
    instrument_info

    Parameters
    ----------
    par_filename: cube parameter reference table
    instrument: Either MIRI or NIRSPEC
    all_channel: all the channels contained in input data
    all_subchannel: all subchannels contained in input data
    all_grating: all the gratings contained in the input data
    all_filter: all the filters contained in the input data
    instrument_info holds the defaults spatial scales, spectral scales, roi size,
    weighting parametes, and min and max wavelengths for each for each band

    Returns
    -------
    The correct elements of instrument_info are filled in

    """
    if instrument == 'MIRI':
        ptab = datamodels.MiriIFUCubeParsModel(par_filename)
        number_bands = len(all_channel)
        # pull out the channels and subcahnnels that cover the data making up the cube
        for i in range(number_bands):
            this_channel = all_channel[i]
            #compare_channel = 'CH'+this_channel
            this_sub = all_subchannel[i]
            # find the table entries for this combination
            for tabdata in ptab.ifucubepars_table:
                table_channel = tabdata['channel']
                table_band = tabdata['band'].lower()
                table_spaxelsize = tabdata['SPAXELSIZE']
                table_spectralstep = tabdata['SPECTRALSTEP']
                table_wavemin = tabdata['WAVEMIN']
                table_wavemax = tabdata['WAVEMAX']
                #match on this_channel and this_sub
                if(this_channel == table_channel and this_sub == table_band):
                    instrument_info.SetSpatialSize(table_spaxelsize, this_channel, this_sub)
                    instrument_info.SetSpectralStep(table_spectralstep, this_channel, this_sub)
                    instrument_info.SetWaveMin(table_wavemin, this_channel, this_sub)
                    instrument_info.SetWaveMax(table_wavemax, this_channel, this_sub)

            for tabdata in ptab.ifucubepars_msn_table:
                table_channel = tabdata['channel']
                table_band = tabdata['band'].lower()
                table_sroi = tabdata['ROISPATIAL']
                table_wroi = tabdata['ROISPECTRAL']
                table_power = tabdata['POWER']
                table_softrad = tabdata['SOFTRAD']
                #match on this_channel and this_sub
                if(this_channel == table_channel and this_sub == table_band):
                    instrument_info.SetSpatialROI(table_sroi, this_channel, this_sub)
                    instrument_info.SetWaveROI(table_wroi, this_channel, this_sub)
                    instrument_info.SetMSMPower(table_power, this_channel, this_sub)
                    instrument_info.SetSoftRad(table_softrad, this_channel, this_sub)
        for tabdata in ptab.ifucubepars_multichannel_wavetable:
            table_wave = tabdata['WAVELENGTH']
            table_sroi = tabdata['ROISPATIAL']
            table_wroi = tabdata['ROISPECTRAL']
            table_power = tabdata['POWER']
            table_softrad = tabdata['SOFTRAD']
            instrument_info.SetMultiChannelTable(table_wave, table_sroi,
                                              table_wroi, table_power,
                                              table_softrad)

#        print('Done reading cubepar reference file')
    elif instrument == 'NIRSPEC':
        ptab = datamodels.NirspecIFUCubeParsModel(par_filename)
        number_gratings = len(all_grating)

        for i in range(number_gratings):
            this_gwa = all_grating[i]
            this_filter = all_filter[i]
            for tabdata in ptab.ifucubepars_table:
                table_grating = tabdata['DISPERSER']
                table_filter = tabdata['FILTER']
                table_spaxelsize = tabdata['SPAXELSIZE']
                table_spectralstep = tabdata['SPECTRALSTEP']
                table_wavemin = tabdata['WAVEMIN']
                table_wavemax = tabdata['WAVEMAX']
#                print(table_grating,table_filter,table_spaxelsize,table_spectralstep,
#                      table_wavemin,table_wavemax)

                if(this_gwa == table_grating and this_filter == table_filter):
                    instrument_info.SetSpatialSize(table_spaxelsize, this_gwa, this_filter)
                    instrument_info.SetSpectralStep(table_spectralstep, this_gwa, this_filter)
                    instrument_info.SetWaveMin(table_wavemin, this_gwa, this_filter)
                    instrument_info.SetWaveMax(table_wavemax, this_gwa, this_filter)

            for tabdata in ptab.ifucubepars_msn_table:
                table_grating = tabdata['DISPERSER']
                table_filter = tabdata['FILTER']
                table_sroi = tabdata['ROISPATIAL']
                table_wroi = tabdata['ROISPECTRAL']
                table_power = tabdata['POWER']
                table_softrad = tabdata['SOFTRAD']

                if(this_gwa == table_grating and this_filter == table_filter):
                    instrument_info.SetSpatialROI(table_sroi, this_gwa, this_filter)
                    instrument_info.SetWaveROI(table_wroi, this_gwa, this_filter)
                    instrument_info.SetMSMPower(table_power, this_gwa, this_filter)
                    instrument_info.SetSoftRad(table_softrad, this_gwa, this_filter)

        for tabdata in ptab.ifucubepars_prism_wavetable:
            table_wave = tabdata['WAVELENGTH']
            table_sroi = tabdata['ROISPATIAL']
            table_wroi = tabdata['ROISPECTRAL']
            table_power = tabdata['POWER']
            table_softrad = tabdata['SOFTRAD']
            instrument_info.SetPrismTable(table_wave, table_sroi,
                                          table_wroi, table_power,
                                          table_softrad)

        for tabdata in ptab.ifucubepars_med_wavetable:
            table_wave = tabdata['WAVELENGTH']
            table_sroi = tabdata['ROISPATIAL']
            table_wroi = tabdata['ROISPECTRAL']
            table_power = tabdata['POWER']
            table_softrad = tabdata['SOFTRAD']
            instrument_info.SetMedTable(table_wave, table_sroi,
                                          table_wroi, table_power,
                                          table_softrad)

        for tabdata in ptab.ifucubepars_high_wavetable:
            table_wave = tabdata['WAVELENGTH']
            table_sroi = tabdata['ROISPATIAL']
            table_wroi = tabdata['ROISPECTRAL']
            table_power = tabdata['POWER']
            table_softrad = tabdata['SOFTRAD']
            instrument_info.SetHighTable(table_wave, table_sroi,
                                          table_wroi, table_power,
                                          table_softrad)
#        print('Done reading cubepar reference file')
#_______________________________________________________________________
# Read MIRI Resolution reference file
#********************************************************************************
def read_resolution_file(resol_filename,
                         channel,
                         all_channel,
                         all_subchannel,
                         instrument_info):
    """
    Short Summary
    -------------
    If this is MIRI data and the weighting is miripsf then read in the MIRI resolition
    file. Read weighting parameters based on which channel and subchannel we are working
    with. Fill in the appropriated values into instrument_info

    Parameters
    ----------
    resol_filename: MIRI resolution reference table
    channel: channels working with
    all_channel: all the channels contained in input data
    all_subchannel: all subchannels contained in input data
    instrument_info holds the  MIRI psf weighting parameters

    Returns
    -------
    The correct elements of instrument_info are filled in

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

    number_bands = len(channel)

        # pull out the channels and subcahnnels that cover the data making up the cube
    for i in range(number_bands):
        this_channel = all_channel[i]
        this_sub = all_subchannel[i]
        compare_band = this_channel + this_sub
        for tabdata in ptab.resolving_power_table:
            table_sub_band = tabdata['SUB_BAND']
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
            #match on this_channel and this_sub
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
#********************************************************************************
