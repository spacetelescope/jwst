import logging
import multiprocessing
from typing import Protocol, Union
import numpy as np

from stdatamodels.jwst import datamodels
from astropy.table import Table
import copy

from .observations import Observation
from .sens1d import get_photom_data

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


def determine_multiprocessing_ncores(max_cores: Union[str, int], num_cores) -> int:

    """Determine the number of cores to use for multiprocessing.

    Parameters
    ----------
    max_cores : string or int
        Number of cores to use for multiprocessing. If set to 'none'
        (the default), then no multiprocessing will be done. The other
        allowable string values are 'quarter', 'half', and 'all', which indicate
        the fraction of cores to use for multi-proc. The total number of
        cores includes the SMT cores (Hyper Threading for Intel).
        If an integer is provided, it will be the exact number of cores used.
    num_cores : int
        Number of cores available on the machine

    Returns
    -------
    ncpus : int
        Number of cores to use for multiprocessing
    """
    match max_cores:
        case 'none':
            return 1
        case None:
            return 1
        case 'quarter':
            return num_cores // 4 or 1
        case 'half':
            return num_cores // 2 or 1
        case 'all':
            return num_cores
        case int():
            if max_cores <= num_cores and max_cores > 0:
                return max_cores
            log.warning(f"Requested {max_cores} cores exceeds the number of cores available on this machine ({num_cores}). Using all available cores.")
            return max_cores
        case _:
            raise ValueError(f"Invalid value for max_cores: {max_cores}")


class UnmatchedSlitIDError(Exception):
    pass


def _find_matching_simul_slit(slit: datamodels.SlitModel,
                              simul_slit_sids: list[int],
                              simul_slit_orders: list[int],
                              ) -> int:
    """
    Parameters
    ----------
    slit : `~jwst.datamodels.SlitModel`
        Source slit model
    simul_slit_sids : list
        List of source IDs for simulated slits
    simul_slit_orders : list
        List of spectral orders for simulated slits

    Returns
    -------
    good_idx : int
        Index of the matching simulated slit in the list of simulated slits
    """

        # Retrieve simulated slit for this source only
    sid = slit.source_id
    order = slit.meta.wcsinfo.spectral_order
    good = (np.array(simul_slit_sids) == sid) * (np.array(simul_slit_orders) == order)
    if not any(good):
        raise UnmatchedSlitIDError(f"Source ID {sid} order {order} requested by input slit model \
                but not found in simulated slits. Setting contamination correction to zero for that slit.")
    return np.where(good)[0][0]


def _cut_frame_to_match_slit(contam: np.ndarray, slit: datamodels.SlitModel) -> np.ndarray:
    
    """Cut out the contamination image to match the extent of the source slit.

    Parameters
    ----------
    contam : 2D array
        Contamination image for the full grism exposure
    slit : `~jwst.datamodels.SlitModel`
        Source slit model

    Returns
    -------
    cutout : 2D array
        Contamination image cutout that matches the extent of the source slit

    """
    x1 = slit.xstart
    y1 = slit.ystart
    cutout = contam[y1:y1 + slit.ysize, x1:x1 + slit.xsize]

    return cutout


class SlitOverlapError(Exception):
    pass

class CommonSlit(Protocol):
    '''
    class protocol for two slits that represent the same source and order, e.g. data and model
    '''
    slit0: datamodels.SlitModel
    slit1: datamodels.SlitModel
    
    def match_backplane(self) -> tuple[datamodels.SlitModel, datamodels.SlitModel]:
        ...


class CommonSlitPreferFirst(CommonSlit):
    '''
    Treat slit0 as the reference slit, and match attributes of slit1 to it
    '''
    def __init__(self, slit0: datamodels.SlitModel, slit1: datamodels.SlitModel):
        self.slit0 = slit0
        self.slit1 = slit1

    def match_backplane(self) -> tuple[datamodels.SlitModel, datamodels.SlitModel]:

        data0 = self.slit0.data
        data1 = self.slit1.data

        x1 = self.slit1.xstart - self.slit0.xstart
        y1 = self.slit1.ystart - self.slit0.ystart
        backplane1 = np.zeros_like(data0)
        
        i0 = max([y1,0])
        i1 = min([y1+data1.shape[0], data0.shape[0], data1.shape[0]])
        j0 = max([x1,0])
        j1 = min([x1+data1.shape[1], data0.shape[1], data1.shape[1]])
        if i0 >= i1 or j0 >= j1:
            raise SlitOverlapError(f"No overlap region between data and model for slit {self.slit0.source_id}, \
                                    order {self.slit0.meta.wcsinfo.spectral_order}. \
                                    Setting contamination correction to zero for that slit.")

        backplane1[i0:i1, j0:j1] = data1[i0:i1, j0:j1]

        self.slit1.data = backplane1
        self.slit1.xstart = self.slit0.xstart
        self.slit1.ystart = self.slit0.ystart
        self.slit1.xsize = self.slit0.xsize
        self.slit1.ysize = self.slit0.ysize
        
        return self.slit0, self.slit1


class CommonSlitEncompass(CommonSlit):
    '''
    Encompass the data from both slits in a common backplane
    '''
    def __init__(self, slit0: datamodels.SlitModel, slit1: datamodels.SlitModel):
        self.slit0 = slit0
        self.slit1 = slit1
    
    def match_backplane(self) -> tuple[datamodels.SlitModel, datamodels.SlitModel]:
        '''
        put data from the two slits into a common backplane
        so outputs have the same dimensions
        and alignment is based on slit.xstart, slit.ystart

        Parameters
        ----------
        slit0 : SlitModel
            First slit model
        slit1 : SlitModel
            Second slit model

        Returns
        -------
        slit0 : SlitModel
            First slit model with data updated to common backplane
        slit1 : SlitModel
            Second slit model with data updated to common backplane
        '''
            
        data0 = self.slit0.data
        data1 = self.slit1.data

        shape = (max(data0.shape[0], data1.shape[0]), max(data0.shape[1], data1.shape[1]))
        xmin = min(self.slit0.xstart, self.slit1.xstart)
        ymin = min(self.slit0.ystart, self.slit1.ystart)
        shape = max(self.slit0.xsize + self.slit0.xstart - xmin, 
                    self.slit1.xsize + self.slit1.xstart - xmin), \
                    max(self.slit0.ysize + self.slit0.ystart - ymin, 
                        self.slit1.ysize + self.slit1.ystart - ymin)
        x0 = self.slit0.xstart - xmin
        y0 = self.slit0.ystart - ymin
        x1 = self.slit1.xstart - xmin
        y1 = self.slit1.ystart - ymin

        backplane0 = np.zeros(shape).T
        backplane0[y0:y0+data0.shape[0], x0:x0+data0.shape[1]] = data0
        backplane1 = np.zeros(shape).T
        backplane1[y1:y1+data1.shape[0], x1:x1+data1.shape[1]] = data1

        self.slit0.data = backplane0
        self.slit1.data = backplane1
        for slit in [self.slit0, self.slit1]:
            slit.xstart = xmin
            slit.ystart = ymin
            slit.xsize = shape[0]
            slit.ysize = shape[1]
        
        return self.slit0, self.slit1


def contam_corr(input_model: datamodels.MultiSlitModel, 
                waverange: datamodels.WavelengthrangeModel, 
                photom: datamodels.NrcWfssPhotomModel | datamodels.NisWfssPhotomModel,
                max_cores: str | int = "none", 
                brightest_n: int = None,
                ) -> tuple[datamodels.MultiSlitModel, datamodels.ImageModel, datamodels.MultiSlitModel, datamodels.MultiSlitModel]:
    """
    The main WFSS contamination correction function

    Parameters
    ----------
    input_model : `~jwst.datamodels.MultiSlitModel`
        Input data model containing 2D spectral cutouts
    waverange : `~jwst.datamodels.WavelengthrangeModel`
        Wavelength range reference file model
    photom : `~jwst.datamodels.NrcWfssPhotomModel` or `~jwst.datamodels.NisWfssPhotomModel`
        Photom (flux cal) reference file model    
    max_cores : string or int
        Number of cores to use for multiprocessing. If set to 'none'
        (the default), then no multiprocessing will be done. The other
        allowable string values are 'quarter', 'half', and 'all', which indicate
        the fraction of cores to use for multi-proc. The total number of
        cores includes the SMT cores (Hyper Threading for Intel).
        If an integer is provided, it will be the exact number of cores used.
    brightest_n : int
        Number of sources to simulate. If None, then all sources in the
        input model will be simulated. Requires loading the source catalog
        file if not None. Note runtime scales non-linearly with this number
        because brightest (and therefore typically largest) sources are
        simulated first.

    Returns
    -------
    output_model : `~jwst.datamodels.MultiSlitModel`
        A copy of the input_model that has been corrected
    simul_model : `~jwst.datamodels.ImageModel`
        Full-frame simulated image of the grism exposure
    contam_model : `~jwst.datamodels.MultiSlitModel`
        Contamination estimate images for each source slit

    """

    num_cores = multiprocessing.cpu_count()
    ncpus = determine_multiprocessing_ncores(max_cores, num_cores)

    # Initialize output model
    output_model = input_model.copy()

    # Get the segmentation map, direct image for this grism exposure
    seg_model = datamodels.open(input_model.meta.segmentation_map)
    direct_file = input_model.meta.direct_image
    image_names = [direct_file]
    log.debug(f"Direct image names={image_names}")

    # Get the grism WCS object and offsets from the first cutout in the input model.
    # This WCS is used to transform from direct image to grism frame for all sources
    # in the segmentation map - the offsets are required so that we can shift
    # each source in the segmentation map to the proper grism image location
    # using this particular wcs, but any cutout's wcs+offsets would work.
    grism_wcs = input_model.slits[0].meta.wcs
    xoffset = input_model.slits[0].xstart - 1
    yoffset = input_model.slits[0].ystart - 1

    # Find out how many spectral orders are defined, based on the
    # array of order values in the Wavelengthrange ref file
    spec_orders = np.asarray(waverange.order)
    spec_orders = spec_orders[spec_orders != 0]  # ignore any order 0 entries
    log.debug(f"Spectral orders defined = {spec_orders}")

    # Get the FILTER and PUPIL wheel positions, for use later
    filter_kwd = input_model.meta.instrument.filter
    pupil_kwd = input_model.meta.instrument.pupil

    # NOTE: The NIRCam WFSS mode uses filters that are in the FILTER wheel
    # with gratings in the PUPIL wheel. NIRISS WFSS mode, however, is just
    # the opposite. It has gratings in the FILTER wheel and filters in the
    # PUPIL wheel. So when processing NIRISS grism exposures the name of
    # filter needs to come from the PUPIL keyword value.
    if input_model.meta.instrument.name == 'NIRISS':
        filter_name = pupil_kwd
    else:
        filter_name = filter_kwd

    # select a subset of the brightest sources using source catalog
    if brightest_n is not None:
        source_catalog = Table.read(input_model.meta.source_catalog, format='ascii.ecsv')
        source_catalog.sort("isophotal_abmag", reverse=False) #magnitudes in ascending order, since brighter is smaller mag number
        selected_IDs = list(source_catalog["label"])[:brightest_n]
    else:
        selected_IDs = None

    obs = Observation(image_names, seg_model, grism_wcs, filter_name,
                      boundaries=[0, 2047, 0, 2047], offsets=[xoffset, yoffset], max_cpu=ncpus,
                      ID=selected_IDs)
    
    good_slits = [slit for slit in output_model.slits if slit.source_id in obs.IDs]
    output_model = datamodels.MultiSlitModel()
    output_model.update(input_model)
    output_model.slits.extend(good_slits)
    log.info(f"Simulating only the brightest {brightest_n} sources")


    simul_all = None
    for order in spec_orders:

        # Load lists of wavelength ranges and flux cal info
        wavelength_range = waverange.get_wfss_wavelength_range(filter_name, [order])
        wmin = wavelength_range[order][0]
        wmax = wavelength_range[order][1]
        log.debug(f"wmin={wmin}, wmax={wmax} for order {order}")
        sens_waves, sens_response = get_photom_data(photom, filter_kwd, pupil_kwd, order)

        # Create simulated grism image for each order and sum them up
        log.info(f"Creating full simulated grism image for order {order}")
        obs.disperse_all(order, wmin, wmax, sens_waves, sens_response)
        if simul_all is None:
            simul_all = obs.simulated_image
        else:
            simul_all += obs.simulated_image

    # Save the full-frame simulated grism image
    simul_model = datamodels.ImageModel(data=simul_all)
    simul_model.update(input_model, only="PRIMARY")

    simul_slit_sids = np.array(obs.simul_slits_sid)
    simul_slit_orders = np.array(obs.simul_slits_order)

    # Loop over all slits/sources to subtract contaminating spectra
    log.info("Creating contamination image for each individual source")
    contam_model = datamodels.MultiSlitModel()
    contam_model.update(input_model)
    slits = []
    for slit in output_model.slits:

        try:
            good_idx = _find_matching_simul_slit(slit, simul_slit_sids, simul_slit_orders)
            this_simul = obs.simul_slits.slits[good_idx]
            slit, this_simul = CommonSlitPreferFirst(slit, this_simul).match_backplane()
            simul_all_cut = _cut_frame_to_match_slit(simul_all, slit)
            contam_cut = simul_all_cut - this_simul.data

        except (UnmatchedSlitIDError, SlitOverlapError) as e:
            log.warning(e)
            contam_cut = np.zeros_like(slit.data)     
        
        contam_slit = copy.copy(slit)
        contam_slit.data = contam_cut
        slits.append(contam_slit)

        # Subtract the contamination from the source slit
        slit.data -= contam_cut

    # Save the contamination estimates for all slits
    contam_model.slits.extend(slits)

    output_model.meta.cal_step.wfss_contam = 'COMPLETE'

    return output_model, simul_model, contam_model, obs.simul_slits
