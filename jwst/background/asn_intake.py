#! /usr/bin/env python
from pathlib import Path
from collections import defaultdict
import logging
import json

from stdatamodels.jwst import datamodels

log = logging.getLogger(__name__)
log.setLevel(logging.DEBUG)


WFSS_TYPES = ["NIS_WFSS", "NRC_GRISM", "NRC_WFSS"]


def asn_get_data(step_input):
    """
    Check if the input is an asn file and get the targets and catalog.

    Parameters
    ----------
    step_input : ImageModel or IFUImageModel or asn file
        Input target data

    Returns
    -------
    step_input : ImageModel or IFUImageModel
        Input target data model
    bkg_list : list
        File name list of background exposures
    """
    try:
        with Path(step_input).open("r") as asn_data:
            asn = json.load(asn_data)
    except UnicodeDecodeError:
        return

    members_by_type = defaultdict(list)

    if len(asn["products"]) > 1:
        log.warning("Multiple products in input association. Output file name will be ignored.")

    # Get the grism image and the catalog, direct image, and segmentation map
    for exp_product in asn["products"]:
        # Find all the member types in the product
        for member in exp_product["members"]:
            members_by_type[member["exptype"].lower()].append(member["expname"])

        # Get the science member. Technically there should only be one. Even if
        # there are more, we'll just get the first one found.
        science_member = members_by_type["science"]
        if len(science_member) != 1:
            log.warning(
                "Wrong number of science exposures found in {}".format(exp_product["name"])  # noqa: E501
            )
            log.warning("    Using only first one.")
        science_member = science_member[0]
        log.info("Working on input %s ...", science_member)

        # Open the datamodel and update it with the relevant info for the background step
        with datamodels.open(science_member) as sci:
            exp_type = sci.meta.exposure.type
            if exp_type in WFSS_TYPES:
                try:
                    sci.meta.source_catalog = Path(members_by_type["sourcecat"][0]).name
                    log.info(f"Using sourcecat file {sci.meta.source_catalog}")
                    sci.meta.segmentation_map = Path(members_by_type["segmap"][0]).name
                    log.info(f"Using segmentation map {sci.meta.segmentation_map}")
                    sci.meta.direct_image = Path(members_by_type["direct_image"][0]).name
                    log.info(f"Using direct image {sci.meta.direct_image}")
                except IndexError:
                    if sci.meta.source_catalog is None:
                        raise IndexError(
                            "No source catalog specified in association or datamodel."
                        ) from None

    return sci, members_by_type
