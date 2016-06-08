"""
This module contains functions to create NIRCAM reference files.

NIRCAM model:
im.meta.instrument.name :'NIRCAM'
im.meta.instrument.channel : 'SHORT'
im.meta.instrument.module : 'B'
im.meta.instrument.detector : 'NRCB1'
im.meta.exposure.type : 'NRC_IMAGE'

???
im.meta.instrument.pupil : 'FLAT' (would be GRISMR or GRISMC for slitless)

Transform Paths for Imaging mode:

NIRCAMASW_1 --> NIRCAMASW --> OTESKY --> OTESKY_AM
NIRCAMASW_2 --> NIRCAMASW --> OTESKY --> OTESKY_AM
NIRCAMASW_3 --> NIRCAMASW --> OTESKY --> OTESKY_AM
NIRCAMASW_4 --> NIRCAMASW --> OTESKY --> OTESKY_AM
NIRCAMALW --> NIRCAMALW --> OTESKY --> OTESKY_AM
NIRCAMBSW_1 --> NIRCAMBSW --> OTESKY --> OTESKY_AM
NIRCAMBSW_2 --> NIRCAMBSW --> OTESKY --> OTESKY_AM
NIRCAMBSW_3 --> NIRCAMBSW --> OTESKY --> OTESKY_AM
NIRCAMBSW_4 --> NIRCAMBSW --> OTESKY --> OTESKY_AM
NIRCAMBLW --> NIRCAMBLW --> OTESKY --> OTESKY_AM


"""
from asdf import AsdfFile
from astropy.modeling.models import Mapping

import read_siaf_table

def create_nircam_distortion(channel, module, detector, outname):
    """
    Create an asdf reference file with all distortion components for the MIRI imager.

    NOTE: The IDT has not provided any distortion information. The files are constructed
    using ISIM transformations provided/(computed?) by the TEL team which they use to
    create the SIAF file.
    These reference files should be replaced when/if the IDT provides us with distortion.

    Parameters
    ----------
    channel : str
        header['channel'] or ImageModel.meta.imstrument.channel
        ("SHORT" or "LONG")
    module : str
        header['module'] or ImageModel.meta.imstrument.module
        ("A" or "B")
    detector : str
        NRCB1, NRCB2, NRCB3, NRCB4, NRCB5, NRCA1, NRCA2, NRCA3, NRCA4, NRCA5
        (NRCA5 and NRCB5 are the LW detectors? or simply NRCA and NRCB)
    outname : str
        Name of output file.

    Examples
    --------

    """
    if channel == "SHORT":
        ch = "SW"
    elif channel == "LONG":
        ch = "LW"
    else:
        raise ValueError("Channel {0} not found".format(ch))
    numdet = detector[-1]
    det = detector[-2]
    if numdet == '5':
        numdet = "1"
    from_system = "NIRCAM" + det + ch + "_" + numdet
    to_system = "NIRCAM" + det + ch
    amodel, bmodel, startunit, endunit = read_siaf_table.get_siaf_transform(
        from_system, to_system, 1, 1)
    to_otesky_x, to_otesky_y, startunit, endunit = read_siaf_table.get_siaf_transform(
        to_system, "OTESKY", 5, 5)
    ote2nrcx, ote2nrcy, startunit, endunit = read_siaf_table.get_siaf_transform(
        "OTESKY", to_system, 5, 5)
    nrc2detx,  nrc2dety, startunit, endunit = read_siaf_table.get_siaf_transform(
        to_system, from_system, 1, 1)
    model = Mapping([0, 1, 0, 1]) | amodel & bmodel | Mapping([0, 1, 0, 1]) | \
          to_otesky_x & to_otesky_y
    model_inv = Mapping([0, 1, 0, 1]) | ote2nrcx & ote2nrcy | Mapping((0, 1, 0, 1)) | \
              nrc2detx & nrc2dety
    model.inverse = model_inv

    tree = {"title": "NIRCAM Distortion",
            "instrument": "NIRCAM",
            "pedigree": "GROUND",
            "reftype" : "DISTORTION",
            "author": "N. Dencheva",
            "detector": detector,
            "module": module,
            "channel": channel,
            "exp_type": "NRC_IMAGE",
            "model": model
            }

    fasdf = AsdfFile()
    fasdf.tree = tree
    fasdf.write_to(outname)
