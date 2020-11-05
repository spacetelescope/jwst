import glob
from astropy.io import fits
from astropy import wcs
from astropy.wcs import InvalidTransformError
from jwst import datamodels
from jwst.datamodels.schema import walk_schema


wcslib_kw_to_remove = ['LONPOLE', 'LATPOLE', 'MJD-OBS', 'DATE-OBS']


def move_wcs(files, remove_asdf=False):
    for name in files:
        print('Working on file {0}.'.format(name))
        f = fits.open(name, mode='update')
        if f[0].header['telescop'] != 'JWST':
            f.close()
            continue

        new_hdr = _collect_wcs_keywords(f)
        update_sci_extension(f, new_hdr)
        clean_primary_header(f, new_hdr)

        # Remove ASDF extension if present:
        if remove_asdf:
            for i, hdu in enumerate(f):
                try:
                    if hdu.header['EXTNAME'] == 'ASDF':
                        break
                except KeyError:
                    continue
            del f[i]
        f.close()


def _build_schema2fits_dict(schema):
    """
    Utility function to create a dict that maps metadata attributes to their
    FITS keyword and FITS HDU locations (if any).

    Parameters
    ----------
    schema : JSON schema fragment
        The schema in which to search.

    Returns
    -------
    results : dict
        Dictionary with schema metadata path as keys and a tuple of FITS
        keyword and FITS HDU as values.

    """
    def build_schema_dict(subschema, path, combiner, ctx, recurse):
        if len(path) and path[0] == 'extra_fits':
            return True
        kw = subschema.get('fits_keyword')
        hdu = subschema.get('fits_hdu')
        if kw is not None:
            results['.'.join(path)] = (kw, hdu)

    results = {}
    walk_schema(schema, build_schema_dict, results)

    return results


def _collect_wcs_keywords(f):
    # Get keywords to go in SCI header from the datamodels schema
    dm = datamodels.open(f, pass_invalid_values=True)
    stsci_wcs_kw = []
    for k, v in _build_schema2fits_dict(dm.meta._schema).items():
        if v[1] == 'SCI':
            if v[0] != 'WCSAXES':
                stsci_wcs_kw.append(v[0])
    dm.close()

    # Convert the WCS to a header
    try:
        w = wcs.WCS(f[0].header)
    except InvalidTransformError:
        if f[0].header['cunit3'] and f[0].header['cunit3'] == 'micron':
            f[0].header['cunit3'] = 'um'
            w = wcs.WCS(f[0].header)
    new_hdr = w.to_header()

    # Remove kw added by wcslib from the new header
    for kw in wcslib_kw_to_remove:
        try:
            del new_hdr[kw]
        except KeyError:
            continue

    # Add STScI specific kw to new_hdr
    for kw in stsci_wcs_kw:
        try:
            new_hdr[kw] = f[0].header[kw]
        except KeyError:
            pass
    new_hdr = add_default_keywords(new_hdr)
    return new_hdr


def clean_primary_header(f, new_hdr):
    # f is a fits file opened in 'update' mode.
    for kw in new_hdr:
        try:
            del f[0].header[kw]
        except KeyError:
            pass


def update_sci_extension(f, new_hdr):
    # Update the SCI extension header
    f[('SCI', 1)].header.update(new_hdr)
    for kw in new_hdr:
        try:
            f[('SCI', 1)].header.comments[kw] = f[0].header.comments[kw]
        except KeyError:
            pass


def add_default_keywords(new_hdr):
    """
    Add keywords with default values which wcslib does not write out.

    PC, CUNIT, CTYPE
    """
    wcsaxes = new_hdr['WCSAXES']
    if wcsaxes == 3:
        default_pc = {"PC1_1": 1, "PC1_2": 0, "PC1_3": 0,
                      "PC2_1": 0, "PC2_2": 1, "PC2_3": 0,
                      "PC3_1": 0, "PC3_2": 0, "PC3_3": 1
                     }
        default_cunit = {"CUNIT1": "deg", "CUNIT2": "deg", "CUNIT3": "um"}
        default_ctype = {"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN", "CTYPE3": "WAVE"}
    elif wcsaxes == 2:
        default_pc = {"PC1_1": 1, "PC1_2": 0,
                      "PC2_1": 0, "PC2_2": 1,
                     }
        default_cunit = {"CUNIT1": "deg", "CUNIT2": "deg"}
        default_ctype = {"CTYPE1": "RA---TAN", "CTYPE2": "DEC--TAN"}

    if "PC1_1" not in new_hdr:
        new_hdr.update(default_pc)
    if "CUNIT1" not in new_hdr:
        new_hdr.update(default_cunit)
    if "CTYPE1" not in new_hdr:
        new_hdr.update(default_ctype)

    return new_hdr


if __name__ == '__main__':
    import argparse
    import os
    parser = argparse.ArgumentParser(description="Move WCS keywords "
                                     "from Primary to SCI extension.")
    parser.add_argument('files', help="A list of FITS filenames or a path "
                        "to a directory with FITS files.")
    res = parser.parse_args()
    files = res.files

    if isinstance(files, str):
        files = os.path.abspath(files)
        if os.path.isdir(files):
            files = glob.glob("files/*.fits")
        else:
            files = glob.glob(files)

    move_wcs(files)
