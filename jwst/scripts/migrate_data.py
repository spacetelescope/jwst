#!/usr/bin/env python

"""
Migrate .fits files whose format has changed between jwst package versions.
"""
import argparse
from datetime import datetime
import os
import re
import traceback
import warnings

import asdf
from astropy.io import fits
from astropy.time import Time
import numpy as np
from packaging.specifiers import SpecifierSet

from stdatamodels.jwst import datamodels

import jwst


def parse_args():
    parser = argparse.ArgumentParser('migrate_data', 'migrate .fits files whose format has changed between jwst package versions')

    parser.add_argument('files', nargs='+', help='one or more .fits files')

    output_group = parser.add_mutually_exclusive_group(required=True)
    output_group.add_argument('--output-dir', help='write modified files to an output directory')
    output_group.add_argument('--in-place', help='modify files in-place', action='store_true')

    return parser.parse_args()


# If there get to be many of these we may want to move
# them to jwst.datamodels somewhere:

def migrate_mt_table_1_2_2(hdul):
    """moving_target.schema has been filled out with actual data.
    """
    schema = asdf.schema.load_schema('http://stsci.edu/schemas/jwst_datamodel/moving_target.schema')
    dtype = asdf.tags.core.ndarray.asdf_datatype_to_numpy_dtype(schema['properties']['moving_target']['datatype'])
    renamed_columns = {
        'moving_target_Dec': 'mt_apparent_Dec',
        'moving_target_RA': 'mt_apparent_RA',
        'moving_target_x': 'mt_sci_x',
        'moving_target_y': 'mt_sci_y',
        'mt_x_helio': 'mt_apparent_x_helio',
        'mt_y_helio': 'mt_apparent_y_helio',
        'mt_z_helio': 'mt_apparent_z_helio',
        'mt_x_jwst': 'mt_apparent_x_jwst',
        'mt_y_jwst': 'mt_apparent_y_jwst',
        'mt_z_jwst': 'mt_apparent_z_jwst',
        'mt_jwst_distance': 'mt_apparent_jwst_distance',
        'mt_sun_distance': 'mt_apparent_sun_distance',
        'phase_angle': 'apparent_phase_angle',
    }

    for hdu in hdul:
        if hdu.name == 'MOVING_TARGET_POSITION':
            new_data = np.zeros(hdu.data.shape, dtype=dtype)
            for column_name in hdu.data.dtype.names:
                new_data[renamed_columns.get(column_name, column_name)] = hdu.data[column_name]

            # Convert from MJD to ISO
            time_data = Time(hdu.data['time'], format='mjd', scale='utc')
            new_data['time'] = [t.isot for t in time_data]
            hdu.data = new_data


def migrate_mt_table_1_4_0(hdul):
    """moving_target.schema has been filled out with actual data.
    """
    schema = asdf.schema.load_schema('http://stsci.edu/schemas/jwst_datamodel/moving_target.schema')
    dtype = asdf.tags.core.ndarray.asdf_datatype_to_numpy_dtype(schema['properties']['moving_target']['datatype'])
    renamed_columns = {
        'mt_detector_x': 'mt_sci_x',
        'mt_detector_y': 'mt_sci_y',
    }

    for hdu in hdul:
        if hdu.name == 'MOVING_TARGET_POSITION':
            new_data = np.zeros(hdu.data.shape, dtype=dtype)
            for column_name in hdu.data.dtype.names:
                new_data[renamed_columns.get(column_name, column_name)] = hdu.data[column_name]

            hdu.data = new_data


def migrate_spec_table_1_1_0(hdul):
    """
    spectable.schema added additional columns and renamed
    two columns.
    """
    schema = asdf.schema.load_schema('http://stsci.edu/schemas/jwst_datamodel/spectable.schema')
    dtype = asdf.tags.core.ndarray.asdf_datatype_to_numpy_dtype(schema['datatype'])
    renamed_columns = {
        'ERROR': 'FLUX_ERROR',
        'BERROR': 'BKGD_ERROR',
    }

    for hdu in hdul:
        if hdu.name == 'EXTRACT1D':
            new_data = np.zeros(hdu.data.shape, dtype=dtype)
            for column_name in hdu.data.dtype.names:
                new_data[renamed_columns.get(column_name, column_name)] = hdu.data[column_name]
            hdu.data = new_data


# The first key is a model class name, the second
# a jwst package version specifier.  The value
# is a method that accepts an HDUList and modifies
# it in-place.
_MIGRATE_METHODS = {
    'Level1bModel': {
        '> 1.2.1, <= 1.3.3': migrate_mt_table_1_4_0,
    },
    'MultiSpecModel': {
        '> 0.13.1, <= 1.1.0': migrate_spec_table_1_1_0,
    },
    'SpecModel': {
        '> 0.13.1, <= 1.1.0': migrate_spec_table_1_1_0,
    },
}


def migrate_file(filename, args):
    if args.in_place:
        mode = 'update'
    else:
        mode = 'readonly'

    with fits.open(filename, memmap=False, mode=mode) as hdul:
        model_type = hdul[0].header.get('DATAMODL')
        jwst_version = hdul[0].header.get('CAL_VER')

        if not (model_type and jwst_version):
            print(f'Unable to migrate {filename}: DATAMODL and CAL_VER keywords are required')
            return

        match = re.match(r'^[0-9]+\.[0-9]+\.[0-9]+', jwst_version)
        if match is None:
            print(f'Unable to migrate {filename}: CAL_VER not understood')
            return
        jwst_version = match.group(0)

        if model_type not in _MIGRATE_METHODS:
            print(f'Migration for {filename} DATAMODL {model_type} not implemented')
            return

        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            exception_raised = False
            try:
                getattr(datamodels, model_type)(hdul, strict_validation=True)
            except Exception:
                exception_raised = True
            if not exception_raised:
                print(f'{filename} is already valid')
                return

        migrate_method = next((m for s, m in _MIGRATE_METHODS[model_type].items() if jwst_version in SpecifierSet(s)), None)
        if migrate_method is None:
            print(f'Migration for {filename} CAL_VER {jwst_version} not implemented')
            return

        migrate_method(hdul)
        hdul[0].header['HISTORY'] = f'Migrated with jwst {jwst.__version__} migrate_data script {datetime.utcnow().isoformat()}'

        try:
            getattr(datamodels, model_type)(hdul, strict_validation=True)
        except Exception:
            print(f'Migration for {filename} failed to produce a valid model:\n')
            traceback.print_exc()
            return

        if args.in_place:
            hdul.flush()
        else:
            output_filename = os.path.join(args.output_dir, os.path.basename(filename))
            hdul.writeto(output_filename, checksum=True, overwrite=True)


def main():
    args = parse_args()

    if args.output_dir:
        os.makedirs(args.output_dir, exist_ok=True)

    for file in args.files:
        try:
            migrate_file(file, args)
        except Exception:
            print(f'Error migrating {file}:\n')
            traceback.print_exc()


if __name__ == '__main__':
    main()
