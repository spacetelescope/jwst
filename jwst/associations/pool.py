"""Association Pool, Assocation Rule, and Data product generator.

For JWST, associations are groups of exposures that are in one way
or another related scientifically or calibration-wise. Examples are:

- Dithers of the same target
- Time-separated visits of the same area
- Filter-separated visits of the same area
- Overlapping observations
- Calibration observations that may be shared among different exposures

The ultimate goal is to be able to say: For this target, deliver a
single product of calibrated data, potentially based on restricting criteria.

To this end, "association pools" will be created of general potential
exposures that may be relatable. The "association generator" will take
as primary input a target and the association pool that target is in,
and the association generator will produce, the list of exposures to
be calibrated and, if possible, combined to create a single product of
the science for the target.

The initial goal will be limited to a common proposal and common RA/DEC.
The main goal of development will be to produce a framework in which
different types of associations can be later built, possibly by
trained users, to define further filters to create associations.

Routine Listings
----------------

Notes
-----

`Project home <https://trac.stsci.edu/trac/DMS/wiki/WebbDMSDataProcessing/Associations>`_

"""
from astropy.io.ascii import convert_numpy

from astropy.table import Table


class AssociationPool(Table):
    """Association Pool

    Parameters
    ----------
    See superclass for initial paramters.

    Notes
    -----
     `Project home <https://trac.stsci.edu/trac/DMS/wiki/WebbDMSDataProcessing/Associations>`_
    """

    @classmethod
    def read(cls, filename, delimiter='|', format='ascii', **kwargs):
        table = Table.read(filename, delimiter=delimiter,
                           format=format,
                           converters=_ConvertToStr(), **kwargs)
        table.meta['pool_file'] = filename
        return table


class _ConvertToStr(dict):
    def __getitem__(self, k):
        return [convert_numpy(str)]

    def get(self, k, default=None):
        return self.__getitem__(k)
