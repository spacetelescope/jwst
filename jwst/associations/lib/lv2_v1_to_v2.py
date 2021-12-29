"""Convert old-style Level2 associations to new style"""
import argparse
from glob import glob
import json
import jsonschema
import logging
import os.path as path
import re

from ...lib.suffix import remove_suffix

logger = logging.getLogger(__name__)
handler = logging.StreamHandler()
logger.addHandler(handler)
handler.setLevel(logging.DEBUG)


def lv2_v1_to_v2(asn_v1):
    """Convert old-style Level2 associations to new style

    Parameters
    ----------
    asn_v1 : dict
        An 'old-style' association

    Returns
    -------
    asn_v2 : dict
        A 'new-style' association
    """

    asn_v2 = dict(asn_v1)
    products = []
    members_to_add = []
    for member in asn_v1['members']:
        product = {
            'name': product_name(member['expname'])
        }
        new_members = []
        new_members.append({
            'expname': member['expname'],
            'exptype': member['exptype']
        })
        for bkg in member.get('bkgexps', []):
            new_members.append({
                'expname': bkg['expname'],
                'exptype': 'BACKGROUND'
            })
            members_to_add.append({
                'expname': bkg['expname'],
                'exptype': 'SCIENCE'
            })
        product['members'] = new_members
        products.append(product)

    for new_member in members_to_add:
        product = {
            'name': product_name(new_member['expname'])
        }
        product['members'] = [{
            'expname': new_member['expname'],
            'exptype': 'SCIENCE'
        }]
        products.append(product)

    asn_v2['products'] = products
    del asn_v2['members']
    return asn_v2


def product_name(expname):
    name = path.splitext(
        path.basename(expname.lower())
    )[0]
    name_nosuffix = remove_suffix(name)
    return name_nosuffix


# ################
# Command line API
# ################
if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Convert old style Level2 associations to new style'
    )

    parser.add_argument(
        'old_asns',
        help='Old style Level2 Association files to convert. May be a glob specification',
        nargs='+'
    )
    parser.add_argument(
        '-p', '--prefix', default='v2_',
        help='Prefix to use for the new files. Default: %(default)s'
    )
    parser.add_argument(
        '-l', '--log-level', default='info',
        help='Logging level'
    )

    args = parser.parse_args()

    numeric_log_level = getattr(logging, args.log_level.upper())
    logger.setLevel(numeric_log_level)

    for fname in args.old_asns:
        logger.info('Working {}'.format(fname))
        with open(fname) as fp:
            asn_v1 = json.load(fp)
        asn_v2 = lv2_v1_to_v2(asn_v1)
        asn_v2_fname = args.prefix + fname
        logger.info('\tWriting to {}'.format(asn_v2_fname))
        with open(asn_v2_fname, 'w') as fp:
            json.dump(asn_v2, fp, indent=4, separators=(',', ': '))
