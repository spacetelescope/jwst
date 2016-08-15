import sys
import argparse

from ..datamodels import MultiSlitModel
from .exp_to_source import exp_to_source


class Main(object):
    """Convert exposure-based slits data to source-based slit data
    Docs from the source.
    """

    def __init__(self, args=None):

        if args is None:
            args = sys.argv[1:]
        if isinstance(args, str):
            args = args.split(' ')

        parser = argparse.ArgumentParser(
            description='Convert exposure-based data to source-based data',
            usage='python -m jwst_pipeline.exp_to_source.exp_to_source files'
        )
        parser.add_argument(
            'files',
            type=str,
            nargs='+',
            help='Files to convert')

        parsed = parser.parse_args(args=args)

        models = [MultiSlitModel(f) for f in parsed.files]
        results = exp_to_source(models)
        for source in results:
            results[source].save('.'.join([source, 'fits']))

if __name__ == '__main__':
    Main()
