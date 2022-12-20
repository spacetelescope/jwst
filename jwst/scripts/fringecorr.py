#!/usr/bin/env python

# ExampleStep

# Import the custom step
from jwst.flatfield import flat_field_step

# Import stpipe.cmdline
from stpipe import cmdline


def main():
    # Pass the step class to cmdline.step_script
    cmdline.step_script(flat_field_step)


if __name__ == '__main__':
    main()
