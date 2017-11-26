# Copyright (C) 2010 Association of Universities for Research in Astronomy(AURA)
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:
#
#     1. Redistributions of source code must retain the above copyright
#       notice, this list of conditions and the following disclaimer.
#
#     2. Redistributions in binary form must reproduce the above
#       copyright notice, this list of conditions and the following
#       disclaimer in the documentation and/or other materials provided
#       with the distribution.
#
#     3. The name of AURA and its representatives may not be used to
#       endorse or promote products derived from this software without
#       specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY AURA ``AS IS'' AND ANY EXPRESS OR IMPLIED
# WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
# MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL AURA BE LIABLE FOR ANY DIRECT, INDIRECT,
# INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
# BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
# OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
# ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
# TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
# USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
# DAMAGE.
"""
Utilities
"""
import inspect
import os
import sys


def import_class(full_name, subclassof=object, config_file=None):
    """
    Import the Python class `full_name` given in full Python package format,
    e.g.::

        package.another_package.class_name

    Return the imported class. Optionally, if `subclassof` is not None
    and is a Python class, make sure that the imported class is a
    subclass of `subclassof`.
    """
    # Understand which class we need to instantiate. The class name is given in
    # full Python package notation, e.g.
    #   package.subPackage.subsubpackage.className
    # in the input parameter `full_name`. This means that
    #   1. We HAVE to be able to say
    #       from package.subPackage.subsubpackage import className
    #   2. If `subclassof` is defined, the newly imported Python class MUST be a
    #      subclass of `subclassof`, which HAS to be a Python class.

    if config_file is not None:
        sys.path.insert(0, os.path.dirname(config_file))

    try:
        full_name = full_name.strip()
        package_name, sep, class_name = full_name.rpartition('.')
        if not package_name:
            raise ImportError("{0} is not a Python class".format(full_name))
        imported = __import__(
            package_name, globals(), locals(), [class_name, ], level=0)

        step_class = getattr(imported, class_name)

        if not isinstance(step_class, type):
            raise TypeError(
                'Object {0} from package {1} is not a class'.format(
                    class_name, package_name))
        elif not issubclass(step_class, subclassof):
            raise TypeError(
                'Class {0} from package {1} is not a subclass of {2}'.format(
                    class_name, package_name, subclassof.__name__))
    finally:
        if config_file is not None:
            del sys.path[0]

    return step_class


def get_spec_file_path(step_class):
    """
    Given a Step (sub)class, divine and return the full path to the
    corresponding spec file. Use the fact that by convention, the spec
    file is in the same directory as the `step_class` source file. It
    has the name of the Step (sub)class and extension .spec.
    """
    try:
        step_source_file = inspect.getfile(step_class)
    except TypeError:
        return None
    step_source_file = os.path.abspath(step_source_file)

    # Since `step_class` could be defined in a file called whatever,
    # we need the source file basedir and the class name.
    dir = os.path.dirname(step_source_file)
    return os.path.join(dir, step_class.__name__ + '.spec')


def find_spec_file(step_class):
    """
    Return the full path of the given Step subclass `step_class`, if
    it exists or None if it does not.
    """
    spec_file = get_spec_file_path(step_class)
    if spec_file is not None and os.path.exists(spec_file):
        return spec_file
    return None


def islist_tuple(obj):
    """
    Return True if `obj` is either a list or a tuple. False otherwise.
    """
    return isinstance(obj, tuple) or isinstance(obj, list)
