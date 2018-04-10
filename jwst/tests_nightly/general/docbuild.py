"""
#
# DISABLED
#

from __future__ import print_function

import pandokia.helpers.pycode as pycode

# The nightly documentation build leaves behind a file specifically for
# us to use for this test.  That is, this test evaluates the performance
# of the most recently nightly documentation build.  It does not perform
# a build of its own.
#
# We report a test result for each html document that is built, and for
# each pdf document that is built.  We build both out of the same source.
#

# The file is in csv format, but it is easier to just read the line
# and split it on commas.
# column 1: directory where we found the document
# column 2: html status
# column 3: pdf status
#   the statuses are '' for pass or 'error' for fail.  This test assumes
#   any non-blank string to be fail.

#f=open("/eng/ssb/websites/ssbpublic/doc/jwst_dev_old/stat_summary.csv","r")
f=open("/eng/ssb/websites/ssbpublic/doc/jwst_git/stat_summary.csv","r")

# Here is the first part of the URL where the report of the
# documentation build lives.
urlbase  = 'http://ssb.stsci.edu/doc/jwst_dev_old'

for x in f :
    name, html, pdf = tuple(x.split(','))
    name = name.strip()
    html = html.strip()
    pdf  = pdf.strip()

    # dotname is the directory name with / changed to '.'  The log
    # file name is based on this.
    dotname = name.replace('/','.')

    with pycode.test( 'pdf/' + name ) :
        print(urlbase + '/' + dotname + '/pdf.txt')
        if pdf != '' :
            print("pdf build reports error")
            assert False, 'The document build failed last night'

    with pycode.test( 'html/' + name ) :
        print(urlbase + '/' + dotname + '/html.txt')
        if html != '' :
            print("html build reports error")
            assert False, 'The document build failed last night'
"""
