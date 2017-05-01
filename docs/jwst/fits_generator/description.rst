
Description
===========

Overview
--------

The FITS generator is used to convert data from several different
types of ground test data to DMS Level1b format data.  This format is
described in the document ``DMS Level 1 and 2 Data Product Design -
JWST-STScI-002111`` by Daryl Swade.  The code uses a collection of
templates that govern the population of Level 1b header keyword values
from the data in the input file headers, with different templates for
different file types.  The FITS generator will transform the input data
(in detector coordinates) to the DMS coordinate system, where all of the
imaging data has the same parity as the sky and very similar orientations.

Input details
-------------

To run the FITS generator, a 'proposal' file is required.  There
should be only one proposal file per directory, and it should have a
name like

  ddddd.prop

where d stands for a decimal digit.  This file gives the names of each
input FITS datafile, whether a subarray needs to be extracted from it
and the exposure type (EXP_TYPE), as well as the relationship between
the files from an operational viewpoint (i.e. Observation, Visit,
ParallelSequenceID, Activity, Exposure, Detector).  The file has a
structure similar to XML with nested groups:

::

    <Proposal title="MIRI FM IMG_OPT_01_FOV">
      <Observation>
        <Visit>
          <VisitGroup>
            <ParallelSequenceID>
              <Activity>
                <Exposure>
                  <Detector>
                    <base>MIRFM1T00012942_1_493_SE_2011-07-13T10h45m00.fits</base>
                    <subarray></subarray>
                    <exp_type>MIR_IMAGE</exp_type>
                  </Detector>
                </Exposure>
              </Activity>
            </ParallelSequenceID>
          </VisitGroup>
        </Visit>
      </Observation>
    </Proposal>

Each nest can be repeated as needed.  The <Detector></Detector> tags
contain the information for each input/output file, with the input
file name inside the `<base></base>` tags, the name of the subarray to
be extracted within the `<subarray></subarray>` tag, and the exposure
type within the `<exp_type></exp_type>` tag.

The files within the `<base></base>` tag should be in the same directory
as the proposal file.

The input FITS files can be from any of several different sources:

1. MIRI VM2 testing
2. MIRI FM testing
3. NIRSPEC FM testing
4. NIRSPEC IPS Simulator
5. NIRCAM NCONT testing (detector only)
6. NIRCAM FM testing
7. NIRISS CV testing
8. FGS CV testing

Most data that has been taken using the FITSWriter tool can be
successfully converted to Level 1b format.
