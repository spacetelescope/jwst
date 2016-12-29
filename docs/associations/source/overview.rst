.. _overview:

.. _pipeline modules: http://ssb.stsci.edu/doc/jwst_git/docs/stpipe/html/

Overview
********

What are Associations?
======================

Associations are basically just lists of things, mostly exposures,
that are somehow related. With respect to JWST and the Data Management
System (DMS), associations have the following characteristics:

- Relationships between multiple exposures are captured in an association.
- An association is a means of identifying a set of exposures that belong together and may be dependent upon one another.
- The association concept permits exposures to be calibrated, archived, retrieved, and reprocessed as a set rather than as individual objects.
-  For each association, DMS will generate the most combined and least combined data products. 

Associations are used as the primary input to various Level2 and
Level3 JWST `pipeline modules`_.

In DMS, associations are created by the :ref:`association generator
<asn-generate>`. The association generator is basically a classifier.
The generator takes, as input, a table and one or more :ref:`rules`.
Based on the rules, the generator takes each row of the table and
classifies it, placing that row into one or more associations. These
relationships are show in the figure below.

.. figure:: graphics/overview.png
   :scale: 50%

   Association Generator Overview

The actual structure of the resulting associations is completely
determined by the rules themselves, though it is expected that the
rows from the originating table that created the association are
actually in the association. The generator was initially developed to
support the JWST pipeline, which requires two types of rules. These
predefined rules, called :ref:`Level2 <level2-associations>` and
:ref:`Level3 <level3-associations>` associations, produce very
different association structures.

As the structures of the associations are determined by the rules,
so is the output. For the JWST predefined rules, the output are a series
of `JSON` or `YAML` formatted files for each association created.
The naming of these files are determined by strict rules.
     
Usage
=====

Users should not need to run the generator. Instead, it is expected
that one edits an already existing association that accompanies the
user's JWST data. Or, if need be, an association can be created based
on the :ref:`Level3 example <asn-level3-example>`.

Once an association is in-hand, one can pass it as input to a pipeline
routine. For example::
  
  % strun calwebb_image3.cfg  jw12345_xxxx_asn.json

Generator Usage
---------------

Basic use of the association generator is done through two methods.
From the command-line, the generator is invoked using the command
:ref:`asn_generate <asn-generate>`. From Python, the generator\'s
:ref:`Main` is instantiated.
