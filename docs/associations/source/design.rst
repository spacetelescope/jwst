.. _design:

******************
Association Design
******************

.. _figure-association-generator-overview:

.. figure:: graphics/overview.png
   :scale: 50%

   Association Generator Overview

As introduced in the :ref:`overview`, the figure above shows all the
major players used in generating associations. Since this section will
be describe the code design, the figure below is the overview but
using the class names involved.

.. _figure-class-overivew:

.. figure:: graphics/overview_classes.png
   :scale: 50%

   Association Class Relationship overview

Starting with the most straightforward component, the input table.
The generator uses astropy data tables. Hence anything file that the
table I/O interface can read, the generator can read. 
