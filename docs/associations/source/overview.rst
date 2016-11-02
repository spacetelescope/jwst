.. _overview:

########
Overview
########

*************
What is this?
*************

The association generator is basically a classifier. The generator
takes, as input, a table and one or more :ref:`rulesets <ruleset>`. Based on the rules found
in the rulesets, the generator takes each row of the table and
classifies it, placing that row into one or more
:ref:`associations <association>`.

The actual structure of the resulting associations is completely
determined by the rulesets themselves. The generator was initially
developed to support the JWST pipeline, which requires two types of
rulesets. These predefined rulesets, called :ref:`Level2 <level2-associations>` and
:ref:`Level3 <level3-associations>` associations,
produce very different association structures.

As the structures of the associations are determined by the rulesets,
so is the output. For the predefined rulesets, the output are a series
of `JSON` or `YAML` formatted files for each association created.
The naming of these files are determined by strict rules.
     
*****
Usage
*****

Basic use of the association generator is done through two methods.
From the command-line, the generator is invoked using the command
:ref:`asn_generate <asn_generate>`. From Python, the generator\'s
:ref:`Main` class is instantiated.

********
Examples
********
