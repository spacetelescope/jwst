Level3 Associations
````````````````````

Level 3 associations associate which level2b files to process to
create higher ordered Level 3 product3. Associations are used as the
primary input to the Level 3 processing tasks. The associations
themselves also define the primary output names to be used by the
level 3 processing tasks when creating the output files.

For Level 3, there is a one-to-one correspondence between associations
and the Level 3 products. However, there are cases where one
association can define multiple Level 3 outputs. Also, Level 3 tasks
are free to create as many auxiliary products as seen fit. Regardless,
any outputs created must use the product names, as defined in the
association, as a template for name creation.

Level 3 associations are created running the `asn_generate` task on an
`association pool` using the default `Level 3 Association Rules`.

The structure of an association is defined by a schema. Level 3
associations follow the Level 3 Association schema. Associations are
currently in JSON format.


Contents:

.. toctree::
   :maxdepth: 2

   level3_asn_technical
   level3_asn_rules
