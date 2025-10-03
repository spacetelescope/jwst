.. _design-registry:

Association Registry
====================

The :class:`~jwst.associations.AssociationRegistry` is the
rule organizer. An `~jwst.associations.AssociationRegistry` is instantiated with the
files containing the desired rules. The
:meth:`~jwst.associations.AssociationRegistry.match` method
is used to find associations that a member belongs to.

`~jwst.associations.AssociationRegistry` is a subclass of :py:obj:`dict` and supports all of
its methods. In particular, multiple `~jwst.associations.AssociationRegistry` instances can be
combined using the :py:meth:`~dict.update` method.
