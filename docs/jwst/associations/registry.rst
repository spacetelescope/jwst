.. _design-registry:

Association Registry
====================

The :class:`~jwst.associations.registry.AssociationRegistry` is the
rule organizer. An `~jwst.associations.registry.AssociationRegistry` is instantiated with the
files containing the desired rules. The
:meth:`~jwst.associations.registry.AssociationRegistry.match` method
is used to find associations that a member belongs to.

`~jwst.associations.registry.AssociationRegistry` is a subclass of :py:obj:`dict` and supports all of
its methods. In particular, multiple `~jwst.associations.registry.AssociationRegistry` instances can be
combined using the :py:meth:`~dict.update` method.
