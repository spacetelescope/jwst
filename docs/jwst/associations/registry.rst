.. _design-registry:

Association Registry
====================

The :class:`~jwst.associations.registry.AssociationRegistry` is the
rule organizer. An `AssociationRegistry` is instantiated with the
files containing the desired rules. The
:meth:`~jwst.associations.registry.AssociationRegistry.match` method
is used to find associations that a member belongs to.

``AssociationRegistry`` is a subclass of :class:`py3:dict` and supports all of
its methods. In particular, multiple ``AssociationRegistry``'s can be
combined using the :meth:`~py3:dict.update` method.
