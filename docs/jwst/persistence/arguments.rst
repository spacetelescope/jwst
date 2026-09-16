.. _persistence_step_args:

Step Arguments
==============
The ``persistence`` step has the following optional arguments.

``--save_persistence`` (str, default=None)
  This is a path for saving the persistence array to an ASDF file. If no
  path is given, no persistence array is saved.

``--persistence_time`` (int, default=None)
  This is the length of time in seconds for the flagging window. For example,
  if 200, then once a flagging window has been determined, the flagging window for
  that pixel will end 200 seconds later. The default is None. A positive integer
  is needed for the persistence flagging to be done.

``--dn_threshold`` (float, default=None)
  Any group in the science data above this threshold will
  open a persistence window for flagging, if there is no existing flagging window
  already open.

``--persistence_array_file`` (str, default=None)
  This is a path to an ASDF file with timing data for the end of a timing window for each
  pixel. The time is given in epoch time.

``--persistence_dnu`` (bool, default=False)
  This flag determines if the ``DO_NOT_USE`` flag will get set when the ``PERSISTENCE``
  flag gets set.
