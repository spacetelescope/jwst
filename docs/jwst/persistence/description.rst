Description
===========

:Class: `jwst.persistence.PersistenceStep`
:Alias: persistence

This step computes persistence flagging based on the ``SATURATED`` flag.
When ``--persistence_time`` is a positive integer, each pixel is searched
for the first group flagged as ``SATURATED``. The exposure time of this
group is then computed and a ``PERSISTENCE`` flagging window is computed
to be this start time as the beginning of the window, with the end time
of the window being this start time plus the number of seconds inputted
for ``--persistence_time``. For each group for that pixel that has an
exposure time that happens before the end time of the frame will then
be flagged ``PERSISTENCE``. This flagging window will persist across
integration times.

An optional ``--persistence_array_file`` can be passed to this step. If
so, this file is expected to be an ASDF file with a 2-D array called 
"persistence_data". This array will contain timing data that will allow
persistence flagging to persist across exposures, not just integrations.
The entries into this array will 0.0, indicating no current timing window
for that pixel, or the time of the end of the current persistence window
in epoch time. If there is no file indicated with this parameter, internally
this array is created with all 0.0 entries. If the ``--save_persistence``
is option is selected, this persistence array will be save as an ASDF file
with the file name used for this parameter.

Input
=====
The input science file is a RampModel.

Output
======
The output science file is a RampModel with the ``PERSISTENCE`` flag set
for identified pixels.

If the user specifies ``save_persistence`` with a file path, a third output
file will be written to that path. This is an ASDF file with ``persistence_data``
attribute. This attribute contains the end of the computed persistence window.
If a pixel has a 0.0 entry, there is no persistence window that persisted past
the end of the processing time for the previous exposure processed.
