
Description
===========

``jwst.assign_mtwcs`` is run in the beginning of the level 3 JWST pipeline.
It assigns a WCS to a ``Moving Target`` exposure. The table below defines the
keywords specific to exposures with moving targets:

+--------------+----------------------+----------+------------------------------------------------+
| FITS keyword | Data model attribute | Type     | Description                                    |
|              |                      | (Value)  |                                                |
+==============+======================+==========+================================================+
| TARGTYPE     | meta.target.type     | string   | Type of target                                 |
|              |                      | (moving) |                                                |
+--------------+----------------------+----------+------------------------------------------------+
| MT_RA,       | meta.wcsinfo.mt_ra,  | number,  | Average right ascension and declination        |
| MT_DEC       | meta.wcsinfo.mt_dec  | number   | of the moving target during exposure [deg]     |
+--------------+----------------------+----------+------------------------------------------------+
| MT_AVRA,     | meta.wcsinfo.mt_avra,| number,  | Right ascension and declination of the         |
| MT_AVDEC     | meta.wcsinfo.mt_avdec| number   | moving target averaged between exposures [deg] |
+--------------+----------------------+----------+------------------------------------------------+

The step takes the original science WCS pipeline and adds another step to it such
that the ``output_frame`` of the final WCS is centered at the average location of
the moving target specified by ``(MT_AVRA, MT_AVDEC)``.

The transform of original WCS associated with the science aperture pointing can
be accessed by executing::

  sci_transform = model.meta.wcs.get_transform('detector', 'world')
