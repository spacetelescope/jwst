This package provides the tools for computing both the heliocentric 
and barycentric times for JWST datasets.

It requires freetds to be installed with support for kerberos and
also requires pymssql to be installed.

For the jpl ephemeris, it requires that jplephem (v 2.5 or higher)

It currently only works with the DMS database that contains the 
predicted JWST position. These particular predictions only provide
one position per day, so cubic interpolation is used to provide
the necessary accuracy (which has been determined to be about
30 ms which translates to an approximately 9000 km accuracy of
JWST relative to the barycenter or the heliocenter).

The following environmental variables are expected:

JPL_EPHEMERIS # to name the directory containing the ephemeris file
METRICS_SERVER # the database server to use
METRICS_DB # name of the database containing the JWST ephemeris info

Notes on the accuracy of interpolating over a day using a cubic
interpolation. I took every other day samples, and interpolated 
on the half way point to compare to the removed values. 
The resulting maximum distance difference (sqrt of the sum of the
squares of the x, y, z components) was about 7.1 km, well
under the 9000 km requirement (and that is overstating the error
since two day spacing is being used.) If I had to guess the 
intrinsic maximum error is at least half that and perhaps well
less than that.

The computations were checked against two calculators on the web:

Two cases were tested, both on UTC time 2009-02-01T00:00:00
for two different target coordinates:
ra, dec
(0, 0)
(150, -37) (degrees, of course)

These require specifying a (0, 0, 0) position for JWST to get 
the corrections for a geocentric frame of reference.

Barycentric:
http://astroutils.astronomy.ohio-state.edu/time/utc2bjd.html

This site does include conversion to tt time from UTC. Results
agreed to about 13 places (the site uses full jd, not mjd)

Heliocentric:
http://britastro.org/computing/applets_dt.html

This site only computes light time delay so tt time was not used
in the call, but the original mjd. Both cases agreed to the precision
given.

The current implementation puts the onus of updating the 
table of time corrections on the operations staff. It may be
possible to automate this in the script, but that raises 
issues of whether it is desired for the script to access 
the internet and update configuration files automatically.