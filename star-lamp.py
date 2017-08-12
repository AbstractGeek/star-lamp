#!/usr/bin/env python
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from datetime import datetime
from astropy import units as u
from os.path import abspath
from math import pow as power
from solid import *
from solid.utils import *  # Not required, but the utils module is useful


def process_ybsc(filename):
    """
    Extract rise ascension and declination.

    Part of the planisphere function from in-the-sky.org (by Dominic Ford)
    """
    stars = {}

    for line in open(filename):
        if len(line) < 100:
            continue
        try:
            hd = int(line[25:31])
            ra_hrs = float(line[75:77])
            ra_min = float(line[77:79])
            ra_sec = float(line[79:82])
            dec_neg = (line[83] == '-')
            dec_deg = float(line[84:86])
            dec_min = float(line[86:88])
            dec_sec = float(line[88:90])
            mag = float(line[102:107])
        except ValueError:
            continue

        RA = (ra_hrs + ra_min / 60 + ra_sec / 3600) / 24 * 360
        DEC = (dec_deg + dec_min / 60 + dec_sec / 3600)
        if dec_neg:
            DEC = -DEC

        stars[hd] = [RA, DEC, mag]

    keys = sorted(stars.keys())

    return stars, keys


def get_bright_stars(stars, keys, magnitude=4.5, radius=5.0):
    """Obtain star radius based on magnitude."""
    new_keys = [k for k in keys if stars[k][2] < magnitude]
    # star_brightness = [power(2.512, -1 * stars[k][2]) for k in new_keys]
    # star_radius = [s_b * radius for s_b in star_brightness]
    bright_stars = [(k, stars[k][0], stars[k][1],
                     (5 - stars[k][2]) / 2) for k in new_keys]

    return bright_stars


# Process yale bright star catalog
stars, keys = process_ybsc(abspath('./Raw-Data/bsc5.dat'))
bright_stars = get_bright_stars(stars, keys)

# import argparse
pos = EarthLocation(lat=12.97, lon=77.59, height=926)
t = Time(datetime.utcnow(), scale='utc', location=pos)
c = SkyCoord(ra=10.625 * u.degree, dec=41.2 * u.degree, frame='icrs')
c_all = SkyCoord(ra=[ra for _, ra, _, _ in bright_stars] * u.degree,
                 dec=[dec for _, _, dec, _ in bright_stars] * u.degree,
                 frame='icrs')
c_azalt = c_all.transform_to(AltAz(obstime=t, location=pos))

# print(c_azalt)

# Convert to scad
SEGMENTS = 50
r = 50
c = difference()(
    sphere(r=r),
    sphere(r=r - 10)
)

for i, az_alt in enumerate(c_azalt):
    # print(az_alt)
    c = difference()(
        c,
        rotate([az_alt.az.deg, az_alt.alt.deg, 0])(
            cylinder(bright_stars[i][3], h=r))
    )

# Cut off bottom Part
c = difference()(
    c,
    down(50)(
        cube(50, 200)
    )
)

scad_render_to_file(c, 'sky-lamp.scad',
                    file_header='$fn = %s;' % SEGMENTS)
