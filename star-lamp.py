#!/usr/bin/env python
import argparse
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from datetime import datetime
from astropy import units as u
from os.path import abspath
import csv
import solid
import solid.utils as sutil

# CONSTANTS
SEGMENTS = 50


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


def visible_bright_stars(stars, keys, obstime, obspos,
                         altitude_cutoff, magnitude, radius):
    """Obtain star radius based on magnitude.

    This function extracts azimuth and elevation of visible bright stars
    at the time and the place given.
    """
    # Find bright stars
    new_keys = [k for k in keys if stars[k][2] < magnitude]
    # Convert to altitude and azimuth
    c_icrs = SkyCoord(ra=[stars[k][0] for k in new_keys] * u.degree,
                      dec=[stars[k][1] for k in new_keys] * u.degree,
                      frame='icrs')
    c_azalt = c_icrs.transform_to(AltAz(obstime=obstime, location=obspos))
    azimuth = list(c_azalt.az.deg)
    altitude = list(c_azalt.alt.deg)
    mag = [stars[k][2] for k in new_keys]
    rad = [(magnitude + 1 - m) * radius / magnitude for m in mag]

    # Save csvs
    with open(abspath("./Processed-Data/bright-stars-raw-data.csv"),
              "w") as csvfile:
        for i, az in enumerate(azimuth):
            csvfile.write("%6d,%17.12f,%17.12f,%17.12f,%17.12f\n" %
                          (new_keys[i], az, altitude[i], mag[i], rad[i]))
    # Obtain only the visible hemisphere
    bright_stars = [(new_keys[i], azimuth[i], alt, mag[i], rad[i])
                    for i, alt in enumerate(altitude)
                    if alt > altitude_cutoff]
    # Save csvs
    with open(abspath("./Processed-Data/bright-stars.csv"), "w") as csvfile:
        for star in bright_stars:
            csvfile.write("%6d,%17.12f,%17.12f,%17.12f,%17.12f\n" % star)
    # Return brigh stars
    return bright_stars


def make_lamp_scad(filename, bright_stars, radius, thickness):
    # Make a hollow sphere
    c = solid.difference()(
        solid.sphere(r=radius),
        solid.sphere(r=radius - thickness)
    )
    # Cut of a part of it
    c = solid.difference()(
        c,
        sutil.down(radius)(
            solid.cube(2 * radius, center=True)
        )
    )
    # Add stars
    for _, az, alt, _, rad in bright_stars:
        # print(az_alt)
        c = solid.difference()(
            c,
            solid.rotate(a=[0, 90 - alt, az])(
                solid.cylinder(rad, h=radius))
        )
    # Render to file
    solid.scad_render_to_file(c, 'star-lamp.scad',
                              file_header='$fn = %s;' % SEGMENTS)


def main():
    # parse input
    parser = argparse.ArgumentParser(
        description=(
            'Creates a scad file with a night sky globe with time'
            ' and positional arguments provided.'))
    # Main arguments
    parser.add_argument("-t", "--time",
                        default=datetime.utcnow(),
                        help="Time in observer location")
    parser.add_argument("-l", "--location",
                        default="12.97:77.59:926",
                        help=("Location of the observer location"
                              "(lat: lon or lat: lon: height)."
                              "Default location is Bangalore, India"))
    # Optional arguments
    parser.add_argument(
        "-m", "--magnitude", default=5.0, type=float,
        help="Minimum brightness magnitude")
    parser.add_argument(
        "-s", "--size", default=100.0, type=float,
        help="Size(radius) of the night sky globe")
    parser.add_argument(
        "-r", "--radius-ratio", default=0.02, type=float,
        help="Radius ratio of the size of the globe and the star size")
    parser.add_argument(
        "-c", "--altitude-cutoff", default=0, type=float,
        help="Altitude cutoff (max visibility angle > 0 default)")
    parser.add_argument(
        "-d", "--thickness-ratio", default=0.05, type=float,
        help="Thickness of the globe/sphere")

    args = parser.parse_args()
    print(args)
    obspos = args.location.split(":")
    # Obtain earthlocation and time
    pos = EarthLocation(lat=float(obspos[0]), lon=float(obspos[1]),
                        height=float(0 if len(obspos) == 2 else obspos[2]))
    t = Time(args.time, scale='utc', location=pos)
    # print double check stuff
    # print(pos)
    # print(t)

    # Process yale bright star catalog
    stars, keys = process_ybsc(abspath('./Raw-Data/bsc5.dat'))
    bright_stars = visible_bright_stars(
        stars, keys, t, pos, args.altitude_cutoff, args.magnitude,
        args.size * args.radius_ratio)

    # Make the scad
    make_lamp_scad(abspath('./sky-lamp.scad'),
                   bright_stars, args.size, args.thickness_ratio * args.size)

    # Done!


if __name__ == '__main__':
    main()
