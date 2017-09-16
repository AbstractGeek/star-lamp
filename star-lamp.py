#!/usr/bin/env python
import argparse
from astropy.time import Time
from astropy.coordinates import EarthLocation
from astropy.coordinates import SkyCoord
from astropy.coordinates import AltAz
from datetime import datetime
from astropy import units as u
from os.path import abspath
from re import search as re_search
from collections import defaultdict
from pprint import pprint
import numpy as np
import solid
import solid.utils as sutil

# CONSTANTS
SEGMENTS = 50


def sph2cart(radius, azimuth, elevation):
    """Convert spherical coordinates to cartesian coordinates."""
    x = radius * np.cos(elevation * np.pi / 180) * \
        np.cos(azimuth * np.pi / 180)
    y = radius * np.cos(elevation * np.pi / 180) * \
        np.sin(azimuth * np.pi / 180)
    z = radius * np.sin(elevation * np.pi / 180)
    return x, y, z


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


def constellation_stick_figures(filename, obstime, obspos, altitude_cutoff):
    """Obtain local coordinates for constellation stick figures."""
    const_name = []
    const_location = []
    for line in open(filename):
        if (line[0] == "#") or (len(line) < 10):
            continue
        try:
            cname, ra1, dec1, ra2, dec2 = re_search(
                "^\s*(\w+)\s+([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s+"
                "([-]?\d+\.\d+)\s+([-]?\d+\.\d+)\s*$",
                line).groups()

            const_name.append(cname)
            const_location.append((float(ra1), float(dec1),
                                   float(ra2),  float(dec2)))
        except ValueError:
            continue
    # Save extracted csvs
    with open(abspath("./Processed-Data/constellations-raw-data.csv"),
              "w") as csvfile:
        for i, cname in enumerate(const_name):
            csvfile.write("%s,%17.12f,%17.12f,%17.12f,%17.12f\n" %
                          (cname, const_location[i][0], const_location[i][1],
                           const_location[i][2], const_location[i][3]))
    # Transform coordinates to azimuth and altitude
    c_start_icrs = SkyCoord(
        ra=[const[0] for const in const_location] * u.degree,
        dec=[const[1] for const in const_location] * u.degree,
        frame='icrs')
    c_start_azalt = c_start_icrs.transform_to(
        AltAz(obstime=obstime, location=obspos))
    c_stop_icrs = SkyCoord(
        ra=[const[2] for const in const_location] * u.degree,
        dec=[const[3] for const in const_location] * u.degree,
        frame='icrs')
    c_stop_azalt = c_stop_icrs.transform_to(
        AltAz(obstime=obstime, location=obspos))
    # Extract transformed coordinates
    az1 = list(c_start_azalt.az.deg)
    alt1 = list(c_start_azalt.alt.deg)
    az2 = list(c_stop_azalt.az.deg)
    alt2 = list(c_stop_azalt.alt.deg)
    stick_figures_all = defaultdict(list)
    # Save transformed coordinates as a defaultdict
    for i, cname in enumerate(const_name):
        stick_figures_all[cname].append((az1[i], alt1[i], az2[i], alt2[i]))
    # Backup data as a csv
    with open(abspath("./Processed-Data/constellations-all.csv"),
              "w") as csvfile:
        for cname in sorted(stick_figures_all.keys()):
            for a1, e1, a2, e2 in stick_figures_all[cname]:
                csvfile.write("%s,%17.12f,%17.12f,%17.12f,%17.12f\n" %
                              (cname, a1, e1, a2, e2))
    # Obtain only the visible hemisphere
    stick_figures = {}
    for cname in sorted(stick_figures_all.keys()):
        save_flag = True
        for a1, e1, a2, e2 in stick_figures_all[cname]:
            if (e1 >= altitude_cutoff) and (e2 >= altitude_cutoff):
                save_flag = (True and save_flag)
            else:
                save_flag = False
                break
        # Check save flag and append to dictionary
        if save_flag:
            stick_figures[cname] = stick_figures_all[cname]
    # Save filtered coordinates as csv
    with open(abspath("./Processed-Data/constellations-visible.csv"),
              "w") as csvfile:
        for cname in sorted(stick_figures.keys()):
            for a1, e1, a2, e2 in stick_figures[cname]:
                csvfile.write("%s,%17.12f,%17.12f,%17.12f,%17.12f\n" %
                              (cname, a1, e1, a2, e2))
    # Return stick_figures data
    return stick_figures


def make_sticks(a1, e1, a2, e2, r, b):
    x1, y1, z1 = sph2cart(r, a1, e1)
    x2, y2, z2 = sph2cart(r, a2, e2)
    # Obtain slope (2D)
    m = -(x2 - x1) / (y2 - y1)
    # Obtain intercept of the perpendicular line
    c1 = y1 - m * x1
    c2 = y2 - m * x2
    # Obtain first point set
    x11 = x1 - np.sqrt(np.square(b) / (np.square(m) + 1))
    x12 = x1 + np.sqrt(np.square(b) / (np.square(m) + 1))
    y11 = m * x11 + c1
    y12 = m * x12 + c1
    # Obtain second point set
    x21 = x2 - np.sqrt(np.square(b) / (np.square(m) + 1))
    x22 = x2 + np.sqrt(np.square(b) / (np.square(m) + 1))
    y21 = m * x21 + c2
    y22 = m * x22 + c2
    # Create sticks tuple and return
    sticks = ((x11, y11), (x12, y12), (x22, y22), (x21, y21))
    return sticks


def make_stick_figures(shell, stick_figures, radius, breadth):
    """Adds stick figure inscription onto the shell."""
    # Create inscribed shell
    i_shell = solid.difference()
    i_shell.add(shell)
    # Add individual polygons to the inscribed shell
    for cname in sorted(stick_figures.keys()):
        for a1, e1, a2, e2 in stick_figures[cname]:
            i_shell.add(solid.linear_extrude(height=radius)(
                        solid.polygon(
                            make_sticks(a1, e1, a2, e2, radius, breadth))))
    return i_shell


def make_lamp_scad(filename, bright_stars, stick_figures, radius, thickness):
    """Use bright stars and input arguements to create the scad lamp."""
    # Create a main shell to be inscribed on
    shell = solid.difference()(
        solid.sphere(r=radius),
        solid.sphere(r=radius - thickness / 2),
        sutil.down(radius)(
            solid.cube(2 * radius, center=True)
        )
    )
    # Add stick figures
    i_shell = make_stick_figures(shell, stick_figures, radius, 0.5)
    # i_shell = shell
    # Add another layer of shell (without inscription)
    c_shell = solid.union()(
        i_shell,
        solid.difference()(
            solid.sphere(r=radius - thickness / 2),
            solid.sphere(r=radius - thickness),
            sutil.down(radius)(
                solid.cube(2 * radius, center=True)
            )
        )
    )

    # Add stars
    c = solid.difference()
    c.add(c_shell)
    for _, az, alt, _, rad in bright_stars:
        c.add(solid.rotate(a=[0, -1 * (90 - alt), az])(
            solid.cylinder(rad, h=radius)))

    # Render to file
    solid.scad_render_to_file(c, filename,
                              file_header='$fn = %s;' % SEGMENTS)


def main():
    """Obtain command line arguments and create the star lamp."""
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
        "-d", "--thickness-ratio", default=0.1, type=float,
        help="Thickness of the globe/sphere")
    parser.add_argument(
        "-o", "--output-file", default="star-lamp.scad", type=str,
        help="Output filename (default: star-lamp.scad)")

    args = parser.parse_args()
    # print(args)       # To debug
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
    # Obtain constellation stick figures
    stick_figures = constellation_stick_figures(
        abspath('./Raw-Data/ConstellationStickFigures.dat'),
        t, pos, args.altitude_cutoff)

    # Make the scad
    make_lamp_scad(abspath('./' + args.output_file),
                   bright_stars, stick_figures, args.size,
                   args.thickness_ratio * args.size)

    # Done!


if __name__ == '__main__':
    main()
