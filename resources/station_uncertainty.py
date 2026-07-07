"""
station_uncertainty.py
The goal of this file is to provide helper functions for determining 
a given stations's positional uncertainty.

WHY THREE NUMBERS PER STATION, NOT ONE
---------------------------------------
(Table 7, https://deepspace.jpl.nasa.gov/dsndocs/810-005/301/301M.pdf)
The CSV gives each station three separate error numbers instead of one:
spin radius, longitude, and z. These come from the standard way DSN
station positions are surveyed and published -- as a point's location
relative to Earth's rotation axis, using ordinary cylindrical
coordinates:

    spin radius : straight-line distance from the point to Earth's spin
                  (rotation) axis -- NOT distance from Earth's center.
    longitude   : angle around that axis (same idea as ordinary
                  longitude).
    z           : how far along the spin axis the point sits, measured
                  from the equatorial plane (north is positive).

Each of the three has its own, separately-estimated error, so they are 
treated as three independent 1-sigma numbers rather than one blended value.
To turn "3 independent errors in a station-specific direction" into "one
3D nudge to the station's XYZ position," this script builds the three
unit-length direction vectors (radial-out, along-longitude, along-axis)
that are natural to THAT station's location, then does:

    offset = (random draw x spin-radius sigma) x radial_direction
           + (random draw x longitude sigma)   x longitude_direction
           + (random draw x z sigma)           x axis_direction

That gives one 3D position error vector per Monte Carlo draw.

REQUIRED SPICE KERNELS (must be furnsh'd before called)
------------------------------------------------------------------------
    naif0012.tls                      (leap seconds)
    de440s.bsp                        (planet/Moon ephemeris)
    pck00011.tpc                      (generic body orientation)
    earth_200101_990827_predict.bpc + its frame kernel (Earth orientation,
                                       needed for the ITRF93 <-> J2000
                                       rotation used below)
    earthstns_itrf93_*.bsp            (DSN station positions)
    earth_topo_*.tf                   (DSS station name/ID assignments)
    dss_17_prelim_itrf93_190814.bsp   (DSS-17 positions)

Station names: try "DSS-14" style first. If SPICE doesn't recognize it,
use the integer ID instead: DSN station IDs are 399000 + the 2-digit
station number, e.g. DSS-14 -> 399014.
"""

import csv

import numpy as np
import spiceypy as spice


# ---------------------------------------------------------------------------
# Read the station-uncertainty CSV.
#
# Expected columns: Name, Spin Radius, Longitude, z (all in meters).
# ---------------------------------------------------------------------------

def read_station_sigma_csv(csv_path):
    """
    Reads the DSN station-uncertainty CSV and returns a plain dict:

        { 14: {"spin_radius_m": 0.024, "longitude_m": 0.035, "z_m": 0.030},
          34: {"spin_radius_m": 0.030, "longitude_m": 0.030, "z_m": 0.030},
          ... }

    keyed by the bare DSS station number (an int), not the raw name
    string -- that way "DSS 14", "DSS-14", and "DSS  14" (extra space)
    all end up under the same key, 14.
    """
    station_sigmas = {}

    with open(csv_path, "r", newline="") as f:
        reader = csv.reader(f)
        header = next(reader)  # skip the "Name, Spin Radius, Longitude, z" row

        for row in reader:
            if not row or not row[0].strip():
                continue  # skip blank lines

            if len(row) < 4:
                print(f"Skipping incomplete row in CSV: {row}")
                continue

            name, spin_radius_m, longitude_m, z_m = row[0], row[1], row[2], row[3]
            dss_number = dss_number_from_name(name)

            station_sigmas[dss_number] = {
                "spin_radius_m": float(spin_radius_m),
                "longitude_m": float(longitude_m),
                "z_m": float(z_m),
            }

    return station_sigmas


def dss_number_from_name(station):
    """
    Pulls the bare station number out of a name like "DSS 14", "DSS-14",
    "DSS  14" (extra space), or a NAIF ID like 399014 -- and returns it
    as a plain int (14).
    """
    text = str(station).strip().upper()

    if "DSS" in text:
        digits = "".join(ch for ch in text if ch.isdigit())
        return int(digits)

    # otherwise, assume it's a NAIF ID: 399000 + 2-digit station number,
    # e.g. 399014 -> 14. The station number is just the last 3 digits.
    digits = "".join(ch for ch in text if ch.isdigit())
    return int(digits[-3:])


# ---------------------------------------------------------------------------
# Figure out each station's own local "cylindrical" directions.
#
# Three per-axis sigmas (spin radius, longitude, z) into one 3D offset.
# ---------------------------------------------------------------------------

def station_direction_vectors(station, et):
    """
    Returns three unit vectors (each a length-3 numpy array), in the
    Earth body-fixed frame ITRF93, for THIS station's local:
        radial_direction     -- straight out from Earth's spin axis
        longitude_direction  -- sideways, the direction of increasing
                                 longitude
        axis_direction       -- straight along Earth's spin axis

    These are found from the station's own nominal position (via SPICE),
    so nothing about any specific station's coordinates is hardcoded.
    """
    station_position_km, _ = spice.spkpos(str(station), et, "ITRF93", "NONE", "EARTH")
    x, y, z = station_position_km
    longitude_rad = np.arctan2(y, x)

    radial_direction = np.array([np.cos(longitude_rad), np.sin(longitude_rad), 0.0])
    longitude_direction = np.array([-np.sin(longitude_rad), np.cos(longitude_rad), 0.0])
    axis_direction = np.array([0.0, 0.0, 1.0])

    return radial_direction, longitude_direction, axis_direction


def random_station_offset_km(station, et, sigmas, rng):
    """
    One random draw of "how far off is this station's true position from
    where SPICE says it is," as a 3D vector in J2000 km.

    sigmas : this station's {"spin_radius_m", "longitude_m", "z_m"} dict,
             e.g. station_sigmas[14] from read_station_sigma_csv().
    rng    : a numpy random Generator, e.g. np.random.default_rng().
    """
    radial_direction, longitude_direction, axis_direction = station_direction_vectors(station, et)

    offset_earth_fixed_km = (
        rng.normal(0.0, sigmas["spin_radius_m"]) * radial_direction
        + rng.normal(0.0, sigmas["longitude_m"]) * longitude_direction
        + rng.normal(0.0, sigmas["z_m"]) * axis_direction
    ) / 1000.0  # meters -> km

    # The offset above is built in ITRF93 (Earth body-fixed, rotating).
    # Everything else in this script works in J2000 (inertial), so rotate
    # it into that frame before using it.
    itrf93_to_j2000 = spice.pxform("ITRF93", "J2000", et)
    return itrf93_to_j2000 @ offset_earth_fixed_km

