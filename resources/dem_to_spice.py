"""
dem_to_spice.py

Bridges LDEM pipeline to SPICE: given a planetocentric lat/lon, pull the
interpolated elevation straight out of an LDEM raster, then package that as
a SPICE-ready Cartesian vector in a body-fixed lunar frame -- optionally
rotated into an inertial frame (e.g. J2000) at a given epoch, which is
needed for the DSN light-time work.

Frame choice
------------
Default body-fixed frame is 'MOON_ME' (mean-Earth/polar-axis), NOT 'MOON_PA'
or 'IAU_MOON'. Reasoning:

  - LOLA/LRO cartographic products define lat/lon in the mean-Earth/polar-axis 
    (ME) frame. That's the cartographic convention, not the dynamical principal-axis
    (PA) frame. If the Cartesian vector is built in MOON_PA using ME-based
    lat/lon, this introduces frame-definition offset (arcsecond-to-sub-100 m
    class error depending on location).

  - 'MOON_ME' here refers to the *generic* SPICE alias, which resolves to
    MOON_ME_DE440_ME421, i.e. matched to de440s.bsp ephemeris. If the generic 
    pck00011.tpc kernel is the only one loaded, 'MOON_ME'/'MOON_PA' aren't defined, 
    causing a fall back to 'IAU_MOON' (lower-fidelity, decoupled
    from DE440) -- pass body_frame="IAU_MOON" explicitly if that's what you
    want.

You can override body_frame on every call below if this default doesn't fit
your setup.

Required kernels (furnsh before calling anything that needs et or rotation):
    naif0012.tls
    de440s.bsp
    pck00011.tpc
    moon_pa_de440_200625.bpc
    moon_de440_220930.tf
"""

from __future__ import annotations

import numpy as np
import rasterio
from pyproj import CRS, Transformer
import spiceypy as spice

# Mean lunar radius -- matches your existing pipeline's `offset` convention
# and the spherical datum LOLA polar-stereographic products are built on.
R_MOON_M = 1737400.0


class LunarDEM:
    """
    Thin wrapper around a single LDEM raster (e.g. LDEM_80S_20M.tif) that
    gives bilinearly-interpolated elevation lookups by lat/lon instead
    of row/col. The lat/lon -> raster-CRS -> pixel chain is handled,
    using whatever CRS is actually baked into the file (so this works
    for north- or south-polar tiles, or anything else, without hardcoding
    a projection like the existing scripts do).

    Mirrors the scaling_factor/offset convention defined in the .LBL file:
        stored_value = raw_pixel * scaling_factor + offset

    Set offset_is_radius=True (default, matches your existing calls with
    offset=1737400.0) if `offset` already represents the Moon's mean radius,
    i.e. stored_value is the absolute radius from the Moon's center. Set it
    to False if `offset` is 0 (or some local datum) and stored_value is a
    height above/below the mean radius instead.
    """

    def __init__(self, path, scaling_factor=0.5, offset=R_MOON_M,
                 offset_is_radius=True, band=1):
        self._ds = rasterio.open(path)
        self._transform = self._ds.transform
        self._inv_transform = ~self._transform
        self.offset_is_radius = offset_is_radius

        arr = self._ds.read(band).astype(np.float32)
        self._array = arr * scaling_factor + offset
        self._nrows, self._ncols = self._array.shape

        # Spherical lunar lat/lon (what your lat/lon columns already are) ->
        # whatever CRS this specific raster was written in.
        moon_geographic = CRS.from_proj4(
            f"+proj=longlat +a={R_MOON_M} +b={R_MOON_M} +no_defs +type=crs"
        )
        self._to_raster_crs = Transformer.from_crs(
            moon_geographic, self._ds.crs, always_xy=True
        )

    def _pixel_coords(self, lat_deg, lon_deg):
        x, y = self._to_raster_crs.transform(lon_deg, lat_deg)
        col, row = self._inv_transform * (x, y)
        return row, col

    def _bilinear(self, row, col):
        r0, c0 = int(np.floor(row)), int(np.floor(col))
        r1, c1 = r0 + 1, c0 + 1
        if r0 < 0 or c0 < 0 or r1 >= self._nrows or c1 >= self._ncols:
            raise ValueError(
                f"lat/lon maps to pixel (row={row:.1f}, col={col:.1f}), "
                f"outside this DEM tile ({self._nrows}x{self._ncols}) -- "
                "wrong tile for this coordinate, or point is off-tile."
            )
        dr, dc = row - r0, col - c0
        h00, h01 = self._array[r0, c0], self._array[r0, c1]
        h10, h11 = self._array[r1, c0], self._array[r1, c1]
        top = h00 + dc * (h01 - h00)
        bot = h10 + dc * (h11 - h10)
        return float(top + dr * (bot - top))

    def radius(self, lat_deg, lon_deg):
        """Distance from the Moon's center at this lat/lon, in meters."""
        row, col = self._pixel_coords(lat_deg, lon_deg)
        value = self._bilinear(row, col)
        return value if self.offset_is_radius else value + R_MOON_M

    def height(self, lat_deg, lon_deg):
        """Height above/below the mean lunar radius, in meters."""
        row, col = self._pixel_coords(lat_deg, lon_deg)
        value = self._bilinear(row, col)
        return value - R_MOON_M if self.offset_is_radius else value

    def close(self):
        self._ds.close()

    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()


def surface_point_to_spice(lat_deg, lon_deg, radius_m,
                            body_frame="MOON_ME",
                            et=None, output_frame=None):
    """
    Package a (lat, lon, radius) triple as a SPICE Cartesian vector.

    body_frame : the body-fixed frame your lat/lon/radius are defined in.
        Default 'MOON_ME' -- see module docstring. Common overrides:
        'MOON_PA' (dynamical principal-axis frame, DE440-matched) or
        'IAU_MOON' (generic IAU rotation model, no extra kernels needed).
    et : ephemeris time (TDB seconds past J2000). Only required if
        output_frame is given.
    output_frame : if set (e.g. 'J2000'), rotates the body-fixed vector
        into this frame at time `et` via spice.pxform. This is what you
        want before differencing against a DSN station's J2000 position
        for a light-time calc.

    Returns a numpy array, position in meters, in body_frame (if
    output_frame is None) or in output_frame (if given).
    """
    lat_rad = np.radians(lat_deg)
    lon_rad = np.radians(lon_deg)
    pos_km = np.array(spice.latrec(radius_m / 1000.0, lon_rad, lat_rad))

    if output_frame is None:
        return pos_km * 1000.0

    if et is None:
        raise ValueError("et is required when output_frame is given")

    rot = spice.pxform(body_frame, output_frame, et)
    return (rot @ pos_km) * 1000.0


def dem_point_to_spice(dem: LunarDEM, lat_deg, lon_deg,
                        body_frame="MOON_ME",
                        et=None, output_frame=None):
    """
    Convenience wrapper: LDEM lookup + surface_point_to_spice in one call.

    Example
    -------
        spice.furnsh("your_metakernel.tm")
        dem = LunarDEM("LDEM_80S_20M.tif", scaling_factor=0.5, offset=R_MOON_M)
        et = spice.str2et("2026-07-01T00:00:00")
        pos_j2000_m = dem_point_to_spice(
            dem, lat_deg=-88.5425, lon_deg=5.4526,
            body_frame="MOON_ME", et=et, output_frame="J2000"
        )
    """
    r = dem.radius(lat_deg, lon_deg)
    return surface_point_to_spice(
        lat_deg, lon_deg, r,
        body_frame=body_frame, et=et, output_frame=output_frame,
    )

