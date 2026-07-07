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
 
import warnings

import numpy as np
import rasterio
from rasterio.windows import Window
from pyproj import CRS, Transformer
import spiceypy as spice

# Mean lunar radius -- matches your existing pipeline's `offset` convention
# and the spherical datum LOLA polar-stereographic products are built on.
R_MOON_M = 1737400.0

def _resolve_path(path):
    """Local file path unchanged; http(s) URL auto-prefixed with GDAL's
    /vsicurl/ so rasterio streams it (byte-range requests) instead of
    requiring you to download it first. Matches the convention used in
    view_ldem.py.
    """
    if isinstance(path, str) and (path.startswith("http://") or path.startswith("https://")):
        return "/vsicurl/" + path
    return path

def _crs_equivalent(crs_a, crs_b):
    """True if two rasterio CRS objects describe the same projection, even
    if their WKT differs cosmetically (different tool/authority naming --
    e.g. an ArcGIS-style "Moon (2015) - Sphere / Ocentric / South Polar"
    WKT vs a GDAL-style "Moon2000_spole" WKT for the identical stereographic
    projection on the identical sphere). Plain `crs_a == crs_b` is too
    strict for this -- it compares WKT/name details that can legitimately
    differ between a DEM and its companion error raster when they were
    produced by different pipelines, even though the actual projection
    (proj type, radius, pole, meridian, false easting/northing) is
    identical. Falls back to comparing the PROJ4 definition, which is
    what the coordinate math in this module actually relies on.
    """
    if crs_a == crs_b:
        return True
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", UserWarning)
            return (
                CRS.from_wkt(crs_a.to_wkt()).to_proj4()
                == CRS.from_wkt(crs_b.to_wkt()).to_proj4()
            )
    except Exception:
        return False

# def _derive_err_path(path):
#     """LDEM_<LAT>_<PIX>MPP_ADJ.TIF -> LDEM_<LAT>_<PIX>MPP_ADJ_ERR.TIF.
#     Matches the naming convention used across the PGDA product page and
#     view_ldem.py's derive_err_url. Preserves the original .tif/.TIF case
#     of the extension rather than forcing uppercase, so this works for
#     real PGDA files (always uppercase .TIF) and locally-downloaded/
#     renamed files (which may be lowercase .tif) alike -- matters on
#     case-sensitive filesystems.
#     """
#     if path.upper().endswith("_ERR.TIF"):
#         return path
#     if path.upper().endswith(".TIF"):
#         ext = path[-4:]  # preserve actual case, e.g. ".tif" or ".TIF"
#         return path[: -len(ext)] + "_ERR" + ext
#     raise ValueError(
#         f"Can't auto-derive an _ERR path from: {path!r}; pass error_path explicitly."
#     )
 
 
def _windowed_bilinear(ds, band, row, col, scaling_factor, offset, nrows, ncols):
    """Bilinear interpolation using only the 2x2 pixel neighborhood around
    (row, col), read via a rasterio windowed read -- shared by LunarDEM's
    height and error lookups so both get the same memory-efficient,
    remote-COG-friendly access pattern.
    """
    r0, c0 = int(np.floor(row)), int(np.floor(col))
    r1, c1 = r0 + 1, c0 + 1
    if r0 < 0 or c0 < 0 or r1 >= nrows or c1 >= ncols:
        raise ValueError(
            f"lat/lon maps to pixel (row={row:.1f}, col={col:.1f}), "
            f"outside this raster ({nrows}x{ncols}) -- wrong tile for "
            "this coordinate, or point is off-tile."
        )
    window = Window(c0, r0, 2, 2)
    raw = ds.read(band, window=window).astype(np.float32)
    vals = raw * scaling_factor + offset
    h00, h01 = vals[0, 0], vals[0, 1]
    h10, h11 = vals[1, 0], vals[1, 1]
    dr, dc = row - r0, col - c0
    top = h00 + dc * (h01 - h00)
    bot = h10 + dc * (h11 - h10)
    return float(top + dr * (bot - top))

class LunarDEM:
    """
    Thin wrapper around a single LDEM raster (e.g. LDEM_80S_20M.tif, or an
    https:// URL to a remote COG) that gives you bilinearly-interpolated
    elevation lookups by lat/lon instead of row/col -- the lat/lon ->
    raster-CRS -> pixel chain is handled for you, using whatever CRS is
    actually baked into the file (so this works for north- or south-polar
    tiles, or anything else, without hardcoding a projection like the
    existing scripts do).
 
    Each radius()/height() call performs a small windowed read (2x2 pixels)
    rather than loading the whole raster -- see the module docstring for
    why that matters for the multi-GB LOLA products.
 
    Mirrors the scaling_factor/offset convention already used in
    slope_analysis_coarse.py and ephemeris_generation.py:
        stored_value = raw_pixel * scaling_factor + offset
 
    Set offset_is_radius=True (default, matches your existing calls with
    offset=1737400.0) if `offset` already represents the Moon's mean radius,
    i.e. stored_value is the absolute radius from the Moon's center. Set it
    to False if `offset` is 0 (or some local datum) and stored_value is a
    height above/below the mean radius instead.
    """
 
    def __init__(self, path, scaling_factor=0.5, offset=R_MOON_M,
                 offset_is_radius=True, band=1,
                 error_path=None, error_scaling_factor=1.0, error_offset=0.0,
                 error_band=1):
        self._orig_path = path
        self._ds = rasterio.open(_resolve_path(path))
        self._transform = self._ds.transform
        self._inv_transform = ~self._transform
        self._band = band
        self._scaling_factor = scaling_factor
        self._offset = offset
        self.offset_is_radius = offset_is_radius
        self._nrows = self._ds.height
        self._ncols = self._ds.width
 
        # Error raster is opened lazily on first .error() call -- see
        # _ensure_error_ds -- so constructing a LunarDEM never does an
        # extra (possibly remote) open you didn't ask for.
        self._err_ds = None
        self._error_path = error_path
        self._error_scaling_factor = error_scaling_factor
        self._error_offset = error_offset
        self._error_band = error_band
 
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
        return _windowed_bilinear(
            self._ds, self._band, row, col,
            self._scaling_factor, self._offset, self._nrows, self._ncols,
        )
 
    def _ensure_error_ds(self):
        if self._err_ds is None:
            path = self._error_path #or _derive_err_path(self._orig_path)
            ds = rasterio.open(_resolve_path(path))
 
            same_shape = (ds.width, ds.height) == (self._ncols, self._nrows)
            same_transform = all(
                abs(a - b) < 1e-6 for a, b in zip(ds.transform, self._transform)
            )
            same_crs = _crs_equivalent(ds.crs, self._ds.crs)
 
            if not (same_shape and same_transform and same_crs):
                ds.close()
                raise ValueError(
                    f"Error raster grid ({ds.width}x{ds.height}, transform={ds.transform}, "
                    f"crs={ds.crs}) doesn't match height raster grid "
                    f"({self._ncols}x{self._nrows}, transform={self._transform}, "
                    f"crs={self._ds.crs}). Pass error_path explicitly if the auto-derived "
                    f"path ({path!r}) is wrong for this tile."
                )
            self._err_ds = ds
        return self._err_ds
 
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
 
    def error(self, lat_deg, lon_deg):
        """Interpolated height error/uncertainty (in meters) at this
        lat/lon, read from the companion _ERR.TIF raster. Uses the same
        bilinear scheme as height() -- both because it's cheap/consistent,
        and because it's equivalent to formal error propagation under the
        (realistic, for adjacent grid cells) assumption that the four
        corner errors are correlated rather than independent.
        """
        err_ds = self._ensure_error_ds()
        # Same grid as the height raster (checked in _ensure_error_ds), so
        # the row/col computed from the height transform applies directly.
        row, col = self._pixel_coords(lat_deg, lon_deg)
        return _windowed_bilinear(
            err_ds, self._error_band, row, col,
            self._error_scaling_factor, self._error_offset,
            self._nrows, self._ncols,
        )
 
    def radius_and_error(self, lat_deg, lon_deg):
        """Convenience: (radius_m, error_m) in one call -- one pixel-coord
        computation shared between the two windowed reads.
        """
        row, col = self._pixel_coords(lat_deg, lon_deg)
        value = self._bilinear(row, col)
        r = value if self.offset_is_radius else value + R_MOON_M
        err_ds = self._ensure_error_ds()
        err = _windowed_bilinear(
            err_ds, self._error_band, row, col,
            self._error_scaling_factor, self._error_offset,
            self._nrows, self._ncols,
        )
        return r, err
 
    def close(self):
        self._ds.close()
        if self._err_ds is not None:
            self._err_ds.close()
 
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
 
    Also works straight off a remote COG, no download step:
        dem = LunarDEM(
            "https://pgda.gsfc.nasa.gov/data/LOLA_20mpp/LDEM_80S_20MPP_ADJ.TIF",
            scaling_factor=0.5, offset=R_MOON_M,
        )
    """
    r = dem.radius(lat_deg, lon_deg)
    return surface_point_to_spice(
        lat_deg, lon_deg, r,
        body_frame=body_frame, et=et, output_frame=output_frame,
    )