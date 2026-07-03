def delay(
    dss_number: int,
    ground_time: str,
    lems_location: tuple[float, float],
) -> tuple[float, float]:
    """
    Calculate the propagation delay from the given DSS-## antenna to the given location on the
    moon for a signal sent at the given time on the ground.

    `dss_number` is the DSN antenna.
    EXAMPLE: `17`

    `ground_time` is an ISO8601 string with UTC time zone information and 4 subsecond decimal points
    EXAMPLE: `"2026-06-30T12:34:56.7890+00:00"`

    `lems_location` is the latitude/longitude of the LEMS landing zone.
    EXAMPLE: `(-88.123456, 123.456789)`

    Return value is a tuple of `(flight_time, error)`
    `flight_time` is the speed-of-light transmission time from the antenna to the LEMS location
    `error` is an estimate of the maximum total accumulated uncertainty
    """
    # TODO: Implement this

    from dateutil import parser
    from datetime import timezone
    from pathlib import Path

    import spiceypy as spice
    import numpy as np
    from resources.dem_to_spice import LunarDEM, surface_point_to_spice, dem_point_to_spice, R_MOON_M

    dem_scaling_factor = 0.5
    dem_offset = 1737400
    max_iter = 50
    tol_s = 1e-9
    curr_folder_location = Path(__file__).resolve().parent
    metakernel_location = curr_folder_location / "resources" / "kernels" / "metakernel.tm"
    dem_location = curr_folder_location / "resources" / "LDEM_80S_20M.JP2"
    error_dem_location = curr_folder_location / "resources" / "LDEM_80S_20MPP_ADJ_ERR.tiff"
    spice_station_id = str(399000 + dss_number) # Ex: DSN-14 -> 14 -> 39014
    body_frame = "MOON_ME" # this matches the convention of the lat/lon

    spice.furnsh(str(metakernel_location)) # only takes strings; loads all the kernels
    c_km_s = spice.clight()

    dt = parser.parse(ground_time).astimezone(timezone.utc)
    dt_iso_str = dt.strftime("%Y-%m-%dT%H:%M:%S.%f") 
    et = spice.str2et(dt_iso_str)

    dem = LunarDEM(dem_location, scaling_factor=dem_scaling_factor, offset=dem_offset)
    r_km = dem.radius(lems_location[0], lems_location[1]) / 1000.0
    lat_rad, lon_rad = np.radians(lems_location[0]), np.radians(lems_location[1])
    body_fixed_km = np.array(spice.latrec(r_km, lon_rad, lat_rad)) # vector for requested location
 
    # error_dem = LunarDEM(error_dem_location, scaling_factor=1, offset=0)
    # height_error_km = error_dem.radius(lems_location[0], lems_location[1]) / 1000.0

    lt = 0.0
    for _ in range(max_iter):
        et_arrival = et + lt  # uplink: signal reaches the Moon after lt
 
        rot = spice.pxform(body_frame, "J2000", et_arrival)
        surface_offset_j2000 = rot @ body_fixed_km
        moon_center, _ = spice.spkpos("MOON", et_arrival, "J2000", "NONE", spice_station_id)
        vec = np.array(moon_center) + surface_offset_j2000
 
        distance_km = float(np.linalg.norm(vec))
        lt_new = distance_km / c_km_s
 
        if abs(lt_new - lt) < tol_s:
            lt = lt_new
            break
        lt = lt_new
    else:
        raise RuntimeError(
            f"Light time iteration did not converge within {max_iter} steps "
            f"(last delta was {abs(lt_new - lt):.3e} s) -- check kernel coverage"
        )
 
    return (lt, 0.0)



if __name__ == "__main__":
    DSS_STATION = 17
    GROUND_TIME = "2027-01-01T00:00:00.0000+00:00"
    LEMS_LOCATION = (-88.0, 0.0)

    time_in_flight, time_in_flight_error = delay(
        dss_number=DSS_STATION,
        ground_time=GROUND_TIME,
        lems_location=LEMS_LOCATION,
    )

    print(time_in_flight)
