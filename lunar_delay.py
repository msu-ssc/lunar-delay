def delay(
    dss_number: int, ground_time: str, lems_location: tuple[float, float]
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
    return (0.0, 0.0)
