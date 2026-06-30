from lunar_delay import delay

DSS_STATION = 17
GROUND_TIME = "2027-01-01T00:00:00.0000+00:00"
LEMS_LOCATION = (-87.0, 0.0)

time_in_flight, time_in_flight_error = delay(
    dss_number=DSS_STATION,
    ground_time=GROUND_TIME,
    lems_location=LEMS_LOCATION,
)

print(f"----- GIVEN -----")
print(f"{DSS_STATION=}")
print(f"{GROUND_TIME=}")
print(f"{LEMS_LOCATION=}")
print(f"----- CALCULATED -----")
print(f"{time_in_flight=}")
print(f"{time_in_flight_error=}")
