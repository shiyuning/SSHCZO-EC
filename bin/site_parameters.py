#!/usr/bin/env python3

# Site specific parameters
SITE = 'SH'                         # Site name
CSAT3_AZIMUTH = 199.0               # Direction between CSAT3 y axis and true north direction
AVERAGING_PERIOD_MINUTES = 30.0     # Averaging period in minutes
FREQUENCY_HZ = 10.0                 # Sampling frequency in Hz

# Name of columns in the high-frequency data file
TIME = 'TmStamp'                    # Time
U = 'Ux'                            # Wind speed in x direction
V = 'Uy'                            # Wind speed in y direction
W = 'Uz'                            # Wind speed in z direction
TS = 'Ts'                           # Sonic temperature
H2O = 'h2o'                         # Water vapor
CO2 = 'co2'                         # Carbon dioxide
# Other parameters for the high-frequency data file
SKIP_ROWS = None                    # Number of rows to skip in the data file (None for automatic detection. Do not skip the header line with column names)
COMMENT = '#'                       # Comment line marker


# Name of columns in the pressure data file
PRESSURE_TIME = 'TmStamp'
PRESSURE = 'pressure_irga_mean'
# Other parameters for the pressure data file
PRESSURE_SKIP_ROWS = None                    # Number of rows to skip in the data file (None for automatic detection. Do not skip the header line with column names)
PRESSURE_COMMENT = '#'                       # Comment line marker
