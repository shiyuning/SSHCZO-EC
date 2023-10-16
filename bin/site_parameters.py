#!/usr/bin/env python3

# Site specific parameters
SITE = 'US-SSH'                         # Site name
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

# Data file Diag configuration for US-SSH
#
# 11 | 10 |  9 |  8 |  7 |  6 |  5 |  4 |  3 |  2 |  1 |  0 |
#    CSAT3 flags    |    IRGA flags     |    AGC/6.25       |
# CSAT3:
# 9: lost trigger special case
# 10: no data special case
# 11: wrong CSAT3 embedded code special case
# 12: SDM error special case
# 13: NaN special case
#
# IRGA:
# 1000: chopper
# 0100: detector
# 0010: pll
# 0001: sync
ANEMOMETER_FILTER = lambda x: (x['diag'] & 3840) == 0   # 3840 = 1111 0000 0000
IRGA_FILTER = lambda x: (x['diag'] & 240) == 0          # 240 = 0000 1111 0000

# Name of columns in the pressure data file
PRESSURE_TIME = 'TmStamp'
PRESSURE = 'pressure_irga_mean'
TA = 'T_hmp_mean'
# Other parameters for the pressure data file
PRESSURE_SKIP_ROWS = None           # Number of rows to skip in the data file (None for automatic detection. Do not skip the header line with column names)
PRESSURE_COMMENT = '#'              # Comment line marker
