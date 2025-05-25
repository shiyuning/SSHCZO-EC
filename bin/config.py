#!/usr/bin/env python3

# Site specific parameters
SITE = 'US-SSH'                     # Site name
CSAT3_AZIMUTH = 199.0               # Direction between CSAT3 y axis and true north direction
AVERAGING_PERIOD_MINUTES = 30       # Averaging period in minutes
FREQUENCY_HZ = 10                   # Sampling frequency in Hz

# Name of columns in the high-frequency data file
TIME = 'TmStamp'                    # Time
U = 'Ux'                            # Wind speed in x direction
V = 'Uy'                            # Wind speed in y direction
W = 'Uz'                            # Wind speed in z direction
T_SONIC = 'Ts'                      # Sonic temperature
H2O = 'h2o'                         # Water vapor
CO2 = 'co2'                         # Carbon dioxide
# Other parameters for the high-frequency data file
SKIP_ROWS = None                    # Number of rows to skip in the data file (None for automatic detection. Do not skip the header line with column names)
COMMENT = '#'                       # Comment line marker

# Unit conversions
# To correctly calculate fluxes, we require wind speed in m/s, sonic temperature in Celsius, water vapor in mg/m3, and
# carbon dioxide in g/m3. If the input is different from those units, please use the following functions to convert
# them.
def wind_speed_m_per_s(u): return u

def tsonic_celsius(t): return t

def h2o_mg_per_m3(h2o): return h2o

def co2_g_per_m3(co2): return co2

# Data file diag configuration for US-SSH
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
ANEMOMETER_FLAGS = 3840 # 3840 = 1111 0000 0000
ANEMOMETER_FILTER = lambda x: (x['diag'] & ANEMOMETER_FLAGS) == 0
IRGA_FLAGS = 240    # 240 = 0000 1111 0000
IRGA_FILTER = lambda x: (x['diag'] & IRGA_FLAGS) == 0

# Name of columns in the pressure data file
PRESSURE_TIME = 'TmStamp'           # Time
PRESSURE = 'pressure_irga_mean'     # Air pressure
T_AIR = 'T_hmp_mean'                # Air temperature
# Other parameters for the pressure data file
PRESSURE_SKIP_ROWS = None           # Number of rows to skip in the data file (None for automatic detection. Do not skip the header line with column names)
PRESSURE_COMMENT = '#'              # Comment line marker

# Unit conversions
# To correctly calculate fluxes, we require air pressure in Pascal, and air temperature in Celsius. If the input is
# different from those units, please use the following functions to convert them.
def pressure_pa(p): return p * 1000.0

def tair_celsius(t): return t

# Special value for flagged data in reports
BADVAL = '-9999'

# Thresholds for quality control
QC_THRESHOLDS = {
    'instrument': 0.95,
    'spike': 0.03,                  # Vickers and Mahrt (1997) value 0.01
    'empty_bins': 0.7,              # Vickers and Mahrt (1997) value 0.7
    'dropouts': 0.1,                # Vickers and Mahrt (1997) value 0.1
    'extreme_dropouts': 0.06,       # Vickers and Mahrt (1997) value 0.06
    'skewness': 3.0,                # Vickers and Mahrt (1997) value 2.0
    'kurtosis': 4.5,                # Vickers and Mahrt (1997) value 3.5
    'discontinuities': 3.0,         # Vickers and Mahrt (1997) value 3.0
    'wind_speed_reduction': 0.5,    # Vickers and Mahrt (1997) value 0.9
    'relative_nonstationarity': 3.0 # Vickers and Mahrt (1997) value 0.5
}

DIAG_FILE = lambda resolution, start, end: f'{SITE}_{resolution}_{start.strftime("%Y%m%d%H%M")}_{end.strftime("%Y%m%d%H%M")}_diag.csv'
FLUX_FILE =  lambda resolution, start, end: f'{SITE}_{resolution}_{start.strftime("%Y%m%d%H%M")}_{end.strftime("%Y%m%d%H%M")}.csv'
