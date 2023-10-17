#!/usr/bin/env python3
import numpy as np
from site_parameters import *

BADVAL = '-9999'
VARIABLES = {
    'USTAR': {'format': '{:.3f}', 'flag': ['U_FLAG', 'V_FLAG', 'W_FLAG']},
    'WD': {'format': '{:.3f}', 'flag': ['U_FLAG', 'V_FLAG']},
    'WS': {'format': '{:.3f}', 'flag': ['U_FLAG', 'V_FLAG']},
    'FC': {'format': '{:.3f}', 'flag': ['W_FLAG', 'CO2_FLAG']},
    'H': {'format': '{:.3f}', 'flag': ['W_FLAG', 'TS_FLAG']},
    'LE': {'format': '{:.3f}', 'flag': ['W_FLAG', 'H2O_FLAG']},
    'CO2': {'format': '{:.3f}', 'flag': ['CO2_FLAG']},
    'H2O': {'format': '{:.3f}', 'flag': ['H2O_FLAG']},
    'PA': {'format': '{:.2f}', 'flag': []},
    'T_SONIC': {'format': '{:.2f}', 'flag': ['TS_FLAG']},
    'TA': {'format': '{:.2f}', 'flag': []},
}

FLAGS = ['U_FLAG', 'V_FLAG', 'W_FLAG', 'TS_FLAG', 'H2O_FLAG', 'CO2_FLAG']


def initialize_csv_files(flux_file, diag_file):
    """Initialize the output CSV files
    """
    with open(diag_file, 'w') as f:
        f.write(','.join(['TIMESTAMP_START', 'TIMESTAMP_END', *list(VARIABLES.keys()), *FLAGS]) + '\n')

    with open(flux_file, 'w') as f:
        f.write(','.join(['TIMESTAMP_START', 'TIMESTAMP_END', *list(VARIABLES.keys())]) + '\n')

def write_csv_files(start, end, fluxes, flags, flux_file, diag_file):
    ts_start = start.strftime("%Y%m%d%H%M")
    ts_end = end.strftime("%Y%m%d%H%M")
    integrated_flags = {}

    if not fluxes:
        for var in VARIABLES:
            fluxes[var] = BADVAL
        for flag in FLAGS:
            integrated_flags[flag] = 0
    else:
        for var in VARIABLES:
            fluxes[var] = BADVAL if np.isnan(fluxes[var]) else VARIABLES[var]['format'].format(fluxes[var])

        for var, flag in zip([U, V, W, TS, H2O, CO2], FLAGS):
            integrated_flags[flag] = (
                flags['instrument'] * 1 +
                flags[f'{var}_spikes'] * 2 +
                flags[f'{var}_resolution'] * 4 +
                flags[f'{var}_dropouts'] * 8 +
                flags[f'{var}_absolute_limits'] * 16 +
                flags[f'{var}_higher_moments'] * 32 +
                flags[f'{var}_discontinuities'] * 64 +
                flags['nonstationary'] * 128
            )

    with open(diag_file, 'a') as f:
        f.write(','.join([ts_start, ts_end, *list(fluxes.values()), *[str(integrated_flags[flag]) for flag in integrated_flags]]) + '\n')

    for var in VARIABLES:
        for flag in VARIABLES[var]['flag']:
            if integrated_flags[flag] > 0:
                fluxes[var] = BADVAL
                break

    with open(flux_file, 'a') as f:
        f.write(','.join([ts_start, ts_end, *list(fluxes.values())]) + '\n')
