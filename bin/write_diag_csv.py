#!/usr/bin/env python3
import numpy as np
from quality_control import VARIABLES
from site_parameters import *

BADVAL = '-9999'
OUTPUT_VARIABLES = {
    'USTAR': {'format': '{:.3f}', 'flag': ['u_flag', 'v_flag', 'w_flag']},
    'WD': {'format': '{:.3f}', 'flag': ['u_flag', 'v_flag']},
    'WS': {'format': '{:.3f}', 'flag': ['u_flag', 'v_flag']},
    'FC': {'format': '{:.3f}', 'flag': ['w_flag', 'co2_flag']},
    'H': {'format': '{:.3f}', 'flag': ['w_flag', 'tsonic_flag']},
    'LE': {'format': '{:.3f}', 'flag': ['w_flag', 'h2o_flag']},
    'CO2': {'format': '{:.3f}', 'flag': ['co2_flag']},
    'H2O': {'format': '{:.3f}', 'flag': ['h2o_flag']},
    'PA': {'format': '{:.2f}', 'flag': []},
    'T_SONIC': {'format': '{:.2f}', 'flag': ['tsonic_flag']},
    'TA': {'format': '{:.2f}', 'flag': []},
}
DIAGNOSTICS = [
    'instrument',
    *[f'{str}_spikes' for str in VARIABLES],
    *[f'{str}_resolution' for str in VARIABLES],
    *[f'{str}_dropouts' for str in VARIABLES],
    *[f'{str}_extreme_dropouts' for str in VARIABLES],
    'u_max', 'v_max', 'w_max', 'tsonic_min', 'tsonic_max', 'co2_min', 'co2_max', 'h2o_min', 'h2o_max',
    *[f'{str}_skewness' for str in VARIABLES],
    *[f'{str}_kurtosis' for str in VARIABLES],
    *[f'{str}_haar_mean' for str in VARIABLES],
    *[f'{str}_haar_variance' for str in VARIABLES],
    'wind_speed_reduction', 'rnu', 'rnv',  'rns',
]


def initialize_diag_file(diag_file):
    """Initialize the output CSV files
    """
    with open(diag_file, 'w') as f:
        f.write(','.join(['TIMESTAMP_START', 'TIMESTAMP_END', *list(OUTPUT_VARIABLES.keys()), *DIAGNOSTICS]) + '\n')


def write_diag_file(start, end, fluxes, diagnostics, diag_file):
    ts_start = start.strftime("%Y%m%d%H%M")
    ts_end = end.strftime("%Y%m%d%H%M")

    _diagnostics = {}

    if not fluxes:
        for var in OUTPUT_VARIABLES:
            fluxes[var] = BADVAL
        for d in DIAGNOSTICS:
            _diagnostics[d] = BADVAL
    else:
        for var in OUTPUT_VARIABLES:
            fluxes[var] = BADVAL if np.isnan(fluxes[var]) else OUTPUT_VARIABLES[var]['format'].format(fluxes[var])
        for d in DIAGNOSTICS:
            _diagnostics[d] = '{:.3f}'.format(diagnostics[d])

    with open(diag_file, 'a') as f:
        f.write(','.join([ts_start, ts_end, *list(fluxes.values()), *list(_diagnostics.values())]) + '\n')
