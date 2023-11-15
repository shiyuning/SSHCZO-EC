#!/usr/bin/env python3
import argparse
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta
from quality_control import VARIABLES
from write_diag_csv import BADVAL
from write_diag_csv import OUTPUT_VARIABLES
from write_diag_csv import DIAGNOSTICS
from site_parameters import *

FLAGS = ['u_flag', 'v_flag', 'w_flag', 'tsonic_flag', 'h2o_flag', 'co2_flag']

def calculate_flags(diagnostics):
    integrated_flags = {}

    # Special case with no data
    if diagnostics['irga'] == BADVAL:
        for flag in FLAGS:
            integrated_flags[flag] = int(BADVAL)

        return integrated_flags

    for d in DIAGNOSTICS:
        diagnostics[d] = float(diagnostics[d]) if diagnostics[d] != BADVAL else diagnostics[d]

    flags = {}

    # Instrument flags
    flags['3dua'] = 1 if (diagnostics['3dua'] < QC_THRESHOLDS['instrument'] or diagnostics['3dua'] > 1.0) else 0
    flags['irga'] = 1 if (diagnostics['irga'] < QC_THRESHOLDS['instrument'] or diagnostics['irga'] > 1.0) else 0

    for var in VARIABLES:
        # Spikes
        if diagnostics[f'{var}_spikes'] == BADVAL:
            flags[f'{var}_spikes'] = 0
        else:
            flags[f'{var}_spikes'] = 1 if (diagnostics[f'{var}_spikes'] > QC_THRESHOLDS['spike']) else 0

        # Amplitude resolution
        if diagnostics[f'{var}_resolution'] == BADVAL:
            flags[f'{var}_resolution'] = 0
        else:
            flags[f'{var}_resolution'] = 1 if (diagnostics[f'{var}_resolution'] > QC_THRESHOLDS['empty_bins'])  else 0

        # Dropouts
        if diagnostics[f'{var}_dropouts'] == BADVAL:
            flags[f'{var}_dropouts'] = 0
        else:
            flags[f'{var}_dropouts'] = 1 if (diagnostics[f'{var}_dropouts'] > QC_THRESHOLDS['dropouts'] or diagnostics[f'{var}_extreme_dropouts'] > QC_THRESHOLDS['extreme_dropouts']) else 0

        # Higher moments
        # The range of kurtosis is (1, 8) in Vickers and Mahrt (1997). It is equivalent to having (kurtosis - 4.5) in the
        # range (-3.5, 3.5) in this implementation.
        if diagnostics[f'{var}_skewness'] == BADVAL:
            flags[f'{var}_higher_moments'] = 0
        else:
            flags[f'{var}_higher_moments'] = 1 if (abs(diagnostics[f'{var}_skewness']) > QC_THRESHOLDS['skewness'] or abs(diagnostics[f'{var}_kurtosis'] - 4.5) > QC_THRESHOLDS['kurtosis']) else 0

        # Discontinuities
        if diagnostics[f'{var}_haar_mean'] == BADVAL:
            flags[f'{var}_discontinuities'] = 0
        else:
            flags[f'{var}_discontinuities'] = 1 if (diagnostics[f'{var}_haar_mean'] > QC_THRESHOLDS['discontinuities'] or diagnostics[f'{var}_haar_variance'] > QC_THRESHOLDS['discontinuities']) else 0

    # Absolute limits
    if diagnostics['u_max'] == BADVAL:
        flags['u_absolute_limits'] = 0
    else:
        flags['u_absolute_limits'] = 1 if diagnostics['u_max'] > 30.0 else 0

    if diagnostics['v_max'] == BADVAL:
        flags['v_absolute_limits'] = 0
    else:
        flags['v_absolute_limits'] = 1 if diagnostics['v_max'] > 30.0 else 0

    if diagnostics['w_max'] == BADVAL:
        flags['w_absolute_limits'] = 0
    else:
        flags['w_absolute_limits'] = 1 if diagnostics['w_max'] > 10.0 else 0

    if diagnostics['tsonic_min'] == BADVAL or diagnostics['tsonic_max'] == BADVAL:
        flags['tsonic_absolute_limits'] = 0
    else:
        flags['tsonic_absolute_limits'] = 1 if diagnostics['tsonic_max'] > 60.0 or diagnostics['tsonic_min'] < -50.0 or diagnostics['tsonic_max'] - diagnostics['tsonic_min'] > 10.0 else 0

    if diagnostics['co2_min'] == BADVAL or diagnostics['co2_max'] == BADVAL:
        flags['co2_absolute_limits'] = 0
    else:
        flags['co2_absolute_limits'] = 1 if diagnostics['co2_max'] > 950.0 or diagnostics['co2_min'] < 550.0 or diagnostics['co2_max'] - diagnostics['co2_min'] > 120.0 else 0

    if diagnostics['h2o_min'] == BADVAL or diagnostics['h2o_max'] == BADVAL:
        flags['h2o_absolute_limits'] = 0
    else:
        flags['h2o_absolute_limits'] = 1 if diagnostics['h2o_max'] > 35.0 or diagnostics['h2o_max'] < 2.5 or diagnostics['h2o_max'] - diagnostics[f'h2o_min'] > 8.0 else 0

    # Nonstationarity
    if diagnostics['wind_speed_reduction'] == BADVAL or diagnostics['rnu'] == BADVAL or diagnostics['rnv'] == BADVAL or diagnostics['rns'] == BADVAL:
        flags['nonstationary'] = 0
    else:
        flags['nonstationary'] = 1 if (
            diagnostics['rnu'] > QC_THRESHOLDS['relative_nonstationarity'] or
            diagnostics['rnv'] > QC_THRESHOLDS['relative_nonstationarity'] or
            diagnostics['rns'] > QC_THRESHOLDS['relative_nonstationarity'] or
            diagnostics['wind_speed_reduction'] < QC_THRESHOLDS['wind_speed_reduction']) else 0

    for flag in FLAGS:
        var = flag.split('_')[0]
        integrated_flags[flag] = (
            flags[f'{var}_spikes'] * 2 +
            flags[f'{var}_resolution'] * 4 +
            flags[f'{var}_dropouts'] * 8 +
            flags[f'{var}_absolute_limits'] * 16 +
            flags[f'{var}_higher_moments'] * 32 +
            flags[f'{var}_discontinuities'] * 64
        )
        if var in ['u', 'v', 'w']:
            integrated_flags[flag] += flags['nonstationary'] * 128
        if var in ['u', 'v', 'w', 'tsonic']:
            integrated_flags[flag] += flags['3dua'] * 1
        if var in ['co2', 'h2o']:
            integrated_flags[flag] += flags['irga'] * 1

    return integrated_flags


def apply_flags(row):
    fluxes = {}
    for var in OUTPUT_VARIABLES:
        fluxes[var]  = row[var]

        for flag in OUTPUT_VARIABLES[var]['flag']:
            if (row[flag] != 0): fluxes[var] = BADVAL

    return fluxes


def write_flux_file(diag_file):
    diag_df = pd.read_csv(diag_file, dtype=str)

    diag_df[FLAGS] = diag_df.apply(lambda x: calculate_flags(x), axis=1, result_type='expand')
    diag_df[list(OUTPUT_VARIABLES.keys())] = diag_df.apply(lambda x: apply_flags(x), axis=1, result_type='expand', )

    flux_file = diag_file.replace('_diag.csv', '.csv')
    diag_df[['TIMESTAMP_START', 'TIMESTAMP_END', *list(OUTPUT_VARIABLES.keys())]].to_csv(flux_file, index=False)

    flag_file = diag_file.replace('_diag.csv', '_flag.csv')
    diag_df[['TIMESTAMP_START', 'TIMESTAMP_END', *FLAGS]].to_csv(flag_file, index=False)

def _main():
    parser = argparse.ArgumentParser(description='Write flux data to csv')
    parser.add_argument(
        '-m',
        '--month',
        required=True,
        type=lambda s: datetime.strptime(s, '%Y-%m'),
        help='Month YYYY-MM',
    )

    start_of_month = parser.parse_args().month
    end_of_month = start_of_month + relativedelta(months=1)

    if AVERAGING_PERIOD_MINUTES == 30.0:
        resolution = 'HH'
    diag_file = f'{SITE}_{resolution}_{start_of_month.strftime("%Y%m%d%H%M")}_{end_of_month.strftime("%Y%m%d%H%M")}_diag.csv'

    write_flux_file(diag_file)


if __name__ == '__main__':
    _main()
