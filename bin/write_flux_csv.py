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

def calculate_flags(row):
    integrated_flags = {}
    if row['instrument'] == BADVAL:
        for flag in FLAGS:
            integrated_flags[flag] = int(BADVAL)

        return integrated_flags

    for d in DIAGNOSTICS:
        row[d] = float(row[d])

    flags = {}
    flags['instrument'] = 1 if (row['instrument'] < QC_THRESHOLDS['instrument'] or row['instrument'] > 1.0) else 0
    for var in VARIABLES:
        flags[f'{var}_spikes'] = 1 if (row[f'{var}_spikes'] > QC_THRESHOLDS['spike']) else 0
        flags[f'{var}_resolution'] = 1 if (row[f'{var}_resolution'] > QC_THRESHOLDS['empty_bins'])  else 0
        flags[f'{var}_dropouts'] = 1 if (row[f'{var}_dropouts'] > QC_THRESHOLDS['dropouts'] or row[f'{var}_extreme_dropouts'] > QC_THRESHOLDS['extreme_dropouts']) else 0
        # The range of kurtosis is (1, 8) in Vickers and Mahrt (1997). It is equivalent to having (kurtosis - 4.5) in the
        # range (-3.5, 3.5) in this implementation.
        flags[f'{var}_higher_moments'] = 1 if (abs(row[f'{var}_skewness']) > QC_THRESHOLDS['skewness'] or abs(row[f'{var}_kurtosis'] - 4.5) > QC_THRESHOLDS['kurtosis']) else 0
        flags[f'{var}_discontinuities'] = 1 if (row[f'{var}_haar_mean'] > QC_THRESHOLDS['discontinuities'] or row[f'{var}_haar_variance'] > QC_THRESHOLDS['discontinuities']) else 0

    flags['u_absolute_limits'] = 1 if row['u_max'] > 30.0 else 0
    flags['v_absolute_limits'] = 1 if row['v_max'] > 30.0 else 0
    flags['w_absolute_limits'] = 1 if row['w_max'] > 10.0 else 0
    flags['tsonic_absolute_limits'] = 1 if row['tsonic_max'] > 60.0 or row['tsonic_min'] < -50.0 or row['tsonic_max'] - row['tsonic_min'] > 10.0 else 0
    flags['co2_absolute_limits'] = 1 if row['co2_max'] > 950.0 or row['co2_min'] < 550.0 or row['co2_max'] - row['co2_min'] > 120.0 else 0
    flags['h2o_absolute_limits'] = 1 if row['h2o_max'] > 35.0 or row['h2o_max'] < 2.5 or row['h2o_max'] - row[f'h2o_min'] > 8.0 else 0

    flags['nonstationary'] = 1 if (
        row['rnu'] > QC_THRESHOLDS['relative_nonstationarity'] or
        row['rnv'] > QC_THRESHOLDS['relative_nonstationarity'] or
        row['rns'] > QC_THRESHOLDS['relative_nonstationarity'] or
        row['wind_speed_reduction'] < QC_THRESHOLDS['wind_speed_reduction']) else 0

    for flag in FLAGS:
        var = flag.split('_')[0]
        integrated_flags[flag] = (
            flags['instrument'] * 1 +
            flags[f'{var}_spikes'] * 2 +
            flags[f'{var}_resolution'] * 4 +
            flags[f'{var}_dropouts'] * 8 +
            flags[f'{var}_absolute_limits'] * 16 +
            flags[f'{var}_higher_moments'] * 32 +
            flags[f'{var}_discontinuities'] * 64
        )
        if var in ['u', 'v', 'w']:
            integrated_flags[flag] += flags['nonstationary'] * 128

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
