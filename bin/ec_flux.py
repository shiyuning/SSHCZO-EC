#!/usr/bin/env python3
import argparse
import numpy as np
import pandas as pd
from datetime import datetime, timedelta
from dateutil.relativedelta import relativedelta
from quality_control import quality_control
from unit_vector_k import unit_vector_k
from calculate_fluxes import calculate_fluxes
from write_csv import *
from site_parameters import *


def get_pressure(start, end, df):
    """Get average air pressure for the averaging time period
    """
    sub_df = df[(df[TIME] >= start) & (df[TIME] < end)]

    if len(sub_df) == 0:
        print('  No pressure data available')
        return None

    return sub_df[PRESSURE].mean()


def main(params):
    files = params['files'][0]
    start_of_month = params['month']
    end_of_month = start_of_month + relativedelta(months=1)
    pres_file = params['WPL']

    # Read all files
    df = pd.DataFrame()
    print(f'Reading {len(files)} {"files" if len(files) > 1 else "file"}:')
    for f in files:
        print(f'  {f}')
        _df = pd.read_csv(
            f,
            skiprows=SKIP_ROWS,
            comment=COMMENT,
        )
        df = pd.concat([df, _df], ignore_index=True)
    df[TIME] = pd.to_datetime(df[TIME])
    df = df[(df[TIME] >= start_of_month) & (df[TIME] < end_of_month)]

    if len(df) == 0:
        print('No data available')
        return

    diag_file = f'{SITE}_{start_of_month.strftime("%Y-%m")}_flux_diag.csv'
    flux_file = f'{SITE}_{start_of_month.strftime("%Y-%m")}_flux.csv'

    initialize_csv_files(flux_file, diag_file)

    if pres_file:
        WPL = True
        pres_df = pd.read_csv(
            pres_file,
            skiprows=PRESSURE_SKIP_ROWS,
            comment=PRESSURE_COMMENT,
        )
        pres_df[PRESSURE_TIME] = pd.to_datetime(pres_df[PRESSURE_TIME])
    else:
        pres_df = pd.DataFrame()
        WPL = False

    # Data file Diag configuration
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

    # Determine unit vector k of planar fit coordinate
    # (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 3, Section 3)
    # unit_k is unit vector parallel to new coordinate z axis
    df = df[df['diag'] < 256]  # Remove bad CSAT3 data
    unit_k, b0 = unit_vector_k(df[U], df[V], df[W])

    periods = pd.date_range(start_of_month,end_of_month,freq=f'{AVERAGING_PERIOD_MINUTES}min').to_list()
    for i in range(len(periods) - 1):
        start = periods[i]
        end = periods[i + 1]
        sub_df = df[(df[TIME] >= start) & (df[TIME] < end) & (df['diag'] < 16)]
        if len(sub_df) == 0:
            continue

        print(f'\n{end.strftime("%Y-%m-%d %H:%M")} {len(sub_df)}')

        flags = quality_control(unit_k, sub_df)

        pressure_kpa = get_pressure(start, end, pres_df) if WPL else np.nan
        fluxes = calculate_fluxes(unit_k, pressure_kpa, sub_df)



def _main():
    parser = argparse.ArgumentParser(description='Process EC data')
    parser.add_argument(
        '--month',
        required=True,
        type=lambda s: datetime.strptime(s, '%Y-%m'),
        help='Month YYYY-MM',
    )
    parser.add_argument(
        '-f',
        '--files',
        action='append',
        nargs='+',
    )
    parser.add_argument(
        '--WPL',
        type=str,
    )
    args = parser.parse_args()

    main(vars(args))


if __name__ == '__main__':
    _main()
