#!/usr/bin/env python3
import argparse
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta
from quality_control import quality_control
from unit_vectors import unit_vector_k, unit_vector_ij
from calculate_fluxes import calculate_fluxes
from write_csv import initialize_csv_files, write_csv_files
from site_parameters import *

def get_pressure(start, end, df):
    """Get average air pressure for the averaging time period
    """
    sub_df = df[(df[TIME] >= start) & (df[TIME] < end)]

    if len(sub_df) == 0:
        print('  No pressure data available')
        return None, None

    return sub_df[PRESSURE].mean(), sub_df[TA].mean()


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

    if AVERAGING_PERIOD_MINUTES == 30.0:
        resolution = 'HH'
    diag_file = f'{SITE}_{resolution}_{start_of_month.strftime("%Y%m%d%H%M")}_{end_of_month.strftime("%Y%m%d%H%M")}_diag.csv'
    flux_file = f'{SITE}_{resolution}_{start_of_month.strftime("%Y%m%d%H%M")}_{end_of_month.strftime("%Y%m%d%H%M")}.csv'

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

    # Determine unit vector k of planar fit coordinate
    # (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 3, Section 3)
    # unit_k is unit vector parallel to new coordinate z axis
    unit_k = unit_vector_k(df[ANEMOMETER_FILTER][U], df[ANEMOMETER_FILTER][V], df[ANEMOMETER_FILTER][W])

    periods = pd.date_range(start_of_month,end_of_month,freq=f'{AVERAGING_PERIOD_MINUTES}min').to_list()
    for k in range(len(periods) - 1):
        [start, end] = [periods[k], periods[k + 1]]
        print(f'\n{end.strftime("%Y-%m-%d %H:%M")}')

        sub_df = df[(df[TIME] >= start) & (df[TIME] < end) & ANEMOMETER_FILTER(df) & IRGA_FILTER(df)]
        if len(sub_df) == 0:
            fluxes = {}
            flags = {}
        else:
            # Determine unit vectors i and j of planar fit coordinate
            unit_i, unit_j = unit_vector_ij(unit_k, sub_df[U], sub_df[V], sub_df[W])

            # Quality control following Vickers and Mahrt (1997)
            flags = quality_control(unit_i, unit_j, unit_k, sub_df)

            # Read pressure data for WPL correction
            if WPL:
                pressure_kpa, ta_c = get_pressure(start, end, pres_df)
            else:
                pressure_kpa = ta_c = None

            # Calculate fluxes
            fluxes = calculate_fluxes(unit_i, pressure_kpa, ta_c, sub_df)

        write_csv_files(start, end, fluxes, flags, flux_file, diag_file)


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
