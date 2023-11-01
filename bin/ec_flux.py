#!/usr/bin/env python3
import argparse
import numpy  as np
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta
from quality_control import quality_control
from unit_vectors import unit_vector_k, unit_vector_ij
from calculate_fluxes import calculate_fluxes
from write_diag_csv import initialize_diag_file, write_diag_file
from write_flux_csv import write_flux_file
from site_parameters import *

def get_pressure(start, end, df):
    """Get average air pressure for the averaging time period
    """
    sub_df = df[(df['time'] >= start) & (df['time'] < end)]

    if len(sub_df[sub_df['pressure'].notna()]) > 0:
        pressure = np.nanmean(sub_df['pressure'].values)
    else:
        pressure = None

    if len(sub_df[sub_df['tair'].notna()]) > 0:
        tair = np.nanmean(sub_df['tair'].values)
    else:
        tair = None

    return pressure, tair


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
    # Rename columns to standard names
    df = df.rename(columns={
        TIME: 'time',
        U: 'u',
        V: 'v',
        W: 'w',
        T_SONIC: 'tsonic',
        CO2: 'co2',
        H2O: 'h2o',
    })
    df['time'] = pd.to_datetime(df['time'])

    df = df[(df['time'] >= start_of_month) & (df['time'] < end_of_month)]
    if len(df) == 0:
        print('No data available')
        exit()

    df['u'] = df['u'].map(lambda x: wind_speed_m_per_s(x))
    df['v'] = df['v'].map(lambda x: wind_speed_m_per_s(x))
    df['w'] = df['w'].map(lambda x: wind_speed_m_per_s(x))
    df['tsonic'] = df['tsonic'].map(lambda x: tsonic_celsius(x))
    df['h2o'] = df['h2o'].map(lambda x: h2o_mg_per_m3(x))
    df['co2'] = df['co2'].map(lambda x: co2_g_per_m3(x))

    if AVERAGING_PERIOD_MINUTES == 30.0:
        resolution = 'HH'
    diag_file = f'{SITE}_{resolution}_{start_of_month.strftime("%Y%m%d%H%M")}_{end_of_month.strftime("%Y%m%d%H%M")}_diag.csv'
    flux_file = f'{SITE}_{resolution}_{start_of_month.strftime("%Y%m%d%H%M")}_{end_of_month.strftime("%Y%m%d%H%M")}.csv'

    initialize_diag_file(diag_file)

    if pres_file:
        WPL = True
        pres_df = pd.read_csv(
            pres_file,
            skiprows=PRESSURE_SKIP_ROWS,
            comment=PRESSURE_COMMENT,
        )
        pres_df = pres_df.rename(columns={
            PRESSURE_TIME: 'time',
            PRESSURE: 'pressure',
            T_AIR: 'tair',
        })
        pres_df['pressure'] = pres_df['pressure'].map(lambda x: pressure_pa(x))
        pres_df['tair'] = pres_df['tair'].map(lambda x: tair_celsius(x))
        pres_df['time'] = pd.to_datetime(pres_df['time'])
    else:
        pres_df = pd.DataFrame()
        WPL = False

    # Determine unit vector k of planar fit coordinate
    # (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 3, Section 3)
    # unit_k is unit vector parallel to new coordinate z axis
    if len(df[ANEMOMETER_FILTER]) > 0:
        unit_k = unit_vector_k(df[ANEMOMETER_FILTER]['u'], df[ANEMOMETER_FILTER]['v'], df[ANEMOMETER_FILTER]['w'])
    else:
        unit_k = None

    periods = pd.date_range(start_of_month,end_of_month,freq=f'{AVERAGING_PERIOD_MINUTES}min').to_list()
    for k in range(len(periods) - 1):
        [start, end] = [periods[k], periods[k + 1]]
        print(f'\n{end.strftime("%Y-%m-%d %H:%M")}')

        sub_df = df[(df['time'] >= start) & (df['time'] < end)].copy().reset_index(drop=True)
        if len(sub_df) == 0:
            fluxes = {}
            diagnostics = {}
        else:
            # Determine unit vectors i and j of planar fit coordinate
            if len(sub_df[ANEMOMETER_FILTER]) > 0:
                unit_i, unit_j = unit_vector_ij(unit_k, sub_df[ANEMOMETER_FILTER]['u'], sub_df[ANEMOMETER_FILTER]['v'], sub_df[ANEMOMETER_FILTER]['w'])
            else:
                unit_i = unit_j = None

            # Quality control following Vickers and Mahrt (1997)
            diagnostics = quality_control(unit_i, unit_j, unit_k, sub_df)

            # Read pressure data for WPL correction
            if WPL:
                pressure_kpa, ta_c = get_pressure(start, end, pres_df)
            else:
                pressure_kpa = ta_c = None

            # Calculate fluxes
            fluxes = calculate_fluxes(unit_i, pressure_kpa, ta_c, sub_df, diagnostics)

        write_diag_file(start, end, fluxes, diagnostics, diag_file)

    write_flux_file(diag_file)


def _main():
    parser = argparse.ArgumentParser(description='Process EC data')
    parser.add_argument(
        '-m',
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
