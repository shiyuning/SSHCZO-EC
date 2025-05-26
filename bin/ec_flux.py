import argparse
import numpy  as np
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta
from config import AVERAGING_PERIOD_MINUTES, ANEMOMETER_FILTER, DIAG_FILE
from config import TIME, U, V, W, T_SONIC, CO2, H2O, SKIP_ROWS, COMMENT
from config import PRESSURE_TIME, PRESSURE, T_AIR, PRESSURE_SKIP_ROWS, PRESSURE_COMMENT
from config import wind_speed_m_per_s, tsonic_celsius, h2o_mg_per_m3, co2_g_per_m3
from config import pressure_pa, tair_celsius
from eddy_covariance import EddyCovariance, INSTANTANEOUS_VARIABLES
from unit_vectors import unit_vector_k
from write_flux_csv import write_flux_file

def read_monthly_data(fns: list[str], start_of_month: datetime, end_of_month: datetime) -> pd.DataFrame:
    df = pd.DataFrame()
    print(f'Reading {len(fns)} {"files" if len(fns) > 1 else "file"}:')
    for f in fns:
        print(f'  {f}')
        _df = pd.read_csv(
            f,
            skiprows=SKIP_ROWS,
            comment=COMMENT,
        )
        df = pd.concat([df, _df], ignore_index=True)

    # Rename columns to standard names
    df.rename(
        columns={TIME: 'time', U: 'u', V: 'v', W: 'w', T_SONIC: 'tsonic', CO2: 'co2', H2O: 'h2o'},
        inplace=True,
    )

    # In the new csv format, timestamps have a 'UTC' string at the end, that needs to be removed
    df['time'] = pd.to_datetime(df['time'].map(lambda x: x[:-4] if x.endswith(' UTC') else x))
    df = df[(df['time'] >= start_of_month) & (df['time'] < end_of_month)]
    if df.empty:
        raise ValueError('No data available.')

    # Mask bad (flagged) data with NaN using data record flags
    for v in INSTANTANEOUS_VARIABLES:
        df.loc[~v.filter(df), v.name] = np.nan

    # Convert to standard units
    df['u'] = df['u'].map(lambda x: wind_speed_m_per_s(x))
    df['v'] = df['v'].map(lambda x: wind_speed_m_per_s(x))
    df['w'] = df['w'].map(lambda x: wind_speed_m_per_s(x))
    df['tsonic'] = df['tsonic'].map(lambda x: tsonic_celsius(x))
    df['h2o'] = df['h2o'].map(lambda x: h2o_mg_per_m3(x))
    df['co2'] = df['co2'].map(lambda x: co2_g_per_m3(x))

    return df


def read_pressure_data(fn: str, start_of_month: datetime, end_of_month: datetime) -> pd.DataFrame:
    df = pd.read_csv(
        fn,
        skiprows=PRESSURE_SKIP_ROWS,
        comment=PRESSURE_COMMENT,
    )

    df.rename(
        columns={PRESSURE_TIME: 'time', PRESSURE: 'pressure', T_AIR: 'tair', },
        inplace=True,
    )

    df['time'] = pd.to_datetime(df['time'].map(lambda x: x[:-4] if x.endswith(' UTC') else x))
    df = df[(df['time'] >= start_of_month) & (df['time'] < end_of_month)]

    if df.empty:
        raise ValueError('No pressure data available.')

    # Convert to standard units
    df['pressure'] = df['pressure'].map(lambda x: pressure_pa(x))
    df['tair'] = df['tair'].map(lambda x: tair_celsius(x))

    return df


def main(params):
    files = params['files'][0]
    start_of_month = params['month']
    end_of_month = start_of_month + relativedelta(months=1)
    pressure_file = params['WPL']

    # Read all files
    df = read_monthly_data(files, start_of_month, end_of_month)

    if AVERAGING_PERIOD_MINUTES == 30:
        resolution = 'HH'
    elif AVERAGING_PERIOD_MINUTES == 60:
        resolution = 'HR'
    else:
        raise ValueError('Please use either 30 min or 60 min for averaging period.')

    pressure_df = read_pressure_data(pressure_file, start_of_month, end_of_month) if pressure_file else pd.DataFrame()

    # Determine unit vector k of planar fit coordinate
    # (Lee, L., W. Massman, and B. Law, 2004: Handbook of Micrometeorology, Chapt 3, Section 3)
    # unit_k is unit vector parallel to new coordinate z axis

    unit_k = unit_vector_k(df['u'].values, df['v'].values, df['w'].values) if not df[ANEMOMETER_FILTER].empty else None

    for time_block in pd.date_range(start_of_month, end_of_month, freq=f'{AVERAGING_PERIOD_MINUTES}min').to_list()[:-1]:
        # Create an EddyCovariance class for data processing
        eddy_covariance = EddyCovariance(time=time_block, unit_k=unit_k, data=df, pressure_data=pressure_df)

        if not eddy_covariance.data.empty:
            # Quality control following Vickers and Mahrt (1997)
            eddy_covariance.quality_control()

            # Calculate fluxes
            eddy_covariance.calculate_fluxes()

        # Write to diagnostic output file
        eddy_covariance.write_diag_file(
            first=time_block == start_of_month,
            fn=DIAG_FILE(resolution, start_of_month, end_of_month),
        )

    # Write flux output file
    write_flux_file(DIAG_FILE(resolution, start_of_month, end_of_month))


def _main():
    parser = argparse.ArgumentParser(description='Process eddy covariance data')
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
