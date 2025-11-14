import argparse
import numpy as np
import pandas as pd
from datetime import datetime
from dateutil.relativedelta import relativedelta
from eddy_covariance import INSTANTANEOUS_VARIABLES, OUTPUT_VARIABLES
from config import SITE, AVERAGING_PERIOD_MINUTES, QC_THRESHOLDS, BADVAL, DIAG_FILE

FLAGS = {
    '3dua': lambda x: np.nan if np.isnan(x['3dua']) else x['3dua'] < QC_THRESHOLDS['instrument'] or x['3dua'] > 1.0,
    'irga': lambda x: np.nan if np.isnan(x['irga']) else x['irga'] < QC_THRESHOLDS['instrument'] or x['irga'] > 1.0,
    'spikes': lambda v, x: np.isnan(x[f'{v}_spikes']) or (x[f'{v}_spikes'] > QC_THRESHOLDS['spike']),
    'resolution': lambda v, x: np.isnan(x[f'{v}_resolution']) or (x[f'{v}_resolution'] > QC_THRESHOLDS['empty_bins']),
    'dropouts': lambda v, x: any(np.isnan([x[f'{v}_dropouts'], x[f'{v}_extreme_dropouts']])) or (x[f'{v}_dropouts'] > QC_THRESHOLDS['dropouts'] or x[f'{v}_extreme_dropouts'] > QC_THRESHOLDS['extreme_dropouts']),
    # Higher moments: The range of kurtosis is (1, 8) in Vickers and Mahrt (1997). It is equivalent to having
    # (kurtosis - 4.5) in the range (-3.5, 3.5) in this implementation.
    'higher_moments': lambda v, x: any(np.isnan([x[f'{v}_skewness'], x[f'{v}_kurtosis']])) or (abs(x[f'{v}_skewness']) > QC_THRESHOLDS['skewness'] or abs(x[f'{v}_kurtosis'] - 4.5) > QC_THRESHOLDS['kurtosis']),
    'discontinuities': lambda v, x: any(np.isnan([x[f'{v}_haar_mean'], x[f'{v}_haar_variance']])) or (x[f'{v}_haar_mean'] > QC_THRESHOLDS['discontinuities'] or x[f'{v}_haar_variance'] > QC_THRESHOLDS['discontinuities']),
    'u_absolute_limits': lambda x: np.isnan(x['u_max']) or (x['u_max'] > 30.0),
    'v_absolute_limits': lambda x: np.isnan(x['v_max']) or (x['v_max'] > 30.0),
    'w_absolute_limits': lambda x: np.isnan(x['w_max']) or (x['w_max'] > 10.0),
    'tsonic_absolute_limits': lambda x: any(np.isnan([x['tsonic_max'], x['tsonic_min']])) or (x['tsonic_max'] > 60.0 or x['tsonic_min'] < -50.0 or x['tsonic_max'] - x['tsonic_min'] > 10.0),
    'co2_absolute_limits': lambda x: any(np.isnan([x['co2_max'], x['co2_min']])) or (x['co2_max'] > 950.0 or x['co2_min'] < 550.0 or x['co2_max'] - x['co2_min'] > 120.0),
    'h2o_absolute_limits': lambda x: any(np.isnan([x['h2o_max'], x['h2o_min']])) or (x['h2o_max'] > 35.0 or x['h2o_min'] < 2.5 or x['h2o_max'] - x[f'h2o_min'] > 8.0),
    'nonstationary': lambda x: any(np.isnan([x['rnu'], x['rnv'], x['rns'], x['wind_speed_reduction']])) or (
        x['rnu'] > QC_THRESHOLDS['relative_nonstationarity'] or x['rnv'] > QC_THRESHOLDS['relative_nonstationarity'] or
        x['rns'] > QC_THRESHOLDS['relative_nonstationarity'] or x['wind_speed_reduction'] < QC_THRESHOLDS['wind_speed_reduction']
    ),
}

def calculate_flags(diagnostics: dict) -> dict:
    flags = {}
    for v in INSTANTANEOUS_VARIABLES:
        flags[f'{v.name}_flag'] = FLAGS['spikes'](v.name, diagnostics) * 2
        flags[f'{v.name}_flag'] += FLAGS['resolution'](v.name, diagnostics) * 4
        flags[f'{v.name}_flag'] += FLAGS['dropouts'](v.name, diagnostics) * 8
        flags[f'{v.name}_flag'] += FLAGS[f'{v.name}_absolute_limits'](diagnostics) * 16
        flags[f'{v.name}_flag'] += FLAGS['higher_moments'](v.name, diagnostics) * 32
        flags[f'{v.name}_flag'] += FLAGS['discontinuities'](v.name, diagnostics) * 64

        if v.name in ['u', 'v', 'w']:
            flags[f'{v.name}_flag'] += FLAGS['nonstationary'](diagnostics) * 128
        if v.name in ['u', 'v', 'w', 'tsonic']:
            flags[f'{v.name}_flag'] += FLAGS['3dua'](diagnostics) * 1
        if v.name in ['co2', 'h2o']:
            flags[f'{v.name}_flag'] += FLAGS['irga'](diagnostics) * 1

    return flags


def write_flux_file(diag_file: str) -> None:
    # Read diagnostic file
    df = pd.read_csv(
        diag_file,
        na_values=BADVAL,
    )

    # Calculate diagnostic flags based on QC thresholds
    df[[f'{v.name}_flag' for v in INSTANTANEOUS_VARIABLES]] = df.apply(lambda x: calculate_flags(x), axis=1, result_type='expand')

    for v in OUTPUT_VARIABLES:
        if not v.flag: continue

        df[v.name] = df.apply(
            lambda x: np.nan if x[v.flag].sum(skipna=False)!= 0 else v.format % x[v.name],
            axis=1,
        )

    df.replace(np.nan, str(BADVAL), inplace=True)

    # Write to report file
    flux_file = diag_file.replace('_diag.csv', '.csv')
    df[['TIMESTAMP_START', 'TIMESTAMP_END', *[v.name for v in OUTPUT_VARIABLES]]].to_csv(flux_file, index=False)

    # Write to diagnostic flag file
    flag_file = diag_file.replace('_diag.csv', '_flag.csv')
    df[['TIMESTAMP_START', 'TIMESTAMP_END', *[f'{v.name}_flag' for v in INSTANTANEOUS_VARIABLES]]].astype(int).to_csv(flag_file, index=False)


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

    if AVERAGING_PERIOD_MINUTES == 30:
        resolution = 'HH'
    elif AVERAGING_PERIOD_MINUTES == 60:
        resolution = 'HR'
    else:
        raise ValueError('Please use either 30 min or 60 min for averaging period.')

    write_flux_file(DIAG_FILE(resolution, start_of_month, end_of_month))


if __name__ == '__main__':
    _main()
