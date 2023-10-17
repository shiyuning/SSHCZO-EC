#!/usr/bin/env python3
import math
import numpy as np
import pandas as pd
from scipy.stats import kurtosis, skew
from itertools import groupby
from operator import itemgetter
from site_parameters import *

# Thresholds for quality control
QC_THRESHOLDS = {
    'instrument': 0.95,
    'spike': 0.03,                      # Vickers and Mahrt 1997 value 0.01
    'empty_bins': 0.7,                  # Vickers and Mahrt 1997 value 0.7
    'dropouts': 0.1,                    # Vickers and Mahrt 1997 value 0.1
    'extreme_dropouts': 0.06,           # Vickers and Mahrt 1997 value 0.06
    'skewness': 3.0,                    # Vickers and Mahrt 1997 value 2.0
    'kurtosis': 3.5,                    # Vickers and Mahrt 1997 value 3.5
    'discontinuities': 3.0,             # Vickers and Mahrt 1997 value 3.0
    'wind_speed_reduction': 0.9,        # Vickers and Mahrt 1997 value 0.9
    'relative_nonstationarity': 0.5     # Vickers and Mahrt 1997 value 0.5

}
SECONDS_IN_MINUTE = 60.0

def quality_control(unit_i, unit_j, unit_k, df):
    """Quality control for eddy covariance flux data
    Vickers, D., and L. Mahrt, 1997: Quality control and flux sampling problems for tower and aircraft data.
    J. Atmos. Oceanic tech., 14, 512-526
    """
    flags = {}

    instrument(df, flags)

    print('%-32s' % '  Diagnostics', end='')
    for var in ['u', 'v', 'w', 'Ts', 'CO2', 'H2O']: print('%20s' % var, end='')
    print()

    print('%-32s' % '  Spike fraction', end='')
    for var in [U, V, W, TS, CO2, H2O]: spikes(df, var, flags)
    print('')

    # Rotate to the natural wind coordinate system after de-spiking
    u = df[U] * unit_i[0] + df[V] * unit_i[1] + df[W] * unit_i[2]
    v = df[U] * unit_j[0] + df[V] * unit_j[1] + df[W] * unit_j[2]
    w = df[U] * unit_k[0] + df[V] * unit_k[1] + df[W] * unit_k[2]
    [df[U], df[V], df[W]] = [u, v, w]

    print('%-32s' % '  Empty bin ratio, dropouts', end='')
    for var in [U, V, W, TS, CO2, H2O]: amplitude_resolution_dropouts(df, var, flags)
    print()

    print('%-32s' % '  Range', end='')
    absolute_limits(df, flags)
    print()

    print('%-32s' % '  Skewness, kurtosis', end='')
    for var in [U, V, W, TS, CO2, H2O]:
        detrend(df, var)
        higher_moment_statistics(df, var, flags)
    print()

    print('%-32s' % '  Normalized Harr mean, variance', end='')
    for var in [U, V, W, TS, CO2, H2O]: discontinuities(df, var, flags)
    print('\n')

    nonstationary(df, flags)
    print()

    first = True
    if any(flags.values()):
        print('  Flags: ', end='')
        for flag in flags:
            if flags[flag] == 1:
                print(f'{flag}' if first else f', {flag}', end='')
                first = False
        print('\n')

    return flags


def instrument(df, flags):
    """Instrument diagnostics
    If available records for current time period is less than a threshold, or is more than 18000, a flag is placed
    """

    ratio = float(len(df)) / (AVERAGING_PERIOD_MINUTES * SECONDS_IN_MINUTE * FREQUENCY_HZ)
    print(f'  Available data fraction: {ratio:.3f}\n')
    flags['instrument'] = 1 if (ratio < QC_THRESHOLDS['instrument'] or ratio > 1.0) else 0


def spikes(df, var, flags):
    """Detect and remove spikes
    The method computes the mean and standard deviation for a series of moving windows of length L1. The window moves
    one point at a time through the series. Any point in the window that is more than 3.5 standard deviations from the
    window mean is considered a spike. The point is replaced using linear interpolation between data points. When four
    or more consecutive points are detected, they are not considered spikes and are not replaced. The entire process is
    repeated until no more spikes are detected. During the second pass, when the standard deviations may be smaller if
    spikes were replaced on the previous pass, the threshold for spike detection increases to 3.6 standard deviations
    and a like amount for each subsequent pass. The record is hard flagged when the total number of spikes replaced
    exceeds 1% of the total number of data points.
    """
    CONSECUTIVE_SPIKES = 3
    L1_SECONDS = 300
    window_position = 0
    std_threshold = 4.5
    threshold_increment = 0.1
    pass_spikes = []
    all_spikes = []

    L1 = min(int(L1_SECONDS * FREQUENCY_HZ), len(df))

    while True:
        while window_position + L1 <= len(df):
            sub_df = df.loc[window_position:window_position + L1, :]
            pass_spikes += sub_df.index[abs(sub_df[var] - np.mean(sub_df[var])) > std_threshold * np.std(sub_df[var])].tolist()
            window_position += 1

        pass_spikes = np.sort(np.unique(pass_spikes))

        if len(pass_spikes) > 0:
            ## Filter out consecutive spikes
            consecutive = [list(map(itemgetter(1), g)) for _, g in groupby(enumerate(pass_spikes), lambda x: x[0] - x[1])]
            for c in consecutive:
                if len(c) > CONSECUTIVE_SPIKES:
                    for ind in c: pass_spikes = np.delete(pass_spikes, np.argwhere(pass_spikes == ind))

        if len(pass_spikes) > 0:
            all_spikes += pass_spikes.tolist()

            # Replace spikes with linear interpolation
            spike_filter = np.full(len(df), False)
            spike_filter[pass_spikes] = True
            df.loc[spike_filter, var] = np.interp(df.loc[spike_filter, TIME], df.loc[~spike_filter, TIME], df.loc[~spike_filter, var])

        if len(pass_spikes) == 0:
            spike_ratio = float(len(np.unique(all_spikes))) / float(len(df))
            flags[f'{var}_spikes'] = 1 if (spike_ratio > QC_THRESHOLDS['spike']) else 0
            print('%20.3f' % (spike_ratio), end='')
            return
        else:
            window_position = 0
            std_threshold += threshold_increment
            threshold_increment += 0.1
            pass_spikes = []


def amplitude_resolution_dropouts(df, var, flags):
    """Detect resolution problems and dropouts
    An amplitude resolution problem is detected by computing a series of discrete frequency distributions for half-
    overlapping windows of length 1000 data points. These windows move one-half the window width at a time through the
    series. For each window position, the number of bins is set to 100 and the interval for the distribution is taken as
    the smaller of seven standard deviations and the range. When the number of empty bins in the discrete frequency
    distribution exceeds a critical threshold value, the record is hard flagged as a resolution problem.

    Dropouts are identified using the same window and frequency distributions used for the resolution problem.
    Consecutive points that fall into the same bin of the frequency distribution are tentatively identified as dropouts.
    When the total number of dropouts in the record exceeds a threshold value, the record is flagged for dropouts.
    """
    L1 = 1000
    N_BINS = 100.0

    L1 = min(L1, len(df))
    empty_bins = 0
    window_position = 0
    dropout_ratio = 0
    extreme_dropout_ratio = 0
    window_position = 0

    while (window_position + L1 <= len(df)):
        sub_df = df.loc[window_position:window_position + L1, var]
        window_mean = np.nanmean(sub_df)
        window_std = np.nanstd(sub_df)
        distribution = min(7.0 * window_std, np.ptp(sub_df))
        bin_edges = np.arange(
            window_mean - distribution / 2.0,
            window_mean + distribution / 2.0 + distribution / N_BINS,
            distribution / N_BINS,
        )

        hist, _ = np.histogram(sub_df, bins=bin_edges)
        empty_bins = max(empty_bins, float(np.count_nonzero(hist==0)) / N_BINS)

        #https://stackoverflow.com/questions/6352425/whats-the-most-pythonic-way-to-identify-consecutive-duplicates-in-a-list
        bins = np.digitize(sub_df, bin_edges)
        grouped_bins = [(k, sum(1 for _ in g)) for k, g in groupby(bins)]

        max_dropouts = max(grouped_bins, key=lambda x: x[1])
        dropout_ratio = max(dropout_ratio, float(max_dropouts[1]) / float(len(df)))
        if  (max_dropouts[0] < 10 or max_dropouts[0] > 90):
            extreme_dropout_ratio = max(extreme_dropout_ratio, float(max_dropouts[1]) / float(len(df)))

        window_position += L1 / 2

    print('%20s' %(f'{empty_bins:.3f}, {dropout_ratio:.3f}({extreme_dropout_ratio:.3f})'), end='')
    flags[f'{var}_resolution'] = 1 if (empty_bins > QC_THRESHOLDS['empty_bins'])  else 0
    flags[f'{var}_dropouts'] = 1 if (dropout_ratio > QC_THRESHOLDS['dropouts'] or extreme_dropout_ratio > QC_THRESHOLDS['extreme_dropouts']) else 0


def absolute_limits(df, flags):
    """Detect unrealistic data out of their physical ranges
    Unrealistic data are detected and hard flagged by simply comparing the minimum and maximum value of all points in
    the record to some fixed limits considered unphysical.
    """
    print('%20s' %(f'[{df[U].min():.3f}, {df[U].max():.3f}]'), end='')
    flags[f'{U}_absolute_limits'] = 1 if any(abs(df[U]) > 30.0) else 0

    print('%20s' %(f'[{df[V].min():.3f}, {df[V].max():.3f}]'), end='')
    flags[f'{V}_absolute_limits'] = 1 if any(abs(df[V]) > 30.0) else 0

    print('%20s' %(f'[{df[W].min():.3f}, {df[W].max():.3f}]'), end='')
    flags[f'{W}_absolute_limits'] = 1 if any(abs(df[W]) > 10.0) else 0

    print('%20s' %(f'[{df[TS].min():.2f}, {df[TS].max():.2f}]'), end='')
    flags[f'{TS}_absolute_limits'] = 1 if any(df[TS] > 60.0) or any(df[TS] < -50.0) or np.ptp(df[TS]) > 10.0 else 0

    print('%20s' %(f'[{df[CO2].min():.2f}, {df[CO2].max():.2f}]'), end='')
    flags[f'{CO2}_absolute_limits'] = 1 if any(df[CO2] > 950.0) or any(df[CO2] < 550.0) or np.ptp(df[CO2]) > 120.0 else 0

    print('%20s' %(f'[{df[H2O].min():.2f}, {df[H2O].max():.2f}]'), end='')
    flags[f'{H2O}_absolute_limits'] = 1 if any(df[H2O] > 35.0) or any(df[H2O] < 2.5) or np.ptp(df[H2O]) > 8.0 else 0


def detrend(df, var):
    """Linear de-trend
    """
    time_array = df[TIME].values.astype(float) / 1.0E9 # Convert to seconds
    time_array -= time_array[0]

    s = np.polyfit(time_array, df[var], 1)

    trend = s[0] * time_array + s[1]
    df[f'{var}_fluct'] = df[var] - trend


def higher_moment_statistics(df, var, flags):
    """Detect possible instrument or recording problems and physical but unusual behavior using higher-moment statistics
    The skewness and kurtosis of the fields are computed for the entire record. The record is hard flagged when the
    skewness is outside the range (-2, 2) or the kurtosis is outside the range (1, 8).
    """
    m3 = skew(df[f'{var}_fluct'])
    m4 = kurtosis(df[f'{var}_fluct'], fisher=False)

    print('%20s' %(f'{m3:.3f}, {m4:.3f}'), end='')

    # The range of kurtosis is (1, 8) in Vickers and Mahrt (1997). It is equivalent to having (kurtosis - 4.5) in the
    # range (-3.5, 3.5) in this implementation.
    flags[f'{var}_higher_moments'] = 1 if (abs(m3) > QC_THRESHOLDS['skewness'] or abs(m4 - 4.5) > QC_THRESHOLDS['kurtosis']) else 0


def discontinuities(df, var, flags):
    """Detect discontinuities in the data using the Haar transform
    The transform is computed for a series of moving windows of width L1 and then normalized by the smaller of the
    standard deviation for the entire record and one-fourth the range for the entire record. The record is hard flagged
    if the absolute value of any single normalized transform exceeds 3 and soft flagged at 2. To identify coherent
    changes over the window width L1 in the intensity of the fluctuations, we compute the variance for each half-window
    and then compute the difference normalized by the variance over the entire record. The record is hard flagged if the
    absolute value of any single normalized transform exceeds 3 and soft flagged at 2.
    """
    L1_SECONDS = 300
    L1 = min(int(L1_SECONDS * FREQUENCY_HZ), len(df))

    record_std = np.std(df[var])
    record_range = np.ptp(df[var])
    normal = min(record_std, record_range / 4.0)

    window_position = 0
    haar_mean_max = 0
    haar_variance_max = 0

    while (window_position + L1 <= len(df)):
        half_point = int(L1 / 2)

        half1_mean = np.mean(df.loc[0:half_point, var])
        half2_mean = np.mean(df.loc[half_point:L1, var])

        half1_variance = np.var(df.loc[0:half_point, var])
        half2_variance = np.var(df.loc[half_point:L1, var])

        haar_mean_max = max(haar_mean_max, abs((half2_mean - half1_mean) / normal))
        haar_variance_max = max(haar_variance_max, abs((half2_variance - half1_variance) / (record_std * record_std)))

        window_position += L1 / 4

    print('%20s' % (f'{haar_mean_max:.3f}, {haar_variance_max:.3f}'), end='')
    flags[f'{var}_discontinuities'] = 1 if (haar_mean_max > QC_THRESHOLDS['discontinuities'] or haar_variance_max > QC_THRESHOLDS['discontinuities']) else 0


def nonstationary(df, flags):
    """Identify nonstaionarity of horizontal wind
    The wind speed reduction is defined as the ratio of the speed of the vector averaged wind to the averaged
    instantaneous speed. When this ratio falls below 0.9, there is some cancellation in the vector average of the wind
    components and a soft flag is raised.
    The alongwind relative nonstationarity is calculated using linear regression to estimate the difference in the
    alongwind component between the beginning and end of the record. This difference normalized by the record mean of
    the alongwind component is used to compute the relative nonstationarity.
    The crosswind relative nonstationarity is computed from the difference based on the regression of the crosswind
    component.
    The flow is classified as nonstationary if RNu, RNv, or RNS > 0.50
    """
    # Wind speed reduction
    vector_average = math.sqrt(df[U].mean() * df[U].mean() + df[V].mean() * df[V].mean())
    instant_average = np.mean(np.sqrt(df[U] * df[U] + df[V] * df[V]))
    wind_speed_reduction = vector_average / instant_average

    print(f'  Wind speed reduction: {wind_speed_reduction:.3f}')

    # Relative nonstationarity
    time_array = df[TIME].values.astype(float) / 1.0E9 # Convert to seconds
    time_array -= time_array[0]
    su = np.polyfit(time_array, df[U], 1)
    du = su[0] * (time_array[-1] - time_array[0])
    sv = np.polyfit(time_array, df[V], 1)
    dv = sv[0] * (time_array[-1] - time_array[0])

    rnu = du / np.mean(df[U])
    rnv = dv / np.mean(df[U])
    rns = math.sqrt(du * du + dv * dv) / np.mean(df[U])

    print(f'  Alongwind relative nonstationarity: {rnu:.3f}')
    print(f'  Crosswind relative nonstationarity: {rnv:.3f}')
    print(f'  Relative nonstationarity: {rns:.3f}')

    flags['nonstationary'] = 1 if (
        rnu > QC_THRESHOLDS['relative_nonstationarity'] or
        rnv > QC_THRESHOLDS['relative_nonstationarity'] or
        rns > QC_THRESHOLDS['relative_nonstationarity'] or
        wind_speed_reduction < QC_THRESHOLDS['wind_speed_reduction']) else 0
